#include <OpenSim/OpenSim.h>
#include "FunctionBasedPathConversionTool.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdio>
#include <cstddef>
#include <sstream>

using namespace OpenSim;

OpenSim::FunctionBasedPathConversionTool::FunctionBasedPathConversionTool() :
    _modelPath{},
    _newModelName{},
    _verbose{false}
{
}

OpenSim::FunctionBasedPathConversionTool::FunctionBasedPathConversionTool(
        const std::string& modelPath,
        const std::string& newModelName) :

    _modelPath{modelPath},
    _newModelName{newModelName},
    _verbose{false}
{
}

const std::string& OpenSim::FunctionBasedPathConversionTool::getModelPath() const
{
     return _modelPath;
}

void OpenSim::FunctionBasedPathConversionTool::setModelPath(const std::string& modelPath)
{
    _modelPath = modelPath;
}

const std::string& OpenSim::FunctionBasedPathConversionTool::getNewModelName() const
{
    return _newModelName;
}

void OpenSim::FunctionBasedPathConversionTool::setNewModelName(const std::string& newModelName)
{
    _newModelName = newModelName;
}

bool OpenSim::FunctionBasedPathConversionTool::run()
{
    Model sourceModel{_modelPath};
    sourceModel.finalizeConnections();
    sourceModel.finalizeFromProperties();
    sourceModel.initSystem();

    Model outputModel{sourceModel};
    outputModel.finalizeConnections();
    outputModel.finalizeFromProperties();

    // print source model structure when runnning in verbose mode
    if (_verbose) {
        sourceModel.printSubcomponentInfo();
    }

    // struct that holds how a PBP in the source maps onto an actuator in the
    // destination
    struct PBPtoActuatorMapping final {
        PointBasedPath& sourcePBP;
        PathActuator& outputActuator;

        PBPtoActuatorMapping(PointBasedPath& sourcePBP_,
                             PathActuator& outputActuator_) :
            sourcePBP{sourcePBP_},
            outputActuator{outputActuator_} {
        }
    };

    // find PBPs in the source and figure out how they map onto `PathActuator`s
    // in the destination
    //
    // (this is because `PathActuator`s are the "owners" of `GeometryPath`s in
    //  most models)
    std::vector<PBPtoActuatorMapping> mappings;
    for (PathActuator& pa : sourceModel.updComponentList<PathActuator>()) {

        PointBasedPath* pbp = dynamic_cast<PointBasedPath*>(&pa.updGeometryPath());

        // if the actuator doesn't use a PBP, ignore it
        if (!pbp) {
            continue;
        }

        // otherwise, find the equivalent path in the destination
        Component* c = const_cast<Component*>(outputModel.findComponent(pa.getAbsolutePath()));

        if (!c) {
            std::stringstream err;
            err << "could not find '" << pa.getAbsolutePathString() << "' in the output model: this is a programming error";
            throw std::runtime_error{move(err).str()};
        }

        PathActuator* paDest = dynamic_cast<PathActuator*>(c);

        if (!paDest) {
            std::stringstream err;
            err << "the component '" << pa.getAbsolutePathString() << "' has a class of '" << pa.getConcreteClassName() << "' in the source model but a class of '" << c->getConcreteClassName() << "' in the destination model: this shouldn't happen";
            throw std::runtime_error{move(err).str()};
        }

        mappings.emplace_back(*pbp, *paDest);
    }

    // for each `PathActuator` that uses a PBP, create an equivalent
    // `FunctionBasedPath` (FBP) by fitting a function against the PBP and
    // replace the PBP owned by the destination's `PathActuator` with the FBP
    for (const PBPtoActuatorMapping& mapping : mappings) {
        // create an FBP in-memory
        FunctionBasedPath::FittingParams p;
        std::unique_ptr<FunctionBasedPath> maybeFbp = FunctionBasedPath::fromPointBasedPath(sourceModel, mapping.sourcePBP, p);

        if (maybeFbp) {
            // assign the FBP over the destination's PBP
            mapping.outputActuator.updProperty_GeometryPath().setValue(*maybeFbp);
        }
    }

    // the output model is now the same as the source model, but each PBP in
    // its `PathActuator`s has been replaced with an FBP. Perform any final
    // model-level fixups and save the output model.
    outputModel.finalizeFromProperties();
    outputModel.finalizeConnections();
    outputModel.initSystem();
    outputModel.print(std::string{_newModelName} + ".osim");

    if (_verbose) {
        std::cerr << "--- interpolation complete ---\n\n"
                  << "model before:\n";
        sourceModel.printSubcomponentInfo();
        std::cerr << "\nmodel after:\n";
        outputModel.printSubcomponentInfo();
    }

    return true;
}
