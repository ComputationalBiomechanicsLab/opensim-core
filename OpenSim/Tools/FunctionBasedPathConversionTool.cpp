#include "FunctionBasedPathConversionTool.h"

#include <OpenSim/Common/Logger.h>
#include <OpenSim/Common/Component.h>
#include <OpenSim/Simulation/Model/Model.h>
#include <OpenSim/Simulation/Model/PathActuator.h>
#include <OpenSim/Simulation/Model/PointBasedPath.h>

#include <memory>
#include <stdexcept>
#include <string>
#include <sstream>
#include <vector>

using namespace OpenSim;

OpenSim::FunctionBasedPathConversionTool::FunctionBasedPathConversionTool() :
    _modelPath{},
    _newModelName{},
    _params{},
    _verbose{false}
{
}

OpenSim::FunctionBasedPathConversionTool::FunctionBasedPathConversionTool(
        const std::string& modelPath,
        const std::string& newModelName) :

    _modelPath{modelPath},
    _newModelName{newModelName},
    _params{},
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

const FunctionBasedPath::FittingParams& OpenSim::FunctionBasedPathConversionTool::getFittingParams() const {
    return _params;
}

void OpenSim::FunctionBasedPathConversionTool::setFittingParams(const FunctionBasedPath::FittingParams& p) {
    _params = p;
}

bool OpenSim::FunctionBasedPathConversionTool::getVerbose() const {
    return _verbose;
}

void OpenSim::FunctionBasedPathConversionTool::setVerbose(bool v) {
    _verbose = v;
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
        OpenSim::log_info("printing source model info");
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
    if (_verbose) {
        OpenSim::log_info("using fitting params: maxCoordsThatCanAffectPath = {}, minProbingMomentArmChange = {}, numDiscretizationStepsPerDimension = {}, numProbingDiscretizations = {}",
                          _params.maxCoordsThatCanAffectPath,
                          _params.minProbingMomentArmChange,
                          _params.numDiscretizationStepsPerDimension,
                          _params.numProbingDiscretizations);
    }
    for (const PBPtoActuatorMapping& mapping : mappings) {
        // create an FBP in-memory

        if (_verbose) {
            OpenSim::log_info("attempting to convert point-based path '{}' into a FunctionBasedPath", mapping.sourcePBP.getAbsolutePathString());
        }

        std::unique_ptr<FunctionBasedPath> maybeFbp = FunctionBasedPath::fromPointBasedPath(sourceModel, mapping.sourcePBP, _params);

        if (_verbose) {
            if (maybeFbp) {
                OpenSim::log_info("successfully converted '{}' into a FunctionBasedPath", mapping.sourcePBP.getAbsolutePathString());
            } else {
                OpenSim::log_info("failed to convert '{}' into a FunctionBasedPath", mapping.sourcePBP.getAbsolutePathString());
            }
        }

        if (maybeFbp) {
            // assign the FBP over the destination's PBP
            mapping.outputActuator.updProperty_GeometryPath().setValue(*maybeFbp);
        }
    }

    if (_verbose) {
        OpenSim::log_info("finished attempting to fit all point based paths: emitting new output model that contains FunctionBasedPaths");
    }

    // the output model is now the same as the source model, but each PBP in
    // its `PathActuator`s has been replaced with an FBP. Perform any final
    // model-level fixups and save the output model.
    outputModel.finalizeFromProperties();
    outputModel.finalizeConnections();
    outputModel.initSystem();
    outputModel.print(_newModelName);

    if (_verbose) {
        OpenSim::log_info("--- interpolation complete ---");
        OpenSim::log_info("model before:");
        sourceModel.printSubcomponentInfo();
        OpenSim::log_info("model after:");
        outputModel.printSubcomponentInfo();
    }

    return true;
}
