#include "OpenSim/Analyses/ForceReporter.h"
#include "OpenSim/Analyses/StatesReporter.h"
#include "OpenSim/Common/Reporter.h"
#include "OpenSim/ExampleComponents/RegisterTypes_osimExampleComponents.h"
#include "OpenSim/Simulation/Model/AbstractGeometryPath.h"
#include "OpenSim/Simulation/Model/AbstractTool.h"
#include "OpenSim/Simulation/Model/GeometryPath.h"
#include "OpenSim/Simulation/Model/ModelComponent.h"
#include "OpenSim/Simulation/Model/Muscle.h"
#include "OpenSim/Simulation/Model/MyCable.h"
#include "OpenSim/Simulation/Model/PathActuator.h"
#include "OpenSim/Simulation/Model/PhysicalFrame.h"
#include "OpenSim/Simulation/RegisterTypes_osimSimulation.h"
#include "OpenSim/Simulation/StatesTrajectoryReporter.h"
#include "OpenSim/Simulation/Wrap/WrapObject.h"
#include "OpenSim/Tools/CMCTool.h"
#include "OpenSim/Tools/ForwardTool.h"
#include "OpenSim/Tools/Tool.h"
#include <chrono>
#include <exception>
#include <memory>
#include <ostream>
#include <simmath/internal/ContactGeometry.h>
#include <stdexcept>
#include <unordered_map>

#include <SimTKcommon/Scalar.h>
#include <SimTKcommon/internal/ExceptionMacros.h>
#include <SimTKcommon/internal/NTraits.h>
#include <SimTKcommon/internal/Transform.h>
#include <SimTKcommon/internal/UnitVec.h>

#include <simbody/internal/common.h>

#include <OpenSim/OpenSim.h>

using namespace OpenSim;
using namespace SimTK;
using Data = AbstractGeometryPath::Data;

struct Timer {
    Timer() = default;

    void tick() { m_Tick = std::chrono::high_resolution_clock::now(); }

    double tock() {
        std::chrono::time_point<std::chrono::high_resolution_clock> tock =
                std::chrono::high_resolution_clock::now();
        return std::chrono::duration<double>(tock - m_Tick).count();
    }

    double measureDuration(const std::function<void()>& f) {
        tick();
        f();
        double dt = tock();
        m_Tick = std::chrono::time_point<std::chrono::high_resolution_clock>::max();
        return dt;
    }

private:
    std::chrono::time_point<std::chrono::high_resolution_clock> m_Tick;
};

Real runAndTime(
        const std::string& modelFile, bool visualize, SimTK::Real finalTime) {
    Model model(modelFile);
    if (visualize) { model.setUseVisualizer(true); }
    Timer timer;
    SimTK::State& s = model.initSystem();
    s.setTime(0.);
    model.realizeReport(s);

    Manager manager(model);
    /* manager.setIntegratorMaximumStepSize(1e-3); */
    manager.initialize(s);

    const Real dt =
            timer.measureDuration([&]() { manager.integrate(finalTime); });

    std::cout << "sim time = " << s.getTime() << "\n";
    std::cout << "steps    = " << manager.getIntegrator().getNumStepsTaken() << "\n";
    std::cout << "attempted    = " << manager.getIntegrator().getNumStepsAttempted() << "\n";
    std::cout << "dt       = " << dt << "\n";

    return dt;
}

template <typename T>
Real runTool(const std::string& modelFile, const std::string& toolFile,
        bool visualize) {

    Model model(modelFile);
    if (visualize) { model.setUseVisualizer(true); }

    T tool(toolFile);
    tool.setModel(model);

    Timer timer;

    const Real dt = timer.measureDuration([&]() { tool.run(); });

    return dt;
}

Real runToolFWD(const std::string& modelFile, const std::string& toolFile,
        bool visualize) {

    Model model(modelFile);
    if (visualize) { model.setUseVisualizer(true); }

    ForwardTool tool(toolFile);
    tool.setModel(model);
    tool.setMaxDT(1e-3);

    Timer timer;

    const Real dt = timer.measureDuration([&]() { tool.run(); });

    return dt;
}

template <typename T> void printType(const Model& model) {
    std::cout << "LIST TYPE" << T::getClassName() << std::endl;
    for (const T& c : model.getComponentList<T>()) {
        std::cout << c.getClassName() << " : " << c.getName() << " @ "
                  << c.getAbsolutePathString() << std::endl;
    }
}

void convertModel(const std::string& file, const std::string& outFile) {

    std::cout << "------------------------------" << std::endl;
    std::cout << "------STARTING CONVERSION-----\n";
    std::cout << "-----------------------------" << std::endl;

    Model model = MyCable::convert(
            file);

    std::cout << "************** DONE CONVERSION ***********" << std::endl;
    std::cout << "Writing converted model to: " << outFile << std::endl;
    model.print(outFile);

    std::cout << "Verify reading the model from: " << outFile << std::endl;
    Model check(outFile);
    printType<GeometryPath>(check);
    printType<MyCable>(check);
    printType<PathActuator>(check);

    std::cout << "----------------" << std::endl;
    std::cout << "COMPLETED CONVERSION\n";
    std::cout << "----------------" << std::endl;
}

int main(int argc, char* argv[]) {
    std::vector<std::string> args(argv + 1, argv + argc);
    // Argument options.
    std::string cmc;
    std::string fwd;
    std::string model;
    std::string modelOut;
    bool run = false;
    SimTK::Real simTime = SimTK::NaN;
    bool visualize = false;
    // Parse the arguments.
    for (auto it = args.begin(); it < args.end(); it++) {
        if (*it == "--model") { model = *++it; }
        if (*it == "--out") { modelOut = *++it; }
        if (*it == "--cmc") { cmc = *++it; }
        if (*it == "--fwd") { fwd = *++it; }
        run |= *it == "--run";
        if (*it == "--run") { simTime = std::stod(*++it); }
        visualize |= *it == "--viz";
    }
    // Check parsed arguments.
    SimTK_ASSERT_ALWAYS(!model.empty(), "No model selected");

    std::cout << "model = " << model << "\n";
    std::cout << "out   = " << modelOut << "\n";
    std::cout << "cmc   = " << cmc << "\n";
    std::cout << "fwd   = " << fwd << "\n";
    std::cout << "run   = " << run << "\n";
    std::cout << "sim   = " << simTime << "\n";
    std::cout << "viz   = " << visualize << "\n";

    RegisterTypes_osimExampleComponents();
    RegisterTypes_osimSimulation();

    if (run) {
        SimTK::Real dt = runAndTime(model, visualize, simTime);

        std::cout << "RESULTS:\n";
        std::cout << "--> dt = " << dt << "\n";
        return 0;
    }

    if (!cmc.empty()) {
        SimTK::Real dt = runTool<CMCTool>(model, cmc, visualize);

        std::cout << "RESULTS:\n";
        std::cout << "--> dt = " << dt << "\n";
        return 0;
    }

    if (!fwd.empty()) {
        SimTK::Real dt = runToolFWD(model, fwd, visualize);

        std::cout << "RESULTS:\n";
        std::cout << "--> dt = " << dt << "\n";
        return 0;
    }

    if (!modelOut.empty()) {
        convertModel(model, modelOut);
        return 0;
    }

    return 0;
}
