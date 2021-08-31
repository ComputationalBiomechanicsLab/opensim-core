/* -------------------------------------------------------------------------- *
 *                     OpenSim:  testFunctionBasedPathConversion.cpp          *
 * -------------------------------------------------------------------------- *
 * ScapulothoracicJoint is offered as an addition to the OpenSim API          *
 * See http://opensim.stanford.edu and the NOTICE file for more information.  *
 *                                                                            *
 * OpenSim is developed at Stanford University and is supported by:           *
 *                                                                            *
 * - The National Institutes of Health (U54 GM072970, R24 HD065690)           *
 * - DARPA, through the Warrior Web program                                   *
 * - The Chan Zuckerberg Initiative (CZI 2020-218896)                         *
 *                                                                            *
 * Copyright (c) 2005-2021 Stanford University, TU Delft, and the Authors     *
 * Author(s): Joris Verhagen                                                  *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */

#include <OpenSim/Common/Reporter.h>
#include <OpenSim/Simulation/Model/Model.h>
#include <OpenSim/Simulation/Manager/Manager.h>
#include <OpenSim/Tools/FunctionBasedPathConversionTool.h>
#include <SimTKcommon.h>

#include <string>
#include <iostream>
#include <memory>
#include <chrono>
#include <cstdio>

using namespace OpenSim;

void testArmModelConversionAccuracy() {
    std::string inputModelPath = "arm26.osim";
    std::string outputModelName = "arm26_FBP";

    // run the function-based path (FBP) conversion tool to create an FBP-based output
    FunctionBasedPathConversionTool tool{inputModelPath, outputModelName};
    tool.run();

    // load both the input and output into memory
    Model inputModel{inputModelPath};
    Model outputModel{outputModelName + ".osim"};

    // init the input model
    inputModel.finalizeConnections();
    inputModel.finalizeFromProperties();
    inputModel.initSystem();

    std::cout << "----- input model details -----\n";
    inputModel.printSubcomponentInfo();

    // init the output model
    outputModel.finalizeConnections();
    outputModel.finalizeFromProperties();
    outputModel.initSystem();

    std::cout << "----- output model details -----\n";
    outputModel.printSubcomponentInfo();

    // connect reporters to each model
    std::string force = "/forceset/TRIlong";
    std::string output = "length";

    // connect reporter to input model
    {
        auto inputModelReporter = std::unique_ptr<ConsoleReporter>{new ConsoleReporter{}};
        inputModelReporter->setName("point_results");
        inputModelReporter->set_report_time_interval(0.05);
        inputModelReporter->addToReport(inputModel.getComponent(force + "/pointbasedpath").getOutput(output));
        inputModel.addComponent(inputModelReporter.release());
    }

    // connect reporter to output model
    {
        auto outputModelReporter = std::unique_ptr<ConsoleReporter>{new ConsoleReporter{}};
        outputModelReporter ->setName("function_results");
        outputModelReporter ->set_report_time_interval(0.05);
        outputModelReporter ->addToReport(outputModel.getComponent(force + "/functionbasedpath").getOutput(output));
        outputModel.addComponent(outputModelReporter.release());
    }

    // init physical system + states
    SimTK::State& inputModelState = inputModel.initSystem();
    SimTK::State& outputModelState = outputModel.initSystem();

    // these are accumulated by each test run
    struct TestStats {
        std::chrono::high_resolution_clock::duration simTime{0};
        int stepsAttempted = 0;
        int stepsTaken = 0;
        std::string name;

        TestStats(std::string name_) : name{name_} {}
    };

    // function that runs a single test
    auto runTest = [](OpenSim::Model& model, SimTK::State& state, TestStats& stats, double finalSimTime) {
        OpenSim::Manager manager{model};
        manager.initialize(state);

        auto before = std::chrono::high_resolution_clock::now();
        manager.integrate(finalSimTime);
        auto after = std::chrono::high_resolution_clock::now();

        stats.simTime += after - before;
        stats.stepsAttempted += manager.getIntegrator().getNumStepsAttempted();
        stats.stepsTaken += manager.getIntegrator().getNumStepsTaken();
    };

    // setup + run the tests
    TestStats inputStats{"Input Model (point-based paths)"};
    TestStats outputStats{"Output Model (function-based paths)"};
    int numRepeats = 10;
    double finalSimTime = 3.5;

    // exercise input model
    std::cerr << "----- running input model (point based path) tests -----\n";
    for (int i = 0; i < numRepeats; ++i) {
        runTest(inputModel, inputModelState, inputStats, finalSimTime);
    }

    // exercise output model
    std::cerr << "----- running output model (function based path) tests -----\n";
    for (int i = 0; i < numRepeats; ++i) {
        runTest(outputModel, outputModelState, outputStats, finalSimTime);
    }

    // average out and print the stats
    for (auto ts : {inputStats, outputStats}) {
        ts.simTime /= numRepeats;
        ts.stepsAttempted /= numRepeats;
        ts.stepsTaken /= numRepeats;

        long millis = static_cast<long>(std::chrono::duration_cast<std::chrono::milliseconds>(ts.simTime).count());

        printf("%s: \n    avg. time (ms) = %ld \n    integration steps attempted = %i \n    integration steps taken = %i)",
               ts.name.c_str(),
               millis,
               ts.stepsAttempted,
               ts.stepsTaken);
    }
}

int main() {
    try {
        SimTK_START_TEST("testFunctionBasedPathConversion");
            SimTK_SUBTEST(testArmModelConversionAccuracy);
        SimTK_END_TEST();
    } catch (const std::exception& ex) {
        OpenSim::log_error(ex.what());
    }
}
