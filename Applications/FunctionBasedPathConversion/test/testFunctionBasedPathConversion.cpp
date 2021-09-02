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

// output reporter that just appends some double output to an output vector
class VectorReporter final : public AbstractReporter {
public:
    struct Datapoint { double time; double v; };
private:
    const Output<double>& m_Source;
    std::vector<Datapoint>& m_Sink;

public:
    VectorReporter(const Output<double>& source,
                   std::vector<Datapoint>& sink) :
        AbstractReporter{},
        m_Source{source},
        m_Sink{sink}
    {
    }

    void implementReport(const SimTK::State& s) const override {
        double t = s.getTime();
        double v = m_Source.getValue(s);

        m_Sink.push_back(Datapoint{t, v});
    }

    VectorReporter* clone() const override {
        return new VectorReporter{*this};
    }

    const std::string& getConcreteClassName() const override {
        static std::string name = "VectorReporter";
        return name;
    }
};

static void testLengthDifferenceBetweenFBPAndPBPIsSmall() {
    std::string inputModelPath = "arm26.osim";
    std::string outputModelName = "arm26_FBP.osim";
    std::string force = "/forceset/TRIlong";
    std::string output = "length";
    double reportingInterval = 0.05;
    int discretizationPoints = 80;

    FunctionBasedPathConversionTool tool{inputModelPath, outputModelName};
    auto params = tool.getFittingParams();
    params.numDiscretizationStepsPerDimension = discretizationPoints;
    tool.setFittingParams(params);
    tool.setVerbose(true);
    tool.run();

    // load + init both models
    Model inputModel{inputModelPath};
    inputModel.finalizeConnections();
    inputModel.finalizeFromProperties();
    inputModel.initSystem();

    Model outputModel{outputModelName};
    outputModel.finalizeConnections();
    outputModel.finalizeFromProperties();
    outputModel.initSystem();

    // connect input model's force output to a vector that records the data
    std::vector<VectorReporter::Datapoint> pbpDatapoints;
    {
        const AbstractOutput& ao = inputModel.getComponent(force + "/pointbasedpath").getOutput(output);
        const Output<double>* o = dynamic_cast<const Output<double>*>(&ao);
        if (!o) {
            throw std::runtime_error{"Input model output is not of `double` type"};
        }
        auto reporter = std::unique_ptr<VectorReporter>(new VectorReporter{*o, pbpDatapoints});
        reporter->set_report_time_interval(reportingInterval);
        inputModel.addComponent(reporter.release());
    }

    // connect output model's force output to a vector that records the data
    std::vector<VectorReporter::Datapoint> fbpDatapoints;
    {
        const AbstractOutput& ao = outputModel.getComponent(force + "/functionbasedpath").getOutput(output);
        const Output<double>* o = dynamic_cast<const Output<double>*>(&ao);
        if (!o) {
            throw std::runtime_error{"Input model output is not of `double` type"};
        }
        auto reporter = std::unique_ptr<VectorReporter>(new VectorReporter{*o, fbpDatapoints});
        reporter->set_report_time_interval(reportingInterval);
        outputModel.addComponent(reporter.release());
    }

    // init initial system + states
    double finalSimTime = 3.5;

    // run FD sim of PBP
    {
        SimTK::State& inputModelState = inputModel.initSystem();
        OpenSim::Manager manager{inputModel};
        manager.initialize(inputModelState);
        manager.integrate(finalSimTime);
    }

    // run FD sim of FBP
    {
        SimTK::State& outputModelState = outputModel.initSystem();
        OpenSim::Manager manager{outputModel};
        manager.initialize(outputModelState);
        manager.integrate(finalSimTime);
    }

    // validate API assumptions
    size_t expectedSteps = static_cast<size_t>(finalSimTime / reportingInterval);
    {
        OPENSIM_THROW_IF(pbpDatapoints.empty(), OpenSim::Exception, "pbpDatapoints empty: is it connected to the model?");
        OPENSIM_THROW_IF(fbpDatapoints.empty(), OpenSim::Exception, "fbpDatapoints empty: is it connected to the model?");
        OPENSIM_THROW_IF(pbpDatapoints.size() == expectedSteps, OpenSim::Exception, "pbpDatapoints has incorrect number of datapoints");
        OPENSIM_THROW_IF(fbpDatapoints.size() == expectedSteps, OpenSim::Exception, "fbpDatapoints has incorrect number of datapoints");

        for (size_t i = 0; i < expectedSteps; ++i) {
            double pbpT = pbpDatapoints[i].time;
            double fbpT = fbpDatapoints[i].time;
            OPENSIM_THROW_IF(!SimTK::isNumericallyEqual(pbpT, fbpT), OpenSim::Exception, "timepoints in PBP and FBP arrays differ");
        }
    }

    struct RelErrStats {
        double min = std::numeric_limits<double>::max();
        double max = std::numeric_limits<double>::min();
        double avg = 0.0;
        int n = 0;
    };
    RelErrStats stats;

    std::cout << "t    \tpbp  \tfbp  \t%err \n";
    for (size_t i = 0; i < expectedSteps; ++i) {
        double t = pbpDatapoints[i].time;
        double pbpV = pbpDatapoints[i].v;
        double fbpV = fbpDatapoints[i].v;
        double relErr = (fbpV-pbpV)/pbpV;
        double relPctErr = 100.0 * relErr;

        stats.min = std::min(stats.min, relErr);
        stats.max = std::max(stats.max, relErr);
        stats.avg = std::fma(stats.avg, stats.n, relErr) / (stats.n+1);
        stats.n += 1;

        std::cout << std::fixed << std::setprecision(5) << t << '\t' << pbpV << '\t' << fbpV << '\t' << relPctErr << " %\n";
    }

    std::cout << "min = " << 100.0*stats.min << " %, max = " << 100.0*stats.max << " %, avg = " << 100.0*stats.avg << " %\n";
}

static void testArmModelConversionAccuracy() {
    std::string inputModelPath = "arm26.osim";
    std::string outputModelName = "arm26_FBP.osim";

    // run the function-based path (FBP) conversion tool to create an FBP-based output
    FunctionBasedPathConversionTool tool{inputModelPath, outputModelName};
    tool.setVerbose(true);
    tool.run();

    // load both the input and output into memory
    Model inputModel{inputModelPath};
    Model outputModel{outputModelName};

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
    if (true) {
        auto inputModelReporter = std::unique_ptr<ConsoleReporter>{new ConsoleReporter{}};
        inputModelReporter->setName("point_results");
        inputModelReporter->set_report_time_interval(0.05);
        inputModelReporter->addToReport(inputModel.getComponent(force + "/pointbasedpath").getOutput(output));
        inputModel.addComponent(inputModelReporter.release());
    }

    // connect reporter to output model
    if (true) {
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

        printf("%s: \n    avg. time (ms) = %ld \n    integration steps attempted = %i \n    integration steps taken = %i)\n",
               ts.name.c_str(),
               millis,
               ts.stepsAttempted,
               ts.stepsTaken);
    }
}

int main() {
    try {
        SimTK_START_TEST("testFunctionBasedPathConversion");
            SimTK_SUBTEST(testLengthDifferenceBetweenFBPAndPBPIsSmall);
            SimTK_SUBTEST(testArmModelConversionAccuracy);
        SimTK_END_TEST();
    } catch (const std::exception& ex) {
        OpenSim::log_error(ex.what());
    }
}
