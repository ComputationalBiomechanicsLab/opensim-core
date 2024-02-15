#include "OpenSim/Common/SmoothSegmentedCubicMonoSpline.h"
#include <iostream>
#include <memory>

#include <OpenSim/OpenSim.h>

using CtrlPoint = OpenSim::MuscleCurveControlPoint;

void writeDataToCSV(
    const std::vector<OpenSim::MuscleCurveControlPoint>& ctrlPts,
    const OpenSim::CurveShape& shapeCurve,
    const OpenSim::SmoothSegmentedCubicMonoSpline& splineCurve,
    const std::string& filename,
    size_t nSamples = 1000)
{
    static constexpr size_t SPLINE_SHAPE_ID = 0;
    static constexpr size_t SPLINE_PLOT_ID  = 1;
    static constexpr size_t CTRL_PTS_ID     = 2;

    std::ofstream outputFile(filename);

    if (!outputFile.is_open()) {
        std::cerr << "Error opening the file: " << filename << std::endl;
        return;
    }

    // Write header
    outputFile << "Index,X,Y\n";

    // Write data

    // Write data points for each shape segment, and the approximating spline.
    for (const OpenSim::QuadraticBezierCurve& shapeSegment :
         shapeCurve.getSegments()) {

        for (size_t i = 0; i <= nSamples; i++) {
            const double factor =
                static_cast<double>(i) / static_cast<double>(nSamples);
            const auto shapePoint = shapeSegment.calcPoint(factor);
            const double x        = shapePoint.x;
            const double ySpline  = splineCurve.calcValue(x);
            outputFile << SPLINE_PLOT_ID << "," << x << "," << ySpline << "\n";
            outputFile << SPLINE_SHAPE_ID << "," << x << "," << shapePoint.y
                       << "\n";
        }
    }

    // Write the original control points.
    for (const OpenSim::MuscleCurveControlPoint& p : ctrlPts) {
        outputFile << CTRL_PTS_ID << "," << p.x << "," << p.y << "\n";
    }

    outputFile.close();
}

int main()
{

    std::vector<CtrlPoint> ctrlPoints;

    {
        CtrlPoint point;
        point.x    = 0.;
        point.y    = 0.;
        point.dydx = 0.;
        ctrlPoints.push_back(point);
    }

    {
        CtrlPoint point;
        point.x    = 1.;
        point.y    = 1.;
        point.dydx = 5.;
        point.curviness = 0.1;
        ctrlPoints.push_back(point);
    }

    {
        CtrlPoint point;
        point.x    = 3.5;
        point.y    = 2.2;
        point.dydx = 0.1;
        point.curviness = 0.9;
        ctrlPoints.push_back(point);
    }

    /* { */
    /* CtrlPoint point; */
    /*     point.x = 2.; */
    /*     point.y = 1.5; */
    /*     point.dydx = 0.; */
    /*     ctrlPoints.push_back(point); */
    /* } */

    /* { */
    /* CtrlPoint point; */
    /*     point.x = 5.; */
    /*     point.y = 0.1; */
    /*     point.dydx = 0.; */
    /*     ctrlPoints.push_back(point); */
    /* } */

    for (const auto& p : ctrlPoints) {
        std::cout << p << std::endl;
    }

    OpenSim::CurveShape cShape(ctrlPoints);

    OpenSim::SmoothSegmentedCubicMonoSpline sCurve(cShape);

    std::string filename = "customSpline.csv";
    writeDataToCSV(ctrlPoints, cShape, sCurve, filename, 1000);

    return 0;
}

/* int convertAndWriteToCsv( */
/*     const OpenSim::SmoothSegmentedFunction* f, */
/*     const OpenSim::SmoothSpline& spline, */
/*     bool checkInverse = false, */
/*     size_t steps      = 20) */
/* { */
/*     { */
/*         auto eval = [&](double x) -> double { */
/*             double yBezier = f->calcValue(x); */
/*             double ySpline = spline.calcDerivative(x, 0); */
/*             double error   = std::abs(ySpline - yBezier); */

/*             std::cout << "x = " << x << std::endl; */
/*             std::cout << "e = " << error << std::endl; */
/*             std::cout << "yBezier = " << yBezier << std::endl; */
/*             std::cout << "ySpline = " << ySpline << std::endl; */

/*             return error; */
/*         }; */
/*         double x0    = 1.; */
/*         double x1    = 1.01453; */
/*         double error = 0.; */

/*         for (size_t i = 0; i <= steps; i++) { */
/*             const double c = */
/*                 static_cast<double>(i) / static_cast<double>(steps); */
/*             const double x = x0 + (x1 - x0) * c; */
/*             error          = std::max(error, eval(x)); */
/*         } */

/*         std::cout << "max error = " << error << std::endl; */
/*     } */

/*     { */
/*         std::string filename = f->getName() + ".csv"; */
/*         writeDataToCSV(spline.getSplines(), *f, filename, checkInverse); */
/*     } */

/*     { */
/*         std::string filenameInt = f->getName() + "_integral.csv"; */
/*         writeDataIntegralToCSV(spline.getSplines(), filenameInt); */
/*     } */

/*     return 0; */
/* } */

/* int plotAndStore() */
/* { */
/*     { */
/*         OpenSim::SmoothSegmentedFunction* f = OpenSim:: */
/*             SmoothSegmentedFunctionFactory::createTendonForceLengthCurve( */
/*                 0.049, */
/*                 28.1, */
/*                 0.67, */
/*                 0.5, */
/*                 false, */
/*                 "tendoncurve"); */
/*         OpenSim::SmoothSpline spline =
 * OpenSim::SmoothSegmentedFunctionFactory:: */
/*             createTendonForceLengthCurve(0.049, 28.1, 0.67, 0.5,
 * "tendoncurve"); */
/*         convertAndWriteToCsv(f, spline, true); */
/*     } */

/*     { */
/*         OpenSim::SmoothSegmentedFunction* f = OpenSim:: */
/*             SmoothSegmentedFunctionFactory::createFiberForceVelocityCurve(
 */
/*                 1.4, */
/*                 0., */
/*                 0.25, */
/*                 5., */
/*                 0., */
/*                 0.15, */
/*                 0.6, */
/*                 0.9, */
/*                 false, */
/*                 "fiberVelocityCurve"); */
/*         OpenSim::SmoothSpline spline =
 * OpenSim::SmoothSegmentedFunctionFactory:: */
/*             SmoothSegmentedFunctionFactory::createFiberForceVelocityCurve(
 */
/*                 1.4, */
/*                 0., */
/*                 0.25, */
/*                 5., */
/*                 0., */
/*                 0.15, */
/*                 0.6, */
/*                 0.9, */
/*                 "fiberVelocityCurve"); */
/*         convertAndWriteToCsv(f, spline); */
/*     } */

/*     { */
/*         OpenSim::SmoothSegmentedFunction* f = OpenSim:: */
/*             SmoothSegmentedFunctionFactory::createFiberActiveForceLengthCurve(
 */
/*                 0.44, */
/*                 0.73, */
/*                 1.0, */
/*                 1.8123, */
/*                 0.1, */
/*                 0.8616, */
/*                 1.0, */
/*                 false, */
/*                 "ActiveForceLengthCurve"); */
/*         OpenSim::SmoothSpline spline =
 * OpenSim::SmoothSegmentedFunctionFactory:: */
/*             SmoothSegmentedFunctionFactory::createFiberActiveForceLengthCurve(
 */
/*                 0.44, */
/*                 0.73, */
/*                 1.0, */
/*                 1.8123, */
/*                 0.1, */
/*                 0.8616, */
/*                 1.0, */
/*                 "ActiveForceLengthCurve"); */
/*         OPENSIM_ASSERT(f->getCurveDomain()[0] ==
 * spline.getCurveDomain()[0]
 */
/*                 && "curve domain start does not match"); */
/*         OPENSIM_ASSERT(f->getCurveDomain()[1] ==
 * spline.getCurveDomain()[1]
 */
/*                 && "curve domain end does not match"); */
/*         convertAndWriteToCsv(f, spline); */
/*     } */

/*     { */
/*         OpenSim::SmoothSegmentedFunction* f = OpenSim:: */
/*             SmoothSegmentedFunctionFactory::createFiberForceLengthCurve(
 */
/*                 0.0, */
/*                 0.7, */
/*                 0.2, */
/*                 2.86, */
/*                 0.75, */
/*                 false, */
/*                 "fiberForceLengthCurve"); */
/*         OpenSim::SmoothSpline spline =
 * OpenSim::SmoothSegmentedFunctionFactory:: */
/*             SmoothSegmentedFunctionFactory::createFiberForceLengthCurve(
 */
/*                 0.0, */
/*                 0.7, */
/*                 0.2, */
/*                 2.86, */
/*                 0.75, */
/*                 "fiberForceLengthCurve"); */
/*         convertAndWriteToCsv(f, spline); */
/*     } */

/*     return 0; */
/* } */

/* int runSimulation(int argc, char *argv[]) { */
/* 	// Create new model. */
/*     OpenSim::Model model = OpenSim::Model(); */
/* 	model.setName("Blips"); */

/* 	// Create two links, each with a mass of 1 kg, center of mass at the
 * body's
 */
/* 	// origin, and moments and products of inertia of zero. */
/* 	double mass = 1.; */
/* 	auto humerus = new OpenSim::Body("humerus", mass, SimTK::Vec3(0),
 * SimTK::Inertia(0)); */
/* 	auto radius = new OpenSim::Body("radius", mass, SimTK::Vec3(0),
 * SimTK::Inertia(0)); */

/* 	// Connect the bodies with pin joints. Assume each body is 1 m long.
 */
/* 	auto shoulder = new OpenSim::PinJoint("shoulder", */
/* 								 model.getGround(), // Parent body */
/* 								 SimTK::Vec3(0, 2, 0),		// Location in
 * parent
 */
/* 								 SimTK::Vec3(0),			// Orientation in
 * parent
 */
/* 								 *humerus,			// Child body */
/* 								 SimTK::Vec3(0, 1, 0),		// Location in
 * child
 */
/* 								 SimTK::Vec3(0)			// Orientation in
 * child
 */
/* 	); */
/* 	auto elbow = new OpenSim::PinJoint("elbow", *humerus,
 * SimTK::Vec3(0), SimTK::Vec3(0), *radius, */
/* 							  SimTK::Vec3(0, 1, 0), SimTK::Vec3(0)); */

/* 	// Add a muscle that flexes the elbow. */
/* 	double maxIsometricForce = 200;	 // N */
/* 	double optimalFiberLength = 0.6; // m */
/* 	double tendonSlackLength = 0.55; // m */
/* 	double pennationAngle = 0.0;	 // rad */

/* 	std::cout << "Muscle properties:" << std::endl; */
/* 	std::cout << "    maxIsometricForce = " << maxIsometricForce <<
 * std::endl;
 */
/* 	std::cout << "    optimalFiberLength = " << optimalFiberLength <<
 * std::endl;
 */
/* 	std::cout << "    tendonSlackLength = " << tendonSlackLength <<
 * std::endl;
 */
/* 	std::cout << "    pennationAngle = " << pennationAngle << std::endl;
 */
/* 	std::cout << "    ratio tendonSlackLength to optimalFiberLength = "
 * << tendonSlackLength / optimalFiberLength << std::endl; */

/*     OpenSim::Muscle *biceps; */
/*     auto millard = new OpenSim::Millard2012EquilibriumMuscle( */
/*             "biceps", maxIsometricForce, optimalFiberLength,
 * tendonSlackLength, */
/*             pennationAngle); */
/*     millard->set_fiber_damping(0.1); */
/*     biceps = millard; */

/* 	biceps->addNewPathPoint("origin", *humerus, SimTK::Vec3(0, 0.8, 0));
 */
/* 	biceps->addNewPathPoint("insertion", *radius, SimTK::Vec3(0, 0.7,
 * 0)); */

/* 	// Add a controller that specifies the excitation of the muscle. */
/* 	auto brain = new OpenSim::PrescribedController(); */
/* 	brain->setName("brain"); */
/* 	brain->addActuator(*biceps); */

/* 	// Muscle excitation. */
/* 	const double sim_time = 10.; */
/* 	auto ctrlfn = */
/* 		new OpenSim::StepFunction(sim_time * 0.1, sim_time * 0.9, */
/* 						 0.1, 1.); */
/* 	brain->prescribeControlForActuator("biceps", ctrlfn); */

/* 	// Add components to the model. */
/* 	model.addBody(humerus); */
/* 	model.addBody(radius); */
/* 	model.addJoint(elbow); */
/* 	model.addJoint(shoulder); */
/* 	model.addForce(biceps); */
/* 	model.addController(brain); */

/* 	// Add a console reporter to print the muscle fiber force and elbow
 */
/* 	// angle.The output will be written to the log file(out.log) in the
 * current
 */
/* 	// directory. */
/* 	auto reporter = new OpenSim::TableReporter(); */

/* 	reporter->setName(model.getName() + "_results"); */
/* 	reporter->set_report_time_interval(1e-3); */

/* 	auto biceps_o_names = biceps->getOutputNames(); */
/* 	for (auto o_name : biceps_o_names) { */
/* 		reporter->addToReport(biceps->getOutput(o_name)); */
/* 	} */

/* 	model.addComponent(reporter); */

/* 	// Add display geometry. */
/* 	auto bodyGeometry = new OpenSim::Ellipsoid(0.1, 0.5, 0.1); */
/* 	bodyGeometry->setColor(SimTK::Vec3(0.5)); // Gray */

/* 	// Attach an ellipsoid to a frame located at the center of each
 * body. */
/* 	auto humerusCenter = new OpenSim::PhysicalOffsetFrame(); */
/* 	humerusCenter->setName("humerusCenter"); */
/* 	humerusCenter->setParentFrame(*humerus); */
/* 	humerusCenter->setOffsetTransform(SimTK::Transform(SimTK::Vec3(0,
 * 0.5, 0)));
 */
/* 	humerus->addComponent(humerusCenter); */
/* 	humerusCenter->attachGeometry(bodyGeometry->clone()); */

/* 	auto radiusCenter = new OpenSim::PhysicalOffsetFrame(); */
/* 	radiusCenter->setName("radiusCenter"); */
/* 	radiusCenter->setParentFrame(*radius); */
/* 	radiusCenter->setOffsetTransform(SimTK::Transform(SimTK::Vec3(0,
 * 0.5, 0)));
 */
/* 	radius->addComponent(radiusCenter); */
/* 	radiusCenter->attachGeometry(bodyGeometry->clone()); */

/* 	// Visualize. */
/*     /1* model.setUseVisualizer(true); *1/ */

/* 	// Configure the model. */
/* 	SimTK::State &state = model.initSystem(); */

/* 	shoulder->getCoordinate().setLocked(state, true); */
/* 	elbow->updCoordinate(OpenSim::PinJoint::Coord::RotationZ) */
/* 		.setDefaultValue(0.5 * SimTK::Pi); */
/* 	elbow->getCoordinate().setValue(state, 0.5 * SimTK::Pi); */
/* 	model.equilibrateMuscles(state); */

/* 	// Simulate. */
/* 	auto finalTime = sim_time; */
/* 	const double accuracy = 1e-5; */
/* 	auto manager = new OpenSim::Manager(model); */
/* 	manager->setIntegratorAccuracy(accuracy); */
/* 	/1* SimTK::Integrator& integ = manager->getIntegrator(); *1/ */
/* 	/1* integ.setUseInfinityNorm(); *1/ */
/* 	manager->initialize(state); */
/* 	manager->integrate(finalTime); */

/* 	model.print("basicMillard.osim"); */

/*     OpenSim::ForwardTool tool; */
/* 	tool.setName("ForwardIntegration"); */
/* 	tool.setModelFilename("basicMillard.osim"); */
/* 	tool.setFinalTime(finalTime); */
/* 	tool.setResultsDir("results"); */
/* 	tool.setErrorTolerance(accuracy); */
/* 	tool.print("setup_basic_millard.xml"); */

/* 	// Print table with simulation summary. */
/* 	SimTK::Integrator integ = manager->getIntegrator(); */
/* 	std::string table_delim = " | "; */
/* 	std::cout << "SimulationSummary" << table_delim << "name" <<
 * table_delim */
/* 			  << "simtime" << table_delim << "numStepsAttempted" <<
 * table_delim
 */
/* 			  << "numStepsTaken" << table_delim << "getNumRealizations"
 */
/* 			  << table_delim << "getNumIterationsstd::string path" */
/* 			  << table_delim << std::endl; */

/* 	std::cout << "SimulationSummary" << table_delim << argv[0] */
/* 			  << table_delim << finalTime << table_delim */
/* 			  << integ.getNumStepsAttempted() << table_delim */
/* 			  << integ.getNumStepsTaken() << table_delim */
/* 			  << integ.getNumRealizations() << table_delim */
/* 			  << integ.getNumIterations() << table_delim << std::endl;
 */

/* 	const std::string prefix = "ok"; */

/* 	// Plot results (print to stdout). */
/* 	{ */
/* 		auto table = reporter->getTable(); */
/* 		auto matrix = reporter->getTable().getMatrix(); */
/* 		for (int i = 0; i < matrix.nrow(); i++) { */
/* 			double time = table.getIndependentColumn()[i]; */
/* 			for (int j = 0; j < matrix.ncol(); j++) { */
/* 				auto label = table.getColumnLabel(j); */
/* 				double value = matrix.row(i)[j]; */
/* 				std::cout << "RUSTCSVPLOT," << prefix << label << "," << time
 */
/* 						  << "," << value << "," << j << std::endl; */
/* 			} */
/* 		} */
/* 	} */

/* 	{ */
/* 		auto table = manager->getStatesTable(); */
/* 		auto matrix = table.getMatrix(); */
/* 		auto table_rep = reporter->getTable(); */
/* 		auto matrix_rep = table_rep.getMatrix(); */
/* 		double prev_time = 0; */
/* 		int prev_nearest_row = 0; */
/* 		int steps = 0; */

/* 		for (int i = 0; i < matrix.nrow(); i++) { */
/* 			// Log timestep. */
/* 			double time = table.getIndependentColumn()[i]; */
/* 			double dt = time - prev_time; */
/* 			prev_time = time; */
/* 			std::cout << "RUSTCSVPLOT," << prefix << "dt," << time << ","
 * << dt
 */
/* 					  << std::endl; */
/* 			std::cout << "RUSTCSVPLOT," << prefix << "steps," << time <<
 * "," */
/* 					  << ++steps << std::endl; */

/* 			// Log all states. */
/* 			for (int j = 0; j < matrix.ncol(); j++) { */
/* 				auto label = table.getColumnLabel(j); */
/* 				double value = matrix.row(i)[j]; */
/* 				std::cout << "RUSTCSVPLOT," << prefix << label << "," << time
 */
/* 						  << "," << value << std::endl; */
/* 			} */

/* 			// Log reported values at state times. */
/* 			int nearest_row = table_rep.getNearestRowIndexForTime(time);
 */
/* 			if (nearest_row == prev_nearest_row) { */
/* 				continue; */
/* 			} */
/* 			prev_nearest_row = nearest_row; */
/* 			double time_rep =
 * table_rep.getIndependentColumn()[nearest_row]; */
/* 			for (int j = 0; j < matrix_rep.ncol(); j++) { */
/* 				auto label_rep = table_rep.getColumnLabel(j); */
/* 				double value_rep = matrix_rep.row(nearest_row)[j]; */
/* 				std::cout << "RUSTCSVPLOT,nearest-" << prefix << label_rep
 */
/* 						  << "," << time_rep << "," << value_rep <<
 * std::endl;
 */
/* 			} */
/* 		} */
/* 	} */

/* 	return 0; */
/* } */
