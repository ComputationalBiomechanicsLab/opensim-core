#include "OpenSim/Simulation/Model/Geometry.h"
#include <OpenSim/Moco/osimMoco.h>
#include <OpenSim/OpenSim.h>

#include <iostream>
#include <memory>
#include <string>

using namespace OpenSim;
using Vec3 = SimTK::Vec3;
using Inertia = SimTK::Inertia;

namespace {
    const std::string MODEL_NAME = std::string("arm");
    const std::string HAND_NAME = std::string("hand");

    const std::string SHOULDER_JOINT_NAME = std::string("shoulder");
    const std::string SHOULDER_ANGLE_NAME = std::string("angle");
    const std::string SHOULDER_MUSCLE_NAME = std::string("muscle");

    const std::string ELBOW_FAT_NAME = std::string("elbow-fat");
    const std::string ELBOW_JOINT_NAME = std::string("elbow");
    const std::string ELBOW_ANGLE_NAME = std::string("angle");

    const std::string ELBOW_MUSCLE_NAME = std::string("elbow-muscle");
} // namespace

std::unique_ptr<Model> createModel() {
    Model model = Model();
    model.setName(MODEL_NAME);

    // Set gravity.
    model.setGravity(SimTK::Vec3(0, -9.80665, 0));

    // Add shoulder geometry.
    {
        // Make pretty.
        std::unique_ptr<Brick> geometry = std::unique_ptr<Brick>(new Brick({0.01, 0.01,1.}));
        geometry->setColor(SimTK::Vec3(0.1, 0.1, 0.8));
        model.updGround().attachGeometry(geometry.release());
    }

    // Add the hand.
    {
        // Hand parameters.
        double mass = 30.;
        double length = 0.05;
        Vec3 cog{0., 0., 0.};
        Inertia inertia = mass * Inertia::brick(Vec3(length / 2.));

        // Pelvis component.
        Body body(HAND_NAME, mass, cog, inertia);

        // Make pretty.
        {
            Brick geometry = Brick(Vec3{length, 3*length, 0.1*length});
            geometry.setColor(SimTK::Vec3(0.1, 0.1, 0.8));
            body.attachGeometry(new Brick(geometry));
        }

        // Add to model.
        model.addComponent(new Body(body));
    }

    // Add the elbow.
    {
        double mass = 0.1;
        double length = 0.2;
        double radius = 0.02;
        Vec3 cog{0., 0., 0.};
        Inertia inertia = mass * Inertia::brick(Vec3(length / 2.));

        // Pelvis component.
        std::unique_ptr<Body> body(new Body(ELBOW_FAT_NAME, mass, cog, inertia));

        // Make pretty.
        std::unique_ptr<Brick> geometry = std::unique_ptr<Brick>(new
                Brick({radius, radius, length}));
        geometry->setColor(SimTK::Vec3(0.1, 0.8, 0.1));
        body->attachGeometry(geometry.release());

        // Add to model.
        model.addComponent(body.release());
    }

    // Attach the elbow to ground using shoulder.
    {
        // Upper arm length.
        double length = 1.;

        // Add the arm to the shoulder at the origin.
        const Ground& parent = model.getGround();
        const Vec3 locationInParent(0.);
        const Vec3 orientationInParent(0.);

        // And add the arm to the elbow.
        const Body& child = model.getComponent<Body>(ELBOW_FAT_NAME);
        const Vec3 locationInChild{0., length, 0.};
        const Vec3 orientationInChild(0.);

        std::unique_ptr<PinJoint> joint(
                new PinJoint(SHOULDER_JOINT_NAME,
                    parent,
                    locationInParent,
                    orientationInParent,
                    child,
                    locationInChild,
                    orientationInChild));

        auto& coord = joint->updCoordinate(PinJoint::Coord::RotationZ);
        coord.setName(SHOULDER_ANGLE_NAME);

        // Make pretty.
        double radius = 0.01;
        std::unique_ptr<Cylinder> geometry = std::unique_ptr<Cylinder>(new
                Cylinder(radius, length / 2.));
        geometry->setColor(SimTK::Vec3(0.8, 0.1, 0.1));
        model.updComponent<Body>(ELBOW_FAT_NAME).attachGeometry(geometry.release());

        model.addComponent(joint.release());
    }

    // Attach the hand to elbow.
    {
        double length = 1.;

        // Add the arm to the elbow at the origin.
        const Body& parent = model.getComponent<Body>(ELBOW_FAT_NAME);
        const Vec3 locationInParent(0.);
        const Vec3 orientationInParent(0.);

        // And add the arm to the hand at arm's length.
        const Body& child = model.getComponent<Body>(HAND_NAME);
        const Vec3 locationInChild{0., length, 0.};
        const Vec3 orientationInChild(0.);

        std::unique_ptr<PinJoint> joint(
                new PinJoint(ELBOW_JOINT_NAME,
                    parent,
                    locationInParent,
                    orientationInParent,
                    child,
                    locationInChild,
                    orientationInChild));

        auto& coord = joint->updCoordinate(PinJoint::Coord::RotationZ);
        coord.setName(ELBOW_ANGLE_NAME);

        // Make pretty.
        double radius = 0.01;
        std::unique_ptr<Cylinder> geometry = std::unique_ptr<Cylinder>(new
                Cylinder(radius, length / 2.));
        geometry->setColor(SimTK::Vec3(0.8, 0.1, 0.1));
        model.updComponent<Body>(HAND_NAME).attachGeometry(geometry.release());

        model.addComponent(joint.release());
    }

    // Connect actuator to shoulder.
    {
        std::unique_ptr<CoordinateActuator> muscle = std::unique_ptr<CoordinateActuator>(new CoordinateActuator());

        OpenSim::Joint& joint = model.updComponent<Joint>(SHOULDER_JOINT_NAME);
        OpenSim::Coordinate& coord = joint.updCoordinate();

        muscle->setName(SHOULDER_MUSCLE_NAME);
        muscle->setCoordinate(&coord);
        muscle->setOptimalForce(1);

        model.addComponent(muscle.release());
    }

    // Connect actuator to elbow.
    {
        std::unique_ptr<CoordinateActuator> muscle = std::unique_ptr<CoordinateActuator>(new CoordinateActuator());

        OpenSim::Joint& joint = model.updComponent<Joint>(ELBOW_JOINT_NAME);
        OpenSim::Coordinate& coord = joint.updCoordinate();

        muscle->setName(ELBOW_MUSCLE_NAME);
        muscle->setCoordinate(&coord);
        muscle->setOptimalForce(1);

        model.addComponent(muscle.release());
    }

    model.finalizeConnections();
    model.initSystem();
    model.printSubcomponentInfo();
    model.printOutputInfo(true);
    auto names = model.getStateVariableNames();
    for (int i = 0; i < names.size(); ++i) {
        std::cout << "name = " << names[i] << std::endl;
    }

    /* std::cout << "HAND abs path = " << model.getComponent<Body>(HAND_NAME).getAbsolutePath().toString() << std::endl; */

    return std::unique_ptr<Model>(new Model(model));
}

int main() {

    MocoStudy study;
    study.setName("FooBar");

    // Define the optimal control problem.
    // ===================================
    MocoProblem& problem = study.updProblem();

    // Model (dynamics).
    // -----------------
    problem.setModel(createModel());

    // Cost.
    // -----
    problem.addGoal<MocoFinalTimeGoal>();

    // Bounds.
    // -------
    double finalTime = 10.0;
    problem.setTimeBounds(0, finalTime);

    // Position must be within [-5, 5] throughout the motion.
    // Initial position must be 0, final position must be 1.
    problem.setStateInfo("/shoulder/angle/value", MocoBounds(-1, 1),
            MocoInitialBounds(0), MocoFinalBounds(1));
    // Speed must be within [-50, 50] throughout the motion.
    // Initial and final speed must be 0. Use compact syntax.
    problem.setStateInfo("/shoulder/angle/speed", {-50, 50}, 0, 0);

    // Applied force must be between -50 and 50.
    problem.setControlInfo("/muscle", MocoBounds(-50, 50));
    problem.setControlInfo("/elbow-muscle", MocoBounds(-50, 50));

    // Cost.
    // -----
    auto* tracking = problem.addGoal<MocoOutputTrackingGoal>("tracking");
    /* TimeSeriesTable ref; */
    /* ref.addTableMetaData<std::string>("inDegrees", "no"); */
    /* ref.setColumnLabels({"/hand/position"}); */
    /* // We supply a reference whose time range is a superset of the problem's */
    /* // time bounds: Moco performs finite differences internally, which may */
    /* // require sampling outside the problem's time bounds. */
    /* for (double time = -0.05; time < finalTime + 0.05; time += 0.01) { */
    /*     ref.appendRow(time, { */
    /*             0.0 * SimTK::Pi * time / finalTime */
    /*     }); */
    /* } */

    tracking->setOutputPath("/hand|position");
    /* tracking->setReference(ref); */

    // A low-weighted cost term minimizing controls helps with convergence.
    problem.addGoal<MocoControlGoal>("effort", 0.001);

    // Configure the solver.
    // =====================
    MocoCasADiSolver& solver = study.initCasADiSolver();
    solver.set_num_mesh_intervals(50);

    // Now that we've finished setting up the tool, print it to a file.
    study.print("exampleMocoArm.omoco");

    // Solve the problem.
    // ==================
    MocoSolution solution = study.solve();

    solution.write("exampleMocoArm_solution.sto");

    // Visualize.
    // ==========
    study.visualize(solution);

    return EXIT_SUCCESS;
}
