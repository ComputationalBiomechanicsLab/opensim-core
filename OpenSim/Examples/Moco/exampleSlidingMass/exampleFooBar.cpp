#include <iostream>
#include <memory>
#include <string>

#include <OpenSim/Moco/osimMoco.h>
#include <OpenSim/OpenSim.h>
#include <OpenSim/Simulation/SimbodyEngine/SliderJoint.h>

using namespace OpenSim;
using Vec3 = SimTK::Vec3;
using Inertia = SimTK::Inertia;

namespace {
    const std::string MODEL_NAME = std::string("arm");
    const std::string HAND_NAME = std::string("hand");
    const std::string SHOULDER_JOINT_NAME = std::string("shoulder");
    const std::string SHOULDER_ANGLE_NAME = std::string("angle");
    const std::string SHOULDER_MUSCLE_NAME = std::string("muscle");
} // namespace

std::unique_ptr<Model> createModel() {
    Model model = Model();
    model.setName(MODEL_NAME);

    // Set gravity.
    // model.setGravity(SimTK::Vec3(0, -9.80665, 0));

    // Add the hand.
    {
        // Hand parameters.
        double mass = 30.;
        double length = 0.2;
        Vec3 cog{0., 0., 0.};
        Inertia inertia = mass * Inertia::brick(Vec3(length / 2.));

        // Pelvis component.
        std::unique_ptr<Body> body(new Body(HAND_NAME, mass, cog, inertia));

        // Make pretty.
        std::unique_ptr<Brick> geometry = std::unique_ptr<Brick>(new
                Brick(
                Vec3(length)));

        // Add to model.
        model.addComponent(body.release());
    }

    // Attach the hand to ground with an arm.
    {
        // Arm length.
        double length = 1.;

        // Add the arm to the shoulder at the origin.
        const Ground& parent = model.getGround();
        const Vec3 locationInParent(0.);
        const Vec3 orientationInParent(0.);

        // And add the arm to the hand at arm's length.
        const Body& child = model.getComponent<Body>(HAND_NAME);
        const Vec3 locationInChild{0., 0., length};
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
        model.updComponent<Body>(HAND_NAME).attachGeometry(geometry.release());

        model.addComponent(joint.release());
    }

    // Connect actuator to shoulder.
    {
        std::unique_ptr<CoordinateActuator> muscle = std::unique_ptr<CoordinateActuator>(new CoordinateActuator());

        OpenSim::Joint& joint = model.updComponent<Joint>(SHOULDER_JOINT_NAME);
        OpenSim::Coordinate& coord = joint.updCoordinate();

        muscle->setCoordinate(&coord);
        muscle->setName(SHOULDER_MUSCLE_NAME);
        muscle->setOptimalForce(1);

        model.addComponent(muscle.release());
    }

    model.finalizeConnections();

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

    // Bounds.
    // -------
    // Initial time must be 0, final time can be within [0, 5].
    problem.setTimeBounds(MocoInitialBounds(0), MocoFinalBounds(0, 5));

    // Position must be within [-5, 5] throughout the motion.
    // Initial position must be 0, final position must be 1.
    problem.setStateInfo("/shoulder/angle/value", MocoBounds(-1, 1),
            MocoInitialBounds(0), MocoFinalBounds(1));
    // Speed must be within [-50, 50] throughout the motion.
    // Initial and final speed must be 0. Use compact syntax.
    problem.setStateInfo("/shoulder/angle/speed", {-50, 50}, 0, 0);

    // Applied force must be between -50 and 50.
    problem.setControlInfo("/muscle", MocoBounds(-50, 50));

    // Cost.
    // -----
    problem.addGoal<MocoFinalTimeGoal>();

    // Configure the solver.
    // =====================
    MocoCasADiSolver& solver = study.initCasADiSolver();
    solver.set_num_mesh_intervals(50);

    // Now that we've finished setting up the tool, print it to a file.
    study.print("sliding_mass.omoco");

    // Solve the problem.
    // ==================
    MocoSolution solution = study.solve();

    solution.write("sliding_mass_solution.sto");

    // Visualize.
    // ==========
    study.visualize(solution);

    return EXIT_SUCCESS;
}
