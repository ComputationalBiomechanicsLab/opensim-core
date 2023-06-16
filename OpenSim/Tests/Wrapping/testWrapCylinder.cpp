/* -------------------------------------------------------------------------- *
 *                         OpenSim:  testWrappingAlgorithm.cpp                *
 * -------------------------------------------------------------------------- *
 * The OpenSim API is a toolkit for musculoskeletal modeling and simulation.  *
 * See http://opensim.stanford.edu and the NOTICE file for more information.  *
 * OpenSim is developed at Stanford University and supported by the US        *
 * National Institutes of Health (U54 GM072970, R24 HD065690) and by DARPA    *
 * through the Warrior Web program.                                           *
 *                                                                            *
 * Copyright (c) 2005-2022 Stanford University and the Authors                *
 * Author(s): Ayman Habib                                                     *
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
// INCLUDE
#include <OpenSim/OpenSim.h>
#include <OpenSim/Auxiliary/auxiliaryTestFunctions.h>

#include <set>
#include <string>
#include <iostream>

using namespace OpenSim;
using namespace SimTK;
using namespace std;

// Fill in this form to test WrapCylinder::wrapLine()
//
// Performed tests:
// - Expected wrap-action must match.
// - Path must be tangent to surface at start and end surface points.
// - Projected path to XY-plane must be a circle (for the surface part).
// - Z-coordinate of path must have constant gradient.
struct WrapCylinderTestCase final {
    // Name of this test.
    std::string name;
    bool visualize = false;

    // Wrapping cylinder parameters.
    double radius = 1.;
    double length = 1.;
    std::string quadrant = "+x";
    SimTK::Vec3 orientation = {0., 0., 0.};

    // Path start and end points.
    SimTK::Vec3 start_point;
    SimTK::Vec3 end_point;

    // The expected answers.
    bool expected_no_wrap = false;
    bool expected_positive_wrapping_direction = true;
    SimTK::Vec3 expected_start_point_on_surface;
    SimTK::Vec3 expected_end_point_on_surface;
    double expected_path_length;

    // Assertion bounds:
    double epsilon_surface_gradient = 1e-10;
    double epsilon_point_position = 1e-10;
    double epsilon_on_circle = 1e-10;
    double epsilon_path_length = 1e-10;

    private:
    // Local cylindrical coordinates for points on cylinder surface.
    struct CylindricalPathCoordinates final {
        double angle;
        double position_along_axis;

        CylindricalPathCoordinates() = default;

        CylindricalPathCoordinates(
            SimTK::Rotation cylinder_orientation,
            SimTK::Vec3 point_on_surface
        ) {
            SimTK::Vec3 point_local =
            cylinder_orientation.transpose() * point_on_surface;
            angle = std::atan2(point_local[1], point_local[0]);
            position_along_axis = point_local[2];
        }

        CylindricalPathCoordinates operator-(
            const CylindricalPathCoordinates& rhs
        ) const {
            CylindricalPathCoordinates diff;
            diff.angle = angle - rhs.angle;
            diff.position_along_axis =
                position_along_axis - rhs.position_along_axis;
            return diff;
        }
    };

    public:
    void do_test(
        SimTK::Array_<std::string> failures
    ) const {
        try {
            do_test();
        }
        catch (const std::exception& e) {
            std::cout << "Exception: " << e.what() << std::endl;
            failures.push_back("TestWrapCylinder: case " + name);
        }
    }

    void do_test() const {
        Model model;
        model.setName("testCylinderWrapping");

        WrapCylinder* cylinder = new WrapCylinder();
        cylinder->setName("cylinder");
        cylinder->set_radius(radius);
        cylinder->set_length(length);
        cylinder->set_quadrant(quadrant);
        cylinder->set_xyz_body_rotation(orientation);
        cylinder->setFrame(model.getGround());
        model.updGround().addWrapObject(cylinder);

        PathSpring* spring = new PathSpring("spring", 1., 1., 1.);
        spring->updGeometryPath().appendNewPathPoint("start_point", model.get_ground(), start_point);
        spring->updGeometryPath().appendNewPathPoint("end_point", model.get_ground(), end_point);
        spring->updGeometryPath().addPathWrap(*cylinder);
        model.addComponent(spring);

        model.finalizeConnections();
        model.setUseVisualizer(visualize);

        SimTK::State& state = model.initSystem();
        model.realizeVelocity(state);
        if (visualize) {
            model.getVisualizer().show(state);
        }

        WrapResult wrap_result =
            spring->get_GeometryPath().getWrapSet().get("pathwrap").getPreviousWrap();

        // =================================================================
        // =============== Test: Expected wrapping action. =================
        // =================================================================

        if (expected_no_wrap) {
            ASSERT(wrap_result.wrap_pts.size() == 0);
            ASSERT_EQUAL(wrap_result.wrap_path_length, 0., 1e-10);
            return;
        } else {
            ASSERT(wrap_result.wrap_pts.size() > 0);
        }

        // =================================================================
        // =============== Test: Match expected points =====================
        // =================================================================

        ASSERT_EQUAL<SimTK::Vec3>(
            wrap_result.r1,
            expected_start_point_on_surface,
            1e-10);

        ASSERT_EQUAL<SimTK::Vec3>(
            wrap_result.r2,
            expected_end_point_on_surface,
            1e-10);

        // =================================================================
        // =============== Test: Points on circle. =========================
        // =================================================================

        SimTK::Rotation cylinder_orientation;
        cylinder_orientation.setRotationToBodyFixedXYZ(orientation);

        auto assert_point_on_surface = [&](
            const SimTK::Vec3& point
        ) -> void {
            SimTK::Vec3 point_local =
                cylinder_orientation * point;
            ASSERT_EQUAL<double>(
                SimTK::Vec2(
                    point_local[0],
                    point_local[1]).norm(),
                radius,
                1e-10);
        };

        assert_point_on_surface(wrap_result.r1);
        assert_point_on_surface(wrap_result.r2);

        for (size_t i = 0; i < wrap_result.wrap_pts.size(); ++i ) {
            assert_point_on_surface(wrap_result.wrap_pts[i]);
        }

        // =================================================================
        // =============== Expected wrapping length. =======================
        // =================================================================

        // Compute surface coordinates in local cylindrical coordinates.
        CylindricalPathCoordinates start_surface_coords(
            cylinder_orientation,
            wrap_result.r1);
        CylindricalPathCoordinates end_surface_coords(
            cylinder_orientation,
            wrap_result.r2);

        // Traversed cylinderical coordinates.
        CylindricalPathCoordinates delta_surface_coords =
            end_surface_coords - start_surface_coords;

        // Remove any wrapping to correctly compute the traversed angle.
        if (expected_positive_wrapping_direction &&
            delta_surface_coords.angle < 0)
        {
            delta_surface_coords.angle += 2. * M_PI;
        }
        if (!expected_positive_wrapping_direction &&
            delta_surface_coords.angle > 0)
        {
            delta_surface_coords.angle -= 2. * M_PI;
        }

        // Length supposing that r1, and r2 in wrap_result are correct:
        double path_length_check =
            (start_point - wrap_result.r1).norm() + // Point-a to surface.
            std::sqrt( // Path on cylinder surface.
                std::pow(delta_surface_coords.angle * radius, 2) +
                std::pow(delta_surface_coords.position_along_axis,2)
            ) +
            (end_point - wrap_result.r2).norm(); // Point-b to surface.

        ASSERT_EQUAL<double>(wrap_result.wrap_path_length, path_length_check, 1e-10);
        ASSERT_EQUAL<double>(wrap_result.wrap_path_length, expected_path_length, 1e-10);

        // =================================================================
        // =============== Test: Is path tangent to surface? ===============
        // =================================================================

        // Compute surface gradient at a cylindrical angle coordinate.
        auto ComputeSurfaceGradient = [&](
            double angle
        ) -> SimTK::Vec3 {
            SimTK::Rotation rot_about_axis;
            rot_about_axis.setRotationFromAngleAboutZ(angle);
            return cylinder_orientation *
                rot_about_axis * (
                    SimTK::Vec3(0., 1., 0.) * radius
                    + SimTK::Vec3(0., 0., 1.)
                    * delta_surface_coords.position_along_axis / delta_surface_coords.angle
                );
        };

        // Path gradient direction at start-point must equal that at surface.
        ASSERT_EQUAL<SimTK::UnitVec3>(
            SimTK::UnitVec3(wrap_result.r1 - start_point),
            SimTK::UnitVec3(ComputeSurfaceGradient(start_surface_coords.angle)),
            1e-10
        );

        // Path gradient direction at end-point must equal that at surface.
        ASSERT_EQUAL<SimTK::UnitVec3>(
            SimTK::UnitVec3(ComputeSurfaceGradient(end_surface_coords.angle)),
            SimTK::UnitVec3(end_point - wrap_result.r2),
            1e-10
        );

        // =================================================================
        // =============== Test: Constant gradient along cylinder axis =====
        // =================================================================


        // =================================================================
        // =============== Test: Constant gradient along cylinder axis =====
        // =================================================================
    }
};

int main()
{
    SimTK::Array_<std::string> failures;

    {
        WrapCylinderTestCase testCase;
        testCase.name = "Unconstrained Miss";
        testCase.radius = 0.5;
        testCase.length = 1.;
        testCase.quadrant = "unconstrained";
        testCase.orientation = {0., 0., 0.};
        testCase.start_point = {1., 2., 3.};
        testCase.end_point = {1, -2, -3.};

        testCase.expected_no_wrap = true;
        testCase.visualize = true;

        testCase.do_test(failures);
    }

    {
        WrapCylinderTestCase testCase;
        testCase.name = "Unconstrained Miss";
        testCase.radius = 1.;
        testCase.length = 1.;
        testCase.quadrant = "+x";
        testCase.orientation = {0., 0., 0.};
        testCase.start_point = {1., 2., 3.};
        testCase.end_point = {1, -2, -3.};

        testCase.expected_no_wrap = true;
        testCase.expected_positive_wrapping_direction = true;
        testCase.expected_start_point_on_surface = {1., 2., 3.};
        testCase.expected_end_point_on_surface = {1., 2., 3.};
        testCase.expected_path_length = 1.;
        testCase.visualize = true;

        testCase.do_test(failures);
    }

    if (!failures.empty()) {
        std::cout << "Done, with failure(s): " << failures << std::endl;
        return 1;
    }

    std::cout << "Done" << std::endl;
    return 0;
}