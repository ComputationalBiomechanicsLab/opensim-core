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

#include <string>
#include <iostream>
#include <math.h>

using namespace OpenSim;

template<typename T = SimTK::Vec3>
struct PathPoints final {
    PathPoints() {}

    PathPoints(T start_point, T end_point) :
        start(start_point),
        end(end_point)
    {}

    template<typename X>
    PathPoints(PathPoints<X> other) :
        PathPoints {
            T(other.start),
            T(other.end)
        }
    {}

    T start;
    T end;
};

PathPoints<> operator*(const SimTK::Rotation& rot, const PathPoints<>& pts) {
    return PathPoints<> {
        rot * pts.start,
        rot * pts.end
      };
}

struct WrapInput final {
    double radius;
    double length;
    std::string quadrant;
    SimTK::Vec3 orientation;
    PathPoints<SimTK::Vec3> path;
};

struct WrappingTolerances final {
    double point_error;
    double length_error;
};

struct CylindricalCoordinates final {
    CylindricalCoordinates() {}

    CylindricalCoordinates( //TODO why this, also change the names
        double radius,
        double angle,
        double position_on_axis) :
        radius(radius),
        angle(angle),
        position_on_axis(position_on_axis)
    {}

    // Convert to Cylindrical Coordinates with angle within [0, 2pi].
    CylindricalCoordinates(const SimTK::Vec3& point):
        angle(std::atan2(point[1], point[0])),
        radius(SimTK::Vec2(point[1], point[0]).norm()),
        position_on_axis(point[2])
    {
        wrap_angle();
    }

    void wrap_angle() {
        if (angle < 0.) {
            angle += 2. * M_PI;
        }
    }

    CylindricalCoordinates operator-()
    {
        return CylindricalCoordinates {
            radius,
            -angle,
            -position_on_axis
        };
    }

    double radius;
    double angle;
    double position_on_axis;
};

// Difference between cylindrical coordinates.
//
// Assumes radius of lhs equals that of rhs.
// Returns angle difference in range [0, 2pi].
CylindricalCoordinates operator-(
    const CylindricalCoordinates& lhs,
    const CylindricalCoordinates& rhs
) {
    ASSERT(lhs.radius == rhs.radius);

    CylindricalCoordinates diff (
        lhs.radius,
        lhs.angle - rhs.angle,
        lhs.position_on_axis - rhs.position_on_axis
    );
    diff.wrap_angle();
    return diff;
}

// The direction in which the cylinder is wrapped along its axis.
//
// The possibilities are:
// - Either positive or negative: normal wrap.
// - Neither positive nor negative: if path length is zero.
// - Both positive and negative: if inconsitent path (is bad).
struct WrappingSign final {
    private:
    WrappingSign(bool is_positive, bool is_negative) :
        _is_positive(is_positive),
        _is_negative(is_negative)
    {}

    public:
    WrappingSign() = default;

    // Determine sign from shortest angular distance around z-axis between start-end.
    WrappingSign(const PathPoints<>& wrapped_points) {
        double sin_theta_scaled = SimTK::dot(
                SimTK::cross(wrapped_points.start, wrapped_points.end),
                SimTK::Vec3(0., 0., 1.));
        _is_positive = sin_theta_scaled > SimTK::Eps;
        _is_negative = sin_theta_scaled < -SimTK::Eps;
    }

    static WrappingSign Positive() {
        return WrappingSign(true, false);
    }

    static WrappingSign Negative() {
        return WrappingSign(false, true);
    }

    bool is_positive() const {
        return _is_positive;
    }

    bool is_negative() const {
        return _is_negative;
    }

    bool is_inconsistent() const {
        return _is_positive && _is_negative;
    }

    bool is_undetermined() const {
        return !_is_positive && !_is_negative;
    }

    friend WrappingSign operator||(
        const WrappingSign& lhs,
        const WrappingSign& rhs);

    friend bool operator==(
        const WrappingSign& lhs,
        const WrappingSign& rhs);

    private:
    bool _is_positive = false;
    bool _is_negative = false;
};

// Returns true if all members equal.
bool operator==(
    const WrappingSign& lhs,
    const WrappingSign& rhs) {
    return lhs.is_positive() == rhs.is_positive()
        && lhs.is_negative() == rhs.is_negative();
}

// Memberwise OR.
WrappingSign operator||(
    const WrappingSign& lhs,
    const WrappingSign& rhs
    ) {
    return WrappingSign {
        lhs.is_positive() || rhs.is_positive(),
        lhs.is_negative() || rhs.is_negative()
     };
}

struct WrapTestResult final {
    static WrapTestResult NoWrap() {
        return WrapTestResult();
    }

    // Path points on cylinder surface.
    PathPoints<> path;
    // Direction of wrapping wrt cylinder axis.
    WrappingSign sign;
    // Length of path.
    double length;
    // If there is no wrapping, the other fields are undefined.
    bool no_wrap = true;
};

WrapTestResult solve(const WrapInput& input, bool visualize = false) {
    Model model;
    model.setName("testCylinderWrapping");

    WrapCylinder* cylinder = new WrapCylinder();
    cylinder->setName("cylinder");
    cylinder->set_radius(input.radius);
    cylinder->set_length(input.length);
    cylinder->set_quadrant(input.quadrant);
    cylinder->set_xyz_body_rotation(input.orientation);
    cylinder->setFrame(model.getGround());
    model.updGround().addWrapObject(cylinder);

    PathSpring* spring = new PathSpring("spring", 1., 1., 1.);
    spring->updGeometryPath().appendNewPathPoint("start_point", model.get_ground(), input.path.start);
    spring->updGeometryPath().appendNewPathPoint("end_point", model.get_ground(), input.path.end);
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

    // Convert the WrapResult to a WrapTestResult for convenient assertions.

    if (wrap_result.wrap_pts.size() == 0) {
        return WrapTestResult::NoWrap();
    }

    // Determine the wrapping sign based on the shortest path between consecutive points.
    // This is easier than asserting exact positions of intermediate path
    // points, while still capturing information for testing.
    SimTK::Rotation cylinder_orientation;
    cylinder_orientation.setRotationToBodyFixedXYZ(input.orientation);
    PathPoints consecutive_points =
        PathPoints{
            wrap_result.r1,
            wrap_result.r1
        };
    WrappingSign sign;
    wrap_result.wrap_pts.append(wrap_result.r2);
    for (int i = 0; i < wrap_result.wrap_pts.size(); ++i) {
        consecutive_points.start = consecutive_points.end;
        consecutive_points.end = wrap_result.wrap_pts[i];
        sign = sign || WrappingSign(cylinder_orientation.transpose() * consecutive_points);
    }
    return WrapTestResult {
        PathPoints{wrap_result.r1, wrap_result.r2}, sign
    };
}

explicit bool IsEqual(double lhs, double rhs, double tolerance) {
    return std::abs(lhs - rhs) < tolerance;
}

explicit bool IsEqual(const SimTK::Vec3& lhs, const SimTK::Vec3& rhs, double tolerance) {
    return SimTK::max(SimTK::abs(lhs - rhs)) < tolerance;
}

explicit bool IsEqual(const WrapTestResult& lhs, const WrapTestResult& rhs, double tolerance) {
    if (lhs.no_wrap && rhs.no_wrap) {
        return true;
    }
    return lhs.no_wrap == rhs.no_wrap
        && lhs.sign == rhs.sign
        && IsEqual(lhs.path.start, rhs.path.start, tolerance)
        && IsEqual(lhs.path.end, rhs.path.end, tolerance);
}

bool TestWrapping(
    const WrapInput& input,
    const WrapTestResult expected,
    const WrappingTolerances& tol,
    bool visualize = false)
{

    WrapTestResult result = solve(input, visualize);

    // =================================================================
    // =============== Test: Match expected result =====================
    // =================================================================

    if (!IsEqual(result, expected, tol.point_error)) {
        return false;
    }

    if (expected.no_wrap) {
        return false;
    }

    // =================================================================
    // ========== Test: Wrapped points distance to surface. ============
    // =================================================================

    SimTK::Rotation cylinder_orientation;
    cylinder_orientation.setRotationToBodyFixedXYZ(input.orientation);

    auto TestDistanceToSurface = [&](const PathPoints<>& points, double tol) -> bool {
        PathPoints<CylindricalCoordinates> points_in_cyl_frame(
            cylinder_orientation * points);
        return IsEqual(points_in_cyl_frame.start.radius, input.radius, tol)
            && IsEqual(points_in_cyl_frame.end.radius, input.radius, tol);
    };

    if (TestDistanceToSurface(result.path, tol.point_error))
    {
        return false;
    }

    // Check expected result for good measure.
    if (TestDistanceToSurface(expected.path, 1e-10))
    {
        return false;
    }

    // =================================================================
    // =============== Expected wrapping length. =======================
    // =================================================================

    // Just compute both path lengths.

    auto TraversedCylindricalCoords = [&](
        const PathPoints<>& path,
        const WrappingSign& sign) -> CylindricalCoordinates
    {
        PathPoints<CylindricalCoordinates> path_in_cylinder_frame(
            cylinder_orientation * path);
        return sign.is_positive()?
            path_in_cylinder_frame.end - path_in_cylinder_frame.start:
            -(path_in_cylinder_frame.start - path_in_cylinder_frame.end);
    };

    auto ComputePathLength = [&](
        const PathPoints<>& wrapped_path,
        const WrappingSign& sign) -> double
    {
        CylindricalCoordinates delta = TraversedCylindricalCoords(wrapped_path, sign);
        return
            (input.path.start - wrapped_path.start).norm() +
            std::sqrt(
                std::pow(delta.angle * delta.radius, 2) +
                std::pow(delta.position_on_axis, 2)
            ) +
            (wrapped_path.end - input.path.end).norm();
    };

    if (!IsEqual(ComputePathLength(result.path, result.sign), result.length, 1e-10) ||
        !IsEqual(ComputePathLength(expected.path, expected.sign), expected.length, 1e-10) ||
        !IsEqual(result.length, expected.length, tol.length_error))
    {
        return false;
    }

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

    // The expected answers.
    bool expected_no_wrap = false;
    bool expected_positive_wrapping_direction = true;
    SimTK::Vec3 expected_start_point_on_surface;
    SimTK::Vec3 expected_end_point_on_surface;
    double expected_path_length;

    // Assertion bounds:

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