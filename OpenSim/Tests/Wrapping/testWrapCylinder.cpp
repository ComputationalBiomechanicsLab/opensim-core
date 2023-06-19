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
    double gradient_error;
};

// Angle defined on range [0, 2pi].
struct Angle {
    Angle(SimTK::Vec2 point):
        _value(std::atan2(point[1], point[0]))
    {
        wrap();
    }

    double get_value() const {
        return _value;
    }

    private:
    void wrap() {
        _value = (_value < 0. )? _value + 2. * M_PI: _value;
    }

    double _value;
};

struct AngleDifference {
    double get_value(const WrappingSign& sign) {
        if (sign.is_positive()) {
            return positive_direction;
        }
    }

    double positive_direction;
    double negative_direction;
};

AngleDifference operator-(const Angle& lhs, const Angle& rhs) {
    double delta = lhs.get_value() - rhs.get_value();
    if (delta < 0.) {
        return AngleDifference{
            delta + 2. * M_PI,
            delta
        };
    };
    return AngleDifference {
        delta,
        delta - 2. * M_PI
    };
}

struct CylindricalCoordinates final {
    CylindricalCoordinates(const SimTK::Vec3& point):
        angle(SimTK::Vec2(point[0], point[1])),
        radius(SimTK::Vec2(point[0], point[1]).norm()),
        position_on_axis(point[2])
    {}

    double radius;
    Angle angle;
    double position_on_axis;
};

struct CylindricalCoordinatesDifference final {
    double radial;
    AngleDifference angular;
    double along_axis;
};

CylindricalCoordinatesDifference operator-(
    const CylindricalCoordinates& lhs,
    const CylindricalCoordinates& rhs)
{
    return CylindricalCoordinatesDifference {
        lhs.radius - rhs.radius,
        lhs.angle - rhs.angle,
        lhs.position_on_axis - rhs.position_on_axis
    };
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
        return _is_positive && !_is_negative;
    }

    bool is_negative() const {
        return _is_negative && !_is_positive;
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
    return lhs._is_positive == rhs._is_positive
        && lhs._is_negative == rhs._is_negative;
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

WrapTestResult solve(const WrapInput& input, bool visualize) {
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
    SimTK::Rotation cylinder_orientation;
    cylinder_orientation.setRotationToBodyFixedXYZ(input.orientation);
    PathPoints consecutive_points =
        PathPoints{
            input.path.start,
            input.path.start
        };
    WrappingSign sign;
    wrap_result.wrap_pts.insert(0, wrap_result.r1);
    wrap_result.wrap_pts.append(wrap_result.r2);
    wrap_result.wrap_pts.append(input.path.end);
    for (int i = 0; i < wrap_result.wrap_pts.size(); ++i) {
        consecutive_points.start = consecutive_points.end;
        consecutive_points.end = wrap_result.wrap_pts[i];
        sign = sign || WrappingSign(cylinder_orientation.invert() * consecutive_points);
    }
    return WrapTestResult {
        PathPoints{wrap_result.r1, wrap_result.r2}, sign
    };
}

explicit bool IsEqual(double lhs, double rhs, double tolerance) {
    return std::abs(lhs - rhs) < tolerance;
}

explicit bool IsEqual(const SimTK::UnitVec3& lhs, const SimTK::UnitVec3& rhs, double tolerance) {
    return IsEqual(lhs[0], rhs[0], tolerance)
        && IsEqual(lhs[1], rhs[1], tolerance)
        && IsEqual(lhs[2], rhs[2], tolerance);
}

explicit bool IsEqual(const SimTK::Vec3& lhs, const SimTK::Vec3& rhs, double tolerance) {
    return IsEqual(lhs[0], rhs[0], tolerance)
        && IsEqual(lhs[1], rhs[1], tolerance)
        && IsEqual(lhs[2], rhs[2], tolerance);
}

explicit bool IsEqual(const PathPoints<>& lhs, const PathPoints<>& rhs, double tolerance) {
    return IsEqual(lhs.start, rhs.start, tolerance)
        && IsEqual(lhs.end, rhs.end, tolerance);
}

explicit bool IsEqual(const WrapTestResult& lhs, const WrapTestResult& rhs, double tolerance) {
    if (lhs.no_wrap && rhs.no_wrap) {
        return true;
    }
    return IsEqual(lhs.path, rhs.path, tolerance)
        && IsEqual(lhs.length, rhs.length, tolerance) // TODO: or use Tolerance.length_error
        && lhs.sign == rhs.sign
        && lhs.no_wrap == rhs.no_wrap;
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

    WrappingSign sign = result.sign;

    auto ComputeWrappedPathLength = [&](
        const PathPoints<>& wrapped_path) -> double
    {
        PathPoints<CylindricalCoordinates> path_in_cylinder_frame(
            cylinder_orientation * wrapped_path);

        CylindricalCoordinatesDifference difference =
            path_in_cylinder_frame.end -
            path_in_cylinder_frame.start;

        return std::sqrt(
                std::pow(difference.angular.get_value(sign) * input.radius, 2) +
                std::pow(difference.along_axis, 2)
            );
    };

    auto ComputePathLength = [&](
        const PathPoints<>& wrapped_path) -> double
    {
        return
            (input.path.start - wrapped_path.start).norm() +
            ComputeWrappedPathLength(wrapped_path) +
            (wrapped_path.end - input.path.end).norm();
    };

    if (!IsEqual(ComputePathLength(result.path), result.length, 1e-10) ||
        !IsEqual(ComputePathLength(expected.path), expected.length, 1e-10) ||
        !IsEqual(result.length, expected.length, tol.length_error))
    {
        return false;
    }

    // =================================================================
    // =============== Test: Is path tangent to surface? ===============
    // =================================================================

    // Compute surface gradient along path.
    auto ComputeSurfaceGradient = [&](
        const PathPoints<>& wrapped_path,
        double t // Between [0, 1] = [start, end]
    ) -> SimTK::Vec3 {
        PathPoints<CylindricalCoordinates> path_in_cylinder_frame(
            cylinder_orientation * wrapped_path);

        CylindricalCoordinatesDifference difference =
            path_in_cylinder_frame.end -
            path_in_cylinder_frame.start;

        double angle = path_in_cylinder_frame.start.angle.get_value() +
            difference.angular.get_value(sign) * t;

        SimTK::Rotation rot_about_axis;
        rot_about_axis.setRotationFromAngleAboutZ(angle);

        return cylinder_orientation *
            rot_about_axis * (
                SimTK::Vec3(0., 1., 0.) * input.radius
                * difference.angular.get_value(sign)
                + SimTK::Vec3(0., 0., 1.)
                * difference.along_axis
            );
    };

    auto TestGradientDirection = [&](
        const PathPoints<>& wrapped_path,
        double tolerance
    ) -> bool {
        return
            IsEqual(
                SimTK::UnitVec3(result.path.start - input.path.start),
                SimTK::UnitVec3(ComputeSurfaceGradient(result.path, 0.)),
                tolerance) &&
            IsEqual(
                SimTK::UnitVec3(input.path.end - result.path.end),
                SimTK::UnitVec3(ComputeSurfaceGradient(result.path, 1.)),
                tolerance);
    };

    if (!TestGradientDirection(result.path, tol.gradient_error) ||
        !TestGradientDirection(expected.path, 1e-10))
    {
        return false;
    }

    // =================================================================
    // =============== Test: Constant gradient along cylinder axis =====
    // =================================================================

}

int main()
{
    WrappingTolerances tolerance {
        1e-10,
        1e-10,
        1e-10
    };

    {
        WrapInput input {
            0.5, // radius;
            1., // length;
            "+x", // quadrant;
            SimTK::Vec3 {0., 0., 0.}, // orientation;
            PathPoints<> {
                SimTK::Vec3 {1., 2., 3.},
                SimTK::Vec3 {1., -2., -3.}
            }
        };

        WrapTestResult expected {
            PathPoints<> {
                SimTK::Vec3 {1., 2., 3.},
                SimTK::Vec3 {1., -2., -3.}
            },
            WrappingSign::Positive(),
            0.1,
        };

        ASSERT(TestWrapping(input, expected, tolerance));
    }

    {
        // TODO more cases
        // ...
    }

    // TODO report on failures...

    std::cout << "Done" << std::endl;
    return 0;
}