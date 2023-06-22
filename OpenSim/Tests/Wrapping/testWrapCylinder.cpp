/* -------------------------------------------------------------------------- *
 *                         OpenSim:  testWrapCylinder.cpp                     *
 * -------------------------------------------------------------------------- *
 * The OpenSim API is a toolkit for musculoskeletal modeling and simulation.  *
 * See http://opensim.stanford.edu and the NOTICE file for more information.  *
 * OpenSim is developed at Stanford University and supported by the US        *
 * National Institutes of Health (U54 GM072970, R24 HD065690) and by DARPA    *
 * through the Warrior Web program.                                           *
 *                                                                            *
 * Copyright (c) 2005-2022 Stanford University and the Authors                *
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
#include <cmath>
#include <sstream>
#include <functional>
#include <memory>

using namespace OpenSim;

// Section with geometry related helper functions and structs.
namespace {

    // A path segment determined in terms of the start and end point.
    template<typename T>
    struct PathSegment final {
        PathSegment(T startPoint, T endPoint) :
            start(std::move(startPoint)),
            end(std::move(endPoint))
        {}

        template<typename X>
        PathSegment(PathSegment<X> other) :
            PathSegment {
                T(other.start),
                T(other.end)
            }
        {}

        T start {};
        T end {};
    };

    using PathSegmentVec3 = PathSegment<SimTK::Vec3>;

    PathSegmentVec3 operator*(const SimTK::Rotation& rot, const PathSegmentVec3& path) {
        return PathSegmentVec3 {
            rot * path.start,
            rot * path.end
          };
    }

    // Returns direction of shortest angular distance for a vector aligned with the
    // start point to become aligned with the end point.
    bool DirectionOfShortestAngularDistance(SimTK::Vec2 start, SimTK::Vec2 end) {
        return std::abs(std::atan2(end[1], end[0])
            - std::atan2(start[1], start[0])) <= M_PI;
    }

    bool DirectionOfShortestAngularDistanceAboutZAxis(const PathSegmentVec3& path) {
        return DirectionOfShortestAngularDistance(
            path.start.getSubVec<2>(0),
            path.end.getSubVec<2>(0));
    }

    // Angular distance from start- to end-angle in either positive or negative direction.
    double AngularDistance(double startAngle, double endAngle, bool positiveRotation) {
        const double distance = std::fmod(endAngle, 2. * M_PI)
            - std::fmod(startAngle, 2. * M_PI);
        if ((distance < 0.) && positiveRotation) {
            return distance + 2.*M_PI;
        }
        if ((distance > 0.) && !positiveRotation) {
            return distance - 2.*M_PI;
        }
        return distance;
    }

    // Representation of a point in space using cylindrical coordinates.
    struct CylindricalCoordinates final {

        CylindricalCoordinates(const SimTK::Vec3& point) :
            radius(SimTK::Vec2(point[0], point[1]).norm()),
            angle(std::atan2(point[1], point[0])),
            axialPosition(point[2])
        {}

        double radius;
        double angle;
        double axialPosition;
    };

    // A geodesic path definition on a cylindrical surface.
    //
    // The wrapping problem is not solved here, instead the geodesic is constructed
    // given a wrapping result, and allows for the surface gradient and path length
    // to be evaluated. These values can then be used for asserting correctness of
    // the given wrapping result.
    class CylinderGeodesic final {

        CylinderGeodesic(
            PathSegment<CylindricalCoordinates> wrappedPath,
            bool wrappingDirection
        ) :
            _startPoint(
                wrappedPath.start),
            _axialDistance(
                wrappedPath.end.axialPosition -
                wrappedPath.start.axialPosition),
            _angularDistance(
                AngularDistance(
                    wrappedPath.start.angle,
                    wrappedPath.end.angle,
                    wrappingDirection)),
            _cylinderOrientation(
                SimTK::Rotation())
        {}

        CylinderGeodesic(
            PathSegmentVec3 wrappedPath,
            bool wrapping_direction
        ) :
            CylinderGeodesic(PathSegment<CylindricalCoordinates>(wrappedPath),
                wrapping_direction)
        {}

        public:
        CylinderGeodesic(
            PathSegmentVec3 wrappedPath,
            bool wrapping_direction,
            SimTK::Rotation cylinder_orientation
        ) :
            CylinderGeodesic(
                cylinder_orientation.invert() * wrappedPath,
                wrapping_direction)
        {
            _cylinderOrientation = cylinder_orientation;
        }

        private:

        // Computes the surface gradient along the geodesic path (start: t=0, end: t=1).
        SimTK::Vec3 computeSurfaceGradient(double t) const
        {
            SimTK::Rotation rotZAxis;
            rotZAxis.setRotationFromAngleAboutZ(_angularDistance * t + _startPoint.angle);
            return _cylinderOrientation *
                rotZAxis * (
                    SimTK::Vec3(0., 1., 0.) * _startPoint.radius
                    * _angularDistance
                    + SimTK::Vec3(0., 0., 1.)
                    * _axialDistance
                );
        }

        public:
        // Compute surface gradient direction at start of geodesic path.
        SimTK::UnitVec3 computeGradientDirectionAtStart() const
        {
            return SimTK::UnitVec3(computeSurfaceGradient(0.));
        }

        // Compute surface gradient direction at end of geodesic path.
        SimTK::UnitVec3 computeGradientDirectionAtEnd() const
        {
            return SimTK::UnitVec3(computeSurfaceGradient(1.));
        }

        // Compute geodesic path length.
        double computeLength() const
        {
            return std::sqrt(
                    std::pow(_angularDistance * _startPoint.radius, 2) +
                    std::pow(_axialDistance, 2)
                );
        }

        double getRadius() const {
            return _startPoint.radius;
        }

        private:
        // Local coordinates of starting point of wrapped path.
        CylindricalCoordinates _startPoint;
        // Distance from start to end of wrapped path.
        double _axialDistance;
        double _angularDistance;
        // Cylinder orientation wrt ground frame.
        SimTK::Rotation _cylinderOrientation;
    };

}

// Section with helpers for evaluating the wrapping result in the upcoming tests.
namespace {

    // Struct for holding the result of the wrapping solution.
    struct WrapTestResult final {

        static WrapTestResult NoWrap() {
            WrapTestResult noWrap;
            noWrap.noWrap = true;
            return noWrap;
        }

        // Path points on cylinder surface.
        PathSegmentVec3 path = {{}, {}};
        // Path segment length on cytlinder surface.
        double length = NAN;
        // Direction of wrapping wrt cylinder axis.
        bool positiveDirection = true;
        // If there is no wrapping (the other fields don't matter).
        bool noWrap = false;
    };

    std::ostream& operator<<(std::ostream& os, const WrapTestResult& result) {
        if (result.noWrap) {
            return os << "WrapTestResult { noWrap = true }";
        }
        return os <<
            "WrapTestResult{" <<
            "path.start: " << result.path.start << ", " <<
            "path.end: " << result.path.end << ", " <<
            "length: " << result.length << ", " <<
            "positiveDirection: " << result.positiveDirection << "}";
    }

    // Struct holding the tolerances when asserting the wrapping result.
    struct WrappingTolerances final {
        WrappingTolerances() = default;

        WrappingTolerances(double eps):
            position(eps),
            length(eps),
            gradientDirection(eps)
        {}

        double position = 1e-3;
        double length = 1e-9;
        double gradientDirection = 1e-3;
    };

    std::ostream& operator<<(std::ostream& os, const WrappingTolerances& tol) {
        return os <<
            "WrappingTolerances {" <<
            "position: " << tol.position << ", " <<
            "length: " << tol.length << ", " <<
            "gradientDirection: " << tol.gradientDirection << "}";
    }

    bool IsEqual(double lhs, double rhs, double tolerance) {
        return std::abs(lhs - rhs) < tolerance;
    }

    bool IsEqual(const SimTK::Vec3& lhs, const SimTK::Vec3& rhs, double tolerance) {
        return IsEqual(lhs[0], rhs[0], tolerance)
            && IsEqual(lhs[1], rhs[1], tolerance)
            && IsEqual(lhs[2], rhs[2], tolerance);
    }

    bool IsEqual(const PathSegmentVec3& lhs, const PathSegmentVec3& rhs, double tolerance) {
        return IsEqual(lhs.start, rhs.start, tolerance)
            && IsEqual(lhs.end, rhs.end, tolerance);
    }

    bool IsEqual(const WrapTestResult& lhs, const WrapTestResult& rhs, double tolerance) {
        if (lhs.noWrap && rhs.noWrap) {
            return true;
        }
        return IsEqual(lhs.path, rhs.path, tolerance)
            && IsEqual(lhs.length, rhs.length, tolerance) // TODO: or use Tolerance.length_error
            && lhs.positiveDirection == rhs.positiveDirection
            && lhs.noWrap == rhs.noWrap;
    }

    double ErrorInfinityNorm(const SimTK::UnitVec3& lhs, const SimTK::UnitVec3& rhs) {
        return SimTK::max(( lhs.asVec3() - rhs.asVec3() ).abs());
    }
}

// Section on configuring and simulating a specific wrapping scenario.
namespace {

    // Parameterization of a wrapping scenario.
    // The simulated scene is one with a cylinder at the origin, with a user
    // defined radius, length, active quadrant, and relative orientation.
    // A specific wrapping case can be configured by choosing the start and end
    // points of the total path.
    struct WrapInput final {

        SimTK::Rotation cylinderOrientation() const {
            SimTK::Rotation orientation;
            orientation.setRotationToBodyFixedXYZ(eulerRotations);
            return orientation;
        }

        // Endpoints of the total path:
        PathSegmentVec3 path = {{}, {}};

        // Wrapping cylinder parameters:
        double radius = 1.;
        double cylinderLength = 1.;
        std::string quadrant = "all";
        SimTK::Vec3 eulerRotations = {0., 0., 0.};
    };

    std::ostream& operator<<(std::ostream& os, const WrapInput& input) {
        return os <<
            "WrapInput{" <<
            "path.start: " << input.path.start << ", " <<
            "path.end: " << input.path.end << ", " <<
            "radius: " << input.radius << ", " <<
            "cylinderLength: " << input.cylinderLength << ", " <<
            "quadrant: " << input.quadrant << ", " <<
            "eulerRotations: " << input.eulerRotations << "}";
    }

    // Simulates the wrapping scenario as configured by the input.
    WrapTestResult solve(const WrapInput& input, bool visualize) {
        Model model;
        model.setName("testWrapCylinderModel");

        std::unique_ptr<WrapCylinder> cylinder(new WrapCylinder());
        cylinder->setName("cylinder");
        cylinder->set_radius(input.radius);
        cylinder->set_length(input.cylinderLength);
        cylinder->set_quadrant(input.quadrant);
        cylinder->set_xyz_body_rotation(input.eulerRotations);
        cylinder->setFrame(model.getGround());

        std::unique_ptr<PathSpring> spring (new PathSpring("spring", 1., 1., 1.));
        spring->updGeometryPath().appendNewPathPoint(
            "start_point",
            model.get_ground(),
            input.path.start);
        spring->updGeometryPath().appendNewPathPoint(
            "end_point",
            model.get_ground(),
            input.path.end);
        spring->updGeometryPath().addPathWrap(*cylinder.get());

        model.addComponent(spring.release());
        model.updGround().addWrapObject(cylinder.release());

        model.finalizeConnections();
        model.setUseVisualizer(visualize);

        const SimTK::State& state = model.initSystem();

        if (visualize) {
            model.realizeVelocity(state);
            model.getVisualizer().show(state);
        }

        model.updComponent<PathSpring>("spring").getLength(state);
        WrapResult wrapResult = model.updComponent<PathSpring>("spring")
                .getGeometryPath()
                .getWrapSet()
                .get("pathwrap")
                .getPreviousWrap();

        // Convert the WrapResult to a WrapTestResult for convenient assertions.
        WrapTestResult result;
        result.path = PathSegmentVec3{
            wrapResult.r1,
            wrapResult.r2,
        };
        result.length = wrapResult.wrap_path_length;
        // Determine the wrapping sign based on the first path segment.
        result.positiveDirection = DirectionOfShortestAngularDistanceAboutZAxis(
            input.cylinderOrientation().invert() * PathSegmentVec3{
                input.path.start,
                wrapResult.r1,
            });
        result.noWrap = wrapResult.wrap_pts.size() == 0;

        return result;
    }

}

namespace {

    std::string TestGeodesicProperties(
        const WrapInput& input,
        const WrapTestResult& result,
        const WrappingTolerances& tol,
        const std::string& name)
    {
        std::ostringstream oss;
        std::string delim = "\n        ";

        // =====================================================================
        // ============ Test: Wrapped points distance to surface. ==============
        // =====================================================================

        // We can assert if the wrapped path points lie on the surface of the
        // cylinder, by computing the radial position in cylindrical
        // coordinates.

        CylinderGeodesic geodesicPath(
            result.path,
            result.positiveDirection,
            input.cylinderOrientation()
        );

        double error = std::abs(geodesicPath.getRadius() - input.radius);
        if (error > tol.position) {
            oss << delim << name << ": Distance to surface error = " << error <<
                " exceeds tolerance = " << tol.position;
        }

        // =====================================================================
        // ===================== Test: Wrapping length. ========================
        // =====================================================================

        // We can test the path length by recomputing it from the path points
        // via the CylinderGeodesic.

        error = std::abs(geodesicPath.computeLength() - result.length);
        if (error > tol.length) {
            oss << delim << name << ": Length error = " << error <<
                " exceeds tolerance = " << tol.length;
        }

        // =====================================================================
        // ===================== Test: Tangency to surface =====================
        // =====================================================================

        // The path must at all times be tangent to the cylinder surface. We can
        // test that the straight-line segments, before and after wrapping the
        // surface, are also tangent to the surface.

        // Straight line segment from path start point to surface must be
        // tangent to cylinder surface.
        error = ErrorInfinityNorm(
            geodesicPath.computeGradientDirectionAtStart(),
            SimTK::UnitVec3(result.path.start - input.path.start));
        if (error > tol.gradientDirection) {
            oss << delim << name << ": Start segment gradient error = " << error
                << " exceeds tolerance = " << tol.gradientDirection;
        }

        // Straight line segment from surface to path end point must be tangent
        // to cylinder surface.
        error = ErrorInfinityNorm(
            geodesicPath.computeGradientDirectionAtEnd(),
            SimTK::UnitVec3(input.path.end - result.path.end));
        if (error > tol.gradientDirection) {
            oss << delim << name << ": End segment gradient error = " << error
                << " exceeds tolerance = " << tol.gradientDirection;
        }

        return oss.str();
    }

    // Simulates and tests a wrapping case:
    // - Computes the wrapping result given the input.
    // - Tests if computed result is equal to the expected result.
    // - Tests if geodesic properties of computed result are correct.
    // - Tests if geodesic properties of expected result are correct.
    // - Returns a string with the reason of failure, if any.
    std::string TestWrapping(
        const WrapInput& input,
        const WrapTestResult& expected,
        WrappingTolerances tol,
        std::string testCase,
        bool visualize = false)
    {

        WrapTestResult result = solve(input, visualize);

        std::ostringstream failedBuf;

        if (!IsEqual(result, expected, tol.position)) {
            failedBuf << "\n        ";
            failedBuf << "Expected and simulated result not equal to within tolerance = " << tol.position;
        }

        if (!result.noWrap) {
            failedBuf << TestGeodesicProperties(
                input,
                result,
                tol,
                "Wrapping result"
            );
        }

        // The supplied expected result should also pass the test.
        if (!expected.noWrap) {
            failedBuf << TestGeodesicProperties(
                input,
                expected,
                WrappingTolerances(1e-13),
                "Supplied expected path"
            );
        }

        std::string failedStr = failedBuf.str();
        if (failedStr.empty()) {
            return failedStr;
        }

        // If failed, print info related to this test case.
        std::ostringstream oss;
        oss << "\nFAILED: case = " << testCase;
        oss << "\n    input = " << input;
        oss << "\n    result = " << result;
        oss << "\n    expected = " << expected;
        oss << "\n    caused by:" << failedStr;
        return oss.str();
    }

}

int main()
{
    // All tests should be invariant to rotating the entire scene.
    auto TestWrappingForDifferentOrientations = [](
        const WrapInput& input,
        const WrapTestResult& expected,
        const WrappingTolerances& tolerance,
        const std::string& name
    ) -> std::string {
        std::string log;

        // A list of interesting rotation perturbations.
        std::vector<SimTK::Vec3> eulerRotationsList {
            {0., 0., 0.,}, // The original case.
            {1., 0., 0.,},
            {0., 1., 0.,},
            {0., 0., 1.,},
            {1., 1., 1.,},
            {-1., 0., 0.,},
            {0., -1., 0.,},
            {0., 0., -1.,},
            {-1., -1., -1.,},
        };

        // Rotating both the input and the expected result should not affect
        // passing the test.
        for (ptrdiff_t i = 0; i < eulerRotationsList.size(); ++i) {
            WrapInput inputRotated = input;
            inputRotated.eulerRotations = eulerRotationsList[i];
            inputRotated.path = inputRotated.cylinderOrientation() * input.path;
            WrapTestResult expectedRotated = expected;
            expectedRotated.path = inputRotated.cylinderOrientation() * expected.path;
            log += TestWrapping(
                inputRotated,
                expectedRotated,
                tolerance,
                (i == 0)? name: name + " (rotated)");

            // Break on first failing test.
            if (!log.empty()) {
                break;
            }
        }

        return log;
    };

    // Some tests should pass for all quadrants.
    auto TestWrappingForAllQuadrants = [&TestWrappingForDifferentOrientations](
        WrapInput input,
        const WrapTestResult& expected,
        const WrappingTolerances& tolerance,
        const std::string& name
    ) -> std::string {
        std::string failLog;
        std::vector<std::string> quadrantsList{
            "+x",
            "-x",
            "+y",
            "-y",
            "all",
        };
        for (ptrdiff_t i = 0; i < quadrantsList.size(); ++i) {
            input.quadrant = quadrantsList[i];
            failLog += TestWrappingForDifferentOrientations(
                input,
                expected,
                tolerance,
                name + " (quadrant " + input.quadrant + ")");
        }
        return failLog;
    };

    std::string failLog;
    WrappingTolerances tolerance;

    {
        std::string name = "Perpendicular unconstrained wrap";

        WrapInput input;
        input.path = PathSegmentVec3{
            {5., 0.1, 0.},
            {-5., 0., 0.},
        };

        WrapTestResult expected;
        expected.path = PathSegmentVec3{
            {0.180327868852459, 0.983606557377049, 0},
            {-0.2,              0.979795897113271, 0},
        };
        expected.length = 0.382677695191821;

        failLog += TestWrappingForDifferentOrientations(input, expected, tolerance, name);
    }

    {
        std::string name = "Unconstrained wrap";

        WrapInput input;
        input.path = PathSegmentVec3{
            {5., 0.1, 5},
            {-5., 1., 1},
        };
        input.cylinderLength = 10;

        WrapTestResult expected;
        expected.path = PathSegmentVec3{
            {0.180327868852459,     0.983606557377049, 3.05581010833836},
            {-1.60812264967664e-16, 1,                 2.98386723638943},
        };
        expected.length = 0.195070852290308;

        failLog += TestWrapping(input, expected, tolerance, name);
    }

    {
        std::string name = "Insider cylinder";

        WrapInput input;
        input.path = PathSegmentVec3{
            {1., 0., 0.},
            {0., 2., 0.},
        };

        WrapTestResult expected = WrapTestResult::NoWrap();

        // Both points inside the cylinder.
        input.radius = 3.;
        failLog += TestWrappingForAllQuadrants(input, expected, tolerance, name);

        // First point inside the cylinder.
        input.radius = 1.5;
        failLog += TestWrappingForAllQuadrants(input, expected, tolerance, name);

        // Second point inside the cylinder.
        input.path = PathSegmentVec3{
            {2., 0., 0.},
            {0., 1., 0.},
        };
        failLog += TestWrappingForAllQuadrants(input, expected, tolerance, name);
    }

    {
        std::string name = "No intersecting straight line";

        WrapInput input;
        input.path = PathSegmentVec3{
            {2., 0., 0.},
            {0., 2., 0.},
        };

        WrapTestResult expected = WrapTestResult::NoWrap();

        failLog += TestWrapping(input, expected, tolerance, name);
    }

    {
        std::string name = "Missed cylinder in axial direction";

        WrapInput input;
        input.path = PathSegmentVec3{
            {1., 0.1, 10.},
            {-1., 0., 10.},
        };

        WrapTestResult expected = WrapTestResult::NoWrap();

        failLog += TestWrapping(input, expected, tolerance, name);
    }

    if (!failLog.empty()) {
        std::cerr << failLog << std::endl;
        return 1;
    }

    std::cout << "Wrap Cylinder Test Done." << std::endl;
    return 0;
}