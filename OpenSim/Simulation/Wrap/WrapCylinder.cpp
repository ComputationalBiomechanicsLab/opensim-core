/* -------------------------------------------------------------------------- *
 *                         OpenSim:  WrapCylinder.cpp                         *
 * -------------------------------------------------------------------------- *
 * The OpenSim API is a toolkit for musculoskeletal modeling and simulation.  *
 * See http://opensim.stanford.edu and the NOTICE file for more information.  *
 * OpenSim is developed at Stanford University and supported by the US        *
 * National Institutes of Health (U54 GM072970, R24 HD065690) and by DARPA    *
 * through the Warrior Web program.                                           *
 *                                                                            *
 * Copyright (c) 2005-2017 Stanford University and the Authors                *
 * Author(s): Peter Loan                                                      *
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

//=============================================================================
// INCLUDES
//=============================================================================
#include <SimTKcommon.h>
#include "WrapCylinder.h"
#include "PathWrap.h"
#include "WrapMath.h"
#include "WrapResult.h"
#include <OpenSim/Simulation/Model/Model.h>
#include <OpenSim/Common/ModelDisplayHints.h>
#include <OpenSim/Common/SimmMacros.h>
#include <OpenSim/Common/ScaleSet.h>
#include <functional>


//=============================================================================
// STATICS
//=============================================================================
using namespace std;
using namespace OpenSim;
using SimTK::Vec3;
using SimTK::UnitVec3;

static const char* wrapTypeName = "cylinder";
static Vec3 p0(0.0, 0.0, -1.0);
static Vec3 dn(0.0, 0.0, 1.0);
#define MAX_ITERATIONS    100
#define TANGENCY_THRESHOLD (0.1 * SimTK_DEGREE_TO_RADIAN) // find tangency to within 1 degree

// Geometry related helper structs and functions.
namespace {
    constexpr double c_TAU = 2. * SimTK_PI;

    template<typename T>
    struct PathSegment final {
        PathSegment(T startPoint, T endPoint) :
            start(std::move(startPoint)),
            end(std::move(endPoint))
        {}

        template<typename X, typename ...Args>
        PathSegment(
            const PathSegment<X>& other,
            Args&& ...args) :
            PathSegment{
                T(other.start, std::forward<Args>(args)...),
                T(other.end, std::forward<Args>(args)...)
            }
        {}

        T start{};
        T end{};
    };

    template<typename T>
    T Difference(const PathSegment<T>& path)
    {
        return path.end - path.start;
    }

    using PathSegmentVec2 = PathSegment<SimTK::Vec2>;
    using PathSegmentVec3 = PathSegment<SimTK::Vec3>;

    // Applies atan2 but maps the resulting range to [0...2pi].
    double ATan2OnPositiveRange(const SimTK::Vec2 point)
    {
        const double angle = std::atan2(point[1], point[0]);
        return angle < 0.? angle + c_TAU : angle;
    }

    class Angle final {
        // Construct from angle already on range [0...2pi].
        Angle(double angleChecked): _value(angleChecked)
        {}

        public:
        Angle() = default;

        Angle(const SimTK::Vec2& point) :
            _value(ATan2OnPositiveRange(point))
        {}

        // Returns a value on range [0..2pi], or NaN if uninitialized.
        double getValue() const {
            return _value;
        }

        friend Angle operator+(Angle lhs, Angle rhs);
        friend Angle operator-(Angle lhs, Angle rhs);

        private:
        double _value = SimTK::NaN;
    };

    Angle operator+(Angle lhs, Angle rhs) {
        double angle = lhs.getValue() + rhs.getValue();
        // Map angle to [0...2pi] range. Using fmod would be more robust against
        // repeated calls, but significantly slower. For the current code this
        // robustness would be overkill, so we go for speed:
        angle = angle > c_TAU? angle - c_TAU : angle; // poor man's fmod()
        return Angle{angle};
    }

    Angle operator-(Angle lhs, Angle rhs) {
        double angle = lhs.getValue() - rhs.getValue();
        // Map angle to [0...2pi] range.
        angle = angle < 0.? angle + c_TAU : angle; // poor man's fmod()
        return Angle{angle};
    }

    struct PolarCoordinates final {
        PolarCoordinates(const SimTK::Vec2& point) :
            angle(point),
            magnitude(point.norm())
        {}

        Angle angle{};
        double magnitude = SimTK::NaN;
    };

    enum class RotationDirection {
        Positive,
        Negative,
    };

    // Angular distance from start- to end angle, given a direction of rotation.
    double AngularDistance(
        Angle startAngle,
        Angle endAngle,
        RotationDirection direction)
    {
        const double diff = endAngle.getValue() - startAngle.getValue();
        if (diff < 0 && direction == RotationDirection::Positive) {
            return diff + c_TAU;
        }
        if (diff > 0 && direction == RotationDirection::Negative) {
            return diff - c_TAU;
        }
        return diff;
    }

    double AngularDistance(
        PathSegment<Angle> path,
        RotationDirection direction)
    {
        return AngularDistance(path.start, path.end, direction);
    }

    // Returns the point on the line A-B that is closest to the origin.
    SimTK::Vec2 ClosestPointOnLineAB(
        const SimTK::Vec2& pointA,
        const SimTK::Vec2& pointB)
    {
        // Given the infinite line through A anb B: If points A and B are on the
        // same side of the shortest distance vector to the infinite line, we
        // know that the shortest distance to line A-B is simply the distance to
        // either point A or point B.
        const SimTK::Vec2 vecAB = pointB - pointA;
        const double vTa = SimTK::dot(vecAB, pointA);
        const double vTb = SimTK::dot(vecAB, pointB);
        const bool closestPointIsEitherAOrB =  vTb * vTa > 0.;

        // Return closest of A and B.
        if (closestPointIsEitherAOrB) {
            const double normASquared = SimTK::dot(pointA, pointA);
            const double normBSquared = SimTK::dot(pointB, pointB);
            return normASquared < normBSquared? pointA : pointB;
        }

        // Compute closest point on infinite line through A and B.
        const double vTv = SimTK::dot(vecAB, vecAB);
        return pointA - vecAB * vTa / vTv;
    }

    // Check if either start or end is inside the circle radius.
    bool IsPointInsideCircle(
        const SimTK::Vec2& point,
        double radiusSquared)
    {
        return SimTK::dot(point, point) < radiusSquared;
    }
}

// Construction of 2D Tangent lines to a circle.
namespace {

    // Representation of a line obtained by connecting a point
    // to a circle, such that the line is tangent to the circle.
    struct TangentLine final {
        operator Angle() const
        {
            return angleOfTangentPoint;
        }

        // Angle at which the line touches the circle.
        Angle angleOfTangentPoint{};
        // Length of the line segment.
        double length = SimTK::NaN;
    };

    // Representation of the two possible TangentLines conntecting a point and
    // circle.
    struct BothTangentLines final {
        BothTangentLines(const PolarCoordinates& point, double circleRadius) :
            _length(SimTK::Vec2{point.magnitude, circleRadius}.norm()),
            _angleMid(point.angle),
            _angleDiff(SimTK::Vec2(_length, circleRadius))
        {}

        BothTangentLines(const SimTK::Vec2& point, double circleRadius) :
            BothTangentLines(PolarCoordinates(point), circleRadius)
        {}

        private:
        TangentLine getTangentLineHelper(bool invertSign) const {
            TangentLine tangentLine;
            tangentLine.angleOfTangentPoint = invertSign?
                _angleMid - _angleDiff : _angleMid + _angleDiff,
            tangentLine.length = _length;
            return tangentLine;
        }

        public:
        TangentLine
            getTangentLineStartingAtPoint(RotationDirection direction) const
        {
            return getTangentLineHelper(
                direction == RotationDirection::Positive);
        }

        TangentLine
            getTangentLineEndingAtPoint(RotationDirection direction) const
        {
            return getTangentLineHelper(
                direction == RotationDirection::Negative);
        }

        private:
        // Line segment length, given either direction.
        double _length = SimTK::NaN;
        // Helpers for constructing the angle of the tangent point.
        Angle _angleMid;
        Angle _angleDiff;
    };
}

// Section defining the active quadrant, and some helpers to check if a point
// lies in it.
namespace {

    enum class WrapSign {
        Unconstrained,
        Positive,
        Negative,
    };

    enum class WrapAxis {
        XAxis,
        YAxis,
    };

    struct Quadrant {

        Quadrant(int wrapSign, int wrapAxis)
        {
            switch (wrapAxis) {
                case 0:
                    axis = WrapAxis::XAxis;
                    break;
                case 1:
                    axis = WrapAxis::YAxis;
                    break;
                default: // todo properly throw...
                    throw std::string("Invalid wrapAxis for wrapCylinder");
            }
            switch (wrapSign) {
                case 0:
                    sign = WrapSign::Unconstrained;
                    break;
                case 1:
                    sign = WrapSign::Positive;
                    break;
                case -1:
                    sign = WrapSign::Negative;
                    break;
                default: // TODO properly throw..
                    throw std::string("Invalid wrapSign for wrapCylinder");
            }
        }

        WrapAxis axis;
        WrapSign sign;
    };

    bool IsAngleInPositiveYQuadrant(Angle angle) {
        return angle.getValue() <= SimTK_PI;
    }

    bool IsAngleInPositiveXQuadrant(Angle angle) {
        return angle.getValue() <= SimTK_PI * 0.5
            || angle.getValue() >= SimTK_PI * 1.5;
    }

    bool IsAngleInQuadrant(
        Angle angle,
        const Quadrant& quadrant)
    {
        bool inPositiveQuadrant = quadrant.axis == WrapAxis::XAxis?
            IsAngleInPositiveXQuadrant(angle):
            IsAngleInPositiveYQuadrant(angle);
        return (quadrant.sign == WrapSign::Positive)?
            inPositiveQuadrant:
            !inPositiveQuadrant;
    }

    bool IsPointInActiveQuadrant(
        const SimTK::Vec2& point,
        Quadrant quadrant)
    {
        const double coord = quadrant.axis == WrapAxis::XAxis? point[0] : point[1];
        return coord > 0. && quadrant.sign == WrapSign::Positive;
    }
}

// Section on wrapping logic of a 2D circle.
namespace {

    // Returns if wrapping should not even be attempted.
    //
    // Checks 4 conditions:
    // - 1. If either point lies inside the circle.
    // - 2. If unconstrained, and straight line from start to end does not
    //      intersect the circle.
    // - 3. If both points lie in the active quadrant.
    // - 4. If one point lies in the non-active quadrant, and closest point in
    //      non-active quadrant.
    bool DoNotAttemptWrapping(
        const PathSegmentVec2& path,
        double radius,
        Quadrant quadrant)
    {
        // Condition 1: Points must not lie inside the circle.
        const double radiusSquared = radius * radius;
        bool doNotAttemptWrapping =
            IsPointInsideCircle(path.start, radiusSquared) ||
            IsPointInsideCircle(path.end, radiusSquared);

        // Condition 2: Do not wrap if unconstrained and the straight line connecting
        // start-end does not intersect the circle.
        const SimTK::Vec2 closestPointOnLine = ClosestPointOnLineAB(
            path.start,
            path.end);
        if (quadrant.sign == WrapSign::Unconstrained) {
            const double normClosestPointOnLineSquared = SimTK::dot(
                closestPointOnLine,
                closestPointOnLine);
            doNotAttemptWrapping |= normClosestPointOnLineSquared > radiusSquared;
            return doNotAttemptWrapping;
        }

        // Condition 3: Do not wrap if both points lie in the active quadrant.
        const bool startInActiveQuadrant =
            IsPointInActiveQuadrant(path.start, quadrant);
        const bool endInActiveQuadrant =
            IsPointInActiveQuadrant(path.end, quadrant);
        doNotAttemptWrapping |= startInActiveQuadrant && endInActiveQuadrant;

        // Condition 4: Do not wrap if only one point in non-active, and nearest
        // in active quadrant.
        const bool closestPointInActiveQuadrant = IsPointInActiveQuadrant(
            closestPointOnLine,
            quadrant);
        doNotAttemptWrapping |= startInActiveQuadrant != endInActiveQuadrant
            && closestPointInActiveQuadrant;

        return doNotAttemptWrapping;
    }

    // Given the wrapping path on a circle, determine if it is valid based on
    // the quadrant.
    //
    // Wrapped path is valid if:
    // 1. Both angles in non-active quadrant.
    // 2. Both angles in active-quadrant, while not wrapping the entire
    //    non-active quadrant.
    bool IsValidQuadrantWrapped(
        const PathSegment<Angle>& path,
        double angularDistance,
        const Quadrant& quadrant)
    {
        const bool startOk = IsAngleInQuadrant(path.start, quadrant);
        const bool endOk = IsAngleInQuadrant(path.end, quadrant);

        // Condition 1: Both angles in non-active.
        bool validWrap = !startOk && !endOk;

        // Condition 2: Both in active, and wrapping entire non-active.
        validWrap |= startOk && endOk && std::abs(angularDistance) < SimTK::Pi;
        return validWrap;
    }

}

// Section on wrapping of a 2D circle.
namespace {

    struct CircleWrap final {
        CircleWrap() = default;

        CircleWrap(
            RotationDirection direction,
            PathSegment<TangentLine> path,
            Quadrant quadrant) :
            wrappedAngularDistance(AngularDistance(
                path,
                direction)),
            tangentLineSegments(path),
            wrapIsValid(IsValidQuadrantWrapped(
                path,
                wrappedAngularDistance,
                quadrant))
        {}

        static CircleWrap NoWrap() {
            return CircleWrap{};
        }

        double getWrappedPathStartAngle() const {
            return tangentLineSegments.start.angleOfTangentPoint.getValue();
        }

        double wrappedAngularDistance = SimTK::NaN;
        PathSegment<TangentLine> tangentLineSegments{{}, {}};
        bool wrapIsValid = false;
    };

    double ComputeTotalPathLength(
        const CircleWrap& wrap,
        double circleRadius)
    {
        return wrap.tangentLineSegments.start.length
            + wrap.tangentLineSegments.end.length
            + circleRadius * std::abs(wrap.wrappedAngularDistance);
    }

}

// Section regarding the two possible directions in which you can compute a
// wrapping path on a circle.
namespace {

    // Struct for holding a generic pair, that differ in the wrapping direction.
    template<typename T>
    struct PositiveAndNegativeRotatingPair final {
        PositiveAndNegativeRotatingPair(
            T positiveRotatingVersion,
            T negativeRotatingVersion) :
            positive(std::move(positiveRotatingVersion)),
            negative(std::move(negativeRotatingVersion))
        {}

        PositiveAndNegativeRotatingPair(std::function<T(RotationDirection)>&& f) :
            PositiveAndNegativeRotatingPair
            {
                f(RotationDirection::Positive),
                f(RotationDirection::Negative)
            }
        {}

        template<typename X, typename ...Args>
        PositiveAndNegativeRotatingPair(
            const PositiveAndNegativeRotatingPair<X>& other,
            Args&& ...args) :
            PositiveAndNegativeRotatingPair
            {
                T(RotationDirection::Positive, other.positive, std::forward<Args>(args)...),
                T(RotationDirection::Negative, other.negative, std::forward<Args>(args)...)
            }
        {}

        T positive;
        T negative;
    };

    // Computes for both positive and negative wrapping directions, the tangent
    // line path segments.
    PositiveAndNegativeRotatingPair<PathSegment<TangentLine>>
        ComputePossibleTangentLineSegments(
            const PathSegmentVec2& path,
            double circleRadius)
    {
        // Compute for each path segment (start/end) both possible tangent lines
        // (positive/negative).
        const PathSegment<BothTangentLines> tangentLines(path, circleRadius);

        // Get one of the possible paths by setting the rotation direction.
        auto GetPathSegmentsGivenDirection = [&tangentLines](
                RotationDirection direction)->PathSegment<TangentLine>
            {
                return {
                    tangentLines.start.getTangentLineStartingAtPoint(direction),
                    tangentLines.end.getTangentLineEndingAtPoint(direction),
                };
            };
        return {GetPathSegmentsGivenDirection};
    }

    // Computes for both positive and negative wrapping directions, the wrapping
    // path on a circle.
    PositiveAndNegativeRotatingPair<CircleWrap>
        ComputePossibleCircleWrapPaths(
            const PathSegmentVec2& path,
            double circleRadius,
            Quadrant quadrant)
    {
        const PositiveAndNegativeRotatingPair<PathSegment<TangentLine>>
            possibleTangentSegments = ComputePossibleTangentLineSegments(
                path,
                circleRadius);

        return {possibleTangentSegments, quadrant};
    }

}

// This section wraps up the 2D circle wrapping, by choosing the best path from
// two possible wrapping paths.
namespace {

    CircleWrap TakeBestPath(
        const CircleWrap& wrapA,
        const CircleWrap& wrapB)
    {
        // Neither wrapping path is valid.
        if (!wrapA.wrapIsValid && !wrapB.wrapIsValid) {
            return CircleWrap::NoWrap();
        }
        // Both are valid: Take the shortest path.
        if (wrapA.wrapIsValid && wrapB.wrapIsValid) {
            return std::abs(wrapA.wrappedAngularDistance)
                < std::abs(wrapB.wrappedAngularDistance)?
                wrapA : wrapB;
        }
        // Return the path that is valid, if any.
        return wrapA.wrapIsValid? wrapA : wrapB;
    }

    CircleWrap TakeBestPath(
        PositiveAndNegativeRotatingPair<CircleWrap> possibleWrappings)
    {
        return TakeBestPath(possibleWrappings.positive,
            possibleWrappings.negative);
    }

    // Computes the best wrapping path on a circle.
    // This function includes the logic for checking the active quadrant, and
    // choosing from both possible wrapping directions.
    CircleWrap ComputeBestCircleWrappingPath(
        const PathSegmentVec2& path,
        double circleRadius,
        Quadrant quadrant)
    {
        if (DoNotAttemptWrapping(path, circleRadius, quadrant)) {
            return CircleWrap::NoWrap();
        }

        const PositiveAndNegativeRotatingPair<CircleWrap> possiblePaths =
            ComputePossibleCircleWrapPaths(
                path,
                circleRadius,
                quadrant);

        return TakeBestPath(possiblePaths);
    }

}

// Section regarding the axial component of a cylinder wrapping path: the
// AxialWrap. Axial component of a cylinder wrapping path can be computed given
// the solution of wrapping a 2D circle. The CircleWrap together with the
// AxialWrap gives the CylinderWrap.
namespace {

    PathSegment<double> ComputeAxialWrapPath(
        const CircleWrap& circleWrap,
        const PathSegment<double>& axialPath,
        double circleRadius)
    {
        // Axial-path start and end can be computed by noting that the
        // ratio of axial-distance to tangential-distance must
        // remain constant for each path segment.
        const double totalTangentialDistance = ComputeTotalPathLength(
            circleWrap,
            circleRadius);
        const double totalAxialDistance = Difference(axialPath);
        const double ratioAxialToTangentialDistance =
            totalAxialDistance / totalTangentialDistance;

        return {
            axialPath.start + circleWrap.tangentLineSegments.start.length
                * ratioAxialToTangentialDistance,
            axialPath.end - circleWrap.tangentLineSegments.end.length
                * ratioAxialToTangentialDistance,
        };
    }

    bool IsAxialPathWithinCylinderLength(
        const PathSegment<double> wrappedAxialPath,
        double cylinderLength)
    {
        return std::abs(wrappedAxialPath.start) <= cylinderLength / 2.
            && std::abs(wrappedAxialPath.end) <= cylinderLength / 2.;
    }

    struct AxialWrap {
        AxialWrap() = default;

        AxialWrap(
            const CircleWrap& circleWrap,
            const PathSegment<double>& axialPath,
            double circleRadius,
            double cylinderLength) :
            wrappedPath(ComputeAxialWrapPath(circleWrap, axialPath, circleRadius)),
            wrapIsValid(IsAxialPathWithinCylinderLength(wrappedPath, cylinderLength))
        {}

        static AxialWrap NoWrap()
        {
            return AxialWrap{};
        }

        PathSegment<double> wrappedPath{SimTK::NaN, SimTK::NaN};
        bool wrapIsValid = false;
    };

}

// Section on the cylinder wrapping: Combines the Circle- and Axial wrapping.
namespace {

    struct SpiralCylinderWrap final {

        SpiralCylinderWrap(
            const PathSegmentVec3& path,
            double radius,
            double cylinderLength,
            const Quadrant& quadrant) :
            radius(radius),
            circleWrap(ComputeBestCircleWrappingPath(
                {path.start.getSubVec<2>(0), path.end.getSubVec<2>(0)},
                radius,
                quadrant))
        {
            if (circleWrap.wrapIsValid) {
                axialWrap = AxialWrap(
                    circleWrap,
                    {path.start[2], path.end[2]},
                    radius,
                    cylinderLength);
            }
        }

        double computeWrappedPathLength() const {
            return std::sqrt(
                std::pow(radius * circleWrap.wrappedAngularDistance, 2) +
                std::pow(Difference(axialWrap.wrappedPath), 2)
            );
        }

        // Get the path point on the wrapped surface for anywhere between start
        // (t=0) and end (t=1).
        SimTK::Vec3 computeWrappedPathPoint(double t) const {
            const double angle = t * circleWrap.wrappedAngularDistance +
                circleWrap.getWrappedPathStartAngle();
            const double axialCoord = t * Difference(axialWrap.wrappedPath) +
                axialWrap.wrappedPath.start;
            return SimTK::Vec3(
                std::cos(angle) * radius,
                std::sin(angle) * radius,
                axialCoord);
        }

        int wrapAction() const {
            // TODO not using: wrapped? = successful wrap, but may not be 'best' path
            // TODO not using: insideRadius
            if (circleWrap.wrapIsValid && axialWrap.wrapIsValid) {
                return WrapObject::WrapAction::mandatoryWrap;
            }
            return WrapObject::WrapAction::noWrap;
        }

        double radius = SimTK::NaN;
        CircleWrap circleWrap = CircleWrap::NoWrap();
        AxialWrap axialWrap = AxialWrap::NoWrap();
    };

}

//=============================================================================
// CONSTRUCTOR(S) AND DESTRUCTOR
//=============================================================================
//_____________________________________________________________________________
/*
* Default constructor.
*/
WrapCylinder::WrapCylinder() : WrapObject()
{
    constructProperties();
}

//_____________________________________________________________________________
/*
* Destructor.
*/
WrapCylinder::~WrapCylinder()
{}


//=============================================================================
// CONSTRUCTION METHODS
//=============================================================================
//_____________________________________________________________________________
/**
* Connect properties to local pointers.
*/
void WrapCylinder::constructProperties()
{
    constructProperty_radius(1.0);
    constructProperty_length(1.0);
}

void WrapCylinder::extendScale(const SimTK::State& s, const ScaleSet& scaleSet)
{
    Super::extendScale(s, scaleSet);

    // Get scale factors (if an entry for the Frame's base Body exists).
    const Vec3& scaleFactors = getScaleFactors(scaleSet, getFrame());
    if (scaleFactors == ModelComponent::InvalidScaleFactors)
        return;

    // _pose.x() holds the ellipsoid's X-axis expressed in the body's reference
    // frame, _pose.y() holds the Y-axis, and _pose.z() holds the Z-axis.
    // Multiplying these vectors by the scale factor vector gives
    // localScaleVector[]. The magnitudes of the localScaleVectors gives the
    // amount to scale the cylinder in the X, Y, and Z dimensions. The wrap
    // cylinder is oriented along the Z-axis, so the length is scaled by the Z
    // scale factor, and the radius is scaled by the average of the X and Y
    // scale factors.
    Vec3 localScaleVector[3];

    localScaleVector[0] = _pose.x().elementwiseMultiply(scaleFactors);
    localScaleVector[1] = _pose.y().elementwiseMultiply(scaleFactors);
    localScaleVector[2] = _pose.z().elementwiseMultiply(scaleFactors);

    upd_radius() *= (localScaleVector[0].norm() + localScaleVector[1].norm()) * 0.5;
    upd_length() *= localScaleVector[2].norm();
}

//_____________________________________________________________________________
/*
* Perform some set up functions that happen after the
* object has been deserialized or copied.
*/
void WrapCylinder::extendFinalizeFromProperties()
{
    // Base class
    Super::extendFinalizeFromProperties();

    // maybe set a parent pointer, _body = aBody;
    OPENSIM_THROW_IF_FRMOBJ(
        get_radius() < 0,
        InvalidPropertyValue,
        getProperty_radius().getName(),
        "Radius cannot be less than zero");

    OPENSIM_THROW_IF_FRMOBJ(
        get_length() < 0,
        InvalidPropertyValue,
        getProperty_length().getName(),
        "Length cannot be less than zero");
}

//_____________________________________________________________________________
/**
 * Get the name of the type of wrap object ("cylinder" in this case)
 *
 * @return A string representing the type of wrap object
 */
const char* WrapCylinder::getWrapTypeName() const
{
    return wrapTypeName;
}

//_____________________________________________________________________________
/**
 * Get a string holding the dimensions definition that SIMM would
 * use to describe this object. This is a rather ugly convenience
 * function for outputting SIMM joint files.
 *
 * @return A string containing the dimensions of the wrap object
 */
string WrapCylinder::getDimensionsString() const
{
    stringstream dimensions;
    dimensions << "radius " << get_radius() << "\nheight " << get_length();

    return dimensions.str();
}

//=============================================================================
// WRAPPING
//=============================================================================
//_____________________________________________________________________________
/**
 * Calculate the wrapping of one line segment over the cylinder.
 *
 * @param aPoint1 One end of the line segment
 * @param aPoint2 The other end of the line segment
 * @param aPathWrap An object holding the parameters for this line/cylinder pairing
 * @param aWrapResult The result of the wrapping (tangent points, etc.)
 * @param aFlag A flag for indicating errors, etc.
 * @return The status, as a WrapAction enum
 */
int WrapCylinder::wrapLine(const SimTK::State& s, SimTK::Vec3& aPoint1, SimTK::Vec3& aPoint2,
                                    const PathWrap& aPathWrap, WrapResult& aWrapResult, bool& aFlag) const
{
    SpiralCylinderWrap wrap = SpiralCylinderWrap({aPoint1, aPoint2}, get_radius(), get_length(), {_wrapAxis, _wrapSign});

    const double& _radius = get_radius();
    double dist, p11_dist, p22_dist, t, dot1, dot2, dot3, dot4, d, sin_theta,
        /* *r11, *r22, */alpha, beta, r_squared = _radius * _radius;
    double dist1, dist2;
    double t12, t00;

    Vec3 pp, vv, uu, r1a, r1b, r2a, r2b, apex, sum_musc, 
        r1am, r1bm, r2am, r2bm, p11, p22, r1p, r2p, axispt, near12, 
        vert1, vert2, mpt, apex1, apex2, l1, l2, near00;
    UnitVec3 plane_normal;
    int i, return_code = wrapped;
    bool r1_inter, r2_inter;
    bool constrained   = (bool) (_wrapSign != 0);
    bool far_side_wrap = false, long_wrap = false;

    // In case you need any variables from the previous wrap, copy them from
    // the PathWrap into the WrapResult, re-normalizing the ones that were
    // un-normalized at the end of the previous wrap calculation.
    const WrapResult& previousWrap = aPathWrap.getPreviousWrap();
    aWrapResult.factor = previousWrap.factor;
    // Use Vec3 operators
    aWrapResult.r1 = previousWrap.r1 * previousWrap.factor;
    aWrapResult.r2 = previousWrap.r2 * previousWrap.factor;
    aWrapResult.c1 = previousWrap.c1;
    aWrapResult.sv = previousWrap.sv;

    aFlag = false;
    aWrapResult.wrap_path_length = 0.0;
    aWrapResult.wrap_pts.setSize(0);

    // abort if aPoint1 or aPoint2 is inside the cylinder.
    if (WrapMath::CalcDistanceSquaredPointToLine(aPoint1, p0, dn) < r_squared ||
        WrapMath::CalcDistanceSquaredPointToLine(aPoint2, p0, dn) < r_squared)
    {
        return insideRadius;
    }

    // Find the closest intersection between the muscle line segment and the line
    // segment from one end of the cylinder to the other. This intersection is
    // used in several places further down in the code to check for various
    // wrapping conditions.
    SimTK::Vec3 cylStart, cylEnd;
    cylStart[0] = cylEnd[0] = 0.0;
    cylStart[1] = cylEnd[1] = 0.0;
    cylStart[2] = -0.5 * get_length();
    cylEnd[2] = 0.5 * get_length();

    WrapMath::IntersectLines(aPoint1, aPoint2, cylStart, cylEnd, near12, t12, near00, t00);

    // abort if the cylinder is unconstrained and p1p2 misses the cylinder.
    // Use the return values from the above call to IntersectLines()
    // to perform the check.
    if ( ! constrained)
    {
        if (WrapMath::CalcDistanceSquaredBetweenPoints(near12, near00) < r_squared && t12 > 0.0 && t12 < 1.0)
        {
            return_code = mandatoryWrap;
        }
        else
        {
            return noWrap;
        }
    }

    // find points p11 & p22 on the cylinder axis closest aPoint1 & aPoint2
    WrapMath::GetClosestPointOnLineToPoint(aPoint1, p0, dn, p11, t);
    WrapMath::GetClosestPointOnLineToPoint(aPoint2, p0, dn, p22, t);

    // find preliminary tangent point candidates r1a & r1b
    vv = aPoint1 - p11;

    p11_dist = WrapMath::NormalizeOrZero(vv, vv);

    sin_theta = _radius / p11_dist;

    dist = _radius * sin_theta;

    for (i = 0; i < 3; i++)
        pp[i] = p11[i] + dist * vv[i];

    dist = sqrt(r_squared - dist * dist);

    uu = dn % vv;

    for (i = 0; i < 3; i++)
    {
        r1a[i] = pp[i] + dist * uu[i];
        r1b[i] = pp[i] - dist * uu[i];
    }

    // find preliminary tangent point candidates r2a & r2b
    vv = aPoint2 - p22;

    p22_dist = WrapMath::NormalizeOrZero(vv, vv);

    sin_theta = _radius / p22_dist;

    dist = _radius * sin_theta;

    for (i = 0; i < 3; i++)
        pp[i] = p22[i] + dist * vv[i];

    dist = sqrt(r_squared - dist * dist);

    uu = dn % vv;

    for (i = 0; i < 3; i++)
    {
        r2a[i] = pp[i] + dist * uu[i];
        r2b[i] = pp[i] - dist * uu[i];
    }

    // choose the best preliminary tangent points r1 & r2 from the 4 candidates.
    if (constrained)
    {
        SimTK::Vec3 sum_r;

        if (DSIGN(aPoint1[_wrapAxis]) == _wrapSign || DSIGN(aPoint2[_wrapAxis]) == _wrapSign)
        {
            // If either muscle point is on the constrained side, then check for intersection
            // of the muscle line and the cylinder. If there is an intersection, then
            // you've found a mandatory wrap. If not, then if one point is not on the constrained
            // side and the closest point on the line is not on the constrained side, you've
            // found a potential wrap. Otherwise, there is no wrap.
            // Use the return values from the previous call to IntersectLines()
            // to perform these checks.
            if (WrapMath::CalcDistanceSquaredBetweenPoints(near12, near00) < r_squared && t12 > 0.0 && t12 < 1.0)
            {
                return_code = mandatoryWrap;
            }
            else
            {
                if (DSIGN(aPoint1[_wrapAxis]) != DSIGN(aPoint2[_wrapAxis]) && DSIGN(near12[_wrapAxis]) != _wrapSign)
                {
                    return_code = wrapped;
                }
                else
                {
                    return noWrap;
                }
            }
        }

        r1am = r1a - p11;
        r1bm = r1b - p11;
        r2am = r2a - p22;
        r2bm = r2b - p22;

        alpha = acos((~r1am * r2bm) / (r1am.norm() * r2bm.norm()));
        beta = acos((~r1bm * r2am) / (r1bm.norm() * r2am.norm()));

        // check to see which of the four tangent points should be chosen by seeing which
        // ones are on the 'active' portion of the cylinder. If r1a and r1b are both on or
        // both off the active portion, then use r2a and r2b to decide.
        if (DSIGN(r1a[_wrapAxis]) == _wrapSign && DSIGN(r1b[_wrapAxis]) == _wrapSign)
        {
            if (DSIGN(r2a[_wrapAxis]) == _wrapSign)
            {
                aWrapResult.r1 = r1b;
                aWrapResult.r2 = r2a;
                if (alpha > beta)
                    far_side_wrap = false;
                else
                    far_side_wrap = true;
            }
            else
            {
                aWrapResult.r1 = r1a;
                aWrapResult.r2 = r2b;
                if (alpha > beta)
                    far_side_wrap = true;
                else
                    far_side_wrap = false;
            }
        }
        else if (DSIGN(r1a[_wrapAxis]) == _wrapSign && DSIGN(r1b[_wrapAxis]) != _wrapSign)
        {
            aWrapResult.r1 = r1a;
            aWrapResult.r2 = r2b;
            if (alpha > beta)
                far_side_wrap = true;
            else
                far_side_wrap = false;
        }
        else if (DSIGN(r1a[_wrapAxis]) != _wrapSign && DSIGN(r1b[_wrapAxis]) == _wrapSign)
        {
            aWrapResult.r1 = r1b;
            aWrapResult.r2 = r2a;
            if (alpha > beta)
                far_side_wrap = false;
            else
                far_side_wrap = true;
        }
        else if (DSIGN(r1a[_wrapAxis]) != _wrapSign && DSIGN(r1b[_wrapAxis]) != _wrapSign)
        {
            if (DSIGN(r2a[_wrapAxis]) == _wrapSign)
            {
                aWrapResult.r1 = r1b;
                aWrapResult.r2 = r2a;
                if (alpha > beta)
                    far_side_wrap = false;
                else
                    far_side_wrap = true;
            }
            else if (DSIGN(r2b[_wrapAxis]) == _wrapSign)
            {
                aWrapResult.r1 = r1a;
                aWrapResult.r2 = r2b;
                if (alpha > beta)
                    far_side_wrap = true;
                else
                    far_side_wrap = false;
            }
            else // none of the four tangent points is on the active portion
            {
                if (alpha > beta)
                {
                    aWrapResult.r1 = r1a;
                    aWrapResult.r2 = r2b;
                    far_side_wrap = true;
                }
                else
                {
                    aWrapResult.r1 = r1b;
                    aWrapResult.r2 = r2a;
                    far_side_wrap = true;
                }
            }
        }
        // determine if the resulting tangent points create a short wrap
        // (less than half the cylinder) or a long wrap.
        for (i = 0; i < 3; i++)
        {
            sum_musc[i] = (aWrapResult.r1[i] - aPoint1[i]) + (aWrapResult.r2[i] - aPoint2[i]);
            sum_r[i] = (aWrapResult.r1[i] - p11[i]) + (aWrapResult.r2[i] - p22[i]);
        }

        if ((~sum_r*sum_musc) < 0.0)
            long_wrap = true;
    }
    else
    {
        r1am = r1a - p11;
        r1bm = r1b - p11;
        r2am = r2a - p22;
        r2bm = r2b - p22;

        WrapMath::NormalizeOrZero(r1am, r1am);
        WrapMath::NormalizeOrZero(r1bm, r1bm);
        WrapMath::NormalizeOrZero(r2am, r2am);
        WrapMath::NormalizeOrZero(r2bm, r2bm);

        dot1 = (~r1am*r2am);
        dot2 = (~r1am*r2bm);
        dot3 = (~r1bm*r2am);
        dot4 = (~r1bm*r2bm);

        if (dot1 > dot2 && dot1 > dot3 && dot1 > dot4)
        {
            aWrapResult.r1 = r1a;
            aWrapResult.r2 = r2a;
        }
        else if (dot2 > dot3 && dot2 > dot4)
        {
            aWrapResult.r1 = r1a;
            aWrapResult.r2 = r2b;
        }
        else if (dot3 > dot4)
        {
            aWrapResult.r1 = r1b;
            aWrapResult.r2 = r2a;
        }
        else
        {
            aWrapResult.r1 = r1b;
            aWrapResult.r2 = r2b;
        }
    }

    // bisect angle between r1 & r2 vectors to find the apex edge of the
    // cylinder:
    uu = aWrapResult.r1 - p11;
    vv = aWrapResult.r2 - p22;

    for (i = 0; i < 3; i++)
        vv[i] = uu[i] + vv[i];

    WrapMath::NormalizeOrZero(vv, vv);

    // find the apex point by using a ratio of the lengths of the
    // aPoint1-p11 and aPoint2-p22 segments:
    t = p11_dist / (p11_dist + p22_dist);

    // find point along muscle line according to calculated t value
    for (i = 0; i < 3; i++)
        mpt[i] = aPoint1[i] + t * (aPoint2[i] - aPoint1[i]);

    // find closest point on cylinder axis to mpt
    WrapMath::GetClosestPointOnLineToPoint(mpt, p0, dn, axispt, t);

    // find normal of plane through aPoint1, aPoint2, axispt
    l1 = aPoint1 - axispt;
    l2 = aPoint2 - axispt;

    WrapMath::NormalizeOrZero(l1, l1);
    WrapMath::NormalizeOrZero(l2, l2);

    plane_normal = UnitVec3(l1 % l2);
    WrapMath::NormalizeOrZero(plane_normal, plane_normal);

    // cross plane normal and cylinder axis (each way) to
    // get vectors pointing from axispt towards mpt and
    // away from mpt (you can't tell which is which yet).
    vert1 = dn % plane_normal;
    WrapMath::NormalizeOrZero(vert1, vert1);
    vert2 = plane_normal % dn;
    WrapMath::NormalizeOrZero(vert2, vert2);

    // now find two potential apex points, one along each vector
    for (i = 0; i < 3; i++)
    {
        apex1[i] = axispt[i] + _radius * vert1[i];
        apex2[i] = axispt[i] + _radius * vert2[i];
    }

    // Now use the distance from these points to mpt to
    // pick the right one.
    dist1 = WrapMath::CalcDistanceSquaredBetweenPoints(mpt, apex1);
    dist2 = WrapMath::CalcDistanceSquaredBetweenPoints(mpt, apex2);
    if (far_side_wrap)
    {
        if (dist1 < dist2)
        {
            for (i = 0; i < 3; i++)
                apex[i] = apex2[i];
        }
        else
        {
            for (i = 0; i < 3; i++)
                apex[i] = apex1[i];
        }
    }
    else
    {
        if (dist1 < dist2)
        {
            for (i = 0; i < 3; i++)
                apex[i] = apex1[i];
        }
        else
        {
            for (i = 0; i < 3; i++)
                apex[i] = apex2[i];
        }
    }

    // determine how far to slide the preliminary r1/r2 along their
    // "edge of tangency" with the cylinder by intersecting the aPoint1-ax
    // line with the plane formed by aPoint1, aPoint2, and apex:
    uu = aPoint1 - apex;
    vv = aPoint2 - apex;

    WrapMath::NormalizeOrZero(uu, uu);
    WrapMath::NormalizeOrZero(vv, vv);

    plane_normal = UnitVec3(uu % vv);

    d = - aPoint1[0] * plane_normal[0] - aPoint1[1] * plane_normal[1] - aPoint1[2] * plane_normal[2];

    for (i = 0; i < 3; i++)
    {
        r1a[i] = aWrapResult.r1[i] - 10.0 * dn[i];
        r2a[i] = aWrapResult.r2[i] - 10.0 * dn[i];

        r1b[i] = aWrapResult.r1[i] + 10.0 * dn[i];
        r2b[i] = aWrapResult.r2[i] + 10.0 * dn[i];
    }

    r1_inter = WrapMath::IntersectLineSegPlane(r1a, r1b, plane_normal, d, r1p);
    r2_inter = WrapMath::IntersectLineSegPlane(r2a, r2b, plane_normal, d, r2p);

    if (r1_inter)
    {
        WrapMath::GetClosestPointOnLineToPoint(r1p, p11, p22, r1a, t);

        if (WrapMath::CalcDistanceSquaredBetweenPoints(r1a, p22) < WrapMath::CalcDistanceSquaredBetweenPoints(p11, p22))
            for (i = 0; i < 3; i++)
                aWrapResult.r1[i] = r1p[i];
    }

    if (r2_inter)
    {
        WrapMath::GetClosestPointOnLineToPoint(r2p, p11, p22, r2a, t);

        if (WrapMath::CalcDistanceSquaredBetweenPoints(r2a, p11) < WrapMath::CalcDistanceSquaredBetweenPoints(p22, p11))
            for (i = 0; i < 3; i++)
                aWrapResult.r2[i] = r2p[i];
    }

    // Now that you have r1 and r2, check to see if they are beyond the
    // [display] length of the cylinder. If both are, there should be
    // no wrapping. Since the axis of the cylinder is the Z axis, and
    // the cylinder is centered on Z=0, check the Z components of r1 and r2
    // to see if they are within _length/2.0 of zero.
    double halfL = get_length() / 2.0;
    if ((aWrapResult.r1[2] < -halfL || aWrapResult.r1[2] > halfL) &&
        (aWrapResult.r2[2] < -halfL || aWrapResult.r2[2] > halfL) )
        return noWrap;

    // make the path and calculate the muscle length:
    _make_spiral_path(aPoint1, aPoint2, long_wrap, aWrapResult);

    aFlag = true;

    return return_code;
}

//_____________________________________________________________________________
/**
 * Calculate the wrapping points along a spiral path on the cylinder from aPoint1
 * to aPoint2. This function may slide aPoint1 and aPoint2 axially along the
 * cylinder's surface to achieve tangency to within 1 degree at the two points.
 *
 * @param aPoint1 One end of the spiral path
 * @param aPoint2 The other end of the spiral path
 * @param far_side_wrap Boolean indicating if the wrap is the long way around
 * @param aWrapResult The result of the wrapping (tangent points, etc.)
 */
void WrapCylinder::_make_spiral_path(SimTK::Vec3& aPoint1,
                                                 SimTK::Vec3& aPoint2,
                                                 bool far_side_wrap,
                                                 WrapResult& aWrapResult) const
{
    double x, y, t, axial_dist, theta;
    Vec3 r1a, r2a, uu, vv, ax, axial_vec, wrap_pt;
    double sense = far_side_wrap ? -1.0 : 1.0;
    int i, iterations = 0;
    const double _radius = get_radius();

restart_spiral_wrap:

    aWrapResult.wrap_pts.setSize(0);

    // determine the axial vector

    WrapMath::GetClosestPointOnLineToPoint(aWrapResult.r1, p0, dn, r1a, t);
    WrapMath::GetClosestPointOnLineToPoint(aWrapResult.r2, p0, dn, r2a, t);

    axial_dist = SimTK::dot(r2a - r1a, dn);

    // determine the radial angle
    uu = aWrapResult.r1 - r1a;
    vv = aWrapResult.r2 - r2a;

    theta = std::atan2(
        SimTK::dot(SimTK::cross(uu, vv), dn), // = sin(theta) * radius^2
        SimTK::dot(uu, vv) // = cos(theta) * radius^2
    );

    if (far_side_wrap)
        theta = 2.0 * SimTK_PI - theta;

    // use Pythagoras to calculate the length of the spiral path (imaging
    // a right triangle wrapping around the surface of a cylinder)
    x = _radius * theta;
    y = axial_dist;

    aWrapResult.wrap_path_length = sqrt(x * x + y * y);

    // Determine wrap point by rotating- and translating r1 along axis.
    auto CalcSpiralWrapPoint = [&](double t) -> SimTK::Vec3 // with t = 0 to 1
    {
        SimTK::Mat33 R = SimTK::Rotation(sense * t * theta, dn);
        return R * aWrapResult.r1 + dn * (t * axial_dist);
    };

    // WrapTorus creates a WrapCyl with no connected model, avoid this hack
    if (!_model.empty() && !getModel().getDisplayHints().isVisualizationEnabled() &&
            aWrapResult.singleWrap) {
        // Use one WrapSegment/cord instead of finely sampled list of wrap_pt(s)
        wrap_pt = CalcSpiralWrapPoint(0.);
        aWrapResult.wrap_pts.append(wrap_pt);
        wrap_pt = CalcSpiralWrapPoint(1.);
        aWrapResult.wrap_pts.append(wrap_pt);
    } 
    else {
        // Each muscle segment on the surface of the cylinder should be
        // 0.002 meters long. This assumes the model is in meters, of course.
        int numWrapSegments = (int)(aWrapResult.wrap_path_length / 0.002);
        if (numWrapSegments < 1) numWrapSegments = 1;

        for (i = 0; i < numWrapSegments; i++) {
            double t = (double)i / numWrapSegments;

            wrap_pt = CalcSpiralWrapPoint(t);

            // adjust r1/r2 tangent points if necessary to achieve tangency with
            // the spiral path:
            if (i == 1 && iterations < MAX_ITERATIONS &&
                    !aWrapResult.singleWrap) {
                bool did_adjust_r2 = false;
                bool did_adjust_r1 = _adjust_tangent_point(
                        aPoint1, dn, aWrapResult.r1, wrap_pt);

                SimTK::Vec3 temp_wrap_pt;

                temp_wrap_pt = CalcSpiralWrapPoint(1.0 - t);

                did_adjust_r2 = _adjust_tangent_point(
                        aPoint2, dn, aWrapResult.r2, temp_wrap_pt);

                if (did_adjust_r1 || did_adjust_r2) {
                    iterations++;
                    goto restart_spiral_wrap;
                }
            }

            aWrapResult.wrap_pts.append(wrap_pt);
        }
    }
}

//_____________________________________________________________________________
/**
 * Determine whether the specified tangent point 'r1'
 * needs to be adjusted or not (so that the line from
 * pt1 to r1 is tangent to the cylinder to within 1 degree).
 * If yes, slide it in the appropriate direction and return
 * true, otherwise return false.
 *
 * @param pt1 One point on the line segment being wrapped
 * @param dn Direction vector of cylinder axis
 * @param r1 The tangent point to be adjusted
 * @param w1 A wrapping point (?)
 * @return Whether or not the point was adjusted
 */
bool WrapCylinder::_adjust_tangent_point(SimTK::Vec3& pt1,
                                                      SimTK::Vec3& dn,
                                                      SimTK::Vec3& r1,
                                                      SimTK::Vec3& w1) const
{
    SimTK::Vec3 pr_vec = r1 - pt1;
    SimTK::Vec3 rw_vec = w1 - r1;
    double alpha, omega, t;
    int i;
    bool did_adust = false;

    WrapMath::NormalizeOrZero(pr_vec, pr_vec);
    WrapMath::NormalizeOrZero(rw_vec, rw_vec);

    alpha = acos((~pr_vec*dn));
    omega = acos((~rw_vec*dn));

    if (fabs(alpha - omega) > TANGENCY_THRESHOLD)
    {
        double p1w1_t, p1aw1a_t;
        SimTK::Vec3 save, p1a, w1a, p1w1_int, p1aw1a_int;

        WrapMath::GetClosestPointOnLineToPoint(pt1, r1, dn, p1a, t);
        WrapMath::GetClosestPointOnLineToPoint(w1, r1, dn, w1a, t);

        WrapMath::IntersectLines(pt1, w1, p1a, w1a, p1w1_int, p1w1_t, p1aw1a_int, p1aw1a_t);

        for (i = 0; i < 3; i++)
        {
            save[i] = r1[i];
            r1[i] += 1.5 * (p1aw1a_int[i] - r1[i]);
        }

        did_adust = true;
    }

    return did_adust;
}

// Implement generateDecorations by WrapCylinder to replace the previous out of place implementation 
// in ModelVisualizer
void WrapCylinder::generateDecorations(bool fixed, const ModelDisplayHints& hints, const SimTK::State& state,
    SimTK::Array_<SimTK::DecorativeGeometry>& appendToThis) const 
{

    Super::generateDecorations(fixed, hints, state, appendToThis);
    if (!fixed) return;

    if (hints.get_show_wrap_geometry()) {
        const Appearance& defaultAppearance = get_Appearance();
        if (!defaultAppearance.get_visible()) return;
        const Vec3 color = defaultAppearance.get_color();

        SimTK::Transform ztoy;
        // Make transform that takes z axis to y axis due to different
        // assumptions between DecorativeCylinder aligned with y  and
        // WrapCylinder aligned with z
        ztoy.updR().setRotationFromAngleAboutX(SimTK_PI / 2);

        const auto X_BP = calcWrapGeometryTransformInBaseFrame();
        SimTK::Transform X_BP_ztoy = X_BP*ztoy;
        appendToThis.push_back(
            SimTK::DecorativeCylinder(get_radius(),
                get_length() / 2)
            .setTransform(X_BP_ztoy).setResolution(2.0)
            .setColor(color).setOpacity(defaultAppearance.get_opacity())
            .setScale(1).setRepresentation(defaultAppearance.get_representation())
            .setBodyId(getFrame().getMobilizedBodyIndex()));
    }
}