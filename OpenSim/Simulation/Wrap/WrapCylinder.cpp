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

// Returns the angle on [0, 2pi] range.
double wrap_tau(double angle) {
    return ( angle < 0. )? angle + 2. * M_PI: angle;
}

// Wrapping path of single point over circle and back to the point.
struct SinglePointCircleWrap final {

    SinglePointCircleWrap(
        const SimTK::Vec2& point,
        double radius)
    {
        // The wrapping path is symmetric about the line through the point.
        // The angle between that line and the first point on the circle can be
        // found from the right-angled triangle with sides:
        // - long side = |point|
        // - near side = radius

        // Long triangle side (squared).
        double long_side_sq = SimTK::dot(point, point);
        // Near triangle side (squared).
        double radius_sq = std::pow(radius, 2);
        // Intermezzo: Point must not lie within cylinder.
        _inside_circle = long_side_sq < radius_sq;
        if (_inside_circle) {
            return;
        }
        // Far triangle side.
        double far_side = std::sqrt(long_side_sq - radius_sq);

        double delta_angle = std::atan2(far_side, radius);
        double angle_of_point = std::atan2(point[1], point[0]);
        _start_angle = wrap_tau(angle_of_point + delta_angle);
        _final_angle = wrap_tau(angle_of_point - delta_angle);
        _surface_to_point_distance = far_side;
    }

    // Start angle of wrapping path on circle: within [0, 2pi].
    double _start_angle = NAN;
    // Final angle of wrapping path on circle: within [0, 2pi].
    double _final_angle = NAN;
    // Distance from point to wrapping-point on circle.
    double _surface_to_point_distance = NAN;

    // Error flag: When point lies inside circle path does not exist.
    bool _inside_circle = false;
};

std::ostream& operator<<(std::ostream& os, const SinglePointCircleWrap& x) {
    os << "SinglePointCircleWrap { ";
    if (x._inside_circle) {
        os << "_inside_circle: " << x._inside_circle << "}";
        return os;
    }
    os << "_start_angle: " << x._start_angle;
    os << ", ";
    os << "_final_angle: " << x._final_angle;
    os << ", ";
    os << "_surface_to_point_distance: " << x._surface_to_point_distance;
    os << " }";
    return os;
}

// 2D wrapping path of point-a to point-b over a circle.
struct CircleWrap final {

    CircleWrap(
        const SimTK::Vec2& point_a,
        const SimTK::Vec2& point_b,
        double radius,
        int wrapAxis,
        int wrapSign
    ) {
        SimTK_ERRCHK_ALWAYS(wrapAxis != 2,
            "CircleWrap",
            "invalid quadrant: wrapAxis equals z-axis.");

        _no_wrap |= pre_wrap_check_failed(point_a, point_b, radius, wrapAxis, wrapSign);
        if (_no_wrap) {
            return;
        }

        const SinglePointCircleWrap angles_a(point_a, radius);
        const SinglePointCircleWrap angles_b(point_b, radius);

        std::cout << "angles_a = " << angles_a << std::endl;
        std::cout << "angles_b = " << angles_b << std::endl;

        _inside_circle = angles_a._inside_circle || angles_b._inside_circle;
        if (_inside_circle) {
            return;
        }

        // There are two possible wrapping paths from a to b: positive and negative rotating paths.
        const double delta_angle_pos = wrap_tau(angles_b._final_angle - angles_a._start_angle);
        const double delta_angle_neg = wrap_tau(angles_a._final_angle - angles_b._start_angle);

        // Check validity of each path.
        const bool wrap_x_axis = wrapAxis == 0;
        const bool pos_invalid = potential_wrap_is_invalid(angles_a._start_angle, angles_b._final_angle, delta_angle_pos, wrap_x_axis, wrapSign);
        const bool neg_invalid = potential_wrap_is_invalid(angles_b._start_angle, angles_a._final_angle, delta_angle_neg, wrap_x_axis, wrapSign);

        // Both invalid? -> no wrapping path.
        _no_wrap |= pos_invalid && neg_invalid;
        if (_no_wrap) {
            return;
        }

        // Determine wrapping direction: positive, or negative.
        bool pos_wrap = !pos_invalid;
        if (pos_wrap && !neg_invalid) {
            // Both wrapping directions are valid: take minimal wrap length.
            pos_wrap = std::abs(delta_angle_neg) > std::abs(delta_angle_pos);
        }

        // Set the path angles.
        if (pos_wrap) {
            _angle = angles_a._start_angle;
            _delta_angle = delta_angle_pos;
        } else {
            _angle = angles_a._final_angle;
            _delta_angle = -delta_angle_neg;
        }

        // Set path segment lengths.
        _length_point_a_to_circle = angles_a._surface_to_point_distance;
        _length_on_circle = std::abs(_delta_angle) * radius;
        _length_total = _length_point_a_to_circle + _length_on_circle +
            angles_b._surface_to_point_distance;
	}

    // Returns if wrapping should be attempted.
    //
    // Checks 3 conditions:
    // - 1. If unconstrained and line a-b does not intersect the cylinder.
    // - 2. If both points in active quadrant.
    // - 3. If one point in non-active quadrant, and closest point in non-active quadrant.
    static bool pre_wrap_check_failed(
        const SimTK::Vec2& point_a,
        const SimTK::Vec2& point_b,
        double radius,
        int wrapAxis,
        int wrapSign)
    {
        bool miss = false;
        SimTK::Vec2 vec_ab = point_b - point_a;
        double vTv = SimTK::dot(vec_ab, vec_ab);
        double aTa = SimTK::dot(point_a, point_a);
        double vTa = SimTK::dot(vec_ab, point_a);

        // Shortest distance (squared) to line b-a:
        double distance_sq = aTa - std::pow(vTa, 2) / vTv;

        // In unconstrained case: Miss if line a-b does not intersect circle.
        if ( wrapSign == 0 ) {
            // If shortest distance larger than radius.
            miss |= distance_sq > std::pow(radius, 2);

            // Otherwise if a,b on same side of shortest distance vector -> miss.
            double dTb = SimTK::dot(vec_ab, point_b);
            bool ab_same_side =  ( dTb * vTa ) > 0.;
            miss |=  ab_same_side;

            return miss;
        }

        // Both points in active quadrant -> miss.
        if (DSIGN(point_a[wrapAxis]) == wrapSign && DSIGN(point_b[wrapAxis]) == wrapSign) {
            miss |= distance_sq > std::pow(radius, 2);
        }

        // If one point in non-active, check if closest is also in non-active quadrant.
        if ( DSIGN(point_a[wrapAxis]) != DSIGN(point_b[wrapAxis]) ) {
            double bTb = std::pow(point_b[0],2 ) + std::pow(point_b[1], 2);
            double near = (aTa < bTb)? point_a[wrapAxis]: point_b[wrapAxis];
            miss |= DSIGN(near) != wrapSign;
        }
        return miss;
    }

    // Returns if the wrapping from start- to final-angle is invalid based on the quadrant.
    //
    // Checks 2 conditions:
    // 1. Both angles in non-active quadrant.
    // 2. Both angles in active-quadrant, but wrapping the entire non-active quadrant.
    //
    // Angles assumed to be within [0, 2pi]
    // wrapXAxis: true if X-quadrant, false if y-quadrant.
    static bool potential_wrap_is_invalid(
        double start_angle,
        double final_angle,
        double delta_angle,
        bool wrapXAxis,
        int wrapSign)
    {
        if (wrapSign == 0) {
            return false;
        }
        const bool negativeQuadrant = wrapSign < 0;

        // Check if points on cylinder lie in the active quadrant.
        auto not_in_active_quadrant = [&](double angle) -> bool {
            bool not_in_active;

            // If +X quadrant: angle must lie in [PI/2 , -PI/2]
            if (wrapXAxis) {
                not_in_active = ( angle + 0.5 * M_PI ) < M_PI;
            // If +Y quadrant: angle must lie in [0 , PI]
            } else {
                not_in_active = angle > M_PI;
            }
            // Invert the quadrant if needed.
            if (negativeQuadrant) {
                return !not_in_active;
            }
            return not_in_active;
        };

        const bool pos_not_in_active = not_in_active_quadrant(start_angle);
        const bool neg_not_in_active = not_in_active_quadrant(final_angle);

        // If both angles in non-active quadrant.
        bool miss = pos_not_in_active && neg_not_in_active;

        // If positive and negative in active quadrant, but wrapping the non-active.
        miss |= !pos_not_in_active && !neg_not_in_active && ( delta_angle > M_PI );

        return miss;
    }

    // Starting angle of path on circle.
	double _angle = NAN;
    // Angle difference between start and final path points on circle.
	double _delta_angle = NAN;
    // Path lengths:
    double _length_point_a_to_circle = NAN;
    double _length_on_circle = NAN;
    double _length_total = NAN;

    // Error flags: no wrapping if either is set.
    bool _no_wrap = false;
    bool _inside_circle = false;
};

std::ostream& operator<<(std::ostream& os, const CircleWrap& x) {
    os << "CircleWrap { ";
    if (x._inside_circle) {
        os << "inside_circle };";
        return os;
    }
    if (x._no_wrap) {
        os << "no_wrap };";
        return os;
    }
    os << "_angle: " << x._angle;
    os << ", _delta_angle: " << x._delta_angle;
    os << ", _length_point_a_to_circle: " << x._length_point_a_to_circle;
    os << ", _length_on_circle: " << x._length_on_circle;
    os << ", _length_total: " << x._length_total;
    os << " }";
    return os;
}

struct ScrewMap final {
	ScrewMap(
        const SimTK::Vec3& point_a,
        const SimTK::Vec3& point_b,
		double radius,
		double cylinder_length,
        int wrapAxis,
        int wrapSign
    ): _radius(radius) {

        // XY-coords can be computed by wrapping a circle.
        CircleWrap circle_wrap(
            SimTK::Vec2(point_a[0], point_a[1]),
            SimTK::Vec2(point_b[0], point_b[1]),
            radius,
            wrapAxis,
            wrapSign);
        std::cout << circle_wrap << std::endl;

        _no_wrap = circle_wrap._no_wrap;
        _inside_radius = circle_wrap._inside_circle;
        if (_no_wrap || _inside_radius) {
            return;
        }
        _delta_angle = circle_wrap._delta_angle;
        _angle = circle_wrap._angle;

        // Path gradient in z-direction is constant:
        const double delta_z_per_delta_xy = ( point_b[2] - point_a[2] ) / circle_wrap._length_total;

        _height = point_a[2] + delta_z_per_delta_xy * circle_wrap._length_point_a_to_circle;
        _delta_height = delta_z_per_delta_xy * circle_wrap._length_on_circle;

        _no_wrap = std::abs(_height) > cylinder_length / 2.;
        _no_wrap |= std::abs(_height + _delta_height) > cylinder_length / 2.;

        _length = std::sqrt(
            std::pow(radius * circle_wrap._delta_angle, 2) +
            std::pow(_delta_height, 2)
        );
	}

	SimTK::Vec3 compute_point_on_surface(double t) const {
		const double angle = t * _delta_angle + _angle;
		const double height = t * _delta_height + _height;

		return SimTK::Vec3(std::cos(angle) * _radius, std::sin(angle) * _radius,
						   height);
	}

    int wrap_action() const {
        // TODO not using: wrapped? = successful wrap, but may not be 'best' path
        if (_no_wrap) {
            return WrapObject::WrapAction::noWrap;
        }
        if (_inside_radius) {
            return WrapObject::WrapAction::insideRadius;
        }
        return WrapObject::WrapAction::mandatoryWrap;
    }

    bool wrap_ok() const {
        return !_no_wrap && !_inside_radius;
    }

    double wrap_length() const {
        return _length;
    }

    // Angle of starting path point on cylinder surface.
	double _delta_angle = 0. / 0.;
    // Angle difference between start and final path points on cylinder surface.
	double _angle = 0. / 0.;
    // Height of starting path point on cylinder surface.
	double _height = 0. / 0.;
    // Height difference between start and final path points on cylinder surface.
	double _delta_height = 0. / 0.;
    // Total path length.
    double _length = 0. / 0.;
    // Cylinder radius.
    double _radius = 0. / 0.;

    // Error flags.
    bool _no_wrap = false;
    bool _inside_radius = false;
};

std::ostream& operator<<(std::ostream& os, const ScrewMap& screw) {
    os << "ScrewMap { ";
    if (screw._no_wrap) {
        os << "no_wrap }";
        return os;
    }
    if (screw._inside_radius) {
        os << "inside_radius }";
        return os;
    }
    os << "delta_angle: " << screw._delta_angle;
    os << ", ";
    os << "angle: " << screw._angle;
    os << ", ";
    os << "delta_height: " << screw._delta_height;
    os << ", ";
    os << "height: " << screw._height;
    os << ", ";
    os << "length: " << screw._length;
    os << ", ";
    os << "radius: " << screw._radius;
    os << " }";
    return os;
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
    ScrewMap screw = ScrewMap(aPoint1, aPoint2, get_radius(), get_length(), _wrapAxis, _wrapSign);

    aWrapResult.wrap_path_length = 0.0;
    aWrapResult.wrap_pts.setSize(0);
    aFlag = false;

    if (screw.wrap_ok()) {
        aFlag = true;
        aWrapResult.r1 = screw.compute_point_on_surface(0.);
        aWrapResult.r2 = screw.compute_point_on_surface(1.);

        aWrapResult.wrap_path_length = screw.wrap_length();

        aWrapResult.wrap_pts.setSize(0);

        int numWrapSegments = (int)(aWrapResult.wrap_path_length / 0.002);
        if (numWrapSegments < 1) {numWrapSegments = 1;}
        for (int i = 0; i < numWrapSegments; i++) {
            double t = (double)i / numWrapSegments;
            aWrapResult.wrap_pts.append(screw.compute_point_on_surface(t));
        }
    }

    return screw.return_code();

    std::cout << "screw = " << screw << std::endl;
    std::cout << "\n    ";
    std::cout << "let radius = " << get_radius() <<";";
    std::cout << "// ";
    std::cout << "wrapSign = " << _wrapSign;
    std::cout << "// ";
    std::cout << "wrapAxis = " << _wrapAxis;
    std::cout << "\n    ";
    std::cout << "let p0 = Vector3::new(" << aPoint1 << ");";
    std::cout << "\n    ";
    std::cout << "let p1 = Vector3::new(" << aPoint2 << ");";
    std::cout << "\n    ";

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
        std::cout << "inside radous wrap\n";

        SimTK_ERRCHK_ALWAYS(screw.wrap_action() == insideRadius,
            "mandatorywrap",
            "mandatorywrap not equal"
        );
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
            SimTK_ERRCHK_ALWAYS(screw.wrap_action() == mandatoryWrap,
                "mandatorywrap",
                "mandatorywrap not equal"
            );
            return_code = mandatoryWrap;
        }
        else
        {
            static int miss_count = 0;
            ++miss_count;
            std::cout << "missed! " << miss_count << "\n";
            // std::cout << "aPoint1 = " << aPoint1 << "\n";
            // std::cout << "aPoint2 = " << aPoint2 << "\n";
            SimTK_ERRCHK_ALWAYS(screw.wrap_action() == noWrap,
                "nowrap",
                "nowrap not equal"
            );
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
                    std::cout << "nowrap2\n";
                    std::cout << "wrapAxis = " << _wrapAxis;
                    std::cout << "wrapSign = " << _wrapSign;
                    SimTK_ERRCHK_ALWAYS(screw.wrap_action() == noWrap,
                        "nowrap",
                        "nowrap not equal"
                    );
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
        {
        // std::cout << "dropped off the wrapp\n";
                    SimTK_ERRCHK_ALWAYS(screw.wrap_action() == noWrap,
                        "nowrap",
                        "nowrap not equal"
                    );
        return noWrap;
        }

    // make the path and calculate the muscle length:
    _make_spiral_path(aPoint1, aPoint2, long_wrap, aWrapResult);

    aFlag = true;

    std::cout << "completed wrap: " << return_code << "\n" ;
    SimTK_ERRCHK_ALWAYS(screw.wrap_action() == return_code,
        "return_code",
        "return_code not equal"
    );
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
void WrapCylinder::_make_spiral_path(
    SimTK::Vec3& aPoint1,
    SimTK::Vec3& aPoint2,
    bool far_side_wrap,
    WrapResult& aWrapResult) const
{
    // std::cout << "start spiral path\n";
    const double _radius = get_radius();

    double x, y, t, axial_dist, theta;
    Vec3 r1a, r2a, uu, vv, ax, axial_vec, wrap_pt;
    double sense = far_side_wrap ? -1.0 : 1.0;
    int i, iterations = 0;

restart_spiral_wrap:

    aWrapResult.wrap_pts.setSize(0);

    // determine the axial vector

    WrapMath::GetClosestPointOnLineToPoint(aWrapResult.r1, p0, dn, r1a, t);
    WrapMath::GetClosestPointOnLineToPoint(aWrapResult.r2, p0, dn, r2a, t);

    axial_vec = r2a - r1a;

    axial_dist = axial_vec.norm();

    // determine the radial angle
    uu = aWrapResult.r1 - r1a;
    vv = aWrapResult.r2 - r2a;

    for (i = 0; i < 3; i++)
    {
        uu[i] /= _radius;
        vv[i] /= _radius;
    }

    theta = acos((~uu * vv) / (uu.norm() * vv.norm()));

    if (far_side_wrap)
        theta = 2.0 * SimTK_PI - theta;

    // use Pythagoras to calculate the length of the spiral path (imaging
    // a right triangle wrapping around the surface of a cylinder)
    x = _radius * theta;
    y = axial_dist;

    aWrapResult.wrap_path_length = sqrt(x * x + y * y);

    // build path segments
    ax = uu % vv;
    WrapMath::NormalizeOrZero(ax, ax);

    vv = ax % uu;

    SimTK::Rotation m;
    m.set(0, 0, ax[0]); m.set(0, 1, ax[1]); m.set(0, 2, ax[2]);
    m.set(1, 0, uu[0]); m.set(1, 1, uu[1]); m.set(1, 2, uu[2]);
    m.set(2, 0, vv[0]); m.set(2, 1, vv[1]); m.set(2, 2, vv[2]);
    // WrapTorus creates a WrapCyl with no connected model, avoid this hack
    if (!_model.empty() && !getModel().getDisplayHints().isVisualizationEnabled() &&
            aWrapResult.singleWrap) {
        SimTK_ERRCHK_ALWAYS(false,
            "wraptoruswhatever",
            "%s: no idea");
        
        // Use one WrapSegment/cord instead of finely sampled list of wrap_pt(s)
        _calc_spiral_wrap_point(
                r1a, axial_vec, m, ax, sense, 0, theta, wrap_pt);
        aWrapResult.wrap_pts.append(wrap_pt);
        _calc_spiral_wrap_point(
                r1a, axial_vec, m, ax, sense, 1.0, theta, wrap_pt);
        aWrapResult.wrap_pts.append(wrap_pt);
    } 
    else {
        // Each muscle segment on the surface of the cylinder should be
        // 0.002 meters long. This assumes the model is in meters, of course.
        int numWrapSegments = (int)(aWrapResult.wrap_path_length / 0.002);
        if (numWrapSegments < 1) numWrapSegments = 1;

        for (i = 0; i < numWrapSegments; i++) {
            double t = (double)i / numWrapSegments;
            // std::cout << " t = " << t << ", ";

            _calc_spiral_wrap_point(
                    r1a, axial_vec, m, ax, sense, t, theta, wrap_pt);

            // adjust r1/r2 tangent points if necessary to achieve tangency with
            // the spiral path:
            if (i == 1 && iterations < MAX_ITERATIONS &&
                    !aWrapResult.singleWrap) {
                bool did_adjust_r2 = false;
                bool did_adjust_r1 = _adjust_tangent_point(
                        aPoint1, dn, aWrapResult.r1, wrap_pt);

                SimTK::Vec3 temp_wrap_pt;

                _calc_spiral_wrap_point(r1a, axial_vec, m, ax, sense, 1.0 - t,
                        theta, temp_wrap_pt);

                did_adjust_r2 = _adjust_tangent_point(
                        aPoint2, dn, aWrapResult.r2, temp_wrap_pt);

                if (did_adjust_r1 || did_adjust_r2) {
                    // std::cout << "ever do any iterations" << iterations << "\n";
                    iterations++;
                    goto restart_spiral_wrap;
                }
            }

            aWrapResult.wrap_pts.append(wrap_pt);
        }
    }
    static double maxLengthErr =0.;

    bool mapSense = _wrapSign > 0;

    ScrewMap map0 = ScrewMap(aPoint1, aPoint2, get_radius(), get_length(), _wrapAxis, _wrapSign);

    double lengthErr = std::abs(aWrapResult.wrap_path_length - map0._length);
    std::cout << "lengtErr = " << lengthErr << "\n";

    std::cout << "// Wrapping";
    std::cout << "\n    ";
    std::cout << "// screw = " << map0;
    std::cout << "\n    ";
    std::cout << "let radius = " << _radius <<";";
    std::cout << "\n    ";
    std::cout << "let sense = " << sense << ";";
    std::cout << "// ";
    std::cout << "wrapSign = " << _wrapSign;
    std::cout << ", ";
    std::cout << "far_side_wrap = " << far_side_wrap;
    std::cout << ", ";
    std::cout << "iterations = " << iterations;
    std::cout << "\n    ";
    std::cout << "// map._angle = " << map0._angle;
    std::cout << "\n    ";
    std::cout << "// map.delta_angle = " << map0._delta_angle;
    std::cout << "\n    ";
    std::cout << "// map._height = " << map0._height;
    std::cout << "\n    ";
    std::cout << "// map._delta_height = " << map0._delta_height;
    std::cout << "\n    ";
    std::cout << "let p0 = Vector3::new(" << aPoint1 << ");";
    std::cout << "\n    ";
    std::cout << "let p1 = Vector3::new(" << aPoint2 << ");";
    std::cout << "\n    ";

    Vec3 wrap_pt0;
    double tvals[5] = {0., 0.25, 0.5, 0.75, 1.};
    for (int i =0; i < 5; ++i) {
        _calc_spiral_wrap_point(r1a, axial_vec, m, ax, sense, tvals[i], theta, wrap_pt0);
        std::cout << "let q" << i << " = Vector3::new(";
        std::string delim = "";
        for (int j =0; j < 3; ++j) {
            std::cout << delim << wrap_pt0[j];
            delim = ", ";
        }
        std::cout << ");\n    ";
    }

    for (int i =0; i < 5; ++i) {
        SimTK::Vec3 check_pt0 = map0.compute_point_on_surface(tvals[i]);
        std::cout << "let c" << i << " = Vector3::new(";
        std::string delim = "";
        for (int j =0; j < 3; ++j) {
            std::cout << delim << check_pt0[j];
            delim = ", ";
        }
        std::cout << ");\n    ";
    }

    std::cout << "//len = " << aWrapResult.wrap_path_length;
    std::cout << ", ";
    std::cout << "//check_len = " << map0._length;
    std::cout << ", ";
    std::cout << "//lengthErr = " << lengthErr;
    std::cout << "\n    ";
    std::cout << "//lengthErr = " << lengthErr;
    std::cout << "\n    ";

    maxLengthErr = std::max(maxLengthErr, lengthErr);
    std::cout << ", ";
    std::cout << "maxLengthErr = " << maxLengthErr;
    std::cout << std::endl;
        SimTK_ERRCHK_ALWAYS(maxLengthErr < 0.05,
            "wrap_points",
            "lengthError to huge");
}

//_____________________________________________________________________________
/**
 * Calculate a new point along a spiral wrapping path.
 *
 * @param r1a An existing point on the spiral path
 * @param axial_vec Vector from r1a parallel to cylinder axis
 * @param m A transform matrix used for the radial component
 * @param axis Axis of the cylinder
 * @param sense The direction of the spiral
 * @param t Parameterized amount of angle for this point
 * @param theta The total angle of the spiral on the cylinder
 * @param wrap_pt The new point on the spiral path
 */
void WrapCylinder::_calc_spiral_wrap_point(const SimTK::Vec3& r1a,
                                             const SimTK::Vec3& axial_vec,
                                             const SimTK::Rotation& m,
                                             const SimTK::Vec3& axis,
                                             double sense,
                                             double t,
                                             double theta,
                                             SimTK::Vec3& wrap_pt) const
{
    SimTK::Rotation R;
    double angle = sense * t * theta;
    R.setRotationFromAngleAboutNonUnitVector(angle, axis);

    SimTK::Mat33 n = m * ~R;

    for (int i = 0; i < 3; i++)
    {
        double radial_component = get_radius() * n[1][i];
        double axial_component = t * axial_vec[i];

        wrap_pt[i] = r1a[i] + radial_component + axial_component;
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