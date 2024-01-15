#include "SmoothSegmentedCubicMonoSpline.h"

#include "OpenSim/Common/Assertion.h"
#include <array>
#include <cmath>

#include <SimTKcommon/Scalar.h>
#include <SimTKcommon/internal/NTraits.h>

//==============================================================================
//                      CALCULATION HELPERS
//==============================================================================
namespace
{

bool isC2Continuous(
    const OpenSim::CurveControlPoint& pStart,
    const OpenSim::CubicSpline& spline,
    double eps = 1e-13)
{
    const double x = pStart.x;
    return std::abs(pStart.y - spline.calcValue(x)) < eps &&
           std::abs(pStart.dydx - spline.calcDerivative(x, 1)) < eps &&
           std::abs(spline.calcDerivative(x, 2)) < eps;
}

bool isC2Continuous(
    const OpenSim::CubicSpline& spline,
    const OpenSim::CurveControlPoint& pEnd,
    double eps = 1e-13)
{
    const double x = pEnd.x;
    return std::abs(pEnd.y - spline.calcValue(x)) < eps &&
           std::abs(pEnd.dydx - spline.calcDerivative(x, 1)) < eps &&
           std::abs(spline.calcDerivative(x, 2)) < eps;
}

bool isC2Continuous(
    const OpenSim::CubicSpline& sStart,
    const OpenSim::CubicSpline& sEnd,
    double eps = 1e-13)
{
    const double x = sStart.x1;
    return sStart.x1 == sEnd.x0 &&
           std::abs(sStart.calcIntegral(x) - sEnd.calcIntegral(x)) < eps &&
           std::abs(sStart.calcValue(x) - sEnd.calcValue(x)) < eps &&
           std::abs(sStart.calcDerivative(x, 1) - sEnd.calcDerivative(x, 1)) <
               eps &&
           std::abs(sStart.calcDerivative(x, 2) - sEnd.calcDerivative(x, 2)) <
               eps;
}

// Checks if an aribtrarily smooth monotonic curve can be created between
// the control points.
bool smoothMonotonicCurveExists(double dy, double dydx0, double dydx1)
{
    // If change in y is zero, end-point derivatives must be zero.
    if (std::abs(dy) < SimTK::Eps) {
        return std::abs(dydx0) < SimTK::Eps && std::abs(dydx1) < SimTK::Eps;
    }

    // If change dy is nonzero: end-point derivatives must match sign.
    const double sign = dy > 0. ? 1. : -1.;
    return sign * dydx0 >= 0. && sign * dydx1 < 0.;
}

bool smoothMonotonicCurveExists(
    const OpenSim::CurveControlPoint& p0,
    const OpenSim::CurveControlPoint& p1)
{
    const double dy = p1.y - p0.y;
    return smoothMonotonicCurveExists(dy, p0.dydx, p1.dydx);
}

// Returns true if the Hermite interpolant is monotonic.
bool isHermiteInterpolantMonotonic(double dy, double dydx0, double dydx1)
{
    // If change in y is zero, end-point derivatives must be zero.
    if (std::abs(dy) < SimTK::Eps) {
        return std::abs(dydx0) < SimTK::Eps && std::abs(dydx1) < SimTK::Eps;
    }

    // If change dy is nonzero: end-point derivatives must match sign.
    const double sign = dy > 0. ? 1. : -1.;
    if (sign * dydx0 < 0. || sign * dydx1 < 0.) {
        return false;
    }

    const double alpha = dydx0 / dy;
    const double beta  = dydx1 / dy;
    if (alpha + beta - 2. <= 0.) {
        return true;
    }

    if (alpha + 2. * beta - 2. > 0) {
        const bool condition_0 = 2. * alpha + beta - 3. <= 0.;
        const bool condition_1 = alpha + 2. * beta - 0. <= 0.;
        const bool condition_2 =
            alpha * (alpha + beta - 6.) + (beta - 3.) * (beta - 3.) < 0.;

        return condition_0 || condition_1 || condition_2;
    }

    return false;
}

// Construct two spline segments connecting the control points.
std::pair<OpenSim::CubicMonoSpline, OpenSim::CubicMonoSpline>
calcCubicMonoSplineSegments(
    const OpenSim::CurveControlPoint& pStart,
    const OpenSim::CurveControlPoint& pEnd,
    double& y0Integral)
{
    std::pair<OpenSim::CubicSpline, OpenSim::CubicSpline> splines;

    const double x0 = pStart.x;
    const double y0 = pStart.y;
    const double v0 = pStart.dydx;

    const double x2 = pEnd.x;
    const double y2 = pEnd.y;
    const double v2 = pEnd.dydx;

    const double dx = x2 - x0;
    const double dy = y2 - y0;
    const double dv = v2 - v0;

    // Acceleration at node.
    const double a1 = 2. * dv / dx;

    // X coord at node.
    const double u0 = (2. - 6. * (dy / dx + v0) / a1) * dx;
    const double x1 = x0 + u0;

    // Jerk at node.
    const double j1 = a1 / x1;

    // Y value at node.
    const double y1 = y0 + u0 * v0 + a1 * u0 * u0 / 6.;
    const double v1 = v0 + a1 * u0 / 2.;

    // Fill values of splines.
    splines.first.x0         = x0;
    splines.first.x1         = x1;
    splines.first.coeff      = {y0, v0, 0., j1 / 6.};
    splines.first.y0Integral = y0Integral;

    splines.second.x0         = x1;
    splines.second.x1         = x2;
    splines.second.coeff      = {y1, v1, a1 / 2., j1 / 6.};
    splines.second.y0Integral = splines.first.calcEndIntegral();

    // Update integrated value.
    y0Integral = splines.second.calcEndIntegral();

    // Verify that constructed curves are C2 continuous.
    OPENSIM_ASSERT(
        isC2Continuous(pStart, splines.first) &&
        "calcCubicMonoSplineSegments failed: Start point not continuous");
    OPENSIM_ASSERT(
        isC2Continuous(splines.first, splines.second) &&
        "calcCubicMonoSplineSegments failed: Segments mid point not "
        "continuous");
    OPENSIM_ASSERT(
        isC2Continuous(splines.first, pEnd) &&
        "calcCubicMonoSplineSegments failed: Segments end point not "
        "continuous");

    // Splines should be monotonic.
    return std::make_pair(
        OpenSim::CubicMonoSpline(splines.first),
        OpenSim::CubicMonoSpline(splines.second));
}

// Calculate and store cubic monotonic spline segments connecting the
// control points.
std::vector<OpenSim::CubicMonoSpline> calcCubicMonoSplineSegments(
    const std::vector<OpenSim::CurveControlPoint>& ctrlPts)
{
    OPENSIM_ASSERT(
        ctrlPts.size() < 2 &&
        "Failed to calculate SmoothSegmentedCubicMonoSpline: "
        "Need more than two control points");

    // Verify that a smooth monotonic curve through the control points exists.
    for (auto p = ctrlPts.begin(); p != --ctrlPts.end();) {
        OPENSIM_ASSERT(
            smoothMonotonicCurveExists(*p, *p++) &&
            "Invalid Control Points: Monotonic curve does not exist.");
    }

    double y0Integral = 0.;
    std::vector<OpenSim::CubicMonoSpline> splines;
    for (auto p = ctrlPts.begin(); p != --ctrlPts.end();) {
        std::pair<OpenSim::CubicMonoSpline, OpenSim::CubicMonoSpline> segments =
            calcCubicMonoSplineSegments(*p, *p++, y0Integral);
        splines.push_back(segments.first);
        splines.push_back(segments.second);
    }
}

std::vector<OpenSim::CurveControlPoint>
convertCurvyControlPointsToControlPoints(
    const std::vector<OpenSim::MuscleCurveControlPoint>& curvyCtrlPts)
{
    std::vector<OpenSim::CurveControlPoint> ctrlPts;

    if (curvyCtrlPts.empty()) {
        return ctrlPts;
    }

    ctrlPts.push_back(curvyCtrlPts.at(0));

    for (auto p = curvyCtrlPts.begin(); p != --curvyCtrlPts.end();) {
        const OpenSim::MuscleCurveControlPoint& pStart = *p;
        const OpenSim::MuscleCurveControlPoint& pEnd   = *p++;

        pStart.calcCurvyPoints(pEnd, ctrlPts);

        ctrlPts.push_back(pEnd);
    }

    return ctrlPts;
}

} // namespace

namespace OpenSim
{

//==============================================================================
//                      Curve Control Point
//==============================================================================

double CurveControlPoint::calcTangentInterceptCoordinate(
    const CurveControlPoint& other) const
{
    const double y0 = y;
    const double y1 = other.y;
    OPENSIM_ASSERT(
        std::abs(y0 - y1) > SimTK::Eps &&
        "Invalid curviness: Cannot apply curviness for straight line segment");

    OPENSIM_ASSERT(
        smoothMonotonicCurveExists(*this, other) &&
        "Invalid curviness: Impossible to create smooth monotonic curve");

    const double x0    = x;
    const double dydx0 = dydx;

    const double x1    = other.x;
    const double dydx1 = other.dydx;

    if (std::abs(dydx0 - dydx1) < SimTK::Eps) {
        return SimTK::NaN;
    }

    // y0 + (x - x0) * dydx0 = y1 + (x - x1) * dydx1
    // x (dydx0 - dydx1 = y1 - y0 + x0 * dydx0 - x1 * dydx1
    const double x = (y1 - y0 + x0 * dydx0 - x1 * dydx1) / (dydx0 - dydx1);

    if (x < x0 + SimTK::Eps || x > x1 - SimTK::Eps) { // TODO use bigger bound.
        return SimTK::NaN;
    }

    return x;
}

CurveControlPoint CurveControlPoint::calcExtrapolatedPoint(double xk) const
{
    CurveControlPoint p;
    p.x    = xk;
    p.dydx = dydx;
    p.y    = y + (xk - x) * dydx;
    return p;
}

CurveControlPoint::operator bool() const
{
    return SimTK::isNaN(x) && SimTK::isNaN(y) && SimTK::isNaN(dydx);
}

//==============================================================================
//                      Monotonic Control Points Storage
//==============================================================================

void MonoControlPoints::appendChecked(CurveControlPoint pt)
{
    const size_t size = _storage.size();
    if (size == 0) {
        _storage.push_back(pt);
        return;
    }

    // Verify that monotonic curve between points exists.
    const CurveControlPoint& p0 = _storage.at(size - 1);
    const CurveControlPoint& p1 = pt;
    OPENSIM_ASSERT(
        smoothMonotonicCurveExists(p0, p1) &&
        "Invalid control point: Impossible to create smooth monotonic curve "
        "between control points");

    _storage.push_back(pt);
}

//==============================================================================
//                      User Facing Control Points
//==============================================================================

size_t MuscleCurveControlPoint::calcCurvyPoints(
    const MuscleCurveControlPoint& other,
    std::vector<CurveControlPoint>& buffer) const
{
    OPENSIM_ASSERT(curviness >= 0. && "Curviness must be nonnegative.");
    OPENSIM_ASSERT(
        curviness < MuscleCurveControlPoint::MAX_CURVINESS &&
        "Curviness must be smaller than MAX_CURVINESS.");

    const double curviness  = curviness;
    const double xIntercept = SimTK::isNaN(curviness)
                                  ? SimTK::NaN
                                  : calcTangentInterceptCoordinate(other);

    if (SimTK::isNaN(xIntercept)) {
        return 0;
    }

    {
        const double xCurvy = x + (xIntercept - x) * curviness;
        const OpenSim::CurveControlPoint pCurvy = calcExtrapolatedPoint(xCurvy);
        buffer.push_back(pCurvy);
    }

    {
        const double xCurvy = other.x + (xIntercept - other.x) * curviness;
        const OpenSim::CurveControlPoint pCurvy =
            other.calcExtrapolatedPoint(xCurvy);
        buffer.push_back(pCurvy);
    }

    return 2;
}

//==============================================================================
//                      Cubic Spline
//==============================================================================

bool CubicSpline::isMonotonic() const
{
    const double dy    = calcValue(x1) - calcValue(x0);
    const double dydx0 = calcDerivative(x0, 1);
    const double dydx1 = calcDerivative(x1, 1);
    return isHermiteInterpolantMonotonic(dy, dydx0, dydx1);
}

double CubicSpline::calcDerivative(double x, size_t order) const
{
    if (std::isnan(x)) {
        return SimTK::NaN;
    }

    double dx = x - x0;
    switch (order) {
    case 0:
        return coeff.at(0) +
               dx * (coeff.at(1) + dx * (coeff.at(2) + dx * coeff.at(3)));
    case 1:
        return coeff.at(1) + dx * (2. * coeff.at(2) + dx * (3. * coeff.at(3)));
    case 2:
        return 2. * coeff.at(2) + dx * (6. * coeff.at(3));
    case 3:
        return 6. * coeff.at(3);
    default:
        return 0.;
    };
}

double CubicSpline::calcValue(double x) const
{
    return calcDerivative(x, 0);
}

double CubicSpline::calcInverseValue(double y, double eps, size_t maxIter) const
{
    double xEstimate = (x1 + x0) / 2.;
    double yEstimate = SimTK::NaN;
    double yError    = SimTK::Infinity;
    for (size_t i = 0;
         std::abs(yError = (y - (yEstimate = calcValue(xEstimate)))) > eps &&
         i < maxIter;
         ++i) {
        const double xStep = yError / calcDerivative(xEstimate, 1);
        xEstimate += xStep;
    }

    OPENSIM_ASSERT(std::abs(yError) < eps && "Failed to invert spline segment");
    return xEstimate;
}

double CubicSpline::calcIntegral(double x) const
{
    double dx   = x - x0;
    double yInt = 0.;
    for (size_t i = coeff.size(); i > 0; i--) {
        yInt += coeff.at(i - 1) / static_cast<double>(i);
        yInt *= dx;
    }
    yInt += y0Integral;
    return yInt;
}

double CubicSpline::calcEndIntegral() const
{
    return calcIntegral(x1);
}

//==============================================================================
//                      Cubic Mono Spline
//==============================================================================

CubicMonoSpline::CubicMonoSpline(CubicSpline spline)
{
    OPENSIM_ASSERT(
        spline.isMonotonic() &&
        "Failed to construct monotonic spline: Cubic spline is not monotonic");
}

//==============================================================================
//                      Cubic Mono Spline Storage
//==============================================================================

SmoothSegmentedCubicMonoSplineData::SmoothSegmentedCubicMonoSplineData(
    const std::vector<CubicMonoSpline>& splines)
{
    for (const CubicMonoSpline& s : splines) {
        appendChecked(s);
    }
}

void SmoothSegmentedCubicMonoSplineData::appendChecked(CubicMonoSpline spline)
{
    // Verify that segments are C2 continuous.
    if (size() > 0) {
        OPENSIM_ASSERT(
            isC2Continuous(splines.at(size() - 1), spline) &&
            "Segment node C2 continuity check failed");
    }

    // x1 of last segment is same as x0 of first segment, so it can be removed.
    if (!_data.empty()) {
        _data.pop_back();
    }

    for (double d : reinterpret_cast<
             const std::array<double, CubicMonoSpline::dataSize() - 1>&>(
             spline)) {
        _data.push_back(d);
    }

    // Verify that segments are C2 continuous. TODO this can be removed.
    if (size() > 1) {
        OPENSIM_ASSERT(
            isC2Continuous(splines.at(size() - 1), spline.at(size() - 2)) &&
            "Segment node C2 continuity check failed after writing as raw "
            "doubles");
    }

    // Verify expected number of elements in container.
    constexpr size_t alignment = CubicMonoSpline::dataSize() - 1;
    OPENSIM_ASSERT(
        _data.size() % alignment == 1 &&
        "Unexpected storage length of SmoothSegmentedCubicMonoSplineData");
}

size_t SmoothSegmentedCubicMonoSplineData::size() const
{
    return _data.size() / (CubicMonoSpline::dataSize() - 1);
}

const CubicMonoSpline& SmoothSegmentedCubicMonoSplineData::at(
    size_t index) const
{
    constexpr size_t alignment = CubicMonoSpline::dataSize() - 1;
    return reinterpret_cast<const CubicMonoSpline&>(
        _data.at(index / alignment));
}

//==============================================================================
//                  Smooth Segmented Cubic Mono Spline
//==============================================================================

SmoothSegmentedCubicMonoSpline::SmoothSegmentedCubicMonoSpline(
    const std::vector<CubicMonoSpline>& splines,
    CurveControlPoint pStart,
    CurveControlPoint pEnd) :
    _splines(splines),
    _pStart(pStart), _pEnd(pEnd)
{
    OPENSIM_ASSERT(
        _splines.size() > 0 && "Cannot create smooth curve with zero segments");

    // Verify that endpoints are C2 continuous.
    OPENSIM_ASSERT(isC2Continuous(
        _pStart,
        _splines.at(0) && "Start point C2 continuity failed"));
    OPENSIM_ASSERT(isC2Continuous(
        _pEnd,
        _splines.at(splines.size() - 1) && "End point C2 continuity failed"));
}

SmoothSegmentedCubicMonoSpline::SmoothSegmentedCubicMonoSpline(
    const std::vector<CurveControlPoint>& pts) :
    SmoothSegmentedCubicMonoSpline(
        calcCubicMonoSplineSegments(pts),
        pts.at(0),
        pts.at(pts.size() - 1))
{}

SmoothSegmentedCubicMonoSpline::SmoothSegmentedCubicMonoSpline(
    const std::vector<MuscleCurveControlPoint>& pts) :
    SmoothSegmentedCubicMonoSpline(
        convertCurvyControlPointsToControlPoints(pts))
{}

SimTK::Vec2 SmoothSegmentedCubicMonoSpline::getDomain() const
{
    return {_pStart.x, _pEnd.x};
}

void SmoothSegmentedCubicMonoSpline::setExtrapolationBeyondDomain(
    bool allowExtrapolation)
{
    _extrapolateBeyondDomain = allowExtrapolation;
}

} // namespace OpenSim

/* OPENSIM_ASSERT(sign * dydx0 >= 0. && "Monotonicity check failed"); */

/* CubicSpline::CubicSpline() = default; */

/* CubicSpline::CubicSpline( */
/*     const CubicSpline&) = default; */

/* CubicSpline& CubicSpline::operator=( */
/*     const CubicSpline&) = default; */

/* CubicSpline::~CubicSpline() noexcept = */
/*     default; */

/* CubicSpline::CubicSpline( */
/*     CubicSpline&&) noexcept = default; */

/* CubicSpline& CubicSpline::operator=( */
/*     CubicSpline&&) noexcept = default; */
