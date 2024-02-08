#include "SmoothSegmentedCubicMonoSpline.h"

#include "OpenSim/Common/Assertion.h"
#include <array>
#include <cmath>
#include <stdexcept>

#include <SimTKcommon/Scalar.h>
#include <SimTKcommon/internal/NTraits.h>

namespace
{

void opensim_assert(bool cond, std::string msg)
{
    if (!cond) {
        throw std::runtime_error(msg);
    }
}

//==============================================================================
//                      CALCULATION HELPERS: C1-CONTINUITY CHECKS
//==============================================================================
bool isNumEq(double lhs, double rhs, double eps = 1e-13)
{
    return std::abs(lhs - rhs) < eps;
}

bool isC1Continuous(
    const OpenSim::CurveKnot& left,
    const OpenSim::CurveKnot& right)
{
    return isNumEq(left.x, right.x) && isNumEq(left.y, right.y) &&
           isNumEq(left.dydx, right.dydx);
}

bool isC1Continuous(
    const OpenSim::CurveKnot& left,
    const OpenSim::CubicSpline& right)
{
    return isC1Continuous(left, right.calcKnot(right.x0));
}

bool isC1Continuous(
    const OpenSim::CubicSpline& left,
    const OpenSim::CurveKnot& right)
{
    return isC1Continuous(left.calcKnot(left.x1), right);
}

bool isC1Continuous(
    const OpenSim::CubicSpline& left,
    const OpenSim::CubicSpline& right)
{
    return isC1Continuous(left.calcKnot(left.x1), right);
}

//==============================================================================
//                      CALCULATION HELPERS: C2-CONTINUITY CHECKS
//==============================================================================

bool isC2Continuous(
    const OpenSim::CubicSpline& left,
    const OpenSim::CubicSpline& right)
{
    const double x = left.x1;
    return isC1Continuous(left, right) &&
           isNumEq(left.calcDerivative(x, 2), right.calcDerivative(x, 2));
}

//==============================================================================
//                      CALCULATION HELPERS: MONOTONICITY CHECKS
//==============================================================================

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
    return sign * dydx0 >= 0. && sign * dydx1 >= 0.;
}

bool smoothMonotonicCurveExists(
    const OpenSim::CurveKnot& p0,
    const OpenSim::CurveKnot& p1)
{
    const double dy = p1.y - p0.y;
    return smoothMonotonicCurveExists(dy, p0.dydx, p1.dydx);
}

// Returns true if the Hermite interpolant is monotonic.
bool isHermiteInterpolantMonotonic(
    double dx,
    double dy,
    double dydx0,
    double dydx1, bool strict = false)
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

    const double alpha = dydx0 / dy * dx;
    const double beta  = dydx1 / dy * dx;

    if (strict) {
        return alpha * alpha + beta * beta < 8.75;
    }

    // Conditions from
    // "Monotone cubic spline interpolation for functions with a strong
    // gradient" by Francesc Arandiga.
    bool condition_0 = alpha + beta <= 3.;
    bool condition_1 =
        alpha * (alpha + beta - 6.) + (beta - 3.) * (beta - 3.) < 0.;
    return condition_0 || condition_1;

    // Or From other literature...
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

bool isHermiteInterpolantMonotonic(
    const OpenSim::CurveKnot& left,
    const OpenSim::CurveKnot& right, bool strict = false)
{
    return isHermiteInterpolantMonotonic(
        right.x - left.x,
        right.y - left.y,
        left.dydx,
        right.dydx, strict);
}

//==============================================================================
//                      HERMITE INTERPOLATION
//==============================================================================

OpenSim::CubicSpline::Coefficients calcCubicHermiteSplineCoeffs(
    const OpenSim::CurveKnot& p0,
    const OpenSim::CurveKnot& p1)
{
    const double dx = p1.x - p0.x;
    const double dv = p1.dydx - p0.dydx;
    const double dy = p1.y - p0.y;

    const double c3 =
        -2. * (dy - p0.dydx * dx - 0.5 * dx * dv) / std::pow(dx, 3);
    const double c2 = (dv / dx - 3. * c3 * dx) / 2.;

    return {
        p0.y,
        p0.dydx,
        c2,
        c3,
    };
}

//==============================================================================
//                      MONO HERMITE INTERPOLATION
//==============================================================================

double calcBehavedMiddleDerivative(
    const OpenSim::CurvePoint& left,
    const OpenSim::CurvePoint& mid,
    const OpenSim::CurvePoint& right)
{
    const double leftSecant  = left.calcSecantLine(mid);
    const double rightSecant = right.calcSecantLine(mid);
    /* std::cout << "leftSecant = " << leftSecant << std::endl; */
    /* std::cout << "rightSecant = " << rightSecant << std::endl; */

    const bool oppositeSigns = (leftSecant > 0 && rightSecant < 0) ||
                               (leftSecant < 0 && rightSecant > 0);

    return oppositeSigns ? 0. : (leftSecant + rightSecant) / 2.;
}

//==============================================================================
//                      PROCESSING CURVINESS
//==============================================================================

// Extrapolate y coordinate using the knot derivative.
double calcExtrapolated(const OpenSim::CurveKnot& knot, double x)
{
    return knot.y + knot.dydx * (x - knot.x);
}

// Calculate the point where the knots intersect when extrapolated.
OpenSim::CurvePoint calcInterceptPoint(
    const OpenSim::CurveKnot& left,
    const OpenSim::CurveKnot& right)
{
    const double x =
        std::abs(left.dydx - right.dydx) < SimTK::Eps ||
                std::abs(left.y - right.y) < SimTK::Eps
            // If zero gradient: take x between left and right.
            ? left.x + (right.x - left.x) / 2.
            // Solve for x:
            // y0 + (x - x0) * dydx0 = y1 + (x - x1) * dydx1
            : (right.y - left.y + left.x * left.dydx - right.x * right.dydx) /
                  (left.dydx - right.dydx);

    const double yL = calcExtrapolated(left, x);
    const double yR = calcExtrapolated(right, x);

    opensim_assert(isNumEq(yL, yR), "error computing intercept point");

    return {x, yL};
}

// Calculate the intermediate knots that enforce the curviness.
std::pair<OpenSim::CurveKnot, OpenSim::CurveKnot> calcCurvyEnforcingKnots(
    const OpenSim::CurveKnot& left,
    const OpenSim::CurveKnot& right,
    double curviness)
{
    OpenSim::CurvePoint intercept = calcInterceptPoint(left, right);

    OpenSim::CurveKnot cLeft(
        left.calcInterpolated(intercept, 1. - curviness),
        left.dydx);
    OpenSim::CurveKnot cRight(
        right.calcInterpolated(intercept, 1. - curviness),
        right.dydx);

    return std::make_pair(cLeft, cRight);
}

// Calculate control knots from curvy control points.
// TODO awkward that curviness is between knots... one is ignored.
std::vector<OpenSim::CurveKnot> calcKnotsFromControlPoints(
    const std::vector<OpenSim::MuscleCurveControlPoint>& ctrlPts)
{
    std::vector<OpenSim::CurveKnot> knots;

    for (auto p = ctrlPts.begin(); p != ctrlPts.end(); p++) {
        // Push all control points as knots.
        auto left = p;
        knots.push_back(*left);

        // Check if the curviness between this and the next knot is set.
        auto right = p + 1;
        if (right == ctrlPts.end() || !(*right).isCurvy()) {
            continue;
        }

        // Insert the extra knots that enforce the curviness.
        auto curvyPoints =
            calcCurvyEnforcingKnots(*left, *right, (*right).curviness);
        knots.push_back(curvyPoints.first);
        knots.push_back(curvyPoints.second);
    }

    return knots;
}

//==============================================================================
//                      QUADRATIC BEZIER FITTING
//==============================================================================

OpenSim::CubicSpline calcCubicSplineFromQuadraticBezier(
    double p0,
    double p1,
    double p2)
{
    OpenSim::CubicSpline s;
    s.x0         = 0.;
    s.x1         = 1.;
    s.y0Integral = 0.;
    // B(t) = P1 + (1-t)^2(P0-P1) + t^2(P2-P1)
    s.coeff = {p0, 2. * (p1 - p0), p0 - 2. * p1 + p2, 0.};
    return s;
}

std::pair<OpenSim::CubicSpline, OpenSim::CubicSpline> calcQuadraticBezierSpline(
    const OpenSim::CurveKnot& left,
    const OpenSim::CurveKnot& right)
{
    const SimTK::Vec2 p0{left.x, left.y};

    const OpenSim::CurvePoint mid = calcInterceptPoint(left, right);
    const SimTK::Vec2 p1{mid.x, mid.y};

    const SimTK::Vec2 p2{right.x, right.y};

    opensim_assert(
        left.x < mid.x && mid.x < right.x,
        "Intercept line falls outside of convex area on x-axis");

    opensim_assert(
        std::min(left.y, right.y) - 1e-13 < mid.y &&
            mid.y < std::max(left.y, right.y) + 1e-13,
        "Intercept line falls outside of convex area on y-axis");

    return std::make_pair(
        calcCubicSplineFromQuadraticBezier(p0[0], p1[0], p2[0]),
        calcCubicSplineFromQuadraticBezier(p0[1], p1[1], p2[1]));
}

//==============================================================================
//                      C1 MONO CUBIC SPLINE ALGO
//==============================================================================

OpenSim::CurveKnot calcCurveKnotWithMeanDerivative(
    const OpenSim::CurvePoint& left,
    const OpenSim::CurvePoint& mid,
    const OpenSim::CurvePoint& right)
{
    const double left_dx = mid.x - left.x;
    const double left_dy = mid.y - left.y;

    const double right_dx = mid.x - right.x;
    const double right_dy = mid.y - right.y;

    const double dydx = calcBehavedMiddleDerivative(left, mid, right);
    return {mid, dydx};
}

OpenSim::CurveKnot calcCurveKnotWithMeanDerivative(
    const OpenSim::QuadraticBezierCurve shape,
    const std::vector<double>& segmentsUs,
    const std::vector<double>::iterator iter)
{
    if (iter == segmentsUs.begin()) {
        return shape.startKnot();
    }
    if (iter + 1 == segmentsUs.end()) {
        return shape.endKnot();
    }
    return calcCurveKnotWithMeanDerivative(
        shape.calcPoint(*(iter - 1)),
        shape.calcPoint(*iter),
        shape.calcPoint(*(iter + 1)));
}

std::vector<double> calcC1CubicMonoSplineUAlgo(
    const OpenSim::QuadraticBezierCurve& shape,
    size_t maxNumSegments)
{
    std::vector<double> segmentsUs = {0., 0.25, 0.5, 1.};
    for (auto it = segmentsUs.begin(); it < segmentsUs.end() - 1; ) {
        const double uL = *it;
        const OpenSim::CurveKnot& left =
            calcCurveKnotWithMeanDerivative(shape, segmentsUs, it);

        const double uR = *++it;
        const OpenSim::CurveKnot& right =
            calcCurveKnotWithMeanDerivative(shape, segmentsUs, it);

        bool segmentAccepted = isHermiteInterpolantMonotonic(left, right, true);
        if (segmentAccepted) {
            double y0Integral = 0.;
            OpenSim::CubicSpline s(left, right, y0Integral);
            segmentAccepted &= shape.isAccurateWithinTol(s);
        }

        if (segmentAccepted) {
            continue;
        }

        it = --segmentsUs.insert(it, (uL + uR) / 2.);

        opensim_assert(
            segmentsUs.size() <= maxNumSegments,
            "Failed to create C1 cubic splines");
    }
    std::cout << "Succesfully fitted " << segmentsUs.size()
              << " spline segments to curve shape using u = ";
    for (double u : segmentsUs) {
        std::cout << u << ", ";
    }
    std::cout << " }" << std::endl;
    return segmentsUs;
}

void calcC1CubicMonoSplineAlgo(
    const OpenSim::QuadraticBezierCurve& shape,
    std::vector<OpenSim::CubicMonoSpline>& splines,
    size_t maxNumSegments)
{
    std::vector<double> segmentsUs =
        calcC1CubicMonoSplineUAlgo(shape, maxNumSegments);
    for (auto it = segmentsUs.begin(); it < segmentsUs.end() - 1;) {
        const OpenSim::CurveKnot& left =
            calcCurveKnotWithMeanDerivative(shape, segmentsUs, it);
        const OpenSim::CurveKnot& right =
            calcCurveKnotWithMeanDerivative(shape, segmentsUs, ++it);

        const double y0Integral =
            splines.empty() ? 0. : splines.end()--->calcIntegral(right.x);
        splines.push_back({left, right, y0Integral});
    }
}

} // namespace

namespace OpenSim
{

//==============================================================================
//                      CURVE POINT
//==============================================================================

double CurvePoint::calcSecantLine(const CurvePoint& other) const
{
    return (other.y - y) / (other.x - x);
}

CurvePoint CurvePoint::calcInterpolated(const CurvePoint& other, double u) const
{
    return {
        x + (other.x - x) * u,
        y + (other.y - y) * u,
    };
}

std::ostream& operator<<(std::ostream& os, const CurvePoint& pt)
{
    return os << "CurvePoint{"
              << "x: " << pt.x << ", "
              << "y: " << pt.y << "}";
}

//==============================================================================
//                      CURVE KNOT
//==============================================================================

std::ostream& operator<<(std::ostream& os, const CurveKnot& knot)
{
    return os << "CurveKnot{"
              << "x: " << knot.x << ", "
              << "y: " << knot.y << ", "
              << "dydx: " << knot.dydx << "}";
}

//==============================================================================
//                      MUSCLE CURVE CONTROL POINT
//==============================================================================

bool MuscleCurveControlPoint::isCurvy() const
{
    // TODO throw if oob.
    return MIN_CURVINESS <= curviness
        && curviness <= MAX_CURVINESS;
}

std::ostream& operator<<(
    std::ostream& os,
    const MuscleCurveControlPoint& ctrlPt)
{
    return os << "MuscleCurveControlPoint{"
              << "x: " << ctrlPt.x << ", "
              << "y: " << ctrlPt.y << ", "
              << "dydx: " << ctrlPt.dydx << ", "
              << "curviness: " << ctrlPt.curviness << "}";
}

//==============================================================================
//                      Cubic Spline
//==============================================================================

CubicSpline::CubicSpline(
    const CurveKnot& left,
    const CurveKnot& right,
    double& startIntegralValue) :
    x0(left.x),
    coeff(calcCubicHermiteSplineCoeffs(left, right)),
    y0Integral(startIntegralValue), x1(right.x)
{
    startIntegralValue = calcIntegral(x1);
    opensim_assert(
        isC1Continuous(left, *this),
        "Hermite interpolation failed: Start knot not connected.");
    opensim_assert(
        isC1Continuous(*this, right),
        "Hermite interpolation failed: End knot not connected.");
}

bool CubicSpline::isMonotonic() const
{
    const double dx    = x1 - x0;
    const double dy    = calcValue(x1) - calcValue(x0);
    const double dydx0 = calcDerivative(x0, 1);
    const double dydx1 = calcDerivative(x1, 1);
    const bool mono    = isHermiteInterpolantMonotonic(dx, dy, dydx0, dydx1);

    /* bool checkMono = true; */
    /* for (size_t i = 0; i < 1000; ++i) { */
    /*     checkMono &= */
    /*     calcDerivative(x0 + dx / 1000. * static_cast<double>(i),1) * dy >=
     * 0.; */
    /* } */
    /* opensim_assert(checkMono == mono, "failed to verify Monotonicity"); */

    return mono;
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

CurveKnot CubicSpline::calcKnot(double x) const
{
    return {x, calcValue(x), calcDerivative(x, 1)};
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

    opensim_assert(std::abs(yError) < eps, "Failed to invert spline segment");
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

std::ostream& operator<<(std::ostream& os, const CubicSpline& spline)
{
    return os << "CubicSpline{"
              << "x0: " << spline.x0 << ", "
              << "coeff: {" << spline.coeff.at(0) << ", " << spline.coeff.at(1)
              << ", " << spline.coeff.at(2) << ", " << spline.coeff.at(3)
              << "}, "
              << "y0Integral: " << spline.y0Integral << ", "
              << "x1: " << spline.x1 << "}";
}

//==============================================================================
//                      Cubic Mono Spline
//==============================================================================

CubicMonoSpline::CubicMonoSpline(CubicSpline spline) : CubicSpline(spline)
{
    opensim_assert(
        spline.isMonotonic(),
        "Failed to construct monotonic spline: Cubic spline is not monotonic");
}

std::ostream& operator<<(std::ostream& os, const CubicMonoSpline& spline)
{
    return os << "Monotonic" << static_cast<CubicSpline>(spline);
}

//==============================================================================
//                      Cubic Mono Spline Storage
//==============================================================================

SmoothSegmentedCubicMonoSplineData::SmoothSegmentedCubicMonoSplineData(
    const CurveShape& curveShape,
    size_t maxNumSegments)
{
    std::vector<CubicMonoSpline> splines;
    for (const QuadraticBezierCurve& shape : curveShape.getSegments()) {
        calcC1CubicMonoSplineAlgo(shape, splines, maxNumSegments);
    }
    *this = SmoothSegmentedCubicMonoSplineData(splines);
    if (size() == 0) {
        return;
    }
}

SmoothSegmentedCubicMonoSplineData::SmoothSegmentedCubicMonoSplineData(
    const std::vector<CubicMonoSpline>& splines)
{
    for (const CubicMonoSpline& s : splines) {
        appendChecked(s);
    }
}

void SmoothSegmentedCubicMonoSplineData::appendChecked(CubicMonoSpline spline)
{
    // Verify that segments are C1 continuous.
    if (size() > 0) {
        opensim_assert(
            isC1Continuous(at(size() - 1), spline),
            "Segment node C1 continuity check failed");
    }

    // x1 of last segment is same as x0 of first segment, so it can be removed.
    if (!_data.empty()) {
        _data.pop_back();
    }

    for (double d : reinterpret_cast<
             const std::array<double, CubicMonoSpline::dataSize()>&>(spline)) {
        _data.push_back(d);
    }

    // Verify that segments are C1 continuous. TODO this can be removed.
    if (size() > 1) {
        opensim_assert(
            isC1Continuous(at(size() - 2), at(size() - 1)),
            "Segment node C1 continuity check failed after writing as raw "
            "doubles");
    }

    // Verify expected number of elements in container.
    constexpr size_t alignment = CubicMonoSpline::dataSize() - 1;
    opensim_assert(
        _data.size() % alignment == 1,
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
        _data.at(index * alignment));
}

SimTK::Vec2 SmoothSegmentedCubicMonoSplineData::getDomain() const
{
    return {_data.front(), _data.back()};
}

//==============================================================================
//              QuadraticBezierCurve
//==============================================================================

QuadraticBezierCurve::QuadraticBezierCurve(
    const CurveKnot& left,
    const CurveKnot& right) :
    _start(left),
    _end(right)
{
    auto xySplines = calcQuadraticBezierSpline(left, right);
    _x             = xySplines.first;
    _y             = xySplines.second;
}

CurvePoint QuadraticBezierCurve::calcPoint(double u) const
{
    if (u == 0.) {
        return {_x.calcValue(0.), _y.calcValue(0.)};
    }
    if (u == 1.) {
        return {_x.calcValue(1.), _y.calcValue(1.)};
    }
    const double x0 = _x.calcValue(0.);
    const double x1 = _x.calcValue(1.);
    const double dx = (x1 - x0) * u;
    const double x  = x0 + dx;
    const double s  = _x.calcInverseValue(x);
    const double y  = _y.calcValue(s);
    return {x, y};
}

double QuadraticBezierCurve::calcValue(double x) const
{
    const double u = _x.calcInverseValue(x);
    return _y.calcValue(u);
}

const CurveKnot& QuadraticBezierCurve::startKnot() const
{
    return _start;
}

const CurveKnot& QuadraticBezierCurve::endKnot() const
{
    return _end;
}

bool OpenSim::QuadraticBezierCurve::isAccurateWithinTol(
    const OpenSim::CubicSpline& spline) const
{
    for (size_t i = 0; i < nAccuracySamples; ++i) {
        const double u =
            static_cast<double>(i) / static_cast<double>(nAccuracySamples);
        const double x     = spline.x0 + (spline.x1 - spline.x0) * u;
        const double y     = calcValue(x);
        const double error = y - spline.calcValue(x);
        const double maxError =
            std::max(std::abs(relAccuracy * y), absAccuracy);
        if (std::abs(error) > maxError) {
            return false;
        }
    }
    return true;
}

//==============================================================================
//              Curve Shape
//==============================================================================

CurveShape::CurveShape(std::vector<MuscleCurveControlPoint> ctrlPts) :
    CurveShape(calcKnotsFromControlPoints(ctrlPts))
{}

CurveShape::CurveShape(std::vector<CurveKnot> knots)
{
    for (size_t i = 0; i < knots.size() - 1;) {
        const CurveKnot& k0 = knots[i++];
        const CurveKnot& k1 = knots[i];

        opensim_assert(
            smoothMonotonicCurveExists(k0, k1),
            "Failed to create curve shape: Knots not monotonic");
        _segments.push_back(QuadraticBezierCurve(k0, k1));
    }
}

//==============================================================================
//                  Smooth Segmented Cubic Mono Spline
//==============================================================================

SmoothSegmentedCubicMonoSpline::SmoothSegmentedCubicMonoSpline(
    const CurveShape& shape,
    size_t maxNumSegments) :
    _splines(shape, maxNumSegments)
{
    for (auto& s : shape.getSegments()) {
        const auto pStart = s.calcPoint(0.);
        std::cout << "pStart = " << pStart << ", yStart = " << calcValue(pStart.x) - pStart.y << std::endl;
        opensim_assert(
            isNumEq(calcValue(pStart.x), pStart.y),
            "Failed to construct CubicMonoSpline: Start Knot points not "
            "fitted.");
        const auto pEnd = s.calcPoint(1.);
        std::cout << "pEnd = " << pEnd << ", yEnd = " << calcValue(pEnd.x) - pEnd.y << std::endl;
        opensim_assert(
            isNumEq(calcValue(pEnd.x), pEnd.y),
            "Failed to construct CubicMonoSpline: End Knot points not fitted.");
    }
}

SimTK::Vec2 SmoothSegmentedCubicMonoSpline::getDomain() const
{
    return _splines.getDomain();
}

const CubicMonoSpline& SmoothSegmentedCubicMonoSpline::findSegment(
    double x) const
{
    // TODO use algo
    size_t s = _splines.size();
    for (size_t idx = 0; idx < s; ++idx) {
        if (_splines.at(idx).x1 >= x) {
            return _splines.at(idx);
        }
    }
    return _splines.at(s - 1);
}

double SmoothSegmentedCubicMonoSpline::calcValue(double x) const
{
    return findSegment(x).calcValue(x);
}

} // namespace OpenSim

//==============================================================================
//                  OLD
//==============================================================================

/* opensim_assert(sign * dydx0 >= 0. && "Monotonicity check failed"); */

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

/* namespace */
/* { */
/* bool splineConverger( */
/*     const OpenSim::CurveControlPoint& pStart, */
/*     const OpenSim::CurveControlPoint& pEnd, */
/*     std::vector<OpenSim::CurveControlPoint>& points, */
/*     double eps     = 1e-13, */
/*     size_t maxIter = 100) */
/* { */
/*     // Initialize x coordinates. */
/*     const double dx = pEnd.x - pEnd.x; */
/*     for (size_t i = 0; i < points.size(); ++i) { */
/*         points[i].x = pStart.x + dx * (static_cast<double>(points.size()) +
 * 1) / */
/*                                      (static_cast<double>(points.size()) -
 * 1); */
/*     } */
/*     // Initialize y coordinates. */
/*     const double dy = pEnd.y - pEnd.y; */
/*     for (OpenSim::CurveControlPoint& p : points) { */
/*         p.y    = pStart.y + dy / dx * (p.x - pStart.x); */
/*         p.dydx = dy / dx; */
/*     } */

/*     // Iteratively find control points. */
/*     for (size_t n = 0; n < maxIter; ++n) { */
/*         double err = 0.; */
/*         for (size_t i = 0; i < maxIter; ++i) { */
/*             const OpenSim::CurveControlPoint pLeft = */
/*                 i == 0 ? pStart : points.at(i - 1); */
/*             const OpenSim::CurveControlPoint pRight = */
/*                 i == points.size() - 1 ? pEnd : points.at(i + 1); */
/*             OpenSim::CurveControlPoint pMid = points.at(i); */

/*             const double leftSecantLine = */
/*                 (pMid.y - pLeft.y) / (pMid.x - pLeft.x); */
/*             const double rightSecantLine = */
/*                 (pRight.y - pMid.y) / (pRight.x - pMid.x); */

/*             pMid.dydx = 2. * leftSecantLine * rightSecantLine / */
/*                         (leftSecantLine + rightSecantLine); */

/*             const double yMidNew = */
/*                 pLeft.y + (pLeft.dydx + pMid.dydx) / 2. * (pMid.x - pLeft.x);
 */
/*             err    = std::max(std::abs(yMidNew - pMid.y), err); */
/*             pMid.y = yMidNew; */
/*         } */
/*         if (err < eps) { */
/*             return true; */
/*         } */
/*     } */
/*     return false; */
/* } */
/* } // namespace */
