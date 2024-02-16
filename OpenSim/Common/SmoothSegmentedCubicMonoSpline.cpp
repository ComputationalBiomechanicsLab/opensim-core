#include "SmoothSegmentedCubicMonoSpline.h"

#include "OpenSim/Common/Assertion.h"
#include <array>
#include <cmath>
#include <stdexcept>

#include <SimTKcommon/Scalar.h>
#include <SimTKcommon/internal/NTraits.h>

namespace
{

constexpr double EPS = 1e-13;
// Allowed error when checking acceleration during C2-continuity check.
// It is slightly larger due to larger values of acceleration.
constexpr double EPS_C2         = 1e3 * EPS;
constexpr double MAX_ITER       = 100;
constexpr double MIN_SEGMENT_DX = 1e-6;

// TODO remove this
void opensim_assert(bool cond, std::string msg)
{
    if (!cond) {
        throw std::runtime_error(msg);
    }
}

template <typename Y, typename U, typename... Args>
std::vector<Y> ConstructFromTwoElements(
    const std::vector<U>& elements,
    Args... args)
{
    std::vector<Y> y;
    for (size_t i = 1; i < elements.size(); ++i) {
        const U& left  = elements.at(i - 1);
        const U& right = elements.at(i);
        y.push_back(Y(left, right, args...));
    }
    return y;
}

struct X0Coordinate
{
    X0Coordinate(double xCoordinate) : x(xCoordinate)
    {}

    X0Coordinate(const OpenSim::CubicMonoSpline& spline) : x(spline.x0)
    {}
    double x;
};

struct Y0Coordinate
{
    Y0Coordinate(double yCoordinate) : y(yCoordinate)
    {}

    Y0Coordinate(const OpenSim::CubicMonoSpline& spline) :
        y(spline.coeff.front())
    {}
    double y;
};

//==============================================================================
//                      CALCULATION HELPERS: C1-CONTINUITY CHECKS
//==============================================================================
bool isNumEq(double lhs, double rhs, double eps = EPS)
{
    return std::abs(lhs - rhs) < eps;
}

bool isC1Continuous(
    const OpenSim::CurveKnot& left,
    const OpenSim::CurveKnot& right)
{
    /* std::cout << "ISC1CONTINUOUS" */
    /*           << "\n" */
    /*           << "left = " << left << "\n" */
    /*           << "right = " << right; */
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

template<typename T>
bool isC1Continuous(const std::vector<T>& splines)
{
    for (size_t i = 0; i + 1 < splines.size(); ++i) {
        if (!isC1Continuous(splines[i], splines[i + 1])) {
            return false;
        }
    }
    return true;
}

template<typename T>
bool isC1Continuous(
    const OpenSim::CurveKnot& startKnot,
    const std::vector<T>& splines,
    const OpenSim::CurveKnot& endKnot)
{
    if (splines.empty()) {
        return isC1Continuous(startKnot, endKnot);
    }
    return isC1Continuous(startKnot, splines.at(0)) &&
           isC1Continuous(splines) && isC1Continuous(splines.back(), endKnot);
}

//==============================================================================
//                      CALCULATION HELPERS: C2-CONTINUITY CHECKS
//==============================================================================

bool isC2Continuous(
    const OpenSim::CubicSpline& left,
    const OpenSim::CubicSpline& right)
{
    const double x = left.x1;
    std::cout << "start C2ContinuousCheck: "
        << "x = " << x << ", "
        << "yL = " << left.calcValue(x) << ", "
        << "dyL = " << left.calcDerivative(x, 1) << ", "
        << "ddyL = " << left.calcDerivative(x, 2) << ", "
        << "yR = " << right.calcValue(x) << ", "
        << "dyR = " << right.calcDerivative(x, 1) << ", "
        << "ddyR = " << right.calcDerivative(x, 2) << ", "
        <<std::endl;
    return isC1Continuous(left, right) && isNumEq(
                                              left.calcDerivative(x, 2),
                                              right.calcDerivative(x, 2),
                                              EPS_C2);
}

template <typename T>
bool isC2Continuous(const std::vector<T>& splines)
{
    std::cout << "isC2Continuous : " << splines.size() << std::endl;
    for (size_t i = 0; i + 1 < splines.size(); ++i) {
        std::cout << "s : " << splines[i] << std::endl;
        if (!isC2Continuous(splines[i], splines[i + 1])) {
            return false;
        }
    }
    return true;
}

template <typename T>
bool isC2Continuous(
    const OpenSim::CurveKnot& startKnot,
    const std::vector<T>& splines,
    const OpenSim::CurveKnot& endKnot)
{
    if (splines.empty()) {
        return isC1Continuous(startKnot, endKnot);
    }
    return isC1Continuous(startKnot, splines.at(0)) &&
           isC2Continuous(splines) && isC1Continuous(splines.back(), endKnot);
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
    double dydx1,
    bool strict = false)
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
    const OpenSim::CurveKnot& right,
    bool strict = false)
{
    return isHermiteInterpolantMonotonic(
        right.x - left.x,
        right.y - left.y,
        left.dydx,
        right.dydx,
        strict);
}

bool isInvertible(const std::vector<OpenSim::CubicMonoSpline>& splines)
{

    auto CalcDirection = [](const OpenSim::CubicMonoSpline& left,
                            const OpenSim::CubicMonoSpline& right) {
        const double x = right.x0;
        return left.calcValue(x) < right.calcValue(x);
    };

    bool increasing = false;
    bool decreasing = false;
    for (auto s = splines.begin(); s + 1 != splines.end(); ++s) {
        increasing |= CalcDirection(*s, *(s + 1));
        decreasing |= CalcDirection(*(s + 1), *s);
    }

    return increasing != decreasing;
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

double calcDerivativeFromSplineCoefficients(
    const OpenSim::CubicSpline::Coefficients& coeff,
    double x0,
    double x,
    size_t order)
{
    double dx = x - x0;

    if (std::isnan(dx)) {
        return SimTK::NaN;
    }

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

//==============================================================================
//              Curve Shape
//==============================================================================

std::vector<OpenSim::CurveKnot> calcMonotonicSegmentKnots(
const OpenSim::C2ContinuousSegmentedCurve& curve)
{
    std::vector<double> knotXCoords = curve.calcMonotonicSegmentXValues();
    for (auto it = knotXCoords.begin(); it + 1 != knotXCoords.end(); ++it) {
        std::cout << "xk = " << *it << ", "
        << *(it+1) << std::endl;
        opensim_assert(
                *it < *(it+1), "Curve s' monotonic curve segments x coordinates are not increasing");
    }
    std::vector<OpenSim::CurveKnot> knots;
    for (double x : knotXCoords) {
        knots.push_back({x, curve.calcValue(x), curve.calcDerivative(x,1)});
    }
    for (size_t i = 0; i < knots.size() - 1;) {
        const OpenSim::CurveKnot& k0 = knots[i++];
        const OpenSim::CurveKnot& k1 = knots[i];

        opensim_assert(
            smoothMonotonicCurveExists(k0, k1),
            "Failed to create curve segments: Knots not monotonic");
    }
    return knots;
}

//==============================================================================
//                      MAX DIFFERENCE BETEEN SHAPE AND SPLINE
//==============================================================================

double calcMaxAbsDifference(
    const OpenSim::C2ContinuousSegmentedCurve& curve,
    const OpenSim::CubicSpline& spline,
    double maxFitError,
    double xEstimate,
    size_t maxIter,
    double eps)
{
    std::cout << "maxFitError = " << maxFitError << std::endl;
    std::cout << "eps = " << eps << std::endl;
    opensim_assert(eps < maxFitError, "Invalid error fit parameter");
    const double diffInit =
        std::abs(curve.calcValue(xEstimate) - spline.calcValue(xEstimate));
    return diffInit; // TODO need to fix this

    /* for (size_t i = 0; i < maxIter; ++i) { */
    /*     // Check if x is out of bounds. */
    /*     const bool oob = spline.x0 > xEstimate || spline.x1 < xEstimate; */
    /*     // Clamp x to segment domain. */
    /*     xEstimate = std::min(std::max(xEstimate, spline.x0), spline.x1); */
    /*     // Compute difference between spline and curve. */
    /*     const double diff = */
    /*         curve.calcValue(xEstimate) - spline.calcValue(xEstimate); */
    /*     const double absDiff              = std::abs(diff); */
    /*     const bool largerThanAllowedError = absDiff > maxFitError; */
    /*     // Compute gradient of difference. */
    /*     const double diffDerivative = curve.calcDerivative(xEstimate, 1) - */
    /*                                   spline.calcDerivative(xEstimate, 1); */
    /*     // Stop if larger than allowed error, or out-of-bounds or if converged. */
    /*         std::cout */
    /*             << "x0 = " << spline.x0 << ", " */
    /*             << "xMax = " << xEstimate << ", " */
    /*             << "x1 = " << spline.x1 << ", " */
    /*             << "diffInit = " << diffInit << ", " */
    /*             << "absDiff = " << absDiff << ", " */
    /*             << "x1 = " << spline.x1 */
    /*             << std::endl; */
    /*     if (largerThanAllowedError || oob || absDiff < eps || */
    /*         std::abs(diffDerivative) < SimTK::Eps) { */
    /*         opensim_assert(diffInit < absDiff + eps, "failed to find max error"); */
    /*         return absDiff; */
    /*     } */
    /*     // Take step to maximise the difference. */
    /*     const double step = diff / diffDerivative; */
    /*     xEstimate += step; */
    /* } */

    // TODO use opensim_assert
    throw std::runtime_error(
        "Failed to compute max abs difference between spline and curve");
}

double calcMaxAbsDifference(
    const OpenSim::C2ContinuousSegmentedCurve& curve,
    const OpenSim::CubicSpline& spline,
    double maxFitError,
    size_t maxIter = MAX_ITER,
    double eps     = EPS)
{
    double xEstimate        = spline.x0;
    double maxAbsDiff       = -SimTK::Infinity;
    const size_t searchGrid = 10;
    for (size_t i = 0; i <= searchGrid; ++i) {
        const double xTest = spline.x0 + (spline.x1 - spline.x0) /
                                             static_cast<double>(searchGrid) *
                                             static_cast<double>(i);
        const double absDiff =
            std::abs(spline.calcValue(xTest) - curve.calcValue(xTest));
        if (absDiff > maxAbsDiff) {
            xEstimate  = xTest;
            maxAbsDiff = absDiff;
        }
        if (maxAbsDiff > maxFitError) {
            std::cout << "xTest " << xEstimate << ", "
                      << "ys = " << spline.calcValue(xTest) << ", "
                      << "yc = " << curve.calcValue(xTest) << ", "
                      << "absDiff = " << absDiff << ", "
                      << "maxabsDiff = " << maxAbsDiff << "\n";
            /* std::cout << "maxFitError = " << maxFitError << std::endl; */
            return maxAbsDiff;
        }
    }
    /* std::cout << "calcMaxAbsDifference:" */
    /*           << "maxDiffGrid = " << maxAbsDiff << ", "; */
    maxAbsDiff = calcMaxAbsDifference(
        curve,
        spline,
        maxFitError,
        xEstimate,
        maxIter,
        eps);
    /* std::cout << "maxDiffIter = " << maxAbsDiff << "\n"; */
    return maxAbsDiff;
}

//==============================================================================
//                      NATURAL SPLINE
//==============================================================================

// Overwrites the derivative values at the knots.
std::vector<OpenSim::CurveKnot>& calcNaturalCubicSplineKnotDerivatives(
   const OpenSim::C2ContinuousSegmentedCurve& curve,
    std::vector<OpenSim::CurveKnot>& knots)
{
    std::cout << "start natural spline fitting" << std::endl;
    for (OpenSim::CurveKnot& k: knots) {
        k.y = curve.calcValue(k.x);
        k.dydx = curve.calcDerivative(k.x,1);
    }
    const size_t nPoints = knots.size();
    const size_t n       = nPoints - 1;

    std::vector<double> a;
    a.reserve(n + 1);
    for (size_t i = 0; i < nPoints; ++i)
        a.push_back(knots.at(i).y);

    std::vector<double> b, d, h;
    h.reserve(n);
    b.reserve(n);
    d.reserve(n);
    for (size_t i = 0; i < n; ++i)
        h.push_back(SimTK::NaN);
    for (size_t i = 0; i < n; ++i)
        b.push_back(SimTK::NaN);
    for (size_t i = 0; i < n; ++i)
        d.push_back(SimTK::NaN);
    for (size_t i = 0; i < n; i++)
        h[i] = knots.at(i + 1).x - knots.at(i).x;

    std::vector<double> alpha;
    alpha.reserve(n);
    for (size_t i = 1; i < n; i++)
        alpha[i] =
            3.0 / h[i] * (a[i + 1] - a[i]) - 3.0 / h[i - 1] * (a[i] - a[i - 1]);
    /* std::cout << "what" << std::endl; */

    std::vector<double> c, l, mu, z;
    c.reserve(n + 1);
    for (size_t i = 0; i < n + 1; ++i)
        c.push_back(SimTK::NaN);
    l.reserve(n + 1);
    z.reserve(n + 1);
    mu.reserve(n + 1);
    l[0] = 1, mu[0] = z[0] = 0;
    for (size_t i = 1; i <= n - 1; i++) {
        l[i] =
            2. * (knots.at(i + 1).x - knots.at(i - 1).x) - h[i - 1] * mu[i - 1];
        mu[i] = h[i] / l[i];
        z[i]  = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
    }
    /* std::cout << "right" << std::endl; */
    l[n] = 1.;
    z[n] = 0.;
    c[n] = 0.;
    for (ptrdiff_t j = n - 1; j >= 0; --j) {
        c[j] = z[j] - mu[j] * c[j + 1];
        b[j] = (a[j + 1] - a[j]) / h[j] - (h[j] * (c[j + 1] + 2. * c[j])) / 3.0;
        d[j] = (c[j + 1] - c[j]) / (3.0 * h[j]);
    }
    std::cout << "natural knots:" << std::endl;
    for (size_t i = 0; i < nPoints; ++i) {
        if (i == n) {
            const double dx  = knots.at(i).x - knots.at(i - 1).x;
            knots.at(i).dydx = b.at(i - 1) + 2. * dx * c.at(i - 1) +
                               3. * dx * dx * d.at(i - 1);
        } else {
            knots.at(i).dydx = b.at(i);
        }
        std::cout << "    knot[" << i << "] = " << knots.at(i) << std::endl;
    }
    return knots;
}

//==============================================================================
//                  FINDING THE APPROXOMATING SPLINE SEGMENTS
//==============================================================================

double calcFittingErrorBound(
    const OpenSim::C2ContinuousSegmentedCurve& curve,
    double accuracy)
{
    SimTK::Vec2 domain = curve.getCurveDomain();
    const double y0    = curve.calcValue(domain(0));
    double y0Integral  = 0.; // TODO remove
    OpenSim::CubicSpline refSpline(
        OpenSim::CurveKnot{domain(0), 0., 0.},
        OpenSim::CurveKnot{domain(1), 0., 0.},
        y0Integral);
    return accuracy * calcMaxAbsDifference(curve, refSpline, SimTK::Infinity);
}

bool evaluateSpline(
    const OpenSim::C2ContinuousSegmentedCurve& curve,
    const OpenSim::CurveKnot& left,
    const OpenSim::CurveKnot& right,
    double maxFitError)
{
    const bool isMonotonic = isHermiteInterpolantMonotonic(left, right, true);
    // TODO add constructor without integral value?
    // class IntegrableCubicSpline
    double y0Integral = 0.;
    const OpenSim::CubicSpline spline(left, right, y0Integral);
    const bool isAccurate =
        calcMaxAbsDifference(curve, spline, maxFitError) < maxFitError;
    double error = calcMaxAbsDifference(curve, spline, maxFitError);
    std::cout << "    ACCURACY = " << error << " maxFitError = " << maxFitError
        << " monotonic = " << isMonotonic << " x+eps = " << curve.calcDerivative(left.x + EPS, 1)
              << std::endl;
    return isMonotonic && isAccurate;
}

bool updateGrid(
    std::vector<OpenSim::CurveKnot>& knots,
    const OpenSim::C2ContinuousSegmentedCurve& curve,
    double maxFittingError,
    size_t maxNumSegments)
{
    std::cout << "start grid update: max error = " << maxFittingError
              << std::endl;
    bool wasUpdated = false;
    for (size_t i = 0; i + 1 < knots.size(); ++i) {
        // Copy, dont take reference.
        const OpenSim::CurveKnot left  = knots.at(i);
        const OpenSim::CurveKnot right = knots.at(i + 1);
        const bool isOk =
            evaluateSpline(curve, left, right, maxFittingError);
        if (isOk) {
            // Nothing to do if requirements were met.
            continue;
        }
        // Did not meet the requirements: need to split this segment.
        wasUpdated = true;

        // Add new node to the grid in middle of left and right node:
        const double x                = (left.x + right.x) / 2.;
        const OpenSim::CurveKnot knot = {
            x,
            curve.calcValue(x),
            curve.calcDerivative(x, 1)};
        std::cout << "    added new knot:\n"
                  << "        left = " << left << "\n"
                  << "        knot = " << knot << "\n"
                  << "        right = " << right << "\n"
                  << "        size = " << knots.size() << "\n"
                  << "        maxNumSegments = " << maxNumSegments << "\n";
        knots.insert(knots.begin() + ++i, knot);
        bool segmentSizeTooSmall =
            std::min(knot.x - left.x, right.x - knot.x) < MIN_SEGMENT_DX;
        if (segmentSizeTooSmall) {
            opensim_assert(
                false,
                "Failed to refine grid: Segment size too small.");
        }

        // Throw if we exceed the max number of allowed segments.
        opensim_assert(
            knots.size() <= maxNumSegments,
            "Failed to refine curve-segment-grid: Exceeded mux number "
            "of segments");
        break;
    }
    return wasUpdated;
}

std::vector<OpenSim::CubicMonoSpline> calcSplineApproximationToCurve(
    const OpenSim::C2ContinuousSegmentedCurve& curve,
    double accuracy,
    size_t maxNumSegments)
{
    std::cout << "start spline approximation" << std::endl;
    // Convert the accuracy to a max allowed fitting error.
    const double maxFitError = calcFittingErrorBound(curve, accuracy);
    std::cout << "maxFitError = " << maxFitError << std::endl;

    // Use the curve segments as the initial guess for the spline knots.
    std::cout << "start calcMonotonicSegmentKnots" << std::endl;
    std::vector<OpenSim::CurveKnot> splineKnots =
        calcMonotonicSegmentKnots(curve);

    // Iteratively compute the natural spline, and refine the grid if needed.
    std::cout << "start algorithm" << std::endl;
    while (updateGrid(
        calcNaturalCubicSplineKnotDerivatives(curve, splineKnots),
        curve,
        maxFitError,
        maxNumSegments)) {
    }

    // Compute the Hermite interpolant connecting the knots.
    std::cout << "start make splines" << std::endl;
    double y0Integral = 0.;
    std::vector<OpenSim::CubicMonoSpline> splines =
        ConstructFromTwoElements<OpenSim::CubicMonoSpline>(
            splineKnots,
            y0Integral);

    // Verify that final spline segments are continuous.
    opensim_assert(
        isC2Continuous(splineKnots.front(), splines, splineKnots.back()),
        "C2 continuity check of spline segments failed");
    return splines;
}

//==============================================================================
//                      SOLVING DIFFERENTIABLE EQUATION
//==============================================================================

double solveDiffentiableScalarEquation(
    const std::function<OpenSim::ValueAndDerivative(double)>& lhs,
    const std::function<OpenSim::ValueAndDerivative(double)>& rhs,
    double xEstimate,
    size_t maxIter,
    double eps)
{
    for (size_t i = 0; i < maxIter; ++i) {
        // Compute difference between spline and curve.
        const OpenSim::ValueAndDerivative error =
            lhs(xEstimate) - rhs(xEstimate);
        // Stop if larger than allowed error, or out-of-bounds or if converged.
        if (error.value < eps || error.derivative < eps) {
            return xEstimate;
        }
        // Take step to minimize the difference.
        const double step = -error.value / error.derivative;
        xEstimate += step;
    }

    // TODO use opensim_assert
    throw std::runtime_error(
        "Failed to solve SmoothSegmentedCubicMonoSpline::calcValue(x) = rhs(x) "
        "for x");
}

} // namespace

namespace OpenSim
{

//==============================================================================
//                      VALUE AND DERIVATIVE
//==============================================================================

ValueAndDerivative operator-(
    const ValueAndDerivative& lhs,
    const ValueAndDerivative& rhs)
{
    return ValueAndDerivative{
        lhs.value - rhs.value,
        lhs.derivative - rhs.derivative};
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
    return isHermiteInterpolantMonotonic(calcKnot(x0), calcKnot(x1));
}

double CubicSpline::calcValue(double x) const
{
    return calcDerivativeFromSplineCoefficients(coeff, x0, x, 0);
}

double CubicSpline::calcDerivative(double x, size_t order) const
{
    return calcDerivativeFromSplineCoefficients(coeff, x0, x, order);
}

ValueAndDerivative CubicSpline::calcValueAndDerivative(double x) const
{
    return {
        calcDerivativeFromSplineCoefficients(coeff, x0, x, 0),
        calcDerivativeFromSplineCoefficients(coeff, x0, x, 1)};
}

CurveKnot CubicSpline::calcKnot(double x) const
{
    return {
        x,
        calcDerivativeFromSplineCoefficients(coeff, x0, x, 0),
        calcDerivativeFromSplineCoefficients(coeff, x0, x, 1)};
}

double CubicSpline::calcInverseValue(double y, double eps, size_t maxIter) const
{
    return solveDiffentiableScalarEquation(
        [&](double x) { return calcValueAndDerivative(x); },
        [=](double) -> ValueAndDerivative {
            return {y, 0.};
        },
        (x1 + x0) / 2.,
        maxIter,
        eps);
    /* double xEstimate = (x1 + x0) / 2.; */
    /* double yEstimate = SimTK::NaN; */
    /* double yError    = SimTK::Infinity; */
    /* for (size_t i = 0; */
    /*      std::abs(yError = (y - (yEstimate = calcValue(xEstimate)))) > eps &&
     */
    /*      i < maxIter; */
    /*      ++i) { */
    /*     const double xStep = yError / calcDerivative(xEstimate, 1); */
    /*     xEstimate += xStep; */
    /* } */

    /* opensim_assert(std::abs(yError) < eps, "Failed to invert spline
     * segment"); */
    /* return xEstimate; */
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
//                  Smooth Segmented Cubic Mono Spline
//==============================================================================

SmoothSegmentedCubicMonoSpline::SmoothSegmentedCubicMonoSpline(
    const C2ContinuousSegmentedCurve& curve,
    size_t maxNumSegments,
    double accuracy) :
    _splines(calcSplineApproximationToCurve(curve, accuracy, maxNumSegments)),
    _isInvertible(isInvertible(_splines))
{
    std::cout << "y 0 = " << curve.calcValue(1.) << "\n"
              << "ys0 = " << _splines.front().calcValue(1.) << "\n";
    std::cout << "Constructed SmoothSegmentedCubicMonoSpline using "
              << _splines.size() << " segments\n";
    for (auto s : _splines) {
        std::cout << "s = " << s << "\n";
    }
}

SimTK::Vec2 SmoothSegmentedCubicMonoSpline::getDomain() const
{
    return {_splines.front().x0, _splines.back().x1};
}

double SmoothSegmentedCubicMonoSpline::calcValue(double x) const
{
    /* std::cout << "calling CalcValue" << std::endl; */
    /* auto& s = findSegment(x); */
    /* std::cout << "ok" << std::endl; */
    /* double y = s */
    /*     .calcValue(x); */
    /* std::cout << "y= "<<y << std::endl; */
    return findSegment(x).calcValue(x);
}

ValueAndDerivative SmoothSegmentedCubicMonoSpline::calcValueAndDerivative(
    double x) const
{
    /* std::cout << "calling CalcValueAndDerivative" << std::endl; */
    return findSegment(x).calcValueAndDerivative(x);
}

double SmoothSegmentedCubicMonoSpline::calcInverseValue(double y) const
{
    std::cout << "calling CalcValueInverse" << std::endl;
    opensim_assert(
        _isInvertible,
        "SmoothSegmentedCubicMonoSpline is not invertible");
    return findInverseSegment(y).calcInverseValue(y);
}

const CubicMonoSpline& SmoothSegmentedCubicMonoSpline::findInverseSegment(
    double y) const
{
    Y0Coordinate searchY(y);
    auto spline = std::lower_bound(
        _splines.begin(),
        _splines.end(),
        searchY,
        [](const Y0Coordinate& a, const Y0Coordinate& b) -> bool {
            return a.y < b.y;
        });
    if (spline == _splines.begin()) {
        return *spline;
    }
    if (spline == _splines.end()) {
        return _splines.front();
    }
    return *spline;
}

const CubicMonoSpline& SmoothSegmentedCubicMonoSpline::findSegment(
    double x) const
{
    /* std::cout << "calling findSegment" << std::endl; */
    X0Coordinate searchX(x);
    auto spline = std::lower_bound(
        _splines.begin(),
        _splines.end(),
        searchX,
        [](const X0Coordinate& a, const X0Coordinate& b) -> bool {
            return a.x < b.x;
        });
    if (spline == _splines.begin()) {
        return *spline;
    }
    if (spline == _splines.end()) {
        return _splines.back();
    }
    return *(spline - 1);
}

double SmoothSegmentedCubicMonoSpline::solve(
    const std::function<ValueAndDerivative(double x)>& rhs,
    double xEstimate,
    size_t maxIter,
    double eps) const
{
    return solveDiffentiableScalarEquation(
        [&](double x) -> ValueAndDerivative {
            return calcValueAndDerivative(x);
        },
        rhs,
        xEstimate,
        maxIter,
        eps);
}

} // namespace OpenSim
