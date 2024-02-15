#ifndef OPENSIM_MUSCLECURVEPARAMS_H_
#define OPENSIM_MUSCLECURVEPARAMS_H_

#include "OpenSim/Common/SmoothSegmentedFunction.h"
#include "SimTKmath.h"
#include <array>
#include <cstddef>
#include <functional>
#include <memory>
#include <utility>

namespace OpenSim
{

struct ValueAndDerivative
{
    double value;
    double derivative;
};

ValueAndDerivative operator-(
    const ValueAndDerivative& lhs,
    const ValueAndDerivative& rhs);

//==============================================================================
//                  CURVE KNOT
//==============================================================================

struct CurveKnot
{
    CurveKnot() = default;

    CurveKnot(double xCoord, double yCoord, double derivative) :
        x(xCoord), y(yCoord), dydx(derivative)
    {}

    CurveKnot(double xCoord, ValueAndDerivative valueAndDerivative):
        x(xCoord), y(valueAndDerivative.value), dydx(valueAndDerivative.derivative)
    {}

    double x    = SimTK::NaN;
    double y    = SimTK::NaN;
    double dydx = SimTK::NaN;
};

std::ostream& operator<<(std::ostream& os, const CurveKnot& knot);

//==============================================================================
//              Curve Shape
//==============================================================================

using C2ContinuousSegmentedCurve = SmoothSegmentedFunction;
/* class C2ContinuousSegmentedCurve */
/* { */
/* public: */
/*     //========================================================================== */
/*     //              Curve Shape Requirements */
/*     //========================================================================== */
/*     virtual double calcValue(double x) const                        = 0; */
/*     virtual double calcFirstDerivative(double x) const              = 0; */
/*     virtual std::vector<double> calcMonotonicSegmentXValues() const = 0; */

/*     virtual double calcDomainMax() const; */
/*     virtual double calcDomainMin() const; */

/*     //========================================================================== */
/*     //              Curve Shape Derived */
/*     //========================================================================== */

/*     std::vector<OpenSim::CurveKnot> calcMonotonicSegmentKnots() const; */
/* }; */

//==============================================================================
//              Cubic Spline
//==============================================================================

class CubicSpline
{
public:
    using Coefficients = std::array<double, 4>;

    CubicSpline() = default;

    // Performs Hermite interpolation.
    CubicSpline(
        const CurveKnot& left,
        const CurveKnot& right,
        double& startIntegralValue);

    bool isMonotonic() const;

    double calcValue(double x) const;

    double calcInverseValue(double y, double eps = 1e-13, size_t maxIter = 20)
        const;

    double calcDerivative(double x, size_t order) const;

    ValueAndDerivative calcValueAndDerivative(double x) const;

    double calcIntegral(double x) const;

    CurveKnot calcKnot(double x) const;

    double x0 = SimTK::NaN;
    Coefficients coeff{SimTK::NaN, SimTK::NaN, SimTK::NaN, SimTK::NaN};
    double y0Integral = SimTK::NaN;
    double x1         = SimTK::NaN;
};

std::ostream& operator<<(std::ostream& os, const CubicSpline& spline);

//==============================================================================
//              Cubic Monotonic Spline
//==============================================================================

class CubicMonoSpline final : public CubicSpline
{
public:
    explicit CubicMonoSpline(CubicSpline spline);

    CubicMonoSpline(
        const CurveKnot& left,
        const CurveKnot& right,
        double y0Integral) :
        CubicMonoSpline(CubicSpline(left, right, y0Integral))
    {}

private:
    CubicMonoSpline() = default;
};

std::ostream& operator<<(std::ostream& os, const CubicMonoSpline& spline);

//==============================================================================
//                  SMOOTH SEGMENTED CUBIC MONO SPLINE
//==============================================================================

class SmoothSegmentedCubicMonoSpline
{
public:
    explicit SmoothSegmentedCubicMonoSpline(
        const C2ContinuousSegmentedCurve& curve,
        size_t maxNumSegments = 100,
        double accuracy       = 1e-3);

    SimTK::Vec2 getDomain() const;

    double calcValue(double x) const;

    ValueAndDerivative calcValueAndDerivative(double x) const;

    double calcInverseValue(double y) const;

    // Solves spline(x) = rhs(x) using newton iteratons starting from xEstimate
    // as initial guess.
    double solve(
        const std::function<ValueAndDerivative(double x)>& rhs,
        double xEstimate,
        size_t maxIter = 100,
        double eps     = 1e-13) const;

private:
    const CubicMonoSpline& findInverseSegment(double y) const;

    const CubicMonoSpline& findSegment(double x) const;

    std::vector<CubicMonoSpline> _splines{};
    bool _isInvertible;
};

} // namespace OpenSim

#endif
