#ifndef OPENSIM_MUSCLECURVEPARAMS_H_
#define OPENSIM_MUSCLECURVEPARAMS_H_

#include "SimTKmath.h"
#include <array>

namespace OpenSim
{
class MuscleCurveControlPoint;

struct CurveControlPoint
{
    CurveControlPoint() = default;

    // Returns NaN if there is no such intercept, e.g. curviness = 0.
    // TODO take base as arg
    double calcTangentInterceptCoordinate(const CurveControlPoint& other) const;

    // Extrapolate y coordinate using the set derivative.
    CurveControlPoint calcExtrapolatedPoint(double xk) const;

    // Returns false if members x, y, and dydx are uninitialized (equal to NaN).
    explicit operator bool() const;

    double x    = SimTK::NaN;
    double y    = SimTK::NaN;
    double dydx = SimTK::NaN;
};

//==============================================================================
//                      Monotonic Control Points Storage
//==============================================================================
class MonoControlPoints // TODO this class is perhaps overengineering.
{
public:
    MonoControlPoints() = default;

    void appendChecked(CurveControlPoint pt);

    const std::vector<CurveControlPoint>& getPoints()
    {
        return _storage;
    }

private:
    std::vector<CurveControlPoint> _storage{};
};

//==============================================================================
//              User Facing Control Points for Curve Generation
//==============================================================================
// Control point with curviness parameter.
class MuscleCurveControlPoint : public CurveControlPoint
{
public:
    MuscleCurveControlPoint() = default;

    size_t calcCurvyPoints(
        const MuscleCurveControlPoint& other,
        std::vector<CurveControlPoint>& buffer) const;

    static constexpr double MAX_CURVINESS =
        0.99; // TODO how to limit the curviness.

    // Between 0 and MAX_CURVINESS.
    double curviness =
        SimTK::NaN; // TODO weird: last curviness of last point is invalid.
};

// MuscleCurveParams -> MuscleCurveControlPoints (x, y, dy/dx, curviness)
//
// MuscleCurveControlPoint (x, y, dy/dx, curviness) -> ControlPoints(x, y,
// dy/dx)
//
// MuscleCurveControlPoint (x, y, dy/dx) -> CubicMonoSpline(x0, yInt, coeffs,
// x1) *** PAIR TO PAIR ***
//
// SmoothSegmentedCubicMonoSpline
//
// SmoothSegmentedSplineFunction
//

class SmoothSegmentedCubicMonoSplineData;

//==============================================================================
//              Cubic Spline
//==============================================================================
class CubicSpline
{
public:
    using Coefficients = std::array<double, 4>;

    CubicSpline() = default;

    bool isMonotonic() const;

    double calcValue(double x) const;

    double calcInverseValue(double y, double eps = 1e-13, size_t maxIter = 20)
        const;

    double calcDerivative(double x, size_t order) const;

    double calcIntegral(double x) const;

    double calcEndIntegral() const;

    double x0 = SimTK::NaN;
    Coefficients coeff{SimTK::NaN, SimTK::NaN, SimTK::NaN, SimTK::NaN};
    double y0Integral = SimTK::NaN;
    double x1         = SimTK::NaN;
};

//==============================================================================
//              Cubic Monotonic Spline
//==============================================================================
class CubicMonoSpline final : public CubicSpline
{
public:
    explicit CubicMonoSpline(CubicSpline spline);

private:
    CubicMonoSpline() = default;

    static constexpr size_t dataSize()
    {
        static_assert(
            offsetof(CubicMonoSpline, x0) == 0,
            "x0 must be the first data member");
        static_assert(
            offsetof(CubicMonoSpline, x1) == sizeof(CubicMonoSpline) - 8,
            "x1 must be the last data member");
        return sizeof(CubicMonoSpline) / 8;
    }

    friend SmoothSegmentedCubicMonoSplineData;
};

//==============================================================================
//                  Smooth Segmented Cubic Mono Spline Storage
//==============================================================================

// C2 continuous segmented cubic monotonic spline storage.
class SmoothSegmentedCubicMonoSplineData
{
public:
    SmoothSegmentedCubicMonoSplineData() = default;

    SmoothSegmentedCubicMonoSplineData(
        const std::vector<CubicMonoSpline>& splines);

    // Checks if spline segments are C2 continuous.
    void appendChecked(CubicMonoSpline spline);

    size_t size() const;

    const CubicMonoSpline& at(size_t index) const;

private:
    std::vector<double> _data;
};

class SmoothSegmentedCubicMonoSpline
{
public:
    SmoothSegmentedCubicMonoSpline(
        const std::vector<CubicMonoSpline>& splines,
        CurveControlPoint pStart,
        CurveControlPoint pEnd);

    SmoothSegmentedCubicMonoSpline(const std::vector<CurveControlPoint>& pts);

    SmoothSegmentedCubicMonoSpline(
        const std::vector<MuscleCurveControlPoint>& pts);

    SimTK::Vec2 getDomain() const;

    void setExtrapolationBeyondDomain(bool allowExtrapolation);

private:
    SmoothSegmentedCubicMonoSplineData _splines;
    CurveControlPoint _pStart;
    CurveControlPoint _pEnd;
    bool _extrapolateBeyondDomain = false;
};

} // namespace OpenSim

#endif
