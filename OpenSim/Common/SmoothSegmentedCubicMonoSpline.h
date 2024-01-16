#ifndef OPENSIM_MUSCLECURVEPARAMS_H_
#define OPENSIM_MUSCLECURVEPARAMS_H_

#include "SimTKmath.h"
#include <array>

namespace OpenSim
{

// Forward declaration.
class MuscleCurveControlPoint;

//==============================================================================
//                  CONTROL POINT
//==============================================================================

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

std::ostream& operator<<(std::ostream& os, const CurveControlPoint& ctrlPt);

//==============================================================================
//              CURVY CONTROL POINT
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

std::ostream& operator<<(std::ostream& os, const MuscleCurveControlPoint& ctrlPt);

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

std::ostream& operator<<(std::ostream& os, const CubicSpline& spline);

//==============================================================================
//              Cubic Monotonic Spline
//==============================================================================

// Forward declaration.
class SmoothSegmentedCubicMonoSplineData;

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
//                  SPLINE STORAGE
//==============================================================================

// C2 continuous segmented cubic monotonic spline storage.
class SmoothSegmentedCubicMonoSplineData
{
public:
    SmoothSegmentedCubicMonoSplineData() = default;

    explicit SmoothSegmentedCubicMonoSplineData(
        const std::vector<CubicMonoSpline>& splines);

    // Checks if spline segments are C2 continuous.
    void appendChecked(CubicMonoSpline spline);

    size_t size() const;

    const CubicMonoSpline& at(size_t index) const;

    const CubicMonoSpline& findSegment(double x) const;

private:
    std::vector<double> _data;
};

//==============================================================================
//                  SMOOTH SEGMENTED CUBIC MONO SPLINE
//==============================================================================

class SmoothSegmentedCubicMonoSpline
{
public:
    explicit SmoothSegmentedCubicMonoSpline(
        std::vector<CubicMonoSpline>&& splines,
        CurveControlPoint pStart,
        CurveControlPoint pEnd);

    explicit SmoothSegmentedCubicMonoSpline(std::vector<CurveControlPoint>&& pts);

    explicit SmoothSegmentedCubicMonoSpline(
        const std::vector<MuscleCurveControlPoint>& pts);

    SimTK::Vec2 getDomain() const;

    double calcValue(double x) const;

    const CubicSpline& findInverseSegment(double y) const;

private:
    SmoothSegmentedCubicMonoSplineData _splines;
    CurveControlPoint _pStart;
    CurveControlPoint _pEnd;
    bool _extrapolateBeyondDomain = false;
};

} // namespace OpenSim

#endif
