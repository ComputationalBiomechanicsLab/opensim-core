#ifndef OPENSIM_MUSCLECURVEPARAMS_H_
#define OPENSIM_MUSCLECURVEPARAMS_H_

#include "SimTKmath.h"
#include <array>
#include <cstddef>
#include <functional>
#include <memory>
#include <utility>

namespace OpenSim
{

// Forward declaration.
class MuscleCurveControlPoint;
class CurveControlPoint;

// Procedure:
//
// 1. Convert Curved control points to CurveControlPoints
// 2. Convert CurveControlPoint into CurveSegmentShape
// 3. Use CurveSegmentShape to generate CurveKnots and CurveControlPoints
// 4. Create C1 interpolating spline of segment:
//      4.1 Create 4 curve knots.
//      4.2 Check if knots are hermite-monotonic
//      4.3 split all segments that are not herminte monotonic
//      4.4 Keep going until done.
// 5. Create C2 interpolating spline of segment:
//      5.1 start with grid from C1 spline
//      5.2 fit the C2 spline
//      5.3 split all segments that are not monotonic
//      5.4 Keep going until done
// 6. Collect spline segments in function
//      6.1 Check continuity at all knots
//      6.2 Check Monotonicity of segments
//

//==============================================================================
//                  CURVE POINT
//==============================================================================

class CurvePoint
{
public:
    CurvePoint() = default;

    CurvePoint(double xCoord, double yCoord) : x(xCoord), y(yCoord)
    {}

    double calcSecantLine(const CurvePoint& other) const;

    CurvePoint calcInterpolated(const CurvePoint& other, double u) const;

    double x = SimTK::NaN;
    double y = SimTK::NaN;
};

std::ostream& operator<<(std::ostream& os, const CurvePoint& pt);

//==============================================================================
//                  CURVE KNOT
//==============================================================================

class CurveKnot : public CurvePoint
{
public:
    CurveKnot() = default;

    CurveKnot(double xCoord, double yCoord, double derivative) :
        CurvePoint(xCoord, yCoord), dydx(derivative)
    {}

    CurveKnot(CurvePoint point, double derivative) :
        CurvePoint(point), dydx(derivative)
    {}

    double dydx = SimTK::NaN;
};

std::ostream& operator<<(std::ostream& os, const CurveKnot& knot);

//==============================================================================
//              CURVY CONTROL POINT
//==============================================================================

// Control point with curviness parameter.
class MuscleCurveControlPoint final : public CurveKnot
{
public:
    MuscleCurveControlPoint() = default;

    // For checking if additional knots should be created.
    bool isCurvy() const;

    // TODO weird: last curviness of last point is invalid.
    double curviness = SimTK::NaN;
};

std::ostream& operator<<(
    std::ostream& os,
    const MuscleCurveControlPoint& ctrlPt);

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

// Forward declaration.
class SmoothSegmentedCubicMonoSplineData;

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

    static constexpr size_t dataSize()
    {
        static_assert(
            offsetof(CubicMonoSpline, x0) == 0,
            "x0 must be the first data member");
        static_assert(
            offsetof(CubicMonoSpline, x1) ==
                sizeof(CubicMonoSpline) - sizeof(double),
            "x1 must be the last data member");
        return sizeof(CubicMonoSpline) / sizeof(double);
    }

    friend SmoothSegmentedCubicMonoSplineData;
};

std::ostream& operator<<(std::ostream& os, const CubicMonoSpline& spline);

//==============================================================================
//              QuadraticBezierCurve
//==============================================================================
class QuadraticBezierCurve
{
public:
    QuadraticBezierCurve(const CurveKnot& left, const CurveKnot& right);

    CurvePoint calcPoint(double u) const;

    const CurveKnot& startKnot() const;
    const CurveKnot& endKnot() const;

private:
    CurveKnot _start;
    CurveKnot _end;
    CubicSpline _x;
    CubicSpline _y;
};

//==============================================================================
//              Curve Shape
//==============================================================================
// If the user gives just a few control points of the curve, it is
// anyones guess what the general shape of the curve is.
// The CurveShape is constructed from a few curve knots, and can generate a
// curve point at any x between the knots.
//
// Input knots must be monotonic.
class CurveShape final
{
public:
    explicit CurveShape(std::vector<MuscleCurveControlPoint> ctrlPts);

    explicit CurveShape(std::vector<CurveKnot> knots);

    const std::vector<QuadraticBezierCurve>& getSegments() const
    {
        return _segments;
    }

private:
    std::vector<QuadraticBezierCurve> _segments;
};

//==============================================================================
//                  SPLINE STORAGE
//==============================================================================

class SmoothSegmentedCubicMonoSpline;

// C2 continuous segmented cubic monotonic spline storage.
class SmoothSegmentedCubicMonoSplineData
{
    SmoothSegmentedCubicMonoSplineData() = default;

    explicit SmoothSegmentedCubicMonoSplineData(
        const CurveShape& shape,
        size_t maxNumSegments);

    explicit SmoothSegmentedCubicMonoSplineData(
        const std::vector<CubicMonoSpline>& splines);

    // Checks if spline segments are C2 continuous.
    void appendChecked(CubicMonoSpline spline);

    size_t size() const;

    const CubicMonoSpline& at(size_t index) const;

    SimTK::Vec2 getDomain() const;

    std::vector<double> _data;

    friend SmoothSegmentedCubicMonoSpline;
};

//==============================================================================
//                  SMOOTH SEGMENTED CUBIC MONO SPLINE
//==============================================================================

// inserted gradient: dy/dx = 2 / ( DX0/DY0 + DX1/DY1), or zero if sign changed.

// RULES:
// - and control point with gradient set must be followed by one without
// gradient set
// - any missing griadient will receive inserted gradient
// - control points must not violate mono-spline condition after all gradients
// are set
// - fit splines
// - check if splines are monotonic
// - check if splines are C1 continuous
//
// Resampling version:
// - Sample a curve and using control points.
// - Insert missing derivative values
// - assert Monotonicity
class SmoothSegmentedCubicMonoSpline
{
public:
    explicit SmoothSegmentedCubicMonoSpline(
        const CurveShape& shape,
        size_t maxNumSegments);

    explicit SmoothSegmentedCubicMonoSpline(
        const std::vector<MuscleCurveControlPoint>& pts);

    SimTK::Vec2 getDomain() const;

    double calcValue(double x) const;

private:
    const CubicSpline& findInverseSegment(double y) const;

    const CubicMonoSpline& findSegment(double x) const;

    SmoothSegmentedCubicMonoSplineData _splines;
};

} // namespace OpenSim

#endif
