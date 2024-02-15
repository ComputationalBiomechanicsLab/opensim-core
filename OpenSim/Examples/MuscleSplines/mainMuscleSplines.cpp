#include "OpenSim/Common/SmoothSegmentedCubicMonoSpline.h"
#include <iostream>
#include <memory>

#include <OpenSim/OpenSim.h>

/* using CtrlPoint = OpenSim::MuscleCurveControlPoint; */

/* void writeDataToCSV( */
/*     const std::vector<OpenSim::MuscleCurveControlPoint>& ctrlPts, */
/*     const OpenSim::CurveShape& shapeCurve, */
/*     const OpenSim::SmoothSegmentedCubicMonoSpline& splineCurve, */
/*     const std::string& filename, */
/*     size_t nSamples = 1000) */
/* { */
/*     static constexpr size_t SPLINE_SHAPE_ID = 0; */
/*     static constexpr size_t SPLINE_PLOT_ID  = 1; */
/*     static constexpr size_t CTRL_PTS_ID     = 2; */

/*     std::ofstream outputFile(filename); */

/*     if (!outputFile.is_open()) { */
/*         std::cerr << "Error opening the file: " << filename << std::endl; */
/*         return; */
/*     } */

/*     // Write header */
/*     outputFile << "Index,X,Y\n"; */

/*     // Write data */

/*     // Write data points for each shape segment, and the approximating spline. */
/*     for (const OpenSim::QuadraticBezierCurve& shapeSegment : */
/*          shapeCurve.getSegments()) { */

/*         for (size_t i = 0; i <= nSamples; i++) { */
/*             const double factor = */
/*                 static_cast<double>(i) / static_cast<double>(nSamples); */
/*             const auto shapePoint = shapeSegment.calcPoint(factor); */
/*             const double x        = shapePoint.x; */
/*             const double ySpline  = splineCurve.calcValue(x); */
/*             outputFile << SPLINE_PLOT_ID << "," << x << "," << ySpline << "\n"; */
/*             outputFile << SPLINE_SHAPE_ID << "," << x << "," << shapePoint.y */
/*                        << "\n"; */
/*         } */
/*     } */

/*     // Write the original control points. */
/*     for (const OpenSim::MuscleCurveControlPoint& p : ctrlPts) { */
/*         outputFile << CTRL_PTS_ID << "," << p.x << "," << p.y << "\n"; */
/*     } */

/*     outputFile.close(); */
/* } */

int main()
{

    /* std::vector<CtrlPoint> ctrlPoints; */

    /* { */
    /*     CtrlPoint point; */
    /*     point.x    = 0.; */
    /*     point.y    = 0.; */
    /*     point.dydx = 0.; */
    /*     ctrlPoints.push_back(point); */
    /* } */

    /* { */
    /*     CtrlPoint point; */
    /*     point.x    = 1.; */
    /*     point.y    = 1.; */
    /*     point.dydx = 5.; */
    /*     point.curviness = 0.1; */
    /*     ctrlPoints.push_back(point); */
    /* } */

    /* { */
    /*     CtrlPoint point; */
    /*     point.x    = 3.5; */
    /*     point.y    = 2.2; */
    /*     point.dydx = 0.1; */
    /*     point.curviness = 0.9; */
    /*     ctrlPoints.push_back(point); */
    /* } */

    /* /1* { *1/ */
    /* /1* CtrlPoint point; *1/ */
    /* /1*     point.x = 2.; *1/ */
    /* /1*     point.y = 1.5; *1/ */
    /* /1*     point.dydx = 0.; *1/ */
    /* /1*     ctrlPoints.push_back(point); *1/ */
    /* /1* } *1/ */

    /* /1* { *1/ */
    /* /1* CtrlPoint point; *1/ */
    /* /1*     point.x = 5.; *1/ */
    /* /1*     point.y = 0.1; *1/ */
    /* /1*     point.dydx = 0.; *1/ */
    /* /1*     ctrlPoints.push_back(point); *1/ */
    /* /1* } *1/ */

    /* for (const auto& p : ctrlPoints) { */
    /*     std::cout << p << std::endl; */
    /* } */

    /* OpenSim::CurveShape cShape(ctrlPoints); */

    /* OpenSim::SmoothSegmentedCubicMonoSpline sCurve(cShape); */

    /* std::string filename = "customSpline.csv"; */
    /* writeDataToCSV(ctrlPoints, cShape, sCurve, filename, 1000); */

    return 0;
}
