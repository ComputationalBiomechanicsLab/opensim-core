#include "OpenSim/Common/SmoothSegmentedCubicMonoSpline.h"
#include "OpenSim/Common/SmoothSegmentedFunction.h"
#include <iostream>
#include <memory>

#include <OpenSim/OpenSim.h>

void writeDataToCSV(
    const OpenSim::SmoothSegmentedFunction& curve,
    size_t nSamples = 100)
{
    static constexpr size_t SPLINE_SHAPE_ID = 0;
    static constexpr size_t SPLINE_PLOT_ID  = 1;
    static constexpr size_t CTRL_PTS_ID     = 2;

    const OpenSim::SmoothSegmentedCubicMonoSpline spline(curve);
    std::cout << "spline.x0 = " << spline.getDomain() << "\n";

    std::string filename = curve.getName() + ".csv";

    std::ofstream outputFile(filename);

    if (!outputFile.is_open()) {
        std::cerr << "Error opening the file: " << filename << std::endl;
        return;
    }

    // Write header
    outputFile << "Index,X,Y\n";

    // Write data

    // Write data points for each shape segment, and the approximating spline.
    std::vector<double> knotCoords = curve.calcMonotonicSegmentXValues();
    for (size_t i = 0; i + 1 < knotCoords.size(); ++i) {
        for (size_t j = 0; j <= nSamples; j++) {
            const double factor =
                static_cast<double>(j) / static_cast<double>(nSamples);
            const double dx        = knotCoords.at(i+1) - knotCoords.at(i);
            const double x        = knotCoords.at(i) + dx * factor;
            const double ySpline  = spline.calcValue(x);
            const double yCurve  = curve.calcValue(x);

            outputFile << SPLINE_PLOT_ID << "," << x << "," << ySpline << "\n";
            outputFile << SPLINE_SHAPE_ID << "," << x << "," << yCurve
                       << "\n";
        }
    }

    // Write the original control points.
    for (size_t i = 0; i < knotCoords.size(); ++i) {
            const double x        = knotCoords.at(i);
            std::cout << "samplling"<< "\n"
                << "    x = " << x << "\n"
                << "    y = " << spline.calcValue(x) << "\n"
                << "    c = " << curve.calcValue(x) << "\n"
                << "    d = " << spline.getDomain() << "\n"
                << "    d = " << curve.getCurveDomain() << "\n";
        outputFile << CTRL_PTS_ID << "," << knotCoords.at(i) << "," << curve.calcValue(knotCoords.at(i)) << "\n";
    }

    outputFile.close();
}

/* int convertAndWriteToCsv( */
/*     const OpenSim::SmoothSegmentedFunction& curve, */
/*     const OpenSim::SmoothSegmentedCubicMonoSpline& spline, */
/*     size_t steps      = 20) */
/* { */
/*     { */
/*         auto eval = [&](double x) -> double { */
/*             double yBezier = f->calcValue(x); */
/*             double ySpline = spline.calcValue(x); */
/*             double error   = std::abs(ySpline - yBezier); */

/*             std::cout << "x = " << x << std::endl; */
/*             std::cout << "e = " << error << std::endl; */
/*             std::cout << "yBezier = " << yBezier << std::endl; */
/*             std::cout << "ySpline = " << ySpline << std::endl; */

/*             return error; */
/*         }; */
/*         double x0    = 1.; */
/*         double x1    = 1.01453; */
/*         double error = 0.; */
/*         for (size_t i = 0; i <= steps; i++) { */
/*             const double c = */
/*                 static_cast<double>(i) / static_cast<double>(steps); */
/*             const double x = x0 + (x1 - x0) * c; */
/*             error          = std::max(error, eval(x)); */
/*         } */

/*         std::cout << "max error = " << error << std::endl; */
/*     } */

/*     { */
/*         std::string filename = f->getName() + ".csv"; */
/*         writeDataToCSV(spline.getSplines(), *f, filename, checkInverse); */
/*     } */

/*     { */
/*         std::string filenameInt = f->getName() + "_integral.csv"; */
/*         writeDataIntegralToCSV(spline.getSplines(), filenameInt); */
/*     } */

/*     return 0; */
/* } */

int main()
{
    if (true)
    {
        OpenSim::SmoothSegmentedFunction* f = OpenSim::
            SmoothSegmentedFunctionFactory::createTendonForceLengthCurve(
                0.049,
                28.1,
                0.67,
                0.5,
                false,
                "tendoncurve");
        OpenSim::SmoothSegmentedFunction* spline =
OpenSim::SmoothSegmentedFunctionFactory::
            createTendonForceLengthCurve(0.049, 28.1, 0.67, 0.5, false,
"tendoncurve");
        writeDataToCSV(*f);
    }

    if (true)
    {
        OpenSim::SmoothSegmentedFunction* f = OpenSim::
            SmoothSegmentedFunctionFactory::createFiberForceVelocityCurve(

                1.4,
                0.,
                0.25,
                5.,
                0.,
                0.15,
                0.6,
                0.9,
                false,
                "fiberVelocityCurve");
        writeDataToCSV(*f);
    }

    if (false)
    {
        OpenSim::SmoothSegmentedFunction* f = OpenSim::
            SmoothSegmentedFunctionFactory::createFiberActiveForceLengthCurve(
                0.44,
                0.73,
                1.0,
                1.8123,
                0.1,
                0.8616,
                1.0,
                false,
                "ActiveForceLengthCurve");
        writeDataToCSV(*f);
    }

    if (true)
    {
        OpenSim::SmoothSegmentedFunction* f = OpenSim::
            SmoothSegmentedFunctionFactory::createFiberForceLengthCurve(
                0.0,
                0.7,
                0.2,
                2.86,
                0.75,
                false,
                "fiberForceLengthCurve");
        writeDataToCSV(*f);
    }

    return 0;
}
