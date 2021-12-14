/* -------------------------------------------------------------------------- *
 *                     OpenSim:  fbp-compiler.cpp                             *
 * -------------------------------------------------------------------------- *
 * ScapulothoracicJoint is offered as an addition to the OpenSim API          *
 * See http://opensim.stanford.edu and the NOTICE file for more information.  *
 *                                                                            *
 * OpenSim is developed at Stanford University and is supported by:           *
 *                                                                            *
 * - The National Institutes of Health (U54 GM072970, R24 HD065690)           *
 * - DARPA, through the Warrior Web program                                   *
 * - The Chan Zuckerberg Initiative (CZI 2020-218896)                         *
 *                                                                            *
 * Copyright (c) 2005-2021 Stanford University, TU Delft, and the Authors     *
 * Author(s): Joris Verhagen                                                  *
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

#include <OpenSim/Common/Function.h>
#include <OpenSim/Common/Logger.h>
#include <OpenSim/Tools/FunctionBasedPathConversionTool.h>

#include <cassert>
#include <cstring>
#include <iostream>
#include <exception>

static const char usage[] = "FunctionBasedPathConversion [--help] INPUT_MODEL OUTPUT_NAME";
static const char help_text[] =
R"(
ARGS

    --help           print this help

DESCRIPTION

    This compiler tries to convert an .osim model input file (INPUT_MODEL), which
    may contain `<PointBasedPath>` or `<GeometryPath>`s, into a functionally
    equivalent model output file (${OUTPUT_NAME}.osim) with associated curve
    fitting data (${OUTPUT_NAME}_DATA_${i}.txt). The output model is the same,
    but with each `<PointBasedPath>` or `<GeometryPath>` potentially replaced by a
    `<FunctionBasedPath>`.

    The phrases "functionally equivalent" and "potentially replaced" are
    important. The implementation tries to parameterize each input (point-based)
    path with the coordinates in the input model. The implementation will
    permute through the coordinates, trying to figure out which of them
    ultimately affect the input path. If the implementation determines the
    coordintate-to-path relationship, it will then try and fit the relationship
    to a multidimensional curve. If that curve appears to produce the same
    same outputs as the input (point-based) path, it will replace that path
    with the function-/curve-based equivalent. This replacement is only
    "functionally equivalent" if the model is then used in ways that ensure its
    coordinates roughly stay in the same value range as was used during curve
    fitting. Input paths are only "potentially replaced" because the implementation
    may not be able to parameterize all input paths.

    The primary reason to replace point-based (input) paths with function-based
    (output) equivalents is runtime performance. Point-based path computations
    can be expensive--especially if those paths are "wrapped" around other geometry--
    whereas function-based paths are (effectively) lookups into precomputed curves.
    Those curves are usually smooth, so the resulting path outpus (length, speed,
    etc.) are usually smoother when integrated over multiple states. This smoothness
    can accelerate error-controlled integration schemes (OpenSim's default).

    If performance isn't an issue for you, you should probably keep using point-based
    paths.

    IMPLEMENTATION STEPS (high-level):

      - Reads INPUT_MODEL
      - Finds `GeometryPath`/`PointBasedPath`s in INPUT_MODEL
      - Parameterizes each input path against each of INPUT_MODEL's coordinates
        to produce an n-dimensional Bezier curve fit of those paths
      - Saves the curve data to ${OUTPUT_NAME}_DATA_${i}.txt, where `i` is an
        arbirary ID that links the `FunctionBasedPath` in the output .osim file
        to the Bezier fit's data
      - Updates the source model to contain `FunctionBasedPath`s
      - Writes the updated model to `${OUTPUT_NAME}.osim`

EXAMPLE USAGE

    FunctionBasedPathModelTool RajagopalModel.osim RajagopalModel_Precomputed
)";

int main(int argc, char **argv) {
    // skip exe name
    --argc;
    ++argv;

    // set during parsing
    char const* sourceModelPath = nullptr;
    char const* outputName = nullptr;

    // parse CLI args
    int nUnnamed = 0;
    while (argc) {
        const char* arg = *argv;

        if (arg[0] != '-') {  // handle unnamed args
            switch (nUnnamed) {
            case 0:
                // INPUT_MODEL
                sourceModelPath = arg;
                break;
            case 1:
                // OUTPUT_NAME
                outputName = arg;
                break;
            default:
                std::cerr << "FunctionBasedPathModelTool: error: too many arguments: should only supply 'MODEL' and 'OUTPUT_NAME'\n\nUSAGE: "
                          << usage
                          << std::endl;
                return -1;
            }

            ++nUnnamed;
            --argc;
            ++argv;
            continue;
        } else {  // else: handle named args
            if (std::strcmp(arg, "--help") == 0) {
                std::cout << "usage: " << usage << '\n' << help_text << std::endl;
                return 0;
            } else {
                std::cerr << "FunctionBasedPathModelTool: error: unknown argument '"
                          << arg
                          << "': see --help for usage info";
                return -1;
            }
        }
    }

    // ensure inputs are correct
    if (nUnnamed != 2) {
        std::cerr << "FunctionBasedPathModelTool: error: too few arguments supplied\n\n"
                  << usage
                  << std::endl;
        return -1;
    }

    assert(argc == 0);
    assert(sourceModelPath != nullptr);
    assert(outputName != nullptr);

    // run the tool
    try {
        OpenSim::FunctionBasedPathConversionTool tool(sourceModelPath,outputName);
        return tool.run() ? 0 : 1;
    } catch(const std::exception& ex) {
        OpenSim::log_error("Exception in ID: {}", ex.what());
        return -1;
    }
}
