#ifndef FUNCTIONBASEDPATHDISCRETIZATION_H
#define FUNCTIONBASEDPATHDISCRETIZATION_H

/* -------------------------------------------------------------------------- *
 *              OpenSim: FunctionBasedPathDiscretization.h                    *
 * -------------------------------------------------------------------------- *
 * The OpenSim API is a toolkit for musculoskeletal modeling and simulation.  *
 * See http://opensim.stanford.edu and the NOTICE file for more information.  *
 * OpenSim is developed at Stanford University and supported by the US        *
 * National Institutes of Health (U54 GM072970, R24 HD065690) and by DARPA    *
 * through the Warrior Web program.                                           *
 *                                                                            *
 * Copyright (c) 2005-2021 TU Delft and the Authors                           *
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

#include <OpenSim/Common/Component.h>
#include <OpenSim/Simulation/osimSimulationDLL.h>

#include <string>

namespace OpenSim {
    class OSIMSIMULATION_API FunctionBasedPathDiscretization : public Component {
        OpenSim_DECLARE_CONCRETE_OBJECT(FunctionBasedPathDiscretization, Component);

    public:
        OpenSim_DECLARE_PROPERTY(coordinate_abspath, std::string, "The absolute path, in the model, to the OpenSim::Coordinate that this discretization was produced from");
        OpenSim_DECLARE_PROPERTY(x_begin, double, "The lowest OpenSim::Coordinate value that was used for discretization");
        OpenSim_DECLARE_PROPERTY(x_end, double, "The highest OpenSim:::Coordinate value that was used for discretization");
        OpenSim_DECLARE_PROPERTY(num_points, int, "The number of evenly-spaced OpenSim::Coordinate values between [x_begin, x_end] (inclusive) that were used for discretization of the OpenSim::Coordinate. E.g. [x_begin, 1*(x_begin+((x_end-x_begin)/3)), 2*(x_begin+((x_end-x_begin)/3)), x_end]");

        FunctionBasedPathDiscretization() {
            constructProperty_coordinate_abspath("");
            constructProperty_x_begin(0.0);
            constructProperty_x_end(0.0);
            constructProperty_num_points(0);
        }
    };
}

#endif // FUNCTIONBASEDPATHDISCRETIZATION_H
