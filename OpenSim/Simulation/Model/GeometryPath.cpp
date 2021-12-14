/* -------------------------------------------------------------------------- *
 *                         OpenSim:  GeometryPath.cpp                         *
 * -------------------------------------------------------------------------- *
 * The OpenSim API is a toolkit for musculoskeletal modeling and simulation.  *
 * See http://opensim.stanford.edu and the NOTICE file for more information.  *
 * OpenSim is developed at Stanford University and supported by the US        *
 * National Institutes of Health (U54 GM072970, R24 HD065690) and by DARPA    *
 * through the Warrior Web program.                                           *
 *                                                                            *
 * Copyright (c) 2005-2017 Stanford University and the Authors                *
 * Author(s): Peter Loan, Ajay Seth                                           *
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

//=============================================================================
// INCLUDES
//=============================================================================
#include "GeometryPath.h"
#include "ConditionalPathPoint.h"
#include "MovingPathPoint.h"
#include "PointForceDirection.h"
#include <OpenSim/Simulation/Wrap/PathWrap.h>
#include "Model.h"

//=============================================================================
// STATICS
//=============================================================================
using namespace std;
using namespace OpenSim;
using namespace SimTK;
using SimTK::Vec3;

//=============================================================================
// CONSTRUCTOR(S) AND DESTRUCTOR
//=============================================================================
//_____________________________________________________________________________
/*
 * Default constructor.
 */
GeometryPath::GeometryPath() :
    ModelComponent(),
    _preScaleLength(0.0)
{
    setAuthors("Peter Loan");
    constructProperties();
 }

//_____________________________________________________________________________
/*
* Perform set up functions after model has been deserialized or copied.
*
*/
void GeometryPath::extendFinalizeFromProperties()
{
    Super::extendFinalizeFromProperties();

    for (int i = 0; i < get_PathWrapSet().getSize(); ++i) {
        if (upd_PathWrapSet()[i].getName().empty()) {
            std::stringstream label;
            label << "pathwrap_" << i;
            upd_PathWrapSet()[i].setName(label.str());
        }
    }
}

void GeometryPath::extendConnectToModel(Model& aModel)
{
    Super::extendConnectToModel(aModel);

    OPENSIM_THROW_IF_FRMOBJ(get_PathPointSet().getSize() < 2,
        InvalidPropertyValue,
        getProperty_PathPointSet().getName(),
        "A valid path must be connected to a model by at least two PathPoints.")

    // Name the path points based on the current path
    // (i.e., the set of currently active points is numbered
    // 1, 2, 3, ...).
    namePathPoints(0);
}

//_____________________________________________________________________________
/*
 * Create the SimTK state, discrete and/or cache for this GeometryPath.
 */
 void GeometryPath::extendAddToSystem(SimTK::MultibodySystem& system) const 
{
    Super::extendAddToSystem(system);

    // Allocate cache entries to save the current length and speed(=d/dt length)
    // of the path in the cache. Length depends only on q's so will be valid
    // after Position stage, speed requires u's also so valid at Velocity stage.
    this->_lengthCV = addCacheVariable("length", 0.0, SimTK::Stage::Position);
    this->_speedCV = addCacheVariable("speed", 0.0, SimTK::Stage::Velocity);

    // Cache the set of points currently defining this path.
    this->_currentPathCV = addCacheVariable("current_path", Array<AbstractPathPoint*>{}, SimTK::Stage::Position);

    // We consider this cache entry valid any time after it has been created
    // and first marked valid, and we won't ever invalidate it.
    this->_colorCV = addCacheVariable("color", get_Appearance().get_color(), SimTK::Stage::Topology);
}

 void GeometryPath::extendInitStateFromProperties(SimTK::State& s) const
{
    Super::extendInitStateFromProperties(s);
    markCacheVariableValid(s, _colorCV);  // it is OK at its default value
}



//_____________________________________________________________________________
/*
 * Connect properties to local pointers.
 */
void GeometryPath::constructProperties()
{
    constructProperty_PathPointSet(PathPointSet());

    constructProperty_PathWrapSet(PathWrapSet());
    
    Appearance appearance;
    appearance.set_color(SimTK::Gray);
    constructProperty_Appearance(appearance);
}


//_____________________________________________________________________________
/*
 * Update the geometric representation of the path.
 * The resulting geometry is maintained at the VisibleObject layer.
 * This function should not be made public. It is called internally
 * by compute() only when the path has changed.
 * 
 */
void GeometryPath::updateGeometry(const SimTK::State& s) const
{
    // Check if the current path needs to recomputed.
    computePath(s);
}

void GeometryPath::setColor(const SimTK::State& s, const SimTK::Vec3& color) const
{
    setCacheVariableValue(s, _colorCV, color);
}

Vec3 GeometryPath::getColor(const SimTK::State& s) const
{
    return getCacheVariableValue(s, _colorCV);
}

void GeometryPath::setPreScaleLength( const SimTK::State& s, double length ) {
    _preScaleLength = length;
}
double GeometryPath::getPreScaleLength( const SimTK::State& s) const {
    return _preScaleLength;
}

//==============================================================================
// SCALING
//==============================================================================
void GeometryPath::
extendPreScale(const SimTK::State& s, const ScaleSet& scaleSet)
{
    Super::extendPreScale(s, scaleSet);
    setPreScaleLength(s, getLength(s));
}

void GeometryPath::
extendPostScale(const SimTK::State& s, const ScaleSet& scaleSet)
{
    Super::extendPostScale(s, scaleSet);
    computePath(s);
}

//_____________________________________________________________________________
// Override default implementation by object to intercept and fix the XML node
// underneath the model to match current version.
void GeometryPath::updateFromXMLNode(SimTK::Xml::Element& aNode, int versionNumber)
{
    if (versionNumber < XMLDocument::getLatestVersion()) {
        if (versionNumber < 30516) {
            // Create Appearance node under GeometryPath
            SimTK::Xml::Element appearanceElement("Appearance");
            aNode.appendNode(appearanceElement);
            SimTK::Xml::element_iterator visObjectIter = aNode.element_begin("VisibleObject");
            if (visObjectIter != aNode.element_end()) {
                SimTK::Xml::element_iterator oldPrefIter = visObjectIter->element_begin("display_preference");
                // old display_preference was used only for hide/show other options unsupported
                if (oldPrefIter != visObjectIter->element_end()) {
                    int oldPrefAsInt = 4;
                    oldPrefIter->getValueAs<int>(oldPrefAsInt);
                    if (oldPrefAsInt == 0) { // Hidden => Appearance/Visible
                        SimTK::Xml::Element visibleElement("visible");
                        visibleElement.setValueAs<bool>(false);
                        appearanceElement.insertNodeAfter(appearanceElement.element_end(), visibleElement);
                    }
                }
            }
            // If default_color existed, copy it to Appearance/color
            SimTK::Xml::element_iterator defaultColorIter = aNode.element_begin("default_color");
            if (defaultColorIter != aNode.element_end()) {
                // Move default_color to Appearance/color
                SimTK::Xml::Element colorElement("color");
                const SimTK::String& colorAsString = defaultColorIter->getValue();
                colorElement.setValue(colorAsString);
                appearanceElement.appendNode(colorElement);
            }
        }
    }
    // Call base class now assuming aNode has been corrected for current version
    Super::updateFromXMLNode(aNode, versionNumber);
}
