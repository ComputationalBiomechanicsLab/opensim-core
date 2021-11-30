#ifndef OPENSIM_GEOMETRY_PATH_H_
#define OPENSIM_GEOMETRY_PATH_H_

/* -------------------------------------------------------------------------- *
 *                          OpenSim:  GeometryPath.h                          *
 * -------------------------------------------------------------------------- *
 * The OpenSim API is a toolkit for musculoskeletal modeling and simulation.  *
 * See http://opensim.stanford.edu and the NOTICE file for more information.  *
 * OpenSim is developed at Stanford University and supported by the US        *
 * National Institutes of Health (U54 GM072970, R24 HD065690) and by DARPA    *
 * through the Warrior Web program.                                           *
 *                                                                            *
 * Copyright (c) 2005-2017 Stanford University and the Authors                *
 * Author(s): Peter Loan                                                      *
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

#include <OpenSim/Simulation/osimSimulationDLL.h>
#include <OpenSim/Simulation/Model/Appearance.h>
#include <OpenSim/Simulation/Model/ModelComponent.h>
#include <OpenSim/Simulation/Model/PathPointSet.h>
#include <OpenSim/Simulation/Wrap/PathWrapSet.h>

#ifdef SWIG
    #ifdef OSIMSIMULATION_API
        #undef OSIMSIMULATION_API
        #define OSIMSIMULATION_API
    #endif
#endif

namespace OpenSim {

class AbstractPathPoint;
class Coordinate;
class PhysicalFrame;
class PointForceDirection;
class ScaleSet;
class WrapResult;
class WrapObject;

/**
  * A base class that represents a path that has a computable length and
  * lengthening speed.
  *
  * This class is typically used in places where the model needs to simulate
  * the changes in a path over time. For example, in `OpenSim::Muscle`s,
  * `OpenSim::Ligament`s, etc.
  *
  * This class *only* defines a length and lenghtning speed. Do not assume that
  * an `OpenSim::GeometryPath` is a straight line between two points, or assume
  * that it is many straight lines between `n` points. The derived implementation
  * may define a path using points, or it may define a path using a curve fit.
  * It may also define a path as `length == 17.3f && lengtheningSpeed == 5.0f`.
  * All of those definitions are *logically* valid - even if they aren't
  * *functionally* valid.
  */
class OSIMSIMULATION_API GeometryPath : public ModelComponent {
    OpenSim_DECLARE_ABSTRACT_OBJECT(GeometryPath, ModelComponent);

//=============================================================================
// OUTPUTS
//=============================================================================
    OpenSim_DECLARE_OUTPUT(length,
            double,
            getLength,
            SimTK::Stage::Position);

    OpenSim_DECLARE_OUTPUT(lengthening_speed,
            double,
            getLengtheningSpeed,
            SimTK::Stage::Velocity);

//=============================================================================
// DATA
//=============================================================================
public:
    OpenSim_DECLARE_UNNAMED_PROPERTY(Appearance,
        "Default appearance attributes for this GeometryPath");

//=============================================================================
// METHODS
//=============================================================================
public:
    GeometryPath();

    virtual const PathPointSet& getPathPointSet() const = 0;
    virtual PathPointSet& updPathPointSet() = 0;

    /**
     * Add a new path point, with default location, to the path.
     *
     * @param aIndex The position in the pathPointSet to put the new point in.
     * @param frame The frame to attach the point to.
     * @return Pointer to the newly created path point.
     */
    virtual AbstractPathPoint* addPathPoint(
            const SimTK::State& s,
            int index,
            const PhysicalFrame& frame) = 0;
    virtual AbstractPathPoint* appendNewPathPoint(
            const std::string& proposedName,
            const PhysicalFrame& frame,
            const SimTK::Vec3& locationOnFrame) = 0;
    /**
     * Returns true if a path point can be deleted. All paths must have at
     * least two active path points to define the path.
     *
     * @param aIndex The index of the point to delete.
     * @return Whether or not the point can be deleted.
     */
    virtual bool canDeletePathPoint(int index) = 0;
    /**
     * Delete a path point.
     *
     * @param aIndex The index of the point to delete.
     * @return Whether or not the point was deleted.
     */
    virtual bool deletePathPoint(const SimTK::State& s, int index) = 0;
    /**
     * Replace a path point in the set with another point. The new one is made a
     * member of all the same groups as the old one, and is inserted in the same
     * place the old one occupied.
     *
     *  @param aOldPathPoint Path point to remove.
     *  @param aNewPathPoint Path point to add.
     */
    virtual bool replacePathPoint(
            const SimTK::State& s,
            AbstractPathPoint* oldPathPoint,
            AbstractPathPoint* newPathPoint) = 0;
    /**
     * Get the current path of the path.
     *
     * @return The array of currently active path points.
     */
    virtual const Array<AbstractPathPoint*>& getCurrentPath(const SimTK::State& s) const = 0;


    virtual const PathWrapSet& getWrapSet() const = 0;
    virtual PathWrapSet& updWrapSet() = 0;
    /**
     * Create a new wrap instance and add it to the set.
     *
     * @param aWrapObject The wrap object to use in the new wrap instance.
     */
    virtual void addPathWrap(WrapObject& aWrapObject) = 0;
    /**
     * Move a wrap instance up in the list. Changing the order of wrap instances for
     * a path may affect how the path wraps over the wrap objects.
     *
     * @param aIndex The index of the wrap instance to move up.
     */
    virtual void moveUpPathWrap(const SimTK::State& s, int index) = 0;
    /**
     * Move a wrap instance down in the list. Changing the order of wrap instances
     * for a path may affect how the path wraps over the wrap objects.
     *
     * @param aIndex The index of the wrap instance to move down.
     */
    virtual void moveDownPathWrap(const SimTK::State& s, int index) = 0;
    /**
     * Delete a wrap instance.
     *
     * @param aIndex The index of the wrap instance to delete.
     */
    virtual void deletePathWrap(const SimTK::State& s, int index) = 0;


    /**
     * Get the default color of the path.
     *
     * Returns the color that will be used to initialize the color cache
     * at the next extendAddToSystem() call. The actual color used to draw the
     * path will be taken from the cache variable, which may have changed.
     */
    virtual const SimTK::Vec3& getDefaultColor() const = 0;

    /**
     * Set the default color of the path.
     *
     * If you call this prior to extendAddToSystem() it will be used to initialize
     * the color cache variable. Otherwise %GeometryPath will choose its own
     * default, which varies depending on owner.
     */
    virtual void setDefaultColor(const SimTK::Vec3& color) = 0;

    /**
     * Get the current color of the path.
     *
     * Internally, gets the value of the color cache entry owned by this
     * %GeometryPath object in the given state. You can access this value any
     * time after the state is initialized, at which point it will have been
     * set to the default color value specified in a call to setDefaultColor()
     * earlier, or it will have the default color value chosen by %GeometryPath.
     *
     * @see setDefaultColor()
     */
    virtual SimTK::Vec3 getColor(const SimTK::State& s) const = 0;

    /**
     * Set the current color of the path.
     *
     * Internally, sets the value of the color cache variable owned by this
     * %GeometryPath object, in the cache of the given state. The value of this
     * variable is used as the color when the path is drawn, which occurs with
     * the state realized to Stage::Dynamics. Therefore, you must call this method
     * during realizeDynamics() or earlier in order for it to have any effect.
     */
    virtual void setColor(const SimTK::State& s, const SimTK::Vec3& color) const = 0;

    virtual double getLength(const SimTK::State& s) const = 0;
    virtual void setLength(const SimTK::State& s, double length) const = 0;

    /**
     * Compute the lengthening speed of the path.
     *
     * @return lengthening speed of the path.
     */
    virtual double getLengtheningSpeed(const SimTK::State& s) const = 0;
    virtual void setLengtheningSpeed(const SimTK::State& s, double speed) const = 0;

    virtual double getPreScaleLength(const SimTK::State& s) const = 0;
    virtual void setPreScaleLength(const SimTK::State& s, double preScaleLength) = 0;

    /**
     * Appends PointForceDirections to the output parameter.
     *
     * These can be used to apply tension to bodies the points are connected to.
     *
     * CAUTION: the return pointers are heap allocated: you must delete them yourself!
     */
    virtual void getPointForceDirections(
            const SimTK::State& s,
            OpenSim::Array<PointForceDirection*>* rPFDs) const = 0;

    /**
     *  Add in the equivalent body and generalized forces to be applied to the
     *  multibody system resulting from a tension along the GeometryPath.
     *
     *  @param         state           state used to evaluate forces
     *  @param[in]     tension         scalar (double) of the applied (+ve) tensile force
     *  @param[in,out] bodyForces      Vector of SpatialVec's (torque, force) on bodies
     *  @param[in,out] mobilityForces  Vector of generalized forces, one per mobility
     */
    virtual void addInEquivalentForces(const SimTK::State& state,
                                       double tension,
                                       SimTK::Vector_<SimTK::SpatialVec>& bodyForces,
                                       SimTK::Vector& mobilityForces) const = 0;

    /**
     * Compute the path's moment arms for  specified coordinate.
     *
     * @param aCoord, the coordinate
     */
    virtual double computeMomentArm(const SimTK::State& s,
                                    const Coordinate& aCoord) const = 0;

    /**
     * Updates the geometry attached to the path (location of path points and
     * connecting segments all in global/inertial frame)
     */
    virtual void updateGeometry(const SimTK::State& s) const = 0;
};
}

#endif // OPENSIM_GEOMETRY_PATH_H_


