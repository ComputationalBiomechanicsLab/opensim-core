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
#include <OpenSim/Simulation/MomentArmSolver.h>
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

class Coordinate;
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
  */
class OSIMSIMULATION_API GeometryPath : public ModelComponent {
OpenSim_DECLARE_CONCRETE_OBJECT(GeometryPath, ModelComponent);

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

private:
    OpenSim_DECLARE_UNNAMED_PROPERTY(PathPointSet,
        "The set of points defining the path");

    OpenSim_DECLARE_UNNAMED_PROPERTY(PathWrapSet,
        "The wrap objects that are associated with this path");

    // used for scaling tendon and fiber lengths
    double _preScaleLength;

    // Solver used to compute moment-arms. The GeometryPath owns this object,
    // but we cannot simply use a unique_ptr because we want the pointer to be
    // cleared on copy.
    SimTK::ResetOnCopy<std::unique_ptr<MomentArmSolver> > _maSolver;

    mutable CacheVariable<double> _lengthCV;
    mutable CacheVariable<double> _speedCV;
    mutable CacheVariable<Array<AbstractPathPoint*>> _currentPathCV;
    mutable CacheVariable<SimTK::Vec3> _colorCV;
    
//=============================================================================
// METHODS
//=============================================================================
public:
    GeometryPath();

    const PathPointSet& getPathPointSet() const;
    PathPointSet& updPathPointSet();

    /**
     * Add a new path point, with default location, to the path.
     *
     * @param aIndex The position in the pathPointSet to put the new point in.
     * @param frame The frame to attach the point to.
     * @return Pointer to the newly created path point.
     */
    AbstractPathPoint* addPathPoint(
            const SimTK::State& s,
            int index,
            const PhysicalFrame& frame);
    AbstractPathPoint* appendNewPathPoint(
            const std::string& proposedName,
            const PhysicalFrame& frame,
            const SimTK::Vec3& locationOnFrame);
    /**
     * Returns true if a path point can be deleted. All paths must have at
     * least two active path points to define the path.
     *
     * @param aIndex The index of the point to delete.
     * @return Whether or not the point can be deleted.
     */
    bool canDeletePathPoint(int index);
    /**
     * Delete a path point.
     *
     * @param aIndex The index of the point to delete.
     * @return Whether or not the point was deleted.
     */
    bool deletePathPoint(const SimTK::State& s, int index);
    /**
     * Replace a path point in the set with another point. The new one is made a
     * member of all the same groups as the old one, and is inserted in the same
     * place the old one occupied.
     *
     *  @param aOldPathPoint Path point to remove.
     *  @param aNewPathPoint Path point to add.
     */
    bool replacePathPoint(
            const SimTK::State& s,
            AbstractPathPoint* oldPathPoint,
            AbstractPathPoint* newPathPoint);
    /**
     * Get the current path of the path.
     *
     * @return The array of currently active path points.
     */
    const Array<AbstractPathPoint*>& getCurrentPath(const SimTK::State& s) const;


    const PathWrapSet& getWrapSet() const;
    PathWrapSet& updWrapSet();
    /**
     * Create a new wrap instance and add it to the set.
     *
     * @param aWrapObject The wrap object to use in the new wrap instance.
     */
    void addPathWrap(WrapObject& aWrapObject);
    /**
     * Move a wrap instance up in the list. Changing the order of wrap instances for
     * a path may affect how the path wraps over the wrap objects.
     *
     * @param aIndex The index of the wrap instance to move up.
     */
    void moveUpPathWrap(const SimTK::State& s, int index);
    /**
     * Move a wrap instance down in the list. Changing the order of wrap instances
     * for a path may affect how the path wraps over the wrap objects.
     *
     * @param aIndex The index of the wrap instance to move down.
     */
    void moveDownPathWrap(const SimTK::State& s, int index);
    /**
     * Delete a wrap instance.
     *
     * @param aIndex The index of the wrap instance to delete.
     */
    void deletePathWrap(const SimTK::State& s, int index);


    /**
     * Get the default color of the path.
     *
     * Returns the color that will be used to initialize the color cache
     * at the next extendAddToSystem() call. The actual color used to draw the
     * path will be taken from the cache variable, which may have changed.
     */
    const SimTK::Vec3& getDefaultColor() const;

    /**
     * Set the default color of the path.
     *
     * If you call this prior to extendAddToSystem() it will be used to initialize
     * the color cache variable. Otherwise %GeometryPath will choose its own
     * default, which varies depending on owner.
     */
    void setDefaultColor(const SimTK::Vec3& color);

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
    SimTK::Vec3 getColor(const SimTK::State& s) const;

    /**
     * Set the current color of the path.
     *
     * Internally, sets the value of the color cache variable owned by this
     * %GeometryPath object, in the cache of the given state. The value of this
     * variable is used as the color when the path is drawn, which occurs with
     * the state realized to Stage::Dynamics. Therefore, you must call this method
     * during realizeDynamics() or earlier in order for it to have any effect.
     */
    void setColor(const SimTK::State& s, const SimTK::Vec3& color) const;

    double getLength(const SimTK::State& s) const;
    void setLength(const SimTK::State& s, double length) const;

    /**
     * Compute the lengthening speed of the path.
     *
     * @return lengthening speed of the path.
     */
    double getLengtheningSpeed(const SimTK::State& s) const;
    void setLengtheningSpeed(const SimTK::State& s, double speed) const;

    double getPreScaleLength(const SimTK::State& s) const;
    void setPreScaleLength(const SimTK::State& s, double preScaleLength);

    /**
     * Appends PointForceDirections to the output parameter.
     *
     * These can be used to apply tension to bodies the points are connected to.
     *
     * CAUTION: the return pointers are heap allocated: you must delete them yourself!
     */
    void getPointForceDirections(const SimTK::State& s,
                                 OpenSim::Array<PointForceDirection*>* rPFDs) const;

    /**
     *  Add in the equivalent body and generalized forces to be applied to the
     *  multibody system resulting from a tension along the GeometryPath.
     *
     *  @param         state           state used to evaluate forces
     *  @param[in]     tension         scalar (double) of the applied (+ve) tensile force
     *  @param[in,out] bodyForces      Vector of SpatialVec's (torque, force) on bodies
     *  @param[in,out] mobilityForces  Vector of generalized forces, one per mobility
     */
    void addInEquivalentForces(const SimTK::State& state,
                               double tension,
                               SimTK::Vector_<SimTK::SpatialVec>& bodyForces,
                               SimTK::Vector& mobilityForces) const;

    /**
     * Compute the path's moment arms for  specified coordinate.
     *
     * @param aCoord, the coordinate
     */
    virtual double computeMomentArm(const SimTK::State& s,
                                    const Coordinate& aCoord) const;

    //--------------------------------------------------------------------------
    // SCALING
    //--------------------------------------------------------------------------

    /**
     * Calculate the path length in the current body position and store it for
     * use after the Model has been scaled.
     */
    void extendPreScale(const SimTK::State& s,
                        const ScaleSet& scaleSet) override;

    /**
     * Recalculate the path after the Model has been scaled.
     */
    void extendPostScale(const SimTK::State& s,
                         const ScaleSet& scaleSet) override;

    /**
     * Updates the geometry attached to the path (location of path points and
     * connecting segments all in global/inertial frame)
     */
    virtual void updateGeometry(const SimTK::State& s) const;

protected:
    void extendFinalizeFromProperties() override;
    void extendConnectToModel(Model& aModel) override;
    void extendAddToSystem(SimTK::MultibodySystem& system) const override;
    void extendInitStateFromProperties(SimTK::State& s) const override;

    void generateDecorations(bool fixed,
                             const ModelDisplayHints& hints,
                             const SimTK::State& state,
                             SimTK::Array_<SimTK::DecorativeGeometry>& appendToThis) const override;

private:

    void constructProperties();

    void computePath(const SimTK::State& s) const;
    void computeLengtheningSpeed(const SimTK::State& s) const;

    void applyWrapObjects(const SimTK::State& s, Array<AbstractPathPoint*>& path) const;

    double calcPathLengthChange(const SimTK::State& s,
                                const WrapObject& wo,
                                const WrapResult& wr,
                                const Array<AbstractPathPoint*>& path) const;

    double calcLengthAfterPathComputation(const SimTK::State& s,
                                          const Array<AbstractPathPoint*>& currentPath) const;


    void namePathPoints(int aStartingIndex);
    void placeNewPathPoint(const SimTK::State& s,
                           SimTK::Vec3& aOffset,
                           int index,
                           const PhysicalFrame& frame);

    // Override of the default implementation to account for versioning.
    void updateFromXMLNode(SimTK::Xml::Element& aNode, int versionNumber = -1) override;
};
}

#endif // OPENSIM_GEOMETRY_PATH_H_


