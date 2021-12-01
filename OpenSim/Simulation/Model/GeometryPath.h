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

    class Impl;
private:
    SimTK::ClonePtr<Impl> _impl;

//=============================================================================
// METHODS
//=============================================================================
public:
    GeometryPath();
    GeometryPath(const GeometryPath&);
    GeometryPath(GeometryPath&&) noexcept;
    ~GeometryPath() noexcept;

    GeometryPath& operator=(const GeometryPath&);
    GeometryPath& operator=(GeometryPath&&) noexcept;


    // INTERFACE METHODS
    //
    // Concrete implementations of `GeometryPath` *must* provide these.

    /**
     * Get the current color of the path.
     *
     * This is the runtime, potentially state-dependent, color of the path. It
     * is the color used to display the path in that state (e.g. for UI rendering).
     *
     * This color value is typically initialized with the default color (see:
     * `getDefaultColor`), but the color can change between simulation states
     * because downstream code (e.g. muscles) might call `setColor` to implement
     * state-dependent path coloring.
     */
    virtual SimTK::Vec3 getColor(const SimTK::State& s) const = 0;

    /**
     * Set the current color of the path.
     *
     * Internally, sets the current color value of the path for the provided state
     * (e.g. using cache variables).
     *
     * The value of this variable is used as the color when the path is drawn, which
     * occurs with the state realized to Stage::Dynamics. Therefore, you must call
     * this method during realizeDynamics() or earlier in order for it to have any
     * effect.
     */
    virtual void setColor(const SimTK::State& s, const SimTK::Vec3& color) const = 0;

    /**
     * Get the current length of the path.
     *
     * Internally, this may use a variety of methods to figure out how long the path
     * is, such as using spline-fits, or computing the distance between points in
     * space. It is up to concrete implementations (e.g. `PointBasedPath`) to provide
     * a relevant implementation.
     */
    virtual double getLength(const SimTK::State& s) const = 0;

    /**
     * Get the lengthening speed of the path.
     *
     * Internally, this may use a variety of methods to figure out the lengthening
     * speed. It might use the finite difference between two lengths, or an analytic
     * solution, or always return `0.0`. It is up to concrete implementations (e.g.
     * `PointBasedPath`) to provide a relevant implementation.
     */
    virtual double getLengtheningSpeed(const SimTK::State& s) const = 0;

    /**
     * Appends PointForceDirections to the output parameter.
     *
     * These can be used to apply tension to bodies the points are connected to.
     *
     * CAUTION: the return pointers are heap allocated: you must delete them yourself!
     */
    DEPRECATED_14("Avoid using GeometryPath::getPointForceDirections(...): prefer GeometryPath::addInEquivalentForces(...) instead.")
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
     * Returns the moment arm of the path in the given state with respect to
     * the specified coordinate.
     */
    virtual double computeMomentArm(const SimTK::State& s,
                                    const Coordinate& aCoord) const = 0;


    // DEFAULTED METHODS
    //
    // These are non-deprecated methods that GeometryPath provides a default
    // implementation of.

    /**
     * Get the default color of the path.
     *
     * Returns the color that will be used to initialize the color cache
     * at the next extendAddToSystem() call. Use `getColor` to retrieve the
     * (potentially different) color that will be used to draw the path.
     */
    const SimTK::Vec3& getDefaultColor() const;

    /**
     * Set the default color of the path.
     *
     * Sets the internal, default, color value for the path. This is the color
     * that's used when the simulation is initialized (specifically, during the
     * `extendAddToSystem` call).
     *
     * This color is not necessarily the *current* color of the path. Other code
     * in the system (e.g. muscle implementations) may change the runtime color
     * with `setColor`. Use `getColor`, with a particular simulation state, to
     * get the color of the path in that state.
     */
    void setDefaultColor(const SimTK::Vec3& color);

    /**
     * Get the current length of the path, *before* the last set of scaling operations
     * were applied to it.
     *
     * Internally, the path stores the original length in a `double` during `extendPreScale`.
     * Therefore, be *very* careful with this method, because the recorded length is dependent
     * on the length as computed during `extendPreScale`, which may have been called with a
     * different state.
     */
    double getPreScaleLength(const SimTK::State& s) const;
    void setPreScaleLength(const SimTK::State& s, double preScaleLength);


    // DEPRECATED METHODS
    //
    // These are here for backwards compatability. Most of these methods are
    // here because `GeometryPath` used to always be a `PointBasedPath`.
    // In later versions of OpenSim, `GeometryPath` was generalized to mean
    // "any path with a length and lengthening speed" to support features
    // like spline-based muscle paths.
    //
    // The base implementation of these methods provides a "stub" implementation
    // that is guaranteed to at least be crash-free and return valid points, but
    // not actually affect the path. Downstream implementations (e.g. an actual
    // `PointBasedPath`) can override these to provide a "real" implementation.

    DEPRECATED_14("Avoid using this method on the GeometryPath base class, because the path may not contain points (e.g. if it is function-based). Instead, test whether the path is a PointBasedPath (e.g. by downcasting with dynamic_cast or similar) and then use the same method on PointBasedPath (which will definitely work as intended)")
    virtual const PathPointSet& getPathPointSet() const;

    DEPRECATED_14("Avoid using this method on the GeometryPath base class, because the path may not contain points (e.g. if it is function-based). Instead, test whether the path is a PointBasedPath (e.g. by downcasting with dynamic_cast or similar) and then use the same method on PointBasedPath (which will definitely work as intended)")
    virtual PathPointSet& updPathPointSet();

    /**
     * Add a new path point, with default location, to the path.
     *
     * @param aIndex The position in the pathPointSet to put the new point in.
     * @param frame The frame to attach the point to.
     * @return Pointer to the newly created path point.
     */
    DEPRECATED_14("Avoid using this method on the GeometryPath base class, because the path may not contain points (e.g. if it is function-based). Instead, test whether the path is a PointBasedPath (e.g. by downcasting with dynamic_cast or similar) and then use the same method on PointBasedPath (which will definitely work as intended)")
    virtual AbstractPathPoint* addPathPoint(
            const SimTK::State& s,
            int index,
            const PhysicalFrame& frame);

    DEPRECATED_14("Avoid using this method on the GeometryPath base class, because the path may not contain points (e.g. if it is function-based). Instead, test whether the path is a PointBasedPath (e.g. by downcasting with dynamic_cast or similar) and then use the same method on PointBasedPath (which will definitely work as intended)")
    virtual AbstractPathPoint* appendNewPathPoint(
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
    DEPRECATED_14("Avoid using this method on the GeometryPath base class, because the path may not contain points (e.g. if it is function-based). Instead, test whether the path is a PointBasedPath (e.g. by downcasting with dynamic_cast or similar) and then use the same method on PointBasedPath (which will definitely work as intended)")
    virtual bool canDeletePathPoint(int index);

    /**
     * Delete a path point.
     *
     * @param aIndex The index of the point to delete.
     * @return Whether or not the point was deleted.
     */
    DEPRECATED_14("Avoid using this method on the GeometryPath base class, because the path may not contain points (e.g. if it is function-based). Instead, test whether the path is a PointBasedPath (e.g. by downcasting with dynamic_cast or similar) and then use the same method on PointBasedPath (which will definitely work as intended)")
    virtual bool deletePathPoint(const SimTK::State& s, int index);

    /**
     * Replace a path point in the set with another point. The new one is made a
     * member of all the same groups as the old one, and is inserted in the same
     * place the old one occupied.
     *
     *  @param aOldPathPoint Path point to remove.
     *  @param aNewPathPoint Path point to add.
     */
    DEPRECATED_14("Avoid using this method on the GeometryPath base class, because the path may not contain points (e.g. if it is function-based). Instead, test whether the path is a PointBasedPath (e.g. by downcasting with dynamic_cast or similar) and then use the same method on PointBasedPath (which will definitely work as intended)")
    virtual bool replacePathPoint(
            const SimTK::State& s,
            AbstractPathPoint* oldPathPoint,
            AbstractPathPoint* newPathPoint);

    /**
     * Get the current path of the path.
     *
     * @return The array of currently active path points.
     */
    DEPRECATED_14("Avoid using this method on the GeometryPath base class, because the path may not contain points (e.g. if it is function-based). Instead, test whether the path is a PointBasedPath (e.g. by downcasting with dynamic_cast or similar) and then use the same method on PointBasedPath (which will definitely work as intended)")
    virtual const Array<AbstractPathPoint*>& getCurrentPath(const SimTK::State& s) const;


    DEPRECATED_14("Avoid using this method on the GeometryPath base class, because the path may not contain points (e.g. if it is function-based). Instead, test whether the path is a PointBasedPath (e.g. by downcasting with dynamic_cast or similar) and then use the same method on PointBasedPath (which will definitely work as intended)")
    virtual const PathWrapSet& getWrapSet() const;

    DEPRECATED_14("Avoid using this method on the GeometryPath base class, because the path may not contain points (e.g. if it is function-based). Instead, test whether the path is a PointBasedPath (e.g. by downcasting with dynamic_cast or similar) and then use the same method on PointBasedPath (which will definitely work as intended)")
    virtual PathWrapSet& updWrapSet();

    /**
     * Create a new wrap instance and add it to the set.
     *
     * @param aWrapObject The wrap object to use in the new wrap instance.
     */
    DEPRECATED_14("Avoid using this method on the GeometryPath base class, because the path may not contain points (e.g. if it is function-based). Instead, test whether the path is a PointBasedPath (e.g. by downcasting with dynamic_cast or similar) and then use the same method on PointBasedPath (which will definitely work as intended)")
    virtual void addPathWrap(WrapObject& aWrapObject);

    /**
     * Move a wrap instance up in the list. Changing the order of wrap instances for
     * a path may affect how the path wraps over the wrap objects.
     *
     * @param aIndex The index of the wrap instance to move up.
     */
    DEPRECATED_14("Avoid using this method on the GeometryPath base class, because the path may not contain points (e.g. if it is function-based). Instead, test whether the path is a PointBasedPath (e.g. by downcasting with dynamic_cast or similar) and then use the same method on PointBasedPath (which will definitely work as intended)")
    virtual void moveUpPathWrap(const SimTK::State& s, int index);

    /**
     * Move a wrap instance down in the list. Changing the order of wrap instances
     * for a path may affect how the path wraps over the wrap objects.
     *
     * @param aIndex The index of the wrap instance to move down.
     */
    DEPRECATED_14("Avoid using this method on the GeometryPath base class, because the path may not contain points (e.g. if it is function-based). Instead, test whether the path is a PointBasedPath (e.g. by downcasting with dynamic_cast or similar) and then use the same method on PointBasedPath (which will definitely work as intended)")
    virtual void moveDownPathWrap(const SimTK::State& s, int index);

    /**
     * Delete a wrap instance.
     *
     * @param aIndex The index of the wrap instance to delete.
     */
    DEPRECATED_14("Avoid using this method on the GeometryPath base class, because the path may not contain points (e.g. if it is function-based). Instead, test whether the path is a PointBasedPath (e.g. by downcasting with dynamic_cast or similar) and then use the same method on PointBasedPath (which will definitely work as intended)")
    virtual void deletePathWrap(const SimTK::State& s, int index);


    DEPRECATED_14("Avoid using GeometryPath::setLength(...): it shouldn't be possible to externally set the length of a (potentially, computed) path.")
    virtual void setLength(const SimTK::State& s, double length) const;

    DEPRECATED_14("Avoid using GeometryPath::setLengtheningSpeed(...): it shouldn't be possible to externally set the lengthening speed of a (potentially, computed) path.")
    virtual void setLengtheningSpeed(const SimTK::State& s, double speed) const;

    /**
     * Proactively updates any decorative geometry attached to the path.
     *
     * This method is an advanced optimization method that attempts to
     * proactively populate any internal datastructures that are used to
     * generate path decorations. E.g. the location of path points, the
     * connecting segments/cylinders of those path points, etc.
     */
    DEPRECATED_14("Avoid using GeometryPath::updateGeometry(...): implementations of GeometryPath should handle any decoration updates, caching, etc. when the caller uses `GeometryPath::generateDecorations(...)`")
    virtual void updateGeometry(const SimTK::State& s) const;
};
}

#endif // OPENSIM_GEOMETRY_PATH_H_


