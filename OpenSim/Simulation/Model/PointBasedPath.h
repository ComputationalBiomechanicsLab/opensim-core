#ifndef OPENSIM_POINTBASED_PATH_H_
#define OPENSIM_POINTBASED_PATH_H_

/* -------------------------------------------------------------------------- *
 *                          OpenSim: PointBasedPath.h                         *
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

#include <OpenSim/Simulation/Model/GeometryPath.h>

#include <OpenSim/Simulation/MomentArmSolver.h>
#include <OpenSim/Simulation/Model/PathPointSet.h>
#include <OpenSim/Simulation/Wrap/PathWrapSet.h>

#ifdef SWIG
    #ifdef OSIMSIMULATION_API
        #undef OSIMSIMULATION_API
        #define OSIMSIMULATION_API
    #endif
#endif

namespace OpenSim {
class OSIMSIMULATION_API PointBasedPath final : public GeometryPath {
    OpenSim_DECLARE_CONCRETE_OBJECT(PointBasedPath, GeometryPath);

private:
    OpenSim_DECLARE_UNNAMED_PROPERTY(PathPointSet,
        "The set of points defining the path");

    OpenSim_DECLARE_UNNAMED_PROPERTY(PathWrapSet,
        "The wrap objects that are associated with this path");

    // Solver used to compute moment-arms. The GeometryPath owns this object,
    // but we cannot simply use a unique_ptr because we want the pointer to be
    // cleared on copy.
    SimTK::ResetOnCopy<std::unique_ptr<MomentArmSolver> > _maSolver;

    mutable CacheVariable<double> _lengthCV;
    mutable CacheVariable<double> _speedCV;
    mutable CacheVariable<Array<AbstractPathPoint*>> _currentPathCV;
    mutable CacheVariable<SimTK::Vec3> _colorCV;

public:
    PointBasedPath();

    const PathPointSet& getPathPointSet() const override;
    PathPointSet& updPathPointSet() override;
    AbstractPathPoint* addPathPoint(
            const SimTK::State& s,
            int index,
            const PhysicalFrame& frame) override;
    AbstractPathPoint* appendNewPathPoint(
            const std::string& proposedName,
            const PhysicalFrame& frame,
            const SimTK::Vec3& locationOnFrame) override;
    bool canDeletePathPoint(int index) override;
    bool deletePathPoint(const SimTK::State& s, int index) override;
    bool replacePathPoint(
            const SimTK::State& s,
            AbstractPathPoint* oldPathPoint,
            AbstractPathPoint* newPathPoint) override;
    const Array<AbstractPathPoint*>& getCurrentPath(const SimTK::State& s) const override;


    const PathWrapSet& getWrapSet() const override;
    PathWrapSet& updWrapSet() override;
    void addPathWrap(WrapObject& aWrapObject) override;
    void moveUpPathWrap(const SimTK::State& s, int index) override;
    void moveDownPathWrap(const SimTK::State& s, int index) override;
    void deletePathWrap(const SimTK::State& s, int index) override;

    SimTK::Vec3 getColor(const SimTK::State& s) const override;
    void setColor(const SimTK::State& s, const SimTK::Vec3& color) const override;

    double getLength(const SimTK::State& s) const override;
    void setLength(const SimTK::State& s, double length) const override;

    double getLengtheningSpeed(const SimTK::State& s) const override;
    void setLengtheningSpeed(const SimTK::State& s, double speed) const override;

    void getPointForceDirections(
            const SimTK::State& s,
            OpenSim::Array<PointForceDirection*>* rPFDs) const override;

    void addInEquivalentForces(const SimTK::State& state,
                               double tension,
                               SimTK::Vector_<SimTK::SpatialVec>& bodyForces,
                               SimTK::Vector& mobilityForces) const override;

    double computeMomentArm(const SimTK::State& s,
                            const Coordinate& aCoord) const override;
    void updateGeometry(const SimTK::State& s) const override;

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

#endif
