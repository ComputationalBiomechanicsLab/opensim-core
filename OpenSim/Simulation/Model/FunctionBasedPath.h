#ifndef OPENSIM_FUNCTIONBASED_PATH_H_
#define OPENSIM_FUNCTIONBASED_PATH_H_

/* -------------------------------------------------------------------------- *
 *                          OpenSim: FunctionBasedPath.h                      *
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

#include <OpenSim/Common/Array.h>
#include <OpenSim/Simulation/Model/AbstractPathPoint.h>
#include <OpenSim/Simulation/Model/GeometryPath.h>
#include <OpenSim/Simulation/Model/FunctionBasedPathDiscretizationSet.h>
#include <SimTKcommon/internal/ClonePtr.h>
#include <SimTKcommon/internal/Vec.h>

#include <string>

#ifdef SWIG
    #ifdef OSIMSIMULATION_API
        #undef OSIMSIMULATION_API
        #define OSIMSIMULATION_API
    #endif
#endif

// forward declarations
namespace OpenSim {
class Component;
class Coordinate;
class Model;
class PointBasedPath;
class PointForceDirection;
class PhysicalFrame;
class WrapObject;
class WrapResult;
}

namespace SimTK {
class State;
}

namespace OpenSim {
/**
 * An `OpenSim::GeometryPath` that computes its length/lengtheningSpeed by
 * interpolating a pre-computed Bezier curve.
 *
 * The precomputed curve can be produced from by fitting an existing
 * `OpenSim::PointBasedPath` in a Model. You can use
 * `FunctionBasedPath::fromPointBasedPath` to do that in code, or use the
 * FunctionBasedPathModelTool do convert all `PointBasedPaths` in an osim into
 * `FunctionBasedPath`s.
 *
 * The curve fitting implementation requires access to the whole OpenSim model
 * that the (concrete) `PointBasedPath` is in. This is because the fitting
 * implementation permutes the model's coordinates to fit the curve, and
 * because a curve's paramaterization depends on its place in the Model as
 * a whole.
 */
class OSIMSIMULATION_API FunctionBasedPath : public GeometryPath {
    OpenSim_DECLARE_CONCRETE_OBJECT(FunctionBasedPath, GeometryPath);

public:
    OpenSim_DECLARE_UNNAMED_PROPERTY(FunctionBasedPathDiscretizationSet, "Discretizations that were used for each OpenSim::Coordinate that the path was fitted against");
    OpenSim_DECLARE_LIST_PROPERTY(Evaluations, double, "The evaluated results of each *permutation* of discretizations. The FunctionBasedPathDiscretizationSet property describes how each OpenSim::Coordinate was discretized. These evaluations are the result of permuting through all possible combinations of discretizations. Effectively, this property contains a N-dimensional 'surface' of points, where each dimension of the surface is a Coordinate, and each dimension of each point is one of the evenly-spaced points in the discretization range [x_begin, x_range] for each dimension");

    struct Impl;
private:
    SimTK::ClonePtr<Impl> _impl;

public:
    /**
     * Returns an in-memory representation of a FunctionBasedPath generated
     * by fitting curves against the supplied PointBasedPath.
     */
    static FunctionBasedPath fromPointBasedPath(const Model& model, const PointBasedPath& pbp);

    FunctionBasedPath();
    FunctionBasedPath(const FunctionBasedPath&);
    FunctionBasedPath(FunctionBasedPath&&);
    FunctionBasedPath& operator=(const FunctionBasedPath&);
    FunctionBasedPath& operator=(FunctionBasedPath&&);
    ~FunctionBasedPath() noexcept override;

    void extendFinalizeFromProperties() override;
    void extendFinalizeConnections(Component& root) override;

    double getLength(const SimTK::State& s) const override;
    void setLength(const SimTK::State& s, double length) const override;

    double getLengtheningSpeed(const SimTK::State& s) const override;
    void setLengtheningSpeed(const SimTK::State& s, double speed) const override;

    double computeMomentArm(const SimTK::State& s, const Coordinate& aCoord) const override;

    void addInEquivalentForces(const SimTK::State& state,
                               const double& tension,
                               SimTK::Vector_<SimTK::SpatialVec>& bodyForces,
                               SimTK::Vector& mobilityForces) const override;

    void generateDecorations(bool fixed,
                             const ModelDisplayHints& hints,
                             const SimTK::State& state,
                             SimTK::Array_<SimTK::DecorativeGeometry>& appendToThis) const override;

protected:
    void computePath(const SimTK::State& s ) const override;
    void computeLengtheningSpeed(const SimTK::State& s) const override;

    double calcLengthAfterPathComputation(const SimTK::State& s,
                                          const Array<AbstractPathPoint*>& currentPath) const override;

    double calcPathLengthChange(const SimTK::State& s,
                                const WrapObject& wo,
                                const WrapResult& wr,
                                const Array<AbstractPathPoint*>& path) const override;

private:
    void getPointForceDirections(const SimTK::State& s,
                                 OpenSim::Array<PointForceDirection*> *rPFDs) const override;

    const Array<AbstractPathPoint*>& getCurrentPath(const SimTK::State& s) const override;

    AbstractPathPoint* addPathPoint(const SimTK::State& s,
                                    int index,
                                    const PhysicalFrame& frame) override;

    AbstractPathPoint* appendNewPathPoint(const std::string& proposedName,
                                          const PhysicalFrame& frame,
                                          const SimTK::Vec3& locationOnFrame) override;

    bool canDeletePathPoint(int index) override;

    bool deletePathPoint(const SimTK::State& s, int index) override;

    bool replacePathPoint(const SimTK::State& s,
                          AbstractPathPoint* oldPathPoint,
                          AbstractPathPoint* newPathPoint) override;


    void addPathWrap(WrapObject& aWrapObject) override;

    void deletePathWrap(const SimTK::State& s, int index) override;

    void moveUpPathWrap(const SimTK::State& s,
                        int index) override;

    void moveDownPathWrap(const SimTK::State& s,
                          int index) override;

    void applyWrapObjects(const SimTK::State& s,
                          Array<AbstractPathPoint*>& path) const override;

    void namePathPoints(int aStartingIndex) override;

    void placeNewPathPoint(const SimTK::State& s,
                           SimTK::Vec3& aOffset,
                           int index,
                           const PhysicalFrame& frame) override;
};
}

#endif
