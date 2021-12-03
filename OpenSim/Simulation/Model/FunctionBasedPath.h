#ifndef OPENSIM_FUNCTIONBASED_PATH_H_
#define OPENSIM_FUNCTIONBASED_PATH_H_

#include <OpenSim/Common/Array.h>
#include <OpenSim/Common/Function.h>
#include <OpenSim/Simulation/Model/GeometryPath.h>
#include <OpenSim/Simulation/Model/ModelComponent.h>
#include <SimTKcommon/internal/Vec.h>
#include <SimTKcommon/internal/Vector_.h>

#ifdef SWIG
    #ifdef OSIMSIMULATION_API
        #undef OSIMSIMULATION_API
        #define OSIMSIMULATION_API
    #endif
#endif

namespace SimTK {
class State;
}

namespace OpenSim {

class Coordinate;
class PointForceDirection;

/**
 * An `OpenSim::GeometryPath` that uses `PathFunction`s to compute its state.
 *
 * A `FunctionBasedPath` uses an externally-provided `PathFunction` to compute the
 * actual path, and also handles any `GeometryPath`-specific concerns, such as
 * handling coloring and caching results.
 */
class OSIMSIMULATION_API FunctionBasedPath final : public GeometryPath {
    OpenSim_DECLARE_CONCRETE_OBJECT(FunctionBasedPath, GeometryPath);

    OpenSim_DECLARE_PROPERTY(PathFunction, OpenSim::Function, "The underlying function that is used at simulation-time to evaluate the length and lengthening speed (derivative) of the path. The function's arity (and argument order) is equal to the coordinates that affect the path. At evaluation-time, the function (value + derivative) is called with a sequence of the coordinate values (the values of each coordinate in the current state)");
    OpenSim_DECLARE_LIST_PROPERTY(Coordinates, std::string, "Absolute paths to coordinates in the model. The number of coordinates provided must equal the arity of the function property. At simulation-time, each of these coordinates are looked up to get their values, which are pumped into the path evaluation function");

    mutable CacheVariable<double> _lengthCV;
    mutable CacheVariable<double> _speedCV;
    mutable CacheVariable<SimTK::Vec3> _colorCV;
    mutable SimTK::Vector _functionArgsBuffer;
    mutable std::vector<int> _derivativeOrderBuffer;

public:
    FunctionBasedPath();
    FunctionBasedPath(const FunctionBasedPath&);

    /**
     * Construct the FunctionBasedPath and immediately assign its PathFunction and
     * Coordinates list.
     *
     * The length of the coordinates list must match the artity (getArgumentSize)
     * of the function. The function must be differentiable.
     */
    FunctionBasedPath(const OpenSim::Function&, std::vector<std::string> coordAbsPaths);
    ~FunctionBasedPath() noexcept;

    FunctionBasedPath& operator=(FunctionBasedPath const&);

    SimTK::Vec3 getColor(const SimTK::State& s) const override;
    void setColor(const SimTK::State& s, const SimTK::Vec3& color) const override;

    double getLength(const SimTK::State& s) const override;
    double getLengtheningSpeed(const SimTK::State& s) const override;

    void addInEquivalentForces(const SimTK::State& state,
                               double tension,
                               SimTK::Vector_<SimTK::SpatialVec>& bodyForces,
                               SimTK::Vector& mobilityForces) const override;

    double computeMomentArm(const SimTK::State& s,
                            const Coordinate& aCoord) const override;

    void extendFinalizeFromProperties() override;
    void extendAddToSystem(SimTK::MultibodySystem& system) const override;
    void extendInitStateFromProperties(SimTK::State& s) const override;
    void extendFinalizeConnections(OpenSim::Component&) override;
private:
    const SimTK::Vector& calcFunctionArguments(const SimTK::State&) const;
    int indexOfCoordinate(const Coordinate&) const;
};
}

#endif // FUNCTIONBASEDPATH_H
