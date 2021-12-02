#ifndef OPENSIM_FUNCTIONBASED_PATH_H_
#define OPENSIM_FUNCTIONBASED_PATH_H_

#include <OpenSim/Common/Array.h>
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
 * An interface for an object that can compute the length, lengthening speed,
 * and moment arm (w.r.t. a particular `OpenSim::Coordinate`) of a path at
 * runtime.
 *
 * See `FunctionBasedPath` for a standard implementation of a `GeometryPath`
 * that automatically handles forwarding calls to the `PathFunction`, state
 * caching, etc.
 */
class OSIMSIMULATION_API PathFunction : public OpenSim::ModelComponent {
    OpenSim_DECLARE_ABSTRACT_OBJECT(PathFunction, OpenSim::Component)

public:
    virtual ~PathFunction() noexcept = default;

    virtual double getLength(const SimTK::State&) const = 0;
    virtual double getLengtheningSpeed(const SimTK::State&) const = 0;
    virtual double computeMomentArm(const SimTK::State&, const OpenSim::Coordinate&) const = 0;
    virtual void getPointForceDirections(const SimTK::State& s, OpenSim::Array<PointForceDirection*>* rPFDs) const = 0;
    virtual void addInEquivalentForces(const SimTK::State& state, double tension, SimTK::Vector_<SimTK::SpatialVec>& bodyForces, SimTK::Vector& mobilityForces) const = 0;
};

/**
 * An `OpenSim::GeometryPath` that uses `PathFunction`s to compute its state.
 *
 * A `FunctionBasedPath` uses an externally-provided `PathFunction` to compute the
 * actual path, and also handles any `GeometryPath`-specific concerns, such as
 * handling coloring and caching results.
 */
class OSIMSIMULATION_API FunctionBasedPath final : public GeometryPath {
    OpenSim_DECLARE_CONCRETE_OBJECT(FunctionBasedPath, GeometryPath);

    OpenSim_DECLARE_PROPERTY(PathFunction, PathFunction, "The underlying function that is used at simulation-time to evaluate the length, lengthening speed, and moment arm of the path.");

    mutable CacheVariable<double> _lengthCV;
    mutable CacheVariable<double> _speedCV;
    mutable CacheVariable<SimTK::Vec3> _colorCV;

public:
    FunctionBasedPath();
    FunctionBasedPath(const FunctionBasedPath&);
    FunctionBasedPath(FunctionBasedPath&&) noexcept;
    explicit FunctionBasedPath(PathFunction const&);
    ~FunctionBasedPath() noexcept;

    FunctionBasedPath& operator=(FunctionBasedPath const&);
    FunctionBasedPath& operator=(FunctionBasedPath&&) noexcept;

    SimTK::Vec3 getColor(const SimTK::State& s) const override;
    void setColor(const SimTK::State& s, const SimTK::Vec3& color) const override;

    double getLength(const SimTK::State& s) const override;
    double getLengtheningSpeed(const SimTK::State& s) const override;

    void getPointForceDirections(
            const SimTK::State& s,
            OpenSim::Array<PointForceDirection*>* rPFDs) const override;

    void addInEquivalentForces(const SimTK::State& state,
                               double tension,
                               SimTK::Vector_<SimTK::SpatialVec>& bodyForces,
                               SimTK::Vector& mobilityForces) const override;

    double computeMomentArm(const SimTK::State& s,
                            const Coordinate& aCoord) const override;

    void extendAddToSystem(SimTK::MultibodySystem& system) const override;
    void extendInitStateFromProperties(SimTK::State& s) const override;
};
}

#endif // FUNCTIONBASEDPATH_H
