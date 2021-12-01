#ifndef OPENSIM_FUNCTIONBASED_PATH_H_
#define OPENSIM_FUNCTIONBASED_PATH_H_

#include <OpenSim/Common/Array.h>
#include <OpenSim/Simulation/Model/GeometryPath.h>
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
 * An `OpenSim::GeometryPath` that uses functions to compute its state.
 */
class FunctionBasedPath final : public GeometryPath {
    OpenSim_DECLARE_CONCRETE_OBJECT(FunctionBasedPath, GeometryPath);

public:
    FunctionBasedPath();
    FunctionBasedPath(const FunctionBasedPath&);
    FunctionBasedPath(FunctionBasedPath&&) noexcept;
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
};
}

#endif // FUNCTIONBASEDPATH_H
