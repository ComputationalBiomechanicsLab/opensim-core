#include "FunctionBasedPath.h"

namespace {
    // used as a stub PathFunction implementation that throws if the user tries
    // to actually use it
    //
    // handy for handling the case where we need *some* PathFunction to stand-in
    // before (e.g.) a "proper" path function is assigned to a property or
    // something.
    static char const g_ThrowingPathErrMsg[] = "Cannot call any member on a ThrowingPathFunction. This exception has *probably* been throwng because a PathFunction has not (yet) been assigned for the parent component. You *probably* need to ensure that a FunctionBasedPath actually has PathFunction assigned to it";
    class ThrowingPathFunction final : public OpenSim::PathFunction {
        OpenSim_DECLARE_CONCRETE_OBJECT(ThrowingPathFunction, OpenSim::PathFunction)

        double getLength(const SimTK::State&) override
        {
            OPENSIM_THROW_FRMOBJ(OpenSim::Exception, g_ThrowingPathErrMsg);
        }

        double getLengtheningSpeed(const SimTK::State&) override
        {
            OPENSIM_THROW_FRMOBJ(OpenSim::Exception, g_ThrowingPathErrMsg);
        }

        double computeMomentArm(const SimTK::State&, const OpenSim::Coordinate&) override
        {
            OPENSIM_THROW_FRMOBJ(OpenSim::Exception, g_ThrowingPathErrMsg);
        }
    };
}

OpenSim::FunctionBasedPath::FunctionBasedPath() : GeometryPath{}
{
    constructProperty_PathFunction(ThrowingPathFunction{});
}

OpenSim::FunctionBasedPath::FunctionBasedPath(const FunctionBasedPath&) = default;

OpenSim::FunctionBasedPath::FunctionBasedPath(FunctionBasedPath&&) noexcept = default;

OpenSim::FunctionBasedPath::~FunctionBasedPath() noexcept = default;

OpenSim::FunctionBasedPath& OpenSim::FunctionBasedPath::operator=(const FunctionBasedPath&) = default;

OpenSim::FunctionBasedPath& OpenSim::FunctionBasedPath::operator=(FunctionBasedPath&&) noexcept = default;

SimTK::Vec3 OpenSim::FunctionBasedPath::getColor(const SimTK::State& s) const
{
    // TODO: read from cache variable or similar
    return SimTK::Vec3{};
}

void OpenSim::FunctionBasedPath::setColor(const SimTK::State& s, const SimTK::Vec3& color) const
{
    // TODO: save to a cache variable or similar
}

double OpenSim::FunctionBasedPath::getLength(const SimTK::State& s) const
{
    // TODO: compute via function abstraction
    return 0.0;
}

double OpenSim::FunctionBasedPath::getLengtheningSpeed(const SimTK::State& s) const
{
    // TODO: compute via function abstraction
    return 0.0;
}

void OpenSim::FunctionBasedPath::getPointForceDirections(const SimTK::State& s, OpenSim::Array<PointForceDirection*>* rPFDs) const
{
    // TODO
}

void OpenSim::FunctionBasedPath::addInEquivalentForces(const SimTK::State& state, double tension, SimTK::Vector_<SimTK::SpatialVec>& bodyForces, SimTK::Vector& mobilityForces) const
{
    // TODO
}

double OpenSim::FunctionBasedPath::computeMomentArm(const SimTK::State& s, const Coordinate& aCoord) const
{
    // TODO: compute via function abstraction
    return 0.0;
}

void OpenSim::FunctionBasedPath::extendAddToSystem(SimTK::MultibodySystem& system) const
{
    Super::extendAddToSystem(system);

    // Allocate cache entries to save the current length and speed(=d/dt length)
    // of the path in the cache. Length depends only on q's so will be valid
    // after Position stage, speed requires u's also so valid at Velocity stage.
    this->_lengthCV = addCacheVariable("length", 0.0, SimTK::Stage::Position);
    this->_speedCV = addCacheVariable("speed", 0.0, SimTK::Stage::Velocity);

    // We consider this cache entry valid any time after it has been created
    // and first marked valid, and we won't ever invalidate it.
    this->_colorCV = addCacheVariable("color", get_Appearance().get_color(), SimTK::Stage::Topology);
}
