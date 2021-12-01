#include "FunctionBasedPath.h"

#include <OpenSim/Common/Component.h>
#include <OpenSim/Simulation/Model/Model.h>
#include <OpenSim/Simulation/SimbodyEngine/Coordinate.h>

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

        double getLength(const SimTK::State&) const override
        {
            OPENSIM_THROW_FRMOBJ(OpenSim::Exception, g_ThrowingPathErrMsg);
        }

        double getLengtheningSpeed(const SimTK::State&) const override
        {
            OPENSIM_THROW_FRMOBJ(OpenSim::Exception, g_ThrowingPathErrMsg);
        }

        double computeMomentArm(const SimTK::State&, const OpenSim::Coordinate&) const override
        {
            OPENSIM_THROW_FRMOBJ(OpenSim::Exception, g_ThrowingPathErrMsg);
        }

        void getPointForceDirections(const SimTK::State& s, OpenSim::Array<OpenSim::PointForceDirection*>* rPFDs) const override
        {
            OPENSIM_THROW_FRMOBJ(OpenSim::Exception, g_ThrowingPathErrMsg);
        }

        void addInEquivalentForces(const SimTK::State& state, double tension, SimTK::Vector_<SimTK::SpatialVec>& bodyForces, SimTK::Vector& mobilityForces) const override
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

OpenSim::FunctionBasedPath::FunctionBasedPath(PathFunction const& pathFn) : GeometryPath{}
{
    constructProperty_PathFunction(pathFn);
}

OpenSim::FunctionBasedPath::~FunctionBasedPath() noexcept = default;

OpenSim::FunctionBasedPath& OpenSim::FunctionBasedPath::operator=(const FunctionBasedPath&) = default;

OpenSim::FunctionBasedPath& OpenSim::FunctionBasedPath::operator=(FunctionBasedPath&&) noexcept = default;

SimTK::Vec3 OpenSim::FunctionBasedPath::getColor(const SimTK::State& s) const
{
    return getCacheVariableValue(s, _colorCV);
}

void OpenSim::FunctionBasedPath::setColor(const SimTK::State& s, const SimTK::Vec3& color) const
{
    setCacheVariableValue(s, _colorCV, color);
}

double OpenSim::FunctionBasedPath::getLength(const SimTK::State& s) const
{
    if (isCacheVariableValid(s, _lengthCV)) {
        return getCacheVariableValue(s, _lengthCV);
    }

    double v = getProperty_PathFunction().getValue().getLength(s);
    setCacheVariableValue(s, _lengthCV, v);
    return v;
}

double OpenSim::FunctionBasedPath::getLengtheningSpeed(const SimTK::State& s) const
{
    if (isCacheVariableValid(s, _speedCV)) {
        return getCacheVariableValue(s, _speedCV);
    }

    double v = getProperty_PathFunction().getValue().getLengtheningSpeed(s);
    setCacheVariableValue(s, _speedCV, v);
    return v;
}

void OpenSim::FunctionBasedPath::getPointForceDirections(const SimTK::State& s, OpenSim::Array<PointForceDirection*>* rPFDs) const
{
    getProperty_PathFunction().getValue().getPointForceDirections(s, rPFDs);
}

void OpenSim::FunctionBasedPath::addInEquivalentForces(const SimTK::State& state, double tension, SimTK::Vector_<SimTK::SpatialVec>& bodyForces, SimTK::Vector& mobilityForces) const
{
    getProperty_PathFunction().getValue().addInEquivalentForces(state, tension, bodyForces, mobilityForces);
}

double OpenSim::FunctionBasedPath::computeMomentArm(const SimTK::State& s, const Coordinate& aCoord) const
{
    return getProperty_PathFunction().getValue().computeMomentArm(s, aCoord);
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

void OpenSim::FunctionBasedPath::extendInitStateFromProperties(SimTK::State& s) const
{
    Super::extendInitStateFromProperties(s);
    markCacheVariableValid(s, _colorCV);
}
