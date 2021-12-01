#include "FunctionBasedPath.h"


OpenSim::FunctionBasedPath::FunctionBasedPath()
{
    // TODO: property init etc.
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
