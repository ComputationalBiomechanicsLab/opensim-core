/* -------------------------------------------------------------------------- *
 *                         OpenSim:  GeometryPath.cpp                         *
 * -------------------------------------------------------------------------- *
 * The OpenSim API is a toolkit for musculoskeletal modeling and simulation.  *
 * See http://opensim.stanford.edu and the NOTICE file for more information.  *
 * OpenSim is developed at Stanford University and supported by the US        *
 * National Institutes of Health (U54 GM072970, R24 HD065690) and by DARPA    *
 * through the Warrior Web program.                                           *
 *                                                                            *
 * Copyright (c) 2005-2017 Stanford University and the Authors                *
 * Author(s): Peter Loan, Ajay Seth                                           *
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

#include "GeometryPath.h"

#include <OpenSim/Simulation/Model/Appearance.h>


class OpenSim::GeometryPath::Impl final {
public:
    // used by `(get|set)PreLengthScale`. Used during `extend(Pre|Post)Scale` by
    // downstream users of GeometryPath to cache the length of the path before
    // scaling.
    //
    // Ideally, downstream classes would perform the caching themselves, because
    // the GeometryPath API isn't an ideal place to store this information. This
    // field is mostly here for backwards-compatability with the API.
    double _preScaleLength = 0.0;

private:
    // used whenever something sneaky is happening, like when a mutation is happening
    // underneath a const method. Multiple threads *should* be able to use the `const`
    // API without any concern about races.
    SimTK::ResetOnCopy<std::mutex> constCastMutex;

    // used as a stand-in if the caller uses the (deprecated) point-based API
    SimTK::ClonePtr<OpenSim::PathPointSet> fixupPathPointSet = nullptr;

    // used as a stand-in if the caller uses the (deprecated) point-based API
    SimTK::ClonePtr<OpenSim::PathPoint> fixupPathPoint = nullptr;

    // used as a stand-in if the caller uses the (deprecated) point-based API
    OpenSim::Array<AbstractPathPoint*> fixupCurrentPathArray;

    // used as a stand-in if the caller uses the (deprecated) point-based API
    SimTK::ClonePtr<OpenSim::PathWrapSet> fixupPathWrapSet = nullptr;

public:
    Impl* clone() const {
        return new Impl{*this};
    }

    OpenSim::PathPointSet& getOrUpdPPSCached() const
    {
        if (!fixupPathPointSet) {
            std::lock_guard<std::mutex> guard{const_cast<Impl&>(*this).constCastMutex.updT()};
            const_cast<Impl&>(*this).fixupPathPointSet = new OpenSim::PathPointSet{};
        }
        return *const_cast<Impl&>(*this).fixupPathPointSet;
    }

    OpenSim::AbstractPathPoint* getOrUpdFixupPathPointCached() const
    {
        if (!fixupPathPoint) {
            std::lock_guard<std::mutex> guard{const_cast<Impl&>(*this).constCastMutex.updT()};
            const_cast<Impl&>(*this).fixupPathPoint = new OpenSim::PathPoint{};
        }
        return const_cast<Impl&>(*this).fixupPathPoint.upd();
    }

    OpenSim::Array<AbstractPathPoint*>& getOrUpdCurrentPathArray() const
    {
        return const_cast<Impl&>(*this).fixupCurrentPathArray;
    }

    PathWrapSet& getOrUpdPathWrapSet() const
    {
        if (!fixupPathWrapSet) {
            std::lock_guard<std::mutex> guard{const_cast<Impl&>(*this).constCastMutex.updT()};
            const_cast<Impl&>(*this).fixupPathWrapSet = new OpenSim::PathWrapSet{};
        }

        return *const_cast<Impl&>(*this).fixupPathWrapSet.upd();
    }
};

OpenSim::GeometryPath::GeometryPath() : ModelComponent{}, _impl{new Impl{}}
{
    setAuthors("Peter Loan");

    Appearance appearance;
    appearance.set_color(SimTK::Gray);
    constructProperty_Appearance(appearance);
}

OpenSim::GeometryPath::GeometryPath(OpenSim::GeometryPath const&) = default;

OpenSim::GeometryPath::GeometryPath(OpenSim::GeometryPath&&) noexcept = default;

OpenSim::GeometryPath::~GeometryPath() noexcept = default;

OpenSim::GeometryPath& OpenSim::GeometryPath::operator=(const GeometryPath&) = default;

OpenSim::GeometryPath& OpenSim::GeometryPath::operator=(GeometryPath&&) noexcept = default;


// DEFAULTED METHODS

const SimTK::Vec3& OpenSim::GeometryPath::getDefaultColor() const
{
    return get_Appearance().get_color();
}

void OpenSim::GeometryPath::setDefaultColor(const SimTK::Vec3& color)
{
    updProperty_Appearance().setValueIsDefault(false);
    upd_Appearance().set_color(color);
}

double OpenSim::GeometryPath::getPreScaleLength(const SimTK::State&) const
{
    return _impl->_preScaleLength;
}

void OpenSim::GeometryPath::setPreScaleLength(const SimTK::State&, double preScaleLength)
{
    _impl->_preScaleLength = preScaleLength;
}


// SHIM/DEPRECATED METHODS
//
// we provide a "shim" implementation of the deprecated API. It isn't logically
// valid, or useful, but providing it here means that other downstream classes
// (i.e. ones that do not use point-based implementations) do not have to
// provide their own shims around a now-deprecated API.
//
// The name of the game here is "crash stabiliity", not "validity". This code
// just needs to ensure that calling code that uses the legacy API doesn't
// segfault or explode when that calling code uses the (deprecated) point-based
// API with a not-point-based (e.g. functional) path implementation. The
// implementations here just need to be crash-stable and (ideally) print out
// a runtime warning that informs the user/developer that they're using shimmed
// (and therefore, probably invalid) code.

static bool emitDeprecationWarning(char const* funcName)
{
    OpenSim::log_warn("deprecation warning: OpenSim tried to call '{}', which is now deprecated. This may result in unusual outputs", funcName);
    return true;
}

const OpenSim::PathPointSet& OpenSim::GeometryPath::getPathPointSet() const
{
    static const bool g_ShownDeprecationWarning = emitDeprecationWarning(__func__);
    return _impl->getOrUpdPPSCached();
}

OpenSim::PathPointSet& OpenSim::GeometryPath::updPathPointSet()
{
    static const bool g_ShownDeprecationWarning = emitDeprecationWarning(__func__);
    return _impl->getOrUpdPPSCached();
}

OpenSim::AbstractPathPoint* OpenSim::GeometryPath::addPathPoint(
        const SimTK::State& s,
        int index,
        const PhysicalFrame& frame)
{
    static const bool g_ShownDeprecationWarning = emitDeprecationWarning(__func__);

    // don't actually add the path point, just return a dummy path point to
    // ensure downstream code reliant on the deprecated API doesn't explode
    // with a segfault
    return _impl->getOrUpdFixupPathPointCached();
}

OpenSim::AbstractPathPoint* OpenSim::GeometryPath::appendNewPathPoint(
        const std::string& proposedName,
        const PhysicalFrame& frame,
        const SimTK::Vec3& locationOnFrame)
{
    static const bool g_ShownDeprecationWarning = emitDeprecationWarning(__func__);

    // don't actually add the path point, just return a dummy path point to
    // ensure downstream code reliant on the deprecated API doesn't explode
    // with a segfault
    return _impl->getOrUpdFixupPathPointCached();
}

bool OpenSim::GeometryPath::canDeletePathPoint(int index)
{
    static const bool g_ShownDeprecationWarning = emitDeprecationWarning(__func__);
    return false;
}

bool OpenSim::GeometryPath::deletePathPoint(const SimTK::State& s, int index)
{
    static const bool g_ShownDeprecationWarning = emitDeprecationWarning(__func__);
    return false;
}

bool OpenSim::GeometryPath::replacePathPoint(
        const SimTK::State& s,
        AbstractPathPoint* oldPathPoint,
        AbstractPathPoint* newPathPoint)
{
    static const bool g_ShownDeprecationWarning = emitDeprecationWarning(__func__);
    return false;
}

const OpenSim::Array<OpenSim::AbstractPathPoint*>& OpenSim::GeometryPath::getCurrentPath(const SimTK::State& s) const
{
    static const bool g_ShownDeprecationWarning = emitDeprecationWarning(__func__);
    return _impl->getOrUpdCurrentPathArray();
}

const OpenSim::PathWrapSet& OpenSim::GeometryPath::getWrapSet() const
{
    static const bool g_ShownDeprecationWarning = emitDeprecationWarning(__func__);
    return _impl->getOrUpdPathWrapSet();
}

OpenSim::PathWrapSet& OpenSim::GeometryPath::updWrapSet()
{
    static const bool g_ShownDeprecationWarning = emitDeprecationWarning(__func__);
    return _impl->getOrUpdPathWrapSet();
}

void OpenSim::GeometryPath::addPathWrap(WrapObject &aWrapObject)
{
    static const bool g_ShownDeprecationWarning = emitDeprecationWarning(__func__);
    return;
}

void OpenSim::GeometryPath::moveUpPathWrap(const SimTK::State& s, int index)
{
    static const bool g_ShownDeprecationWarning = emitDeprecationWarning(__func__);
    return;
}

void OpenSim::GeometryPath::moveDownPathWrap(const SimTK::State& s, int index)
{
    static const bool g_ShownDeprecationWarning = emitDeprecationWarning(__func__);
    return;
}

void OpenSim::GeometryPath::deletePathWrap(const SimTK::State& s, int index)
{
    static const bool g_ShownDeprecationWarning = emitDeprecationWarning(__func__);
    return;
}

void OpenSim::GeometryPath::setLength(const SimTK::State& s, double length) const
{
    static const bool g_ShownDeprecationWarning = emitDeprecationWarning(__func__);
    return;
}

void OpenSim::GeometryPath::setLengtheningSpeed(const SimTK::State& s, double speed) const
{
    static const bool g_ShownDeprecationWarning = emitDeprecationWarning(__func__);
    return;
}

void OpenSim::GeometryPath::updateGeometry(const SimTK::State& s) const
{
    static const bool g_ShownDeprecationWarning = emitDeprecationWarning(__func__);
    return;
}
