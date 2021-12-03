#include "FunctionBasedPath.h"

#include <OpenSim/Common/Component.h>
#include <OpenSim/Simulation/Model/Model.h>
#include <OpenSim/Simulation/SimbodyEngine/Coordinate.h>

#include <sstream>

namespace {
    OpenSim::Exception InvalidFunctionError(std::string parentPath, char const* funcName)
    {
        std::stringstream ss;
        ss << parentPath << ':' << funcName << ": cannot call underlying function: it has not been initialized yet: ensure you have assigned a function to this FunctionBasedPath and finalized its properites";
        return OpenSim::Exception{ss.str()};
    }

    // stub function that will throw an exception with an information message if called
    //
    // used as a stand-in for the object state where the caller has allocated a FunctionBasedPath
    // but hasn't set its underlying function yet.
    class ThrowingSimTKFunction final : public SimTK::Function {
        std::string m_ParentPath;
    public:
        explicit ThrowingSimTKFunction(const std::string& parentPath) :
            m_ParentPath{parentPath}
        {
        }

        double calcValue(const SimTK::Vector&) const override
        {
            throw InvalidFunctionError(m_ParentPath, __func__);
        }

        double calcDerivative(const SimTK::Array_<int>&, const SimTK::Vector&) const override
        {
            throw InvalidFunctionError(m_ParentPath, __func__);
        }

        int getArgumentSize() const override
        {
            throw InvalidFunctionError(m_ParentPath, __func__);
        }

        int getMaxDerivativeOrder() const override
        {
            throw InvalidFunctionError(m_ParentPath, __func__);
        }
    };

    // stub function that will throw an exception with an information message if called
    //
    // used as a stand-in for the object state where the caller has allocated a FunctionBasedPath
    // but hasn't set its underlying function yet.
    class ThrowingOpenSimFunction final : public OpenSim::Function {
        OpenSim_DECLARE_CONCRETE_OBJECT(ThrowingOpenSimFunction, OpenSim::Function)

        std::string m_ParentPath;
    public:
        explicit ThrowingOpenSimFunction(const std::string& parentPath) :
          m_ParentPath{parentPath}
        {
        }

        double calcValue(const SimTK::Vector&) const override
        {
            throw InvalidFunctionError(m_ParentPath, __func__);
        }

        double calcDerivative(const std::vector<int>&, const SimTK::Vector&) const override
        {
            throw InvalidFunctionError(m_ParentPath, __func__);
        }

        int getArgumentSize() const override
        {
            throw InvalidFunctionError(m_ParentPath, __func__);
        }

        int getMaxDerivativeOrder() const override
        {
            throw InvalidFunctionError(m_ParentPath, __func__);
        }

        SimTK::Function* createSimTKFunction() const override
        {
            return new ThrowingSimTKFunction{m_ParentPath};
        }
    };

    void ValidatePropertiesOrThrow(const OpenSim::FunctionBasedPath& parent,
                                   const OpenSim::Property<OpenSim::Function>& funcProp,
                                   const OpenSim::Property<std::string>& coordPathsProp)
    {
        if (funcProp.size() == 0) {
            throw OpenSim::Exception(__FILE__, __LINE__, __func__, parent, "A FunctionBasedPath's function property has not been set");
        }

        if (funcProp.getValue().getMaxDerivativeOrder() < 1) {
            throw OpenSim::Exception(__FILE__, __LINE__, __func__, parent, "The function provided to a FunctionBasedPath is not differentiable: it cannot be used as a path function");
        }

        if (funcProp.getValue().getArgumentSize() != coordPathsProp.size()) {
            std::stringstream ss;
            ss << "The number of coordinate paths provided (" << coordPathsProp.size() << ") does not match the number of arguments the function, " << funcProp.getValue().getName() << ", takes (" << funcProp.getValue().getArgumentSize() << ")";
            throw OpenSim::Exception(__FILE__, __LINE__, __func__, parent, ss.str());
        }

        int nCoords = coordPathsProp.size();
        for (int i = 0; i < nCoords; ++i) {
            const std::string& coordPath = coordPathsProp.getValue(i);
            if (coordPath.empty()) {
                throw OpenSim::Exception(__FILE__, __LINE__, __func__, parent, "An empty coordinate string was provided to a FunctionBasedPath: all coordinate strings must be absolute paths to coordinates within the model");
            }
        }
    }
}

OpenSim::FunctionBasedPath::FunctionBasedPath() : GeometryPath{}
{
    constructProperty_PathFunction(ThrowingOpenSimFunction{this->getAbsolutePathString()});
    constructProperty_Coordinates();
}

OpenSim::FunctionBasedPath::FunctionBasedPath(const FunctionBasedPath&) = default;

OpenSim::FunctionBasedPath::FunctionBasedPath(const OpenSim::Function& func, std::vector<std::string> coordAbsPaths)
{
    constructProperty_PathFunction(func);
    constructProperty_Coordinates();
    for (std::string& coordAbsPath : coordAbsPaths) {
        updProperty_Coordinates().appendValue(std::move(coordAbsPath));
    }

    ValidatePropertiesOrThrow(*this, getProperty_PathFunction(), getProperty_Coordinates());
}

OpenSim::FunctionBasedPath::~FunctionBasedPath() noexcept = default;

OpenSim::FunctionBasedPath& OpenSim::FunctionBasedPath::operator=(const FunctionBasedPath&) = default;

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

    const SimTK::Vector& args = calcFunctionArguments(s);
    double v = getProperty_PathFunction().getValue().calcValue(args);
    setCacheVariableValue(s, _lengthCV, v);
    return v;
}

double OpenSim::FunctionBasedPath::getLengtheningSpeed(const SimTK::State& s) const
{
    if (isCacheVariableValid(s, _speedCV)) {
        return getCacheVariableValue(s, _speedCV);
    }

    // the lengthening speed is the sum of: pathDeriv (w.r.t. coord) * coordinateSpeed
    _derivativeOrderBuffer.resize(1);
    const SimTK::Vector& args = calcFunctionArguments(s);
    int nCoords = getProperty_Coordinates().size();
    double acc = 0.0;
    for (int i = 0; i < nCoords; ++i) {
        const std::string& coordAbsPath = get_Coordinates(i);
        const OpenSim::Coordinate& coord = getRoot().getComponent<OpenSim::Coordinate>(coordAbsPath);
        _derivativeOrderBuffer[0] = i;
        double deriv = get_PathFunction().calcDerivative(_derivativeOrderBuffer, args);
        double coordSpeed = coord.getSpeedValue(s);

        acc = acc + deriv*coordSpeed;
    }
    setCacheVariableValue(s, _speedCV, acc);
    return acc;
}

void OpenSim::FunctionBasedPath::addInEquivalentForces(const SimTK::State& state, double tension, SimTK::Vector_<SimTK::SpatialVec>&, SimTK::Vector& mobilityForces) const
{
    const SimTK::SimbodyMatterSubsystem& matter = getModel().getMatterSubsystem();

    const SimTK::Vector args = calcFunctionArguments(state);

    _derivativeOrderBuffer.resize(1);
    int nCoords = getProperty_Coordinates().size();
    for (int i = 0; i < nCoords; ++i) {
        const std::string& coordAbsPath = get_Coordinates(i);
        const OpenSim::Coordinate& coord = getRoot().getComponent<OpenSim::Coordinate>(coordAbsPath);
        _derivativeOrderBuffer[0] = i;
        double momentArm = get_PathFunction().calcDerivative(_derivativeOrderBuffer, args);

        double torque = -tension*momentArm;

        matter.addInMobilityForce(state,
                                  SimTK::MobilizedBodyIndex(coord.getBodyIndex()),
                                  SimTK::MobilizerUIndex(coord.getMobilizerQIndex()),
                                  torque,
                                  mobilityForces);
    }
}

double OpenSim::FunctionBasedPath::computeMomentArm(const SimTK::State& st, const Coordinate& coord) const
{
    // the moment arm of a path with respect to a coordinate is the path's
    // length derivative with respect to the coordinate

    int coordIndex = indexOfCoordinate(coord);
    if (coordIndex == -1) {
        // the provided coordinate does not affect this path, so it
        // has no moment arm w.r.t. it
        return 0.0;
    }

    _derivativeOrderBuffer.resize(1);
    _derivativeOrderBuffer[0] = coordIndex;

    const SimTK::Vector& args = calcFunctionArguments(st);

    return getProperty_PathFunction().getValue().calcDerivative(_derivativeOrderBuffer, args);
}

void OpenSim::FunctionBasedPath::extendFinalizeFromProperties()
{
    ValidatePropertiesOrThrow(*this, getProperty_PathFunction(), getProperty_Coordinates());
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

void OpenSim::FunctionBasedPath::extendFinalizeConnections(OpenSim::Component& root)
{
    // populate pointer-based coordinate lookups
    //
    // the reason this isn't done in `extendFinalizeFromProperties` is because the
    // not-yet-property-finalized Model hasn't necessarily "connected" to the
    // coordinates that the coordinate files refer to, so the implementation
    // can't lookup the `OpenSim::Coordinate*` pointers during that phase

    // Allow (model) component to include its own subcomponents
    // before calling the base method which automatically invokes
    // connect all the subcomponents.
    {
        Model* model = dynamic_cast<Model*>(&root);
        if (model) {
            connectToModel(*model);
        }
    }

    int nCoords = getProperty_Coordinates().size();
    for (int i = 0; i < nCoords; ++i) {
        root.getComponent(get_Coordinates(i));  // should throw if missing
    }
}

const SimTK::Vector& OpenSim::FunctionBasedPath::calcFunctionArguments(const SimTK::State& st) const
{
    int nargs = getProperty_Coordinates().size();
    _functionArgsBuffer.resize(nargs);
    for (int i = 0; i < nargs; ++i) {
        // HACK: this lookup is horrible
        const std::string& coordName = get_Coordinates(i);
        const OpenSim::Coordinate& coord = getRoot().getComponent<OpenSim::Coordinate>(coordName);
        _functionArgsBuffer[i] = coord.getValue(st);
    }
    return _functionArgsBuffer;
}

int OpenSim::FunctionBasedPath::indexOfCoordinate(const Coordinate& c) const
{
    std::string absPath = c.getAbsolutePathString();

    int nCoords = getProperty_Coordinates().size();
    for (int i = 0; i < nCoords; ++i) {
        const std::string& coordAbsPath = get_Coordinates(i);
        if (coordAbsPath == absPath) {
            return i;
        }
    }
    return -1;
}
