#include <OpenSim/Simulation/Model/FunctionBasedPath.h>

#include <OpenSim/Auxiliary/auxiliaryTestFunctions.h>
#include <OpenSim/Common/Function.h>
#include <OpenSim/Simulation/Model/Model.h>
#include <OpenSim/Simulation/SimbodyEngine/PinJoint.h>

#include <chrono>
#include <random>
#include <stdexcept>
#include <vector>

// test support code
//
// enables the OSIM_TEST macro, and other nice-to-haves
namespace {
    struct Test { char const* suiteName; char const* name; void(*testFunc)(void); };

    std::ostream& operator<<(std::ostream& o, Test const& t)
    {
        return o << t.suiteName << '/' << t.name;
    }

    static std::vector<Test>& GetTestList()
    {
        static std::vector<Test> g_TestList;
        return g_TestList;
    }

    static void RegisterTest(char const* suiteName, char const* name, void(*testFunc)(void))
    {
        GetTestList().push_back(Test{suiteName, name, testFunc});
    }

    #define OSIM_TEST( SuiteName , TestName ) \
        static void Test_ ## SuiteName ## TestName(); \
        static bool g_Add_ ## SuiteName ## TestName ## _ToTestList = [](){ RegisterTest(#SuiteName, #TestName, Test_ ## SuiteName ## TestName); return true; }(); \
        static void Test_ ## SuiteName ## TestName()

    static int RunAllTests()
    {
        std::cerr << "[ ====== ] Running " << GetTestList().size() << " tests\n";

        int numTests = 0;
        int numFailures = 0;
        auto timerAllTestsStart = std::chrono::high_resolution_clock::now();
        for (Test const& t : GetTestList()) {
            ++numTests;
            std::cerr << "[ RUN    ] " << t << '\n';
            auto timerTestStart = std::chrono::high_resolution_clock::now();
            try {
                t.testFunc();
                auto timerTestEnd = std::chrono::high_resolution_clock::now();
                auto timerMillis = std::chrono::duration_cast<std::chrono::milliseconds>(timerTestEnd - timerTestStart);
                std::cerr << "[     OK ] " << t << " (" << timerMillis.count() << " ms)\n";
            } catch (std::exception const& ex) {
                auto timerTestEnd = std::chrono::high_resolution_clock::now();
                auto timerMillis = std::chrono::duration_cast<std::chrono::milliseconds>(timerTestEnd - timerTestStart);
                ++numFailures;
                std::cerr << "[ FAILED ] " << t << " (" << timerMillis.count() << " ms)\n";
                std::cerr << "Exception message: " << ex.what() << '\n';
            }
        }
        auto timerAllTestsEnd = std::chrono::high_resolution_clock::now();
        std::chrono::milliseconds allTestsMillis = std::chrono::duration_cast<std::chrono::milliseconds>(timerAllTestsEnd - timerAllTestsStart);

        std::cerr << "[ ====== ] " << numTests << " tests ran (" << allTestsMillis.count() << " ms total)\n";

        if (numFailures == 0) {
            std::cerr << "All tests (" << numTests << ") passed\n";
            return 0;
        } else {
            std::cerr << "There were " << numFailures << " failed tests\n";
            return 1;
        }
    }

    static double GenerateDouble()
    {
        static std::default_random_engine prng{std::random_device{}()};
        return std::uniform_real_distribution<double>{}(prng);
    }

    static SimTK::Vec3 GenerateRandomVector()
    {
        return SimTK::Vec3{GenerateDouble(), GenerateDouble(), GenerateDouble()};
    }

    static std::vector<std::string> GenerateNPrefixedStrings(int n, std::string prefix = "str")
    {
        std::vector<std::string> rv;
        for (int i = 0; i < n; ++i) {
            std::stringstream ss;
            ss << prefix << i;
            rv.push_back(ss.str());
        }
        return rv;
    }
}

// FunctionBasedPath test support code
//
// enables mocking functions, etc.
namespace {
    template<typename T>
    std::vector<T> ToStdVector(const SimTK::Array_<T>& stkArray)
    {
        std::vector<T> rv;
        rv.reserve(stkArray.size());
        for (const T& v : stkArray) {
            rv.push_back(v);
        }
        return rv;
    }

    template<typename T>
    std::vector<T> ToStdVector(const SimTK::Vector_<T>& stkVec)
    {
        std::vector<T> rv;
        rv.reserve(stkVec.size());
        for (int i = 0; i < stkVec.size(); ++i) {
            rv.push_back(stkVec[i]);
        }
        return rv;
    }

    struct SharedFunctionData final {
        double value = 0.0;
        double derivative = 0.0;
        int argumentSize = 0;
        int maxDerivativeOrder = 1;  // these functions must be differentiable at least once
        std::vector<double> lastCalcValueArg = {};
        std::vector<int> lastCalcDerivativeIntArg = {};
        std::vector<double> lastCalcDerivativeValArgs = {};

        int numTimesCalcValueCalled = 0;
        int numTimesCalcDerivativeCalled = 0;
        int numTimesGetArgumentSizeCalled = 0;
        int numTimesGetMaxDerivativeOrderCalled = 0;
    };

    class MockSimTKFunction final : public SimTK::Function {
    public:
        mutable std::shared_ptr<SharedFunctionData> m_SharedData;

        MockSimTKFunction(std::shared_ptr<SharedFunctionData> sharedData) :
            m_SharedData{sharedData}
        {
        }

        double calcValue(const SimTK::Vector& args) const override
        {
            ++m_SharedData->numTimesCalcValueCalled;
            m_SharedData->lastCalcValueArg = ToStdVector(args);
            return m_SharedData->value;
        }

        double calcDerivative(const SimTK::Array_<int>& intArg, const SimTK::Vector& valArgs) const override
        {
            ++m_SharedData->numTimesCalcDerivativeCalled;
            m_SharedData->lastCalcDerivativeIntArg = ToStdVector(intArg);
            m_SharedData->lastCalcDerivativeValArgs = ToStdVector(valArgs);
            return m_SharedData->derivative;
        }

        int getArgumentSize() const override
        {
            ++m_SharedData->numTimesGetArgumentSizeCalled;
            return m_SharedData->argumentSize;
        }

        int getMaxDerivativeOrder() const override
        {
            ++m_SharedData->numTimesGetMaxDerivativeOrderCalled;
            return m_SharedData->maxDerivativeOrder;
        }

        SimTK::Function* clone() const override
        {
            return new MockSimTKFunction{*this};
        }
    };

    // used to test that the FunctionBasedPath is forwarding things correctly
    class MockPathFunction : public OpenSim::Function {
        OpenSim_DECLARE_CONCRETE_OBJECT(MockPathFunction, OpenSim::Function);

    public:
        mutable std::shared_ptr<SharedFunctionData> m_SharedData = std::make_shared<SharedFunctionData>();

        double calcValue(const SimTK::Vector& args) const override
        {
            ++m_SharedData->numTimesCalcValueCalled;
            m_SharedData->lastCalcValueArg = ToStdVector(args);
            return m_SharedData->value;
        }

        double calcDerivative(const std::vector<int>& intArgs, const SimTK::Vector& valArgs) const override
        {
            ++m_SharedData->numTimesCalcDerivativeCalled;
            m_SharedData->lastCalcDerivativeIntArg = intArgs;
            m_SharedData->lastCalcDerivativeValArgs = ToStdVector(valArgs);
            return m_SharedData->derivative;
        }

        int getArgumentSize() const override
        {
            ++m_SharedData->numTimesGetArgumentSizeCalled;
            return m_SharedData->argumentSize;
        }

        int getMaxDerivativeOrder() const override
        {
            ++m_SharedData->numTimesGetMaxDerivativeOrderCalled;
            return m_SharedData->maxDerivativeOrder;
        }

        SimTK::Function* createSimTKFunction() const override
        {
            return new MockSimTKFunction{m_SharedData};
        }
    };

    // helper struct for generating a model with N coordinates
    struct ModelWithNCoordinates {
        OpenSim::Model model;
        std::vector<std::string> coordinateAbsPaths;
    };

    ModelWithNCoordinates GenerateModelWithNCoordinates(int n)
    {
        ModelWithNCoordinates rv;

        for (int i = 0; i < n; ++i) {
            // create a dummy body the joint is joining to ground
            OpenSim::Body* body = new OpenSim::Body{};
            body->setName(std::string{"body_"} + std::to_string(i));
            body->setMass(1.0);
            rv.model.addBody(body);

            // create the joint and connect it to the body + ground
            OpenSim::PinJoint* joint = new OpenSim::PinJoint{};
            joint->setName(std::string{"pinjoint_"} + std::to_string(i));
            joint->connectSocket_parent_frame(rv.model.getGround());
            joint->connectSocket_child_frame(*body);
            rv.model.addComponent(joint);

            // save the coordinate's path (used by the FBP API)
            rv.coordinateAbsPaths.push_back(joint->getCoordinate().getAbsolutePathString());
        }

        return rv;
    }
}


// construction tests

OSIM_TEST(FunctionBasedPath, CanBeDefaultConstructedWithoutThrowing)
{
    OpenSim::FunctionBasedPath fbp;  // shouldn't throw
}

OSIM_TEST(FunctionBasedPath, CanBeCopyConstructedWithoutThrowing)
{
    OpenSim::FunctionBasedPath fbp1;
    OpenSim::FunctionBasedPath fbp2{fbp1};  // shouldn't throw
}

OSIM_TEST(FunctionBasedPath, CanBeMoveConstructedWithoutThrowing)
{
    OpenSim::FunctionBasedPath fbp1;
    OpenSim::FunctionBasedPath fbp2{std::move(fbp1)};  // shouldn't throw
}

OSIM_TEST(FunctionBasedPath, CanBeConstructedWithA0ArgPathFunctionAndEmptyCoordinateListWithoutThrowing)
{
    MockPathFunction pathFn;
    std::vector<std::string> coords;
    OpenSim::FunctionBasedPath fbp{pathFn, coords};  // shouldn't throw
}

OSIM_TEST(FunctionBasedPath, CanBeConstructedWithAnNArityFunctionAndNCoordinates)
{
    // function-based paths have to support being parameterized with a varying
    // number of coordinates - not just 0, 1, or 2 (etc.)

    for (int i = 0; i < 10; ++i) {
        MockPathFunction fn;
        fn.m_SharedData->argumentSize = i;
        std::vector<std::string> coords = GenerateNPrefixedStrings(i);

        OpenSim::FunctionBasedPath fbp{fn, coords};  // shouldn't throw
    }
}

OSIM_TEST(FunctionBasedPath, ConstructingWithDifferentArityAndNumberOfCoordinatesThrows)
{
    // the arity (number of arguments) the function accepts should be equal to the
    // number of coordinates the caller specifies
    //
    // this is because, at simulation time, the FunctionBasedPath pumps the coordinate
    // values into the function as arguments. Effectively, the generic function is
    // parameterized by coordinate values.

    MockPathFunction fn;
    fn.m_SharedData->argumentSize = 3;
    std::vector<std::string> coords = { "onlyone" };

    ASSERT_THROW(OpenSim::Exception, OpenSim::FunctionBasedPath fbp(fn, coords));
}

OSIM_TEST(FunctionBasedPath, ConstructingWithAFunctionThatHasNoDerivativeThrows)
{
    // the provided function must be differentiable
    //
    // this is because the function needs to provide 0-order results (i.e. the
    // length of the path) and 1st-order results (i.e. the lengthening speed of
    // the path)

    MockPathFunction fn;
    fn.m_SharedData->maxDerivativeOrder = 0;  // non-differentiable function

    ASSERT_THROW(OpenSim::Exception, OpenSim::FunctionBasedPath fbp(fn, {}));
}

OSIM_TEST(FunctionBasedPath, ConstructingWithAnEmptyCoordNameThrows)
{
    // the provided coordinate names should be non-empty
    //
    // this is a sanity-check by the implementation. It's entirely feasible
    // that some external process (e.g. a file reader) pumps a blank line
    // or whitespace string into this object's constructor - we try to
    // aggressively find + report it early, rather than during some later
    // runtime lookup

    MockPathFunction fn;
    fn.m_SharedData->argumentSize = 1;

    ASSERT_THROW(OpenSim::Exception, OpenSim::FunctionBasedPath fbp(fn, {""}));
}

OSIM_TEST(FunctionBasedPath, WhenConstructedWithPathFunctionUsesTheFunctionInLengthEvaluation)
{
    // this is partially an implementation detail - you can remove this test if you need
    //
    // the test is providing a 0-arity, 0-coordinate function to the path, but this
    // test is ensuring that this (to be fair, bizzare) edge-case still ends up
    // calling `getLength` with 0 args and ultimately returns some value
    //
    // see later API tests for how `getLength` etc. *should* be used: this is
    // mostly just a sanity check on the "construct with a function" pattern that
    // the test suite uses (prod uses are likely to just set property values on
    // loading an OSIM file or something)

    OpenSim::Model model;

    MockPathFunction pathFn;
    std::vector<std::string> coords;
    OpenSim::FunctionBasedPath* fbp = new OpenSim::FunctionBasedPath{pathFn, coords};

    model.addComponent(fbp);

    SimTK_TEST(pathFn.m_SharedData->numTimesCalcValueCalled == 0);

    SimTK::State& s = model.initSystem();

    SimTK_TEST(pathFn.m_SharedData->numTimesCalcValueCalled == 0);

    SimTK_TEST(fbp->getLength(s) == pathFn.m_SharedData->value);

    SimTK_TEST(pathFn.m_SharedData->numTimesCalcValueCalled == 1);
}


// assignment tests

OSIM_TEST(FunctionBasedPath, CanBeCopyAssignedWithoutThrowing)
{
    OpenSim::FunctionBasedPath fbp1;
    OpenSim::FunctionBasedPath fbp2;
    fbp2 = fbp1;
}

OSIM_TEST(FunctionBasedPath, CanBeMoveAssignedWithoutThrowing)
{
    OpenSim::FunctionBasedPath fbp1;
    OpenSim::FunctionBasedPath fbp2;
    fbp2 = std::move(fbp1);
}


// api tests

OSIM_TEST(FunctionBasedPath, FinalizeFromPropertiesWithoutSettingFunctionThrows)
{
    // calling finalizeFromProperties implies the FBP's properties are correct
    //
    // it is incorrect to leave the Function property not-set

    OpenSim::FunctionBasedPath fbp;

    ASSERT_THROW(OpenSim::Exception, fbp.finalizeFromProperties());
}

OSIM_TEST(FunctionBasedPath, FinalizeFromPropertiesWithNonDifferentiableFunctionThrows)
{
    // calling finalizeFromProperties implies the FBP's properties are correct
    //
    // it is incorrect to set the Function property to a function that cannot
    // be differentiated (an API requirement)

    OpenSim::FunctionBasedPath fbp;

    MockPathFunction fn;
    fn.m_SharedData->maxDerivativeOrder = 0;  // non-differentiable

    auto& prop = dynamic_cast<OpenSim::Property<OpenSim::Function>&>(fbp.updPropertyByName("PathFunction"));
    prop.setValue(fn);

    ASSERT_THROW(OpenSim::Exception, fbp.finalizeFromProperties());
}

OSIM_TEST(FunctionBasedPath, FinalizeFromPropertiesWithInvalidFunctionArityThrows)
{
    // calling finalizeFromProperties implies the FBP's properties are correct
    //
    // it is incorrect to set the Function property to a function that has a
    // different arity (number of arguments) from the number of coordinates
    // defined in the coordinates property

    OpenSim::FunctionBasedPath fbp;

    MockPathFunction fn;
    fn.m_SharedData->argumentSize = 10;

    auto& prop = dynamic_cast<OpenSim::Property<OpenSim::Function>&>(fbp.updPropertyByName("PathFunction"));
    prop.setValue(fn);

    ASSERT_THROW(OpenSim::Exception, fbp.finalizeFromProperties());
}

OSIM_TEST(FunctionBasedPath, FinalizeFromPropertiesWithEmptyCoordinateStringThrows)
{
    // calling finalizeFromProperties implies the FBP's properties are correct
    //
    // it is incorrect to provide an empty coordinate string as a property (every
    // string should be a nonempty absolute path to a coordinate)

    OpenSim::FunctionBasedPath fbp;

    MockPathFunction fn;
    fn.m_SharedData->argumentSize = 1;
    auto& funcProp = dynamic_cast<OpenSim::Property<OpenSim::Function>&>(fbp.updPropertyByName("PathFunction"));
    funcProp.setValue(fn);

    auto& coordProp = dynamic_cast<OpenSim::Property<std::string>&>(fbp.updPropertyByName("Coordinates"));
    coordProp.adoptAndAppendValue(new std::string{""});

    ASSERT_THROW(OpenSim::Exception, fbp.finalizeFromProperties());
}

OSIM_TEST(FunctionBasedPath, FinalizeFromPropertiesWithCorrectPropertiesDoesNotThrow)
{
    // this is a throwaway test that ensures the sane-case doesn't throw

    OpenSim::FunctionBasedPath fbp;

    MockPathFunction fn;
    fn.m_SharedData->argumentSize = 0;
    auto& funcProp = dynamic_cast<OpenSim::Property<OpenSim::Function>&>(fbp.updPropertyByName("PathFunction"));
    funcProp.setValue(fn);

    fbp.finalizeFromProperties();  // shouldn't throw
}

OSIM_TEST(FunctionBasedPath, FinalizeConnectionsWithCorrectPropertiesDoesNotThrow)
{
    // calling `finalizeConnections` should be fine if every coordinate path
    // provided to the FBP constructor is actually a real coordinate

    ModelWithNCoordinates m = GenerateModelWithNCoordinates(5);

    MockPathFunction fn;
    fn.m_SharedData->argumentSize = static_cast<int>(m.coordinateAbsPaths.size());

    OpenSim::FunctionBasedPath* fbp = new OpenSim::FunctionBasedPath{fn, m.coordinateAbsPaths};
    m.model.addComponent(fbp);

    m.model.finalizeFromProperties();
    m.model.finalizeConnections();
}

OSIM_TEST(FunctionBasedPath, FinalizeConnectionsWithIncorrectCoordinatePathThrows)
{
    // calling `finalizeConnections` should throw with a (hopefully) handy error
    // message at conneciton-time if an FBP contains an invalid coordinate path

    ModelWithNCoordinates m = GenerateModelWithNCoordinates(5);

    MockPathFunction fn;
    fn.m_SharedData->argumentSize = static_cast<int>(m.coordinateAbsPaths.size());

    auto paths = m.coordinateAbsPaths;
    paths[2] = "not-a-real-coordinate-path";

    OpenSim::FunctionBasedPath* fbp = new OpenSim::FunctionBasedPath{fn, paths};
    m.model.addComponent(fbp);

    m.model.finalizeFromProperties();
    ASSERT_THROW(OpenSim::Exception, m.model.finalizeConnections());
}

OSIM_TEST(FunctionBasedPath, CanGetColorWithoutThrowing)
{
    OpenSim::Model model;
    OpenSim::FunctionBasedPath* fbp = new OpenSim::FunctionBasedPath{MockPathFunction{}, {}};
    model.addComponent(fbp);

    SimTK::State& s = model.initSystem();

    fbp->getColor(s);
}

OSIM_TEST(FunctionBasedPath, GetColorReturnsDefaultColorIfSetColorIsNotCalled)
{
    SimTK::Vec3 randomColor = GenerateRandomVector();

    OpenSim::Model model;
    OpenSim::FunctionBasedPath* fbp = new OpenSim::FunctionBasedPath{MockPathFunction{}, {}};
    fbp->setDefaultColor(randomColor);
    model.addComponent(fbp);

    SimTK::State& s = model.initSystem();

    SimTK_TEST(fbp->getColor(s) == randomColor);
}

OSIM_TEST(FunctionBasedPath, SetColorSetsTheColorInTheState)
{
    SimTK::Vec3 randomColor = GenerateRandomVector();

    OpenSim::Model model;
    OpenSim::FunctionBasedPath* fbp = new OpenSim::FunctionBasedPath{MockPathFunction{}, {}};
    fbp->setDefaultColor(randomColor);
    model.addComponent(fbp);

    SimTK::State& s = model.initSystem();

    fbp->setColor(s, randomColor);

    SimTK_TEST(fbp->getColor(s) == randomColor);
}

OSIM_TEST(FunctionBasedPath, GetLengthUsesUnderlyingFunctionCalcValue)
{
    // the FunctionBasedPath ultimately gets its result values from the underlying
    // Function object, rather than generating them itself

    ModelWithNCoordinates m = GenerateModelWithNCoordinates(1);
    OpenSim::Coordinate const& coord = m.model.getComponent<OpenSim::Coordinate>(m.coordinateAbsPaths.at(0));

    MockPathFunction fn;
    fn.m_SharedData->argumentSize = 1;
    std::vector<std::string> coordPaths = {coord.getAbsolutePathString()};

    OpenSim::FunctionBasedPath* fbp = new OpenSim::FunctionBasedPath{fn, coordPaths};
    m.model.addComponent(fbp);

    SimTK_TEST(fn.m_SharedData->numTimesCalcValueCalled == 0);

    SimTK::State& s = m.model.initSystem();

    SimTK_TEST(fbp->getLength(s) == fn.m_SharedData->value);
    SimTK_TEST(fn.m_SharedData->numTimesCalcValueCalled == 1);
    SimTK_TEST(fn.m_SharedData->lastCalcValueArg.size() == 1);
    SimTK_TEST(fn.m_SharedData->lastCalcValueArg == std::vector<double>{coord.getValue(s)});
}

OSIM_TEST(FunctionBasedPath, GetLengthIsCached)
{
    // FunctionBasedPath::getLength should be "cached" depending on the simulation
    // stage.
    //
    // this is an implementation detail that's good for implementors to know. Remove
    // this test if it's a bad assumption.

    ModelWithNCoordinates m = GenerateModelWithNCoordinates(1);

    MockPathFunction fn;
    fn.m_SharedData->argumentSize = static_cast<int>(m.coordinateAbsPaths.size());
    OpenSim::FunctionBasedPath* fbp = new OpenSim::FunctionBasedPath{fn, m.coordinateAbsPaths};
    m.model.addComponent(fbp);

    SimTK_TEST(fn.m_SharedData->numTimesCalcValueCalled == 0);

    SimTK::State& s = m.model.initSystem();
    m.model.realizePosition(s);

    SimTK_TEST(fn.m_SharedData->numTimesCalcValueCalled == 0);
    SimTK_TEST(fbp->getLength(s) == fn.m_SharedData->value);
    SimTK_TEST(fn.m_SharedData->numTimesCalcValueCalled == 1);
    fbp->getLength(s);  // should cache
    SimTK_TEST(fn.m_SharedData->numTimesCalcValueCalled == 1);
}

OSIM_TEST(FunctionBasedPath, GetLengtheningSpeedUsesUnderlyingFunctionCalcValue)
{
    // the FunctionBasedPath ultimately calculates its lengthening speed
    // by using the underlying Function object, rather than computing it
    // itself

    ModelWithNCoordinates m = GenerateModelWithNCoordinates(1);

    MockPathFunction fn;
    fn.m_SharedData->argumentSize = static_cast<int>(m.coordinateAbsPaths.size());
    OpenSim::FunctionBasedPath* fbp = new OpenSim::FunctionBasedPath{fn, m.coordinateAbsPaths};
    m.model.addComponent(fbp);

    SimTK_TEST(fn.m_SharedData->numTimesCalcDerivativeCalled == 0);

    SimTK::State& s = m.model.initSystem();
    m.model.realizeAcceleration(s);

    double lengtheningSpeed = fbp->getLengtheningSpeed(s);

    SimTK_TEST(fn.m_SharedData->numTimesCalcDerivativeCalled == 1);
    SimTK_TEST(lengtheningSpeed == fn.m_SharedData->derivative);
    SimTK_TEST(fn.m_SharedData->lastCalcDerivativeIntArg == std::vector<int>{0});
}

OSIM_TEST(FunctionBasedPath, GetLengtheningSpeedIsCached)
{
    // FunctionBasedPath::getLengtheningSpeed should be "cached", depending on
    // the simulation stage.
    //
    // this is an implementation detail that's good for implementors to know. Remove
    // this test if it's a bad assumption.

    ModelWithNCoordinates m = GenerateModelWithNCoordinates(1);

    MockPathFunction fn;
    fn.m_SharedData->argumentSize = static_cast<int>(m.coordinateAbsPaths.size());
    OpenSim::FunctionBasedPath* fbp = new OpenSim::FunctionBasedPath{fn, m.coordinateAbsPaths};
    m.model.addComponent(fbp);

    SimTK_TEST(fn.m_SharedData->numTimesCalcDerivativeCalled == 0);

    SimTK::State& s = m.model.initSystem();
    m.model.realizeAcceleration(s);

    SimTK_TEST(fn.m_SharedData->numTimesCalcDerivativeCalled == 0);
    SimTK_TEST(fbp->getLengtheningSpeed(s) == fn.m_SharedData->value);
    SimTK_TEST(fn.m_SharedData->numTimesCalcDerivativeCalled == 1);
    fbp->getLengtheningSpeed(s);  // should cache
    SimTK_TEST(fn.m_SharedData->numTimesCalcDerivativeCalled == 1);
}

OSIM_TEST(FunctionBasedPath, AddInEquivalentForcesOnValidFBPDoesNotThrow)
{
    // calling `addInEquivalentForces` on a fully-initialized function-based
    // path shouldn't throw, and should probably add *something* to the mobility
    // forces outparam

    ModelWithNCoordinates m = GenerateModelWithNCoordinates(5);
    MockPathFunction fn;
    fn.m_SharedData->argumentSize = static_cast<int>(m.coordinateAbsPaths.size());  // all coords affect this path

    OpenSim::FunctionBasedPath* fbp = new OpenSim::FunctionBasedPath{fn, m.coordinateAbsPaths};
    m.model.addComponent(fbp);

    SimTK::State& s = m.model.initSystem();
    m.model.realizeAcceleration(s);

    double tension = 5.0;
    SimTK::Vector_<SimTK::SpatialVec> bodyForces;
    SimTK::Vector mobilityForces;
    mobilityForces.resize(100); // required because impl doesn't actually add entries (no range checking) it's assigning a force per body index or something

    fbp->addInEquivalentForces(s, tension, bodyForces, mobilityForces);  // shouldn't throw
}

OSIM_TEST(FunctionBasedPath, AddInEquivalentForcesOnEmptyCoordlistWorks)
{
    // edge-case: you can still add in equivalent forces (i.e. nothing) when
    // using an FBP that isn't parameterized by any Coordinates

    ModelWithNCoordinates m = GenerateModelWithNCoordinates(5);
    MockPathFunction fn;
    fn.m_SharedData->argumentSize = 0;  // no coords affect this path

    OpenSim::FunctionBasedPath* fbp = new OpenSim::FunctionBasedPath{fn, {}};
    m.model.addComponent(fbp);

    SimTK::State& s = m.model.initSystem();
    m.model.realizeAcceleration(s);

    double tension = 5.0;
    SimTK::Vector_<SimTK::SpatialVec> bodyForces;
    SimTK::Vector mobilityForces;
    mobilityForces.resize(100); // required because impl doesn't actually add entries (no range checking) it's assigning a force per body index or something

    fbp->addInEquivalentForces(s, tension, bodyForces, mobilityForces);  // shouldn't throw
}

OSIM_TEST(FunctionBasedPath, ComputeMomentArmReturns0IfCalledWithNonAffectingCoord)
{
    // `computeMomentArm` may be called with a coordinate that does not actually
    // affect the given FBP
    //
    // in that case, the returned moment arm shall always be 0.0

    ModelWithNCoordinates m = GenerateModelWithNCoordinates(5);
    MockPathFunction fn;
    fn.m_SharedData->argumentSize = 4;  // the first 4 coords affect the path

    auto affectingCoords = m.coordinateAbsPaths;
    affectingCoords.resize(affectingCoords.size()-1);
    OpenSim::FunctionBasedPath* fbp = new OpenSim::FunctionBasedPath{fn, affectingCoords};
    m.model.addComponent(fbp);

    SimTK::State& s = m.model.initSystem();
    m.model.realizeAcceleration(s);

    const auto& nonAffectingCoord = m.model.getComponent<OpenSim::Coordinate>(m.coordinateAbsPaths.back());

    SimTK_TEST(fbp->computeMomentArm(s, nonAffectingCoord) == 0);
}

OSIM_TEST(FunctionBasedPath, ComputeMomentArmCallsCalcDerivOfFunction)
{
    // this is an implementation detail - remove this test if the detail becomes non-true
    //
    // computing a moment arm should use the path's derivative to compute the
    // moment arm. This test just tries to ensure that that's what's happening,
    // as a sanity-check on the underlying plumbing

    ModelWithNCoordinates m = GenerateModelWithNCoordinates(5);
    MockPathFunction fn;
    fn.m_SharedData->derivative = GenerateDouble();
    fn.m_SharedData->argumentSize = m.coordinateAbsPaths.size();  // all coords affect the path

    OpenSim::FunctionBasedPath* fbp = new OpenSim::FunctionBasedPath{fn, m.coordinateAbsPaths};
    m.model.addComponent(fbp);

    SimTK::State& s = m.model.initSystem();
    m.model.realizeAcceleration(s);

    SimTK_TEST(fn.m_SharedData->numTimesCalcDerivativeCalled == 0);
    int coordIdx = 0;
    double momentArm = fbp->computeMomentArm(s, m.model.getComponent<OpenSim::Coordinate>(m.coordinateAbsPaths[coordIdx]));
    SimTK_TEST(fn.m_SharedData->numTimesCalcDerivativeCalled == 1);
    SimTK_TEST(momentArm == fn.m_SharedData->derivative);
}




// HACK: test Joris's implementation here

#include <OpenSim/Simulation/Model/PointBasedPath.h>


// core runtime algorithm
namespace joris {

    static constexpr size_t g_MaxNumDimensions = 8;  // important: this is an upper limit that's used for stack allocations

    // an `nPoints` evenly-spaced discretization of the range [begin, end] (inclusive)
    //
    // e.g. [0, 10], 0 points == [], step size = 10
    //             , 1 point  == [0], step size = 10
    //             , 2 points == [0, 10], step size = 10
    //             , 3 points == [0, 5, 10], step size = 5
    //             , 4 points == [0, 3.33, 6.66, 10], step size = 3.33
    struct Discretization final {
        double begin;
        double end;
        int nPoints;
    };

    double stepSize(const Discretization& d)
    {
        double diff = d.end - d.begin;
        return d.nPoints <= 2 ? diff : diff/(d.nPoints-1);
    }

    double realIndexOf(const Discretization& d, double v)
    {
        // solve: `x = d.begin + stepSize(d)*n` for `n`
        return (v - d.begin) / stepSize(d);
    }

    std::pair<int, double> splitIntoWholeAndFractional(double v)
    {
        double wholePart;
        double fractionalPart = std::modf(v, &wholePart);
        return {static_cast<int>(wholePart), fractionalPart};
    }

    std::pair<int, double> splitIntoWholeAndFractionalForBetaCalc(const Discretization& d, double x)
    {
        std::pair<int, double> wholeAndFrac = splitIntoWholeAndFractional(realIndexOf(d, x));

        if (wholeAndFrac.first < 1) {
            // edge-case: x is below the second graduation, use the 2nd discretization to ensure the fitting
            // impl has access to all 4 (1 before, 2 after) datapoints
            return {1, 0.0};
        } else if (wholeAndFrac.first > d.nPoints-3) {
            // edge-case: x is greater than the third-to-last graduation, use the third-to-last discretization
            // to ensurethe fitting impl. has asscess to all 4 (1 before, 2 after) datapoints
            return {d.nPoints-3, 0.0};
        } else {
            return wholeAndFrac;
        }
    }

    using Polynomial = std::array<double, 4>;

    Polynomial ComputeBeta(double frac)
    {
        // compute polynomial based on fraction the point is toward the next point

        double frac2 = frac*frac;
        double frac3 = frac2*frac;
        double frac4 = frac3*frac;
        double fracMinusOne = frac - 1;
        double fracMinusOne3 = fracMinusOne*fracMinusOne*fracMinusOne;

        Polynomial p;
        p[0] =  0.5 * fracMinusOne3 * frac*(2.0*frac + 1.0);
        p[1] = -0.5 * fracMinusOne  * (6.0*frac4 - 9.0*frac3 + 2.0*frac + 2.0);
        p[2] =  0.5 * frac          * (6.0*frac4 - 15.0*frac3 + 9.0*frac2 + frac + 1.0);
        p[3] = -0.5 * fracMinusOne  * frac3*(2.0*frac - 3.0);

        return p;
    }

    Polynomial ComputeBetaDerivative(double frac)
    {
        // compute polynomial derivative based on fraction the point is toward the next point

        double frac2 = frac*frac;
        double frac3 = frac2*frac;
        double frac4 = frac3*frac;

        Polynomial p;
        p[0] = 5*frac4 - 10*frac3 + 4.5*frac2 + frac - 0.5;
        p[1] = -15*frac4 + 30*frac3 - 13.5*frac2 - 2*frac;
        p[2] = 15*frac4 - 30*frac3 + 13.5*frac2 + frac + 0.5;
        p[3] = frac2*(-5*frac2 + 10*frac - 4.5);
        return p;
    }

    // computes the interpolated Y value for a given permutation of X values
    //
    // this is the "heart" of the FPB algorithm. It's loosely based on the algorithm
    // described here:
    //
    //     "Two hierarchies of spline interpolations. Practical algorithms for multivariate higher order splines"
    //     https://arxiv.org/abs/0905.3564
    //
    // `xVals` are the current values of each dimension (independent var: e.g. coordinate values - the "x"es)
    double Impl_ComputeValueOrDeriv(
            const std::vector<Discretization>& fittingDiscretizations,
            const std::vector<double>& fittingEvals,
            const SimTK::Vector& xVals,
            int derivIndex = -1)
    {
        SimTK_ASSERT_ALWAYS(!fittingDiscretizations.empty(), "this implementation requires that at least one X (e.g. coordinate) affects the path");
        SimTK_ASSERT_ALWAYS(xVals.size() == static_cast<int>(fittingDiscretizations.size()), "You must call this function with the correct number of (precomputed) x values");
        SimTK_ASSERT_ALWAYS(xVals.size() < static_cast<int>(g_MaxNumDimensions), "too many dimensions in this fit - the implementation cannot handle this many");

        // compute:
        //
        // - the index of the first discretization step *before* the input value
        //
        // - the polynomial of the curve at that step, given its fractional distance
        //   toward the next step
        std::array<int, g_MaxNumDimensions> closestDiscretizationSteps;
        std::array<Polynomial, g_MaxNumDimensions> betas;
        for (int dim = 0; dim < xVals.size(); ++dim) {
            std::pair<int, double> wholeFrac = splitIntoWholeAndFractionalForBetaCalc(fittingDiscretizations[dim], xVals[dim]);
            closestDiscretizationSteps[dim] = wholeFrac.first;
            betas[dim] = dim != derivIndex ? ComputeBeta(wholeFrac.second) : ComputeBetaDerivative(wholeFrac.second);
        }

        // for each dim, permute through 4 (discretized) locations *around* the x location:
        //
        // - A one step before B
        // - B the first discretization step before the input value
        // - C one step after B
        // - D one step after C
        //
        // where:
        //
        // - betas are coefficients that affect each location. beta[0] affects A,
        //   betas[1] affects B, betas[2] affects C, and betas[3] affects D

        // represent permuting through each location around each x as a sequence
        // of integer offsets that can be -1, 0, 1, or 2
        //
        // the algorithm increments this array as it goes through each permutation
        std::array<int, g_MaxNumDimensions> dimIdxOffsets;
        for (int dim = 0; dim < xVals.size(); ++dim) {
            dimIdxOffsets[dim] = -1;
        }

        // permute through all locations around the input value
        //
        // e.g. the location permutations for 3 dims iterate like this for each
        //      crank of the loop
        //
        //     [-1, -1, -1]
        //     [-1, -1,  0]
        //     [-1, -1,  1]
        //     [-1, -1,  2]
        //     [-1,  0, -1]
        //     ...(4^3 steps total)...
        //     [ 2,  2,  1]
        //     [ 2,  2,  2]
        //     [ 3,  0,  0]   (the termination condition)

        double z = 0.0;  // result
        int cnt = 0;  // sanity-check counter

        // seems to be used in original implementation to scale the eval step size in the deriv calculation specifically
        double evalScaler = derivIndex == -1 ? 1.0 : 1.0/stepSize(fittingDiscretizations[derivIndex]);

        while (dimIdxOffsets[0] < 3) {

            // compute `beta` (weighted coefficient per dimension) for this particular
            // permutation's x locations (e.g. -1, 0, 0, 2) and figure out
            // what the closest input value was at the weighted location. Add the
            // result the the output

            double beta = 1.0;
            int strideInFittingEvals = 1;
            int idxInFittingEvals = 0;

            // go backwards, from least-significant dim (highest idx) to figure
            // out where the evaluation is in the fittingEvals array
            //
            // this is so that we can compute the stride as the algorithm runs
            for (int coord = xVals.size()-1; coord >= 0; --coord) {
                int offset = dimIdxOffsets[coord];  // -1, 0, 1, or 2
                int closestStep = closestDiscretizationSteps[coord];
                int step = closestStep + offset;

                beta *= betas[coord][offset+1];  // index into the polynomial
                idxInFittingEvals += strideInFittingEvals * step;
                strideInFittingEvals *= fittingDiscretizations[coord].nPoints;
            }

            // equivalent to z += b*v, but handles rounding errors when the rhs
            // is very small
            z = std::fma(beta, evalScaler * fittingEvals.at(idxInFittingEvals), z);

            // increment the offsets
            //
            // this is effectively the step that permutes [-1, -1, 2] --> [-1,  0, -1]
            {
                int pos = xVals.size()-1;
                ++dimIdxOffsets[pos];  // perform least-significant increment (may overflow)
                while (pos > 0 && dimIdxOffsets[pos] > 2) {  // handle overflows + carry propagation
                    dimIdxOffsets[pos] = -1;  // overflow
                    ++dimIdxOffsets[pos-1];  // carry propagation
                    --pos;
                }
            }

            ++cnt;
        }

        // sanity check: is `z` accumulated from the expected number of iterations?
        {
            int expectedIterations = 1 << (2*xVals.size());
            if (cnt != expectedIterations) {
                std::stringstream msg;
                msg << "invalid number of permutations explored: expected = " << expectedIterations << ", got = " << cnt;
                OPENSIM_THROW(OpenSim::Exception, std::move(msg).str());
            }
        }

        return z;
    }
}


// TODO: test the low-level algorithm by providing it low-level data etc
namespace joris {

}


// wireup of core runtime algorithm to OpenSim
namespace joris {

    // the underlying SimTK::Function that actually calls into the function
    class JorisPathSimTKFunction : public SimTK::Function {
        std::shared_ptr<std::vector<Discretization>> _discretizations;
        std::shared_ptr<std::vector<double>> _evaluations;
    public:
        JorisPathSimTKFunction(
                std::shared_ptr<std::vector<Discretization>> discretizations,
                std::shared_ptr<std::vector<double>> evaluations) :
            _discretizations{std::move(discretizations)},
            _evaluations{std::move(evaluations)}
        {
        }

        double calcValue(const SimTK::Vector& coordVals) const override
        {
            return Impl_ComputeValueOrDeriv(*_discretizations, *_evaluations, coordVals);
        }

        double calcDerivative(const SimTK::Array_<int>& derivComponents, const SimTK::Vector& coordVals) const override
        {
            SimTK_ASSERT_ALWAYS(derivComponents.size() == 1, "Can only find first-order derivative w.r.t. one coord");
            return Impl_ComputeValueOrDeriv(*_discretizations, *_evaluations, coordVals, derivComponents[0]);
        }

        int getArgumentSize() const override
        {
            return static_cast<int>(_discretizations->size());
        }

        int getMaxDerivativeOrder() const override
        {
            return 1;
        }
    };

    using namespace OpenSim;  // required by the property macros...

    class FunctionBasedPathDiscretization : public OpenSim::Component {
        OpenSim_DECLARE_CONCRETE_OBJECT(FunctionBasedPathDiscretization, OpenSim::Component);

    public:
        OpenSim_DECLARE_PROPERTY(coordinate_abspath, std::string, "The absolute path, in the model, to the OpenSim::Coordinate that this discretization was produced from");
        OpenSim_DECLARE_PROPERTY(x_begin, double, "The lowest OpenSim::Coordinate value that was used for discretization");
        OpenSim_DECLARE_PROPERTY(x_end, double, "The highest OpenSim:::Coordinate value that was used for discretization");
        OpenSim_DECLARE_PROPERTY(num_points, int, "The number of evenly-spaced OpenSim::Coordinate values between [x_begin, x_end] (inclusive) that were used for discretization of the OpenSim::Coordinate. E.g. [x_begin, 1*(x_begin+((x_end-x_begin)/3)), 2*(x_begin+((x_end-x_begin)/3)), x_end]");

        FunctionBasedPathDiscretization()
        {
            constructProperty_coordinate_abspath("");
            constructProperty_x_begin(0.0);
            constructProperty_x_end(0.0);
            constructProperty_num_points(0);
        }
    };

    class FunctionBasedPathDiscretizationSet : public OpenSim::Set<FunctionBasedPathDiscretization> {
        OpenSim_DECLARE_CONCRETE_OBJECT(FunctionBasedPathDiscretizationSet, OpenSim::Set<FunctionBasedPathDiscretization>);
    };

    class JorisPathFunction final : public OpenSim::Function {
        OpenSim_DECLARE_CONCRETE_OBJECT(JorisPathFunction, OpenSim::Function);
        // TODO: this needs to have PROPERTYs to store the discretizations + evaluations
    public:
        std::shared_ptr<std::vector<Discretization>> _discretizations;
        std::shared_ptr<std::vector<double>> _evaluations;
    public:
        JorisPathFunction() :
            _discretizations{std::make_shared<std::vector<Discretization>>()},
            _evaluations{std::make_shared<std::vector<double>>()}
        {
        }

        JorisPathFunction(
                std::shared_ptr<std::vector<Discretization>> discretizations,
                std::shared_ptr<std::vector<double>> evaluations) :
            _discretizations{std::move(discretizations)},
            _evaluations{std::move(evaluations)}
        {
        }

        SimTK::Function* createSimTKFunction() const override
        {
            return new JorisPathSimTKFunction{_discretizations, _evaluations};
        }

        // TODO: flash vectors from properties with `extendFinalizeFromProperties`
    };

    /** TODO: handle computing a fresh FBP from a PBP, flashing from props, etc.

    // init underlying implementation data from a `FunctionBasedPath`s (precomputed) properties
    //
    // the properties being set in the FBP usually implies that the FBP has already been built
    // from a PBP at some previous point in time
    void Impl_InitFromFBPProperties(JorisFBP& impl)
    {
        FunctionBasedPathDiscretizationSet const& discSet = impl.getProperty_FunctionBasedPathDiscretizationSet().getValue();

        // set `coords` pointers to null
        //
        // they are lazily looked up in a later phase (after the model is connected up)
        impl.coords.clear();
        impl.coords.resize(discSet.getSize(), nullptr);

        // set `coordAbsPaths` from discretizations property
        impl.coordAbsPaths.clear();
        impl.coordAbsPaths.reserve(discSet.getSize());
        for (int i = 0; i < discSet.getSize(); ++i) {
            impl.coordAbsPaths.push_back(discSet[i].getProperty_coordinate_abspath().getValue());
        }

        // set `discretizations` from discretizations property
        impl.discretizations.clear();
        impl.discretizations.reserve(discSet.getSize());
        for (int i = 0; i < discSet.getSize(); ++i) {
            FunctionBasedPathDiscretization const& disc = discSet[i];
            Discretization d;
            d.begin = disc.getProperty_x_begin().getValue();
            d.end = disc.getProperty_x_end().getValue();
            d.nsteps = disc.getProperty_num_points().getValue();
            impl.discretizations.push_back(d);
        }

        // set `evals` from evaluations property
        auto const& evalsProp = impl.getProperty_Evaluations();
        impl.evals.clear();
        impl.evals.reserve(evalsProp.size());
        for (int i = 0; i < evalsProp.size(); ++i) {
            impl.evals.push_back(evalsProp.getValue(i));
        }
    }
    */
}


// TODO: integration test that quickly ensures the wireup is using Joris's alg
namespace joris {
}


// TODO: code that compiles a new "FunctionBasedPath" from a "PointBasedPath"
namespace joris {
    static constexpr int g_MaxCoordsThatCanAffectPathDefault = static_cast<int>(g_MaxNumDimensions);
    static constexpr int g_NumProbingDiscretizationsDefault = 8;
    static constexpr double g_MinProbingMomentArmChangeDefault = 0.001;
    static constexpr int g_NumDiscretizationStepsPerDimensionDefault = 8;

    // returns `true` if changing the supplied `Coordinate` changes the moment arm
    // of the supplied `PointBasedPath` (PBP)
    bool coordAffectsPBP(
            OpenSim::PointBasedPath const& pbp,
            OpenSim::Coordinate const& c,
            SimTK::State& state,
            int numProbingSteps,
            double minMomentArmChangeRequired)
    {
        bool initialLockedState = c.getLocked(state);
        double initialValue = c.getValue(state);

        c.setLocked(state, false);

        double start = c.getRangeMin();
        double end = c.getRangeMax();
        double step = (end - start) / (numProbingSteps-1);

        bool affectsCoord = false;
        for (double v = start; v <= end; v += step) {
            c.setValue(state, v);
            double ma = pbp.computeMomentArm(state, c);

            if (std::abs(ma) >= minMomentArmChangeRequired) {
                affectsCoord = true;
                break;
            }
        }

        c.setValue(state, initialValue);
        c.setLocked(state, initialLockedState);

        return affectsCoord;
    }

    // returns a sequence of `OpenSim::Coordinate`s that affect the supplied
    // point-based path (PBP)
    //
    // is not guaranteed to find *all* coordinates that affect the supplied PBP,
    // because that may involve extreme probing (which this implementation does not
    // do)
    std::vector<OpenSim::Coordinate const*> coordsThatAffectPBP(
            OpenSim::Model const& model,
            OpenSim::PointBasedPath const& pbp,
            SimTK::State& st,
            int numProbingSteps,
            double minMomentArmChangeRequired)
    {
        std::vector<const OpenSim::Coordinate*> rv;
        for (OpenSim::Coordinate const& c : model.getComponentList<OpenSim::Coordinate>()){
            if (c.getMotionType() == OpenSim::Coordinate::MotionType::Coupled) {
                continue;
            }

            if (!coordAffectsPBP(pbp, c, st, numProbingSteps, minMomentArmChangeRequired)) {
                continue;
            }

            rv.push_back(&c);
        }
        return rv;
    }

    // compute ideal discretization of the given coordinate
    Discretization discretizationForCoord(OpenSim::Coordinate const& c, int numDiscretizationSteps)
    {
        SimTK_ASSERT_ALWAYS(numDiscretizationSteps >= 4, "need to supply more than 4 discretization steps");

        Discretization d;
        //d.begin = -static_cast<double>(SimTK_PI)/2;
        //d.end = static_cast<double>(SimTK_PI)/2;
        d.begin = std::max(c.getRangeMin(), -static_cast<double>(SimTK_PI));
        d.end = std::min(c.getRangeMax(), static_cast<double>(SimTK_PI));
        d.nPoints = numDiscretizationSteps - 3;
        double step = stepSize(d);

        // expand range slightly in either direction to ensure interpolation is
        // clean around the edges
        d.begin -= step;
        d.end += 2.0 * step;
        d.nPoints += 3;

        return d;
    }

    // compute all permutations of the coordinates for the given discretizations
    //
    // these output values are stored "lexographically", with coords being "big endian", so:
    //
    // - [coord[0].begin, coord[1].begin, ..., coord[n-1].begin]
    // - [coord[0].begin, coord[1].begin, ..., (coord[n-1].begin + step)]
    // - [coord[0].begin, coord[1].begin, ..., (coord[n-1].begin + 2*step)]
    // - ...
    // - [coord[0].begin, (coord[1].begin + step), ..., coord[n-1].begin]
    // - [coord[0].begin, (coord[1].begin + step), ..., (coord[n-1].begin + step)]
    // - ...
    // - [(coord[0].begin + step), coord[1].begin, ..., coord[n-1].begin]
    // - [(coord[0].begin + step), coord[1].begin, ..., (coord[n-1].begin + step)]
    std::vector<double> computeEvaluationsFromPBP(
            OpenSim::PointBasedPath const& pbp,
            SimTK::State& state,
            OpenSim::Coordinate const** coords,
            Discretization const* discs,
            size_t ncoords)
    {
        std::vector<double> rv;

        if (ncoords == 0) {  // edge-case: logic below assumes ncoords > 0
            return rv;
        }

        OPENSIM_THROW_IF(ncoords > g_MaxNumDimensions, OpenSim::Exception, "too many coordinates affect this path - the FunctionBasedPath implementation cannot handle this");

        // number of evaluations is the total number of permutations of all dimensions for
        // all discretizations
        int expectedEvals = 1;
        for (size_t i = 0; i < ncoords; ++i) {
            expectedEvals *= discs[i].nPoints;
        }
        rv.reserve(expectedEvals);

        // holds which "step" in each Coordinate's [begin, end] discretization we
        // have evaluated up to
        std::array<int, g_MaxNumDimensions> discStepIdx{};
        while (discStepIdx[0] < discs[0].nPoints) {

            // set all coordinate values for this step
            for (size_t coord = 0; coord < ncoords; ++coord) {
                Discretization const& discr = discs[coord];

                double stepSz = stepSize(discr);
                int step = discStepIdx[coord];
                double val =  discr.begin + step*stepSz;

                coords[coord]->setValue(state, val);
            }

            // eval the length of the PBP for this permutation of coordinate values
            {
                double eval = pbp.getLength(state);
                rv.push_back(eval);
            }

            // update which coordinate steps we're up to for each coordinate
            //
            // always updates "least significant" coordinate first, then performs
            // "carry propagation" to the "more significant" coordinates
            int pos = ncoords - 1;
            discStepIdx[pos]++;
            while (pos > 0 && discStepIdx[pos] >= discs[pos].nPoints) {
                discStepIdx[pos] = 0;  // overflow
                ++discStepIdx[pos-1];  // carry
                --pos;
            }
        }

        SimTK_ASSERT_ALWAYS(discStepIdx[0] == discs[0].nPoints, "should be true, after the final overflow");
        for (size_t i = 1; i < discStepIdx.size(); ++i) {
            SimTK_ASSERT_ALWAYS(discStepIdx[i] == 0, "these less-significant coordinates should all be overflow-n by the end of the alg");
        }
        SimTK_ASSERT_ALWAYS(rv.size() == static_cast<size_t>(expectedEvals), "these two values should match, given the above alg");

        return rv;
    }

    struct FittingParams final {

        // maximum coords that can affect the given PointBasedPath
        //
        // if this is higher, more paths may be eligible for
        // PointBasedPath --> FunctionBasedPath conversion, because some paths may be
        // affected by more coordinates than other paths. However, be careful. Increasing
        // this also *significantly* increases the memory usage of the function-based fit
        //
        // must be 0 < v <= 16, or -1 to mean "use a sensible default"
        int maxCoordsThatCanAffectPath;

        // number of discretization steps to use for each coordinate during the "probing
        // phase"
        //
        // in the "probing phase", each coordinate is set to this number of evenly-spaced
        // values in the range [getRangeMin()..getRangeMax()] (inclusive) to see if changing
        // that coordinate has any affect on the path. The higher this value is, the longer
        // the probing phase takes, but the higher chance it has of spotting a pertubation
        //
        // must be >0, or -1 to mean "use a sensible default"
        int numProbingDiscretizations;

        // minimum amount that the moment arm of the path must change by during the "probing phase"
        // for the coorinate to be classified as affecting the path
        //
        // must be >0, or <0 to mean "use a sensible default"
        double minProbingMomentArmChange;

        // the number of discretization steps for each coordinate that passes the "probing phase" and,
        // therefore, is deemed to affect the input (point-based) path
        //
        // this is effectively "grid granulaity". More discretizations == better fit, but it can increase
        // the memory usage of the fit significantly. Assume the path is parameterized as an n-dimensional
        // surface. E.g. if you discretize 10 points over 10 dimensions then you may end up with
        // 10^10 datapoints (ouch).
        //
        // must be >0, or -1 to mean "use a sensible default"
        int numDiscretizationStepsPerDimension;

        FittingParams() :
            maxCoordsThatCanAffectPath{g_MaxCoordsThatCanAffectPathDefault},
            numProbingDiscretizations{g_NumProbingDiscretizationsDefault},
            minProbingMomentArmChange{g_MinProbingMomentArmChangeDefault},
            numDiscretizationStepsPerDimension{g_NumDiscretizationStepsPerDimensionDefault}
        {
        }
    };

    /* todo
    std::unique_ptr<JorisFBP> fromPointBasedPath(
            const Model& model,
            const PointBasedPath& pbp,
            FittingParams params)
    {
        // sanitize + validate params
        {
            if (params.maxCoordsThatCanAffectPath == -1) {
                params.maxCoordsThatCanAffectPath = g_MaxCoordsThatCanBeInterpolated;
            }

            if (params.numProbingDiscretizations == -1) {
                params.numProbingDiscretizations = g_NumProbingDiscretizationsDefault;
            }

            if (params.minProbingMomentArmChange < 0.0) {
                params.minProbingMomentArmChange = g_MinProbingMomentArmChangeDefault;
            }

            if (params.numDiscretizationStepsPerDimension == -1) {
                params.numDiscretizationStepsPerDimension = g_NumDiscretizationStepsPerDimensionDefault;
            }

            OPENSIM_THROW_IF(params.maxCoordsThatCanAffectPath <= 0, OpenSim::Exception, "maxCoordsThatCanAffectPath must be a positive number that is <=8");
            OPENSIM_THROW_IF(params.maxCoordsThatCanAffectPath > static_cast<int>(g_MaxCoordsThatCanBeInterpolated), OpenSim::Exception, "maxCoordsThatCanAffectPath must be a positive number that is <=8");
            OPENSIM_THROW_IF(params.numProbingDiscretizations <= 0, OpenSim::Exception, "numProbingDiscretizations must be a positive number");
            OPENSIM_THROW_IF(params.minProbingMomentArmChange <= 0, OpenSim::Exception, "minProbingMomentArmChange must be a positive number");
            OPENSIM_THROW_IF(params.numDiscretizationStepsPerDimension <= 0, OpenSim::Exception, "numDiscretizationStepsPerDimension must be a positive number");
        }

        std::unique_ptr<JorisFBP> fbp{new JorisFBP{}};

        JorisFBP& impl = *fbp;

        // compute underlying impl data from the PBP
        if (!Impl_ComputeFromPBP(impl, model, pbp, params)) {
            return nullptr;
        }

        // write impl discretizations into the `Discretizations` property
        FunctionBasedPathDiscretizationSet& set = fbp->updProperty_FunctionBasedPathDiscretizationSet().updValue();
        for (size_t i = 0; i < impl.coords.size(); ++i) {
            auto disc = std::unique_ptr<FunctionBasedPathDiscretization>{new FunctionBasedPathDiscretization{}};
            disc->set_x_begin(impl.discretizations[i].begin);
            disc->set_x_end(impl.discretizations[i].end);
            disc->set_num_points(impl.discretizations[i].nsteps);
            disc->set_coordinate_abspath(impl.coordAbsPaths[i]);
            set.adoptAndAppend(disc.release());
        }

        // write evals into `Evaluations` property
        auto& evalsProp = fbp->updProperty_Evaluations();
        for (double eval : fbp->evals) {
            evalsProp.appendValue(eval);
        }

        return fbp;
    }

    // compute fresh implementation data from an existing PointBasedPath by
    // evaluating it and fitting it to a function-based curve
    //
    // returns false if too many/too little coordinates affect the path
    bool Impl_ComputeFromPBP(
            const OpenSim::Model& model,
            const OpenSim::PointBasedPath& pbp,
            const FittingParams& params,
            std::vector<Discretization>& discretizationsOut,
            std::vector<double>& evalsOut,
            std::vector<std::string>& coordAbsPathsOut)
    {
        // copy model, so we can independently equilibrate + realize + modify the
        // copy without having to touch the source model
        std::unique_ptr<OpenSim::Model> modelClone{model.clone()};
        SimTK::State& initialState = modelClone->initSystem();
        modelClone->equilibrateMuscles(initialState);
        modelClone->realizeVelocity(initialState);

        // set `coords`
        impl.coords = coordsThatAffectPBP(*modelClone, pbp, initialState, params.numProbingDiscretizations, params.minProbingMomentArmChange);
        if (static_cast<int>(impl.coords.size()) > params.maxCoordsThatCanAffectPath || impl.coords.size() == 0) {
            impl.coords.clear();
            return false;
        }

        // set `coordAbsPaths`
        impl.coordAbsPaths.clear();
        impl.coordAbsPaths.reserve(impl.coords.size());
        for (const OpenSim::Coordinate* c : impl.coords) {
            impl.coordAbsPaths.push_back(c->getAbsolutePathString());
        }

        // set `discretizations`
        impl.discretizations.clear();
        impl.discretizations.reserve(impl.coords.size());
        for (const OpenSim::Coordinate* c : impl.coords) {
            impl.discretizations.push_back(discretizationForCoord(*c, params.numDiscretizationStepsPerDimension));
        }

        // set `evals`
        SimTK_ASSERT_ALWAYS(impl.coords.size() == impl.discretizations.size(), "these should be equal by now");
        impl.evals = computeEvaluationsFromPBP(pbp, initialState, impl.coords.data(), impl.discretizations.data(), impl.coords.size());

        return true;
    }
        */
}

// TODO: test fitting some basic `PointBasedPath`s using the funciton-fitting implementation
namespace joris {
}

// TODO: implement conversion tool
namespace joris {
}

// TODO: test conversion tool works as intended
namespace joris {
}

// TODO: implement CLI tool
namespace joris {
}

// TODO: implement CLI tests
namespace joris {
}

int main()
{
    return RunAllTests();
}
