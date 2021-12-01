#include <OpenSim/Simulation/Model/FunctionBasedPath.h>
#include <OpenSim/Simulation/Model/Model.h>
#include <OpenSim/Simulation/SimbodyEngine/PinJoint.h>

#include <chrono>
#include <random>
#include <stdexcept>
#include <vector>

// this is a poor-man's GoogleTest, but I still prefer it to the madness of
// having a massive `main()` containing many `try..catch` blocks
namespace {
    struct Test { char const* suiteName; char const* name; void(*testFunc)(void); };

    std::ostream& operator<<(std::ostream& o, Test const& t)
    {
        return o << t.suiteName << ':' << t.name;
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

    struct MockPathFunctionData {
        int numTimesGetLengthCalled = 0;
        double lengthValue = GenerateDouble();

        int numTimesLengtheningSpeedCalled = 0;
        double lengtheningSpeedValue = GenerateDouble();

        int numTimesComputeMomentArmCalled = 0;
        double computeMomentArmValue = GenerateDouble();

        int numTimesGetPointForceDirectionsCalled = 0;
        int numTimesAddInEquivalentForcesCalled = 0;

        int numTimesExtendConnectToModelCalled = 0;
        int numTimesExtendInitStateFromPropertiesCalled = 0;
        int numTimesExtendAddToSystemCalled = 0;
        int numTimesExtendFinalizeFromPropertiesCalled = 0;
    };

    // used to test that the FunctionBasedPath is forwarding things correctly
    class MockPathFunction : public OpenSim::PathFunction {
        OpenSim_DECLARE_CONCRETE_OBJECT(MockPathFunction, OpenSim::PathFunction);
    public:
        mutable std::shared_ptr<MockPathFunctionData> data = std::make_shared<MockPathFunctionData>();

        double getLength(const SimTK::State&) const override
        {
            ++data->numTimesGetLengthCalled;
            return data->lengthValue;
        }

        double getLengtheningSpeed(const SimTK::State&) const override
        {
            ++data->numTimesLengtheningSpeedCalled;
            return data->lengtheningSpeedValue;
        }

        double computeMomentArm(const SimTK::State&, const OpenSim::Coordinate&) const override
        {
            ++data->numTimesComputeMomentArmCalled;
            return data->computeMomentArmValue;
        }

        void getPointForceDirections(const SimTK::State& s, OpenSim::Array<OpenSim::PointForceDirection*>* rPFDs) const override
        {
            ++data->numTimesGetPointForceDirectionsCalled;
        }

        void addInEquivalentForces(const SimTK::State& state, double tension, SimTK::Vector_<SimTK::SpatialVec>& bodyForces, SimTK::Vector& mobilityForces) const override
        {
            ++data->numTimesAddInEquivalentForcesCalled;
        }

        // these are redundantly checked (the Component/Property tests *should*
        // also test these) because existing implementations do use these methods
        // to hook into the model/system/state at various steps (e.g. to cache
        // coordinates, or whatever they need) and it's handy to redundantly
        // ensure these are hooked up to the PathFunction via the FunctionBasedPath
        // correctly.

        void extendConnectToModel(OpenSim::Model&) override
        {
            ++data->numTimesExtendConnectToModelCalled;
        }

        void extendInitStateFromProperties(SimTK::State&) const override
        {
            ++data->numTimesExtendInitStateFromPropertiesCalled;
        }

        void extendAddToSystem(SimTK::MultibodySystem&) const override
        {
            ++data->numTimesExtendAddToSystemCalled;
        }

        void extendFinalizeFromProperties() override
        {
            ++data->numTimesExtendFinalizeFromPropertiesCalled;
        }
    };
}


OSIM_TEST(FunctionBasedPath, CanBeDefaultConstructedWithoutThrowing)
{
    OpenSim::FunctionBasedPath fbp;  // shouldn't throw
}

OSIM_TEST(FunctionBasedPath, CanBeConstructedWithAPathFunction)
{
    MockPathFunction pathFn;
    OpenSim::FunctionBasedPath fbp{pathFn};
}

OSIM_TEST(FunctionBasedPath, WhenConstructedWithPathFunctionUsesTheFunctionInLengthEvaluation)
{
    OpenSim::Model model;
    MockPathFunction pathFn;
    OpenSim::FunctionBasedPath* fbp = new OpenSim::FunctionBasedPath{pathFn};
    model.addComponent(fbp);

    SimTK_TEST(pathFn.data->numTimesGetLengthCalled == 0);

    SimTK::State& s = model.initSystem();

    SimTK_TEST(pathFn.data->numTimesGetLengthCalled == 0);

    SimTK_TEST(fbp->getLength(s) == pathFn.data->lengthValue);

    SimTK_TEST(pathFn.data->numTimesGetLengthCalled == 1);
}

OSIM_TEST(FunctionBasedPath, CanBeCopyConstructedWithoutThrowing)
{
    OpenSim::FunctionBasedPath fbp1;
    OpenSim::FunctionBasedPath fbp2{fbp1};  // shouldn't throw
}

OSIM_TEST(FunctionBasedPath, CanBeMoveConstructedWithoutThrowing)
{
    OpenSim::FunctionBasedPath fbp1;
    OpenSim::FunctionBasedPath fbp2{std::move(fbp1)};
}

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

OSIM_TEST(FunctionBasedPath, CanGetColorWithoutThrowing)
{
    OpenSim::Model model;
    OpenSim::FunctionBasedPath* fbp = new OpenSim::FunctionBasedPath{};
    model.addComponent(fbp);

    SimTK::State& s = model.initSystem();

    fbp->getColor(s);
}

OSIM_TEST(FunctionBasedPath, GetColorReturnsDefaultColorIfSetColorIsNotCalled)
{
    SimTK::Vec3 randomColor = GenerateRandomVector();

    OpenSim::Model model;
    OpenSim::FunctionBasedPath* fbp = new OpenSim::FunctionBasedPath{};
    fbp->setDefaultColor(randomColor);
    model.addComponent(fbp);

    SimTK::State& s = model.initSystem();

    SimTK_TEST(fbp->getColor(s) == randomColor);
}

OSIM_TEST(FunctionBasedPath, SetColorSetsTheColorInTheState)
{
    SimTK::Vec3 randomColor = GenerateRandomVector();

    OpenSim::Model model;
    OpenSim::FunctionBasedPath* fbp = new OpenSim::FunctionBasedPath{};
    fbp->setDefaultColor(randomColor);
    model.addComponent(fbp);

    SimTK::State& s = model.initSystem();

    fbp->setColor(s, randomColor);

    SimTK_TEST(fbp->getColor(s) == randomColor);
}

OSIM_TEST(FunctionBasedPath, CanCallGetLengthWithoutThrowing)
{
    OpenSim::Model model;
    MockPathFunction pathFn;
    OpenSim::FunctionBasedPath* fbp = new OpenSim::FunctionBasedPath{pathFn};
    model.addComponent(fbp);

    SimTK::State& s = model.initSystem();
    model.realizePosition(s);

    fbp->getLength(s);  // shouldn't throw
}

OSIM_TEST(FunctionBasedPath, GetLengthUsesPathFunctionImpl)
{
    OpenSim::Model model;
    MockPathFunction pathFn;
    OpenSim::FunctionBasedPath* fbp = new OpenSim::FunctionBasedPath{pathFn};
    model.addComponent(fbp);

    SimTK_TEST(pathFn.data->numTimesGetLengthCalled == 0);

    SimTK::State& s = model.initSystem();
    model.realizePosition(s);

    SimTK_TEST(pathFn.data->numTimesGetLengthCalled == 0);

    SimTK_TEST(fbp->getLength(s) == pathFn.data->lengthValue);

    SimTK_TEST(pathFn.data->numTimesGetLengthCalled == 1);
}

OSIM_TEST(FunctionBasedPath, GetLengthIsCached)
{
    OpenSim::Model model;
    MockPathFunction pathFn;
    OpenSim::FunctionBasedPath* fbp = new OpenSim::FunctionBasedPath{pathFn};
    model.addComponent(fbp);

    SimTK_TEST(pathFn.data->numTimesGetLengthCalled == 0);

    SimTK::State& s = model.initSystem();
    model.realizePosition(s);

    SimTK_TEST(pathFn.data->numTimesGetLengthCalled == 0);

    SimTK_TEST(fbp->getLength(s) == pathFn.data->lengthValue);

    SimTK_TEST(pathFn.data->numTimesGetLengthCalled == 1);

    SimTK_TEST(fbp->getLength(s) == pathFn.data->lengthValue);  // should cache

    SimTK_TEST(pathFn.data->numTimesGetLengthCalled == 1);
}

OSIM_TEST(FunctionBasedPath, CanCallGetLengtheningSpeedWithoutThrowing)
{
    OpenSim::Model model;
    MockPathFunction pathFn;
    OpenSim::FunctionBasedPath* fbp = new OpenSim::FunctionBasedPath{pathFn};
    model.addComponent(fbp);

    SimTK::State& s = model.initSystem();
    model.realizeAcceleration(s);

    fbp->getLengtheningSpeed(s);  // shouldn't throw
}

OSIM_TEST(FunctionBasedPath, GetLengtheningSpeedUsesPathFunctionImpl)
{
    OpenSim::Model model;
    MockPathFunction pathFn;
    OpenSim::FunctionBasedPath* fbp = new OpenSim::FunctionBasedPath{pathFn};
    model.addComponent(fbp);

    SimTK_TEST(pathFn.data->numTimesLengtheningSpeedCalled == 0);

    SimTK::State& s = model.initSystem();
    model.realizeAcceleration(s);

    SimTK_TEST(pathFn.data->numTimesLengtheningSpeedCalled == 0);

    SimTK_TEST(fbp->getLengtheningSpeed(s) == pathFn.data->lengtheningSpeedValue);

    SimTK_TEST(pathFn.data->numTimesLengtheningSpeedCalled == 1);
}

OSIM_TEST(FunctionBasedPath, GetLengtheningSpeedIsCached)
{
    OpenSim::Model model;
    MockPathFunction pathFn;
    OpenSim::FunctionBasedPath* fbp = new OpenSim::FunctionBasedPath{pathFn};
    model.addComponent(fbp);

    SimTK_TEST(pathFn.data->numTimesLengtheningSpeedCalled == 0);

    SimTK::State& s = model.initSystem();
    model.realizeAcceleration(s);

    SimTK_TEST(pathFn.data->numTimesLengtheningSpeedCalled == 0);

    SimTK_TEST(fbp->getLengtheningSpeed(s) == pathFn.data->lengtheningSpeedValue);

    SimTK_TEST(pathFn.data->numTimesLengtheningSpeedCalled == 1);

    SimTK_TEST(fbp->getLengtheningSpeed(s) == pathFn.data->lengtheningSpeedValue);  // should cache

    SimTK_TEST(pathFn.data->numTimesLengtheningSpeedCalled == 1);
}

OSIM_TEST(FunctionBasedPath, CanCallGetPointForceDirectionWithoutThrowing)
{
    OpenSim::Model model;
    MockPathFunction pathFn;
    OpenSim::FunctionBasedPath* fbp = new OpenSim::FunctionBasedPath{pathFn};
    model.addComponent(fbp);

    SimTK::State& s = model.initSystem();

    OpenSim::Array<OpenSim::PointForceDirection*> pfds;
    fbp->getPointForceDirections(s, &pfds);  // shouldn't throw
}

OSIM_TEST(FunctionBasedPath, GetPointForceDirectionUsesPathFunctionImpl)
{
    OpenSim::Model model;
    MockPathFunction pathFn;
    OpenSim::FunctionBasedPath* fbp = new OpenSim::FunctionBasedPath{pathFn};
    model.addComponent(fbp);

    SimTK_TEST(pathFn.data->numTimesGetPointForceDirectionsCalled == 0);

    SimTK::State& s = model.initSystem();

    SimTK_TEST(pathFn.data->numTimesGetPointForceDirectionsCalled == 0);

    {
        OpenSim::Array<OpenSim::PointForceDirection*> pfds;
        fbp->getPointForceDirections(s, &pfds);
    }

    SimTK_TEST(pathFn.data->numTimesGetPointForceDirectionsCalled == 1);
}

OSIM_TEST(FunctionBasedPath, GetPointForceDirectionsIsNotCached)
{
    OpenSim::Model model;
    MockPathFunction pathFn;
    OpenSim::FunctionBasedPath* fbp = new OpenSim::FunctionBasedPath{pathFn};
    model.addComponent(fbp);

    SimTK_TEST(pathFn.data->numTimesGetPointForceDirectionsCalled == 0);

    SimTK::State& s = model.initSystem();

    SimTK_TEST(pathFn.data->numTimesGetPointForceDirectionsCalled == 0);

    {
        OpenSim::Array<OpenSim::PointForceDirection*> pfds;
        fbp->getPointForceDirections(s, &pfds);
    }

    SimTK_TEST(pathFn.data->numTimesGetPointForceDirectionsCalled == 1);

    {
        OpenSim::Array<OpenSim::PointForceDirection*> pfds;
        fbp->getPointForceDirections(s, &pfds);  // not cached
    }

    SimTK_TEST(pathFn.data->numTimesGetPointForceDirectionsCalled == 2);
}

OSIM_TEST(FunctionBasedPath, CanCallAddInEquivalentForcesWithoutThrowing)
{
    OpenSim::Model model;
    MockPathFunction pathFn;
    OpenSim::FunctionBasedPath* fbp = new OpenSim::FunctionBasedPath{pathFn};
    model.addComponent(fbp);

    SimTK::State& s = model.initSystem();

    SimTK::Vector_<SimTK::SpatialVec> bodyForces;
    SimTK::Vector mobilityForces;

    fbp->addInEquivalentForces(s, 1.0, bodyForces, mobilityForces);
}

OSIM_TEST(FunctionBasedPath, AddInEquivalentForcesUsesPathFunctionImpl)
{
    OpenSim::Model model;
    MockPathFunction pathFn;
    OpenSim::FunctionBasedPath* fbp = new OpenSim::FunctionBasedPath{pathFn};
    model.addComponent(fbp);

    SimTK_TEST(pathFn.data->numTimesAddInEquivalentForcesCalled == 0);

    SimTK::State& s = model.initSystem();

    SimTK_TEST(pathFn.data->numTimesAddInEquivalentForcesCalled == 0);

    {
        SimTK::Vector_<SimTK::SpatialVec> bodyForces;
        SimTK::Vector mobilityForces;
        fbp->addInEquivalentForces(s, 1.0, bodyForces, mobilityForces);
    }

    SimTK_TEST(pathFn.data->numTimesAddInEquivalentForcesCalled == 1);
}

OSIM_TEST(FunctionBasedPath, AddInEquivalentForcesIsNotCached)
{
    OpenSim::Model model;
    MockPathFunction pathFn;
    OpenSim::FunctionBasedPath* fbp = new OpenSim::FunctionBasedPath{pathFn};
    model.addComponent(fbp);

    SimTK_TEST(pathFn.data->numTimesAddInEquivalentForcesCalled == 0);

    SimTK::State& s = model.initSystem();

    SimTK_TEST(pathFn.data->numTimesAddInEquivalentForcesCalled == 0);

    {
        SimTK::Vector_<SimTK::SpatialVec> bodyForces;
        SimTK::Vector mobilityForces;
        fbp->addInEquivalentForces(s, 1.0, bodyForces, mobilityForces);
    }

    SimTK_TEST(pathFn.data->numTimesAddInEquivalentForcesCalled == 1);

    {
        SimTK::Vector_<SimTK::SpatialVec> bodyForces;
        SimTK::Vector mobilityForces;
        fbp->addInEquivalentForces(s, 1.0, bodyForces, mobilityForces);
    }

    SimTK_TEST(pathFn.data->numTimesAddInEquivalentForcesCalled == 2);
}


OSIM_TEST(FunctionBasedPath, CanCallComputeMomentArmWithoutThrowing)
{
    OpenSim::Model model;
    MockPathFunction pathFn;
    OpenSim::FunctionBasedPath* fbp = new OpenSim::FunctionBasedPath{pathFn};
    model.addComponent(fbp);

    OpenSim::Body* body = new OpenSim::Body{};
    body->setMass(1.0);
    model.addComponent(body);

    OpenSim::PinJoint* joint = new OpenSim::PinJoint{};
    joint->connectSocket_parent_frame(model.getGround());
    joint->connectSocket_child_frame(*body);
    model.addJoint(joint);

    SimTK::State& s = model.initSystem();

    fbp->computeMomentArm(s, joint->getCoordinate());
}

OSIM_TEST(FunctionBasedPath, ComputeMomentArmUsesPathFunctionImpl)
{
    OpenSim::Model model;
    MockPathFunction pathFn;
    OpenSim::FunctionBasedPath* fbp = new OpenSim::FunctionBasedPath{pathFn};
    model.addComponent(fbp);

    OpenSim::Body* body = new OpenSim::Body{};
    body->setMass(1.0);
    model.addComponent(body);

    OpenSim::PinJoint* joint = new OpenSim::PinJoint{};
    joint->connectSocket_parent_frame(model.getGround());
    joint->connectSocket_child_frame(*body);
    model.addJoint(joint);

    SimTK_TEST(pathFn.data->numTimesComputeMomentArmCalled == 0);

    SimTK::State& s = model.initSystem();

    SimTK_TEST(pathFn.data->numTimesComputeMomentArmCalled == 0);

    SimTK_TEST(fbp->computeMomentArm(s, joint->getCoordinate()) == pathFn.data->computeMomentArmValue);

    SimTK_TEST(pathFn.data->numTimesComputeMomentArmCalled == 1);
}

OSIM_TEST(FunctionBasedPath, PathFunctionExtendConnectToModelIsCalledWhenFinalizeConnectionCalledOnTopLevelModel)
{
    OpenSim::Model model;
    MockPathFunction pathFn;
    OpenSim::FunctionBasedPath* fbp = new OpenSim::FunctionBasedPath{pathFn};
    model.addComponent(fbp);

    SimTK_TEST(pathFn.data->numTimesExtendConnectToModelCalled == 0);
    model.finalizeConnections();
    SimTK_TEST(pathFn.data->numTimesExtendConnectToModelCalled == 1);
}

OSIM_TEST(FunctionBasedPath, PathFunctionExtendInitStateFromPropertiesCalledWhenCalledOnTopLevelModel)
{
    OpenSim::Model model;
    MockPathFunction pathFn;
    OpenSim::FunctionBasedPath* fbp = new OpenSim::FunctionBasedPath{pathFn};
    model.addComponent(fbp);

    SimTK_TEST(pathFn.data->numTimesExtendInitStateFromPropertiesCalled == 0);
    model.initSystem();
    SimTK_TEST(pathFn.data->numTimesExtendConnectToModelCalled == 1);
}

OSIM_TEST(FunctionBasedPath, PathFunctionExtendAddToSystemCalledWhenInitSystemCalledOnTopLevelModel)
{
    OpenSim::Model model;
    MockPathFunction pathFn;
    OpenSim::FunctionBasedPath* fbp = new OpenSim::FunctionBasedPath{pathFn};
    model.addComponent(fbp);

    SimTK_TEST(pathFn.data->numTimesExtendAddToSystemCalled == 0);
    model.initSystem();
    SimTK_TEST(pathFn.data->numTimesExtendAddToSystemCalled == 1);
}

OSIM_TEST(FunctionBasedPath, PathFunctionExtendFinalizePropertiesCalledWhenFinalizeFromPropertiesCalledOnTopLevelModel)
{
    OpenSim::Model model;
    MockPathFunction pathFn;
    OpenSim::FunctionBasedPath* fbp = new OpenSim::FunctionBasedPath{pathFn};
    model.addComponent(fbp);  // also calls it?

    SimTK_TEST(pathFn.data->numTimesExtendFinalizeFromPropertiesCalled == 1);
    model.finalizeFromProperties();
    SimTK_TEST(pathFn.data->numTimesExtendFinalizeFromPropertiesCalled == 2);
}

int main()
{
    return RunAllTests();
}
