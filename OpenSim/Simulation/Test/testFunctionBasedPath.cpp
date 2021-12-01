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
}


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
    OpenSim::FunctionBasedPath* fbp = new OpenSim::FunctionBasedPath{};
    model.addComponent(fbp);

    SimTK::State& s = model.initSystem();

    fbp->getLength(s);  // shouldn't throw
}

OSIM_TEST(FunctionBasedPath, CanCallGetLengtheningSpeedWithoutThrowing)
{
    OpenSim::Model model;
    OpenSim::FunctionBasedPath* fbp = new OpenSim::FunctionBasedPath{};
    model.addComponent(fbp);

    SimTK::State& s = model.initSystem();

    fbp->getLengtheningSpeed(s);  // shouldn't throw
}

OSIM_TEST(FunctionBasedPath, CanCallGetPointForceDirectionWithoutThrowing)
{
    OpenSim::Model model;
    OpenSim::FunctionBasedPath* fbp = new OpenSim::FunctionBasedPath{};
    model.addComponent(fbp);

    SimTK::State& s = model.initSystem();

    OpenSim::Array<OpenSim::PointForceDirection*> pfds;
    fbp->getPointForceDirections(s, &pfds);  // shouldn't throw
}

OSIM_TEST(FunctionBasedPath, CanCallAddInEquivalentForcesWithoutThrowing)
{
    OpenSim::Model model;
    OpenSim::FunctionBasedPath* fbp = new OpenSim::FunctionBasedPath{};
    model.addComponent(fbp);

    SimTK::State& s = model.initSystem();

    SimTK::Vector_<SimTK::SpatialVec> bodyForces;
    SimTK::Vector mobilityForces;

    fbp->addInEquivalentForces(s, 1.0, bodyForces, mobilityForces);
}

OSIM_TEST(FunctionBasedPath, CanCallComputeMomentArmWithoutThrowing)
{
    OpenSim::Model model;
    OpenSim::FunctionBasedPath* fbp = new OpenSim::FunctionBasedPath{};
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

int main()
{
    return RunAllTests();
}
