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

    struct MockPathFunctionData {
        int numTimesGetLengthCalled = 0;
        double lengthValue = GenerateDouble();

        int numTimesLengtheningSpeedCalled = 0;
        double lengtheningSpeedValue = GenerateDouble();

        int numTimesComputeMomentArmCalled = 0;
        double computeMomentArmValue = GenerateDouble();

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

        void addInEquivalentForces(const SimTK::State&, double, SimTK::Vector_<SimTK::SpatialVec>&, SimTK::Vector&) const override
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


// HACK: test Joris's implementation here
//
// this is because I cba'd splitting OSIM_TEST into another compilation unit, and because
// it's currently handy to have all this code in one unit while we flesh out the implementation

#include <OpenSim/Simulation/Model/PointBasedPath.h>

namespace joris {
    static constexpr size_t g_MaxCoordsThatCanBeInterpolated = 8;  // important: this is an upper limit that's used for stack allocations
    static constexpr int g_MaxCoordsThatCanAffectPathDefault = static_cast<int>(g_MaxCoordsThatCanBeInterpolated);
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

    // discretization of a particular coordinate
    //
    // assumes `nsteps` evenly-spaced points ranging from [begin, end] (inclusive)
    struct Discretization final {
        double begin;
        double end;
        int nsteps;
    };

    // compute ideal discretization of the given coordinate
    Discretization discretizationForCoord(OpenSim::Coordinate const& c, int numDiscretizationSteps)
    {
        SimTK_ASSERT_ALWAYS(numDiscretizationSteps >= 4, "need to supply more than 4 discretization steps");

        Discretization d;
        //d.begin = -static_cast<double>(SimTK_PI)/2;
        //d.end = static_cast<double>(SimTK_PI)/2;
        d.begin = std::max(c.getRangeMin(), -static_cast<double>(SimTK_PI));
        d.end = std::min(c.getRangeMax(), static_cast<double>(SimTK_PI));
        d.nsteps = numDiscretizationSteps - 3;
        double step = (d.end-d.begin) / (d.nsteps-1);

        // expand range slightly in either direction to ensure interpolation is
        // clean around the edges
        d.begin -= step;
        d.end += 2.0 * step;
        d.nsteps += 3;

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

        OPENSIM_THROW_IF(ncoords > g_MaxCoordsThatCanBeInterpolated, OpenSim::Exception, "too many coordinates affect this path - the FunctionBasedPath implementation cannot handle this");

        // number of evaluations is the total number of permutations of all dimensions for
        // all discretizations
        int expectedEvals = 1;
        for (size_t i = 0; i < ncoords; ++i) {
            expectedEvals *= discs[i].nsteps;
        }
        rv.reserve(expectedEvals);

        // holds which "step" in each Coordinate's [begin, end] discretization we
        // have evaluated up to
        std::array<int, g_MaxCoordsThatCanBeInterpolated> discStepIdx{};
        while (discStepIdx[0] < discs[0].nsteps) {

            // set all coordinate values for this step
            for (size_t coord = 0; coord < ncoords; ++coord) {
                Discretization const& discr = discs[coord];

                double stepSz = (discr.end - discr.begin) / (discr.nsteps - 1);
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
            while (pos > 0 && discStepIdx[pos] >= discs[pos].nsteps) {
                discStepIdx[pos] = 0;  // overflow
                ++discStepIdx[pos-1];  // carry
                --pos;
            }
        }

        SimTK_ASSERT_ALWAYS(discStepIdx[0] == discs[0].nsteps, "should be true, after the final overflow");
        for (size_t i = 1; i < discStepIdx.size(); ++i) {
            SimTK_ASSERT_ALWAYS(discStepIdx[i] == 0, "these less-significant coordinates should all be overflow-n by the end of the alg");
        }
        SimTK_ASSERT_ALWAYS(rv.size() == static_cast<size_t>(expectedEvals), "these two values should match, given the above alg");

        return rv;
    }

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

    class JorisFBP : public OpenSim::PathFunction {
        OpenSim_DECLARE_CONCRETE_OBJECT(JorisFBP, OpenSim::PathFunction);

    public:
        OpenSim_DECLARE_UNNAMED_PROPERTY(FunctionBasedPathDiscretizationSet, "Discretizations that were used for each OpenSim::Coordinate that the path was fitted against");
        OpenSim_DECLARE_LIST_PROPERTY(Evaluations, double, "The evaluated results of each *permutation* of discretizations. The FunctionBasedPathDiscretizationSet property describes how each OpenSim::Coordinate was discretized. These evaluations are the result of permuting through all possible combinations of discretizations. Effectively, this property contains a N-dimensional 'surface' of points, where each dimension of the surface is a Coordinate, and each dimension of each point is one of the evenly-spaced points in the discretization range [x_begin, x_range] for each dimension");

        // direct pointers to each coordinate
        std::vector<OpenSim::Coordinate const*> coords;

        // absolute paths of each coordinate (1:1 with coords)
        std::vector<std::string> coordAbsPaths;

        // discretizations ranges for each coordinate (1:1 with coords)
        std::vector<Discretization> discretizations;

        // evaluations for each permutation of coordinates' discretizations
        std::vector<double> evals;

        JorisFBP();
        double getLength(const SimTK::State&) const override;
        double getLengtheningSpeed(const SimTK::State&) const override;
        double computeMomentArm(const SimTK::State&, const OpenSim::Coordinate&) const override;
        void addInEquivalentForces(const SimTK::State& state, double tension, SimTK::Vector_<SimTK::SpatialVec>& bodyForces, SimTK::Vector& mobilityForces) const override;
        void extendFinalizeFromProperties() override;
        void extendFinalizeConnections(OpenSim::Component&) override;
    };

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

    // compute fresh implementation data from an existing PointBasedPath by
    // evaluating it and fitting it to a function-based curve
    //
    // returns false if too many/too little coordinates affect the path
    bool Impl_ComputeFromPBP(
            JorisFBP& impl,
            const OpenSim::Model& model,
            const OpenSim::PointBasedPath& pbp,
            const FittingParams& params)
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

    // ensure that the OpenSim::Coordinate* pointers held in Impl are up-to-date
    //
    // the pointers are there to reduce runtime path lookups
    static void Impl_SetCoordinatePointersFromCoordinatePaths(JorisFBP& impl,
                                                              OpenSim::Component const& c) {

        for (size_t i = 0; i < impl.coords.size(); ++i) {
            impl.coords[i] = &c.getComponent<OpenSim::Coordinate>(impl.coordAbsPaths[i]);
        }
    }

    // returns interpolated path length for a given permutation of coordinate
    // values
    //
    // this is the "heart" of the FPB algorithm. It's loosely based on the algorithm
    // described here:
    //
    //     "Two hierarchies of spline interpolations. Practical algorithms for multivariate higher order splines"
    //     https://arxiv.org/abs/0905.3564
    //
    // `inputVals` points to a sequence of `nCoords` values that were probably
    // retrieved via `Coordinate::getValue(SimTK::State const&)`. The reason
    // that `inputVals` is provided externally (rather than have this implementation
    // handle calling `getValue`) is because derivative calculations need to fiddle
    // the input values slightly
    double Impl_GetPathLength(JorisFBP const& impl,
                                     double const* inputVals,
                                     int nCoords)
    {
        SimTK_ASSERT_ALWAYS(!impl.coords.empty(), "FBPs require at least one coordinate to affect the path");
        SimTK_ASSERT_ALWAYS(nCoords == static_cast<int>(impl.coords.size()), "You must call this function with the correct number of (precomputed) coordinate values");

        // compute:
        //
        // - the index of the first discretization step *before* the input value
        //
        // - the polynomial of the curve at that step, given its fractional distance
        //   toward the next step
        using Polynomial = std::array<double, 4>;
        std::array<int, g_MaxCoordsThatCanBeInterpolated> closestDiscretizationSteps;
        std::array<Polynomial, g_MaxCoordsThatCanBeInterpolated> betas;
        for (int coord = 0; coord < nCoords; ++coord) {
            double inputVal = inputVals[coord];
            Discretization const& disc = impl.discretizations[coord];
            double step = (disc.end - disc.begin) / (disc.nsteps - 1);

            // compute index of first discretization step *before* the input value and
            // the fraction that the input value is towards the *next* discretization step
            int idx;
            double frac;
            if (inputVal < disc.begin+step) {
                idx = 1;
                frac = 0.0;
            } else if (inputVal > disc.end-2*step) {
                idx = disc.nsteps-3;
                frac = 0.0;
            } else {
                // solve for `n`: inputVal = begin + n*step
                double n = (inputVal - disc.begin) / step;
                double wholePart;
                double fractionalPart = std::modf(n, &wholePart);

                idx = static_cast<int>(wholePart);
                frac = fractionalPart;
            }

            // compute polynomial based on fraction the point is toward the next point
            double frac2 = frac*frac;
            double frac3 = frac2*frac;
            double frac4 = frac3*frac;
            double fracMinusOne = frac - 1;
            double fracMinusOne3 = fracMinusOne*fracMinusOne*fracMinusOne;

            Polynomial p;
            p[0] =  0.5 * fracMinusOne3*frac*(2*frac + 1);
            p[1] = -0.5 * (frac - 1)*(6*frac4 - 9*frac3 + 2*frac + 2);
            p[2] =  0.5 * frac*(6*frac4 - 15*frac3 + 9*frac2 + frac + 1);
            p[3] = -0.5 * (frac - 1)*frac3*(2*frac - 3);

            closestDiscretizationSteps[coord] = idx;
            betas[coord] = p;
        }

        // for each coord, permute through 4 locations *around* the input's location:
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

        // represent permuting through each location around each coordinate as a string
        // of integer offsets that can be -1, 0, 1, or 2
        //
        // the algorithm increments this array as it goes through each permutation
        std::array<int, g_MaxCoordsThatCanBeInterpolated> dimIdxOffsets;
        for (int coord = 0; coord < nCoords; ++coord) {
            dimIdxOffsets[coord] = -1;
        }

        // permute through all locations around the input value
        //
        // e.g. the location permutations for 3 coords iterate like this for each
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

        double z = 0.0;
        int cnt = 0;
        while (dimIdxOffsets[0] < 3) {

            // compute `beta` (weighted coefficient per coord) for this particular
            // permutation's coordinate locations (e.g. -1, 0, 0, 2) and figure out
            // what the closest input value was at the weighted location. Add the
            // result the the output

            double beta = 1.0;
            int evalStride = 1;
            int evalIdx = 0;

            // go backwards, from least-significant coordinate (highest idx)
            //
            // this is so that we can compute the stride as the algorithm runs
            for (int coord = nCoords-1; coord >= 0; --coord) {
                int offset = dimIdxOffsets[coord];  // -1, 0, 1, or 2
                int closestStep = closestDiscretizationSteps[coord];
                int step = closestStep + offset;

                beta *= betas[coord][offset+1];
                evalIdx += evalStride * step;
                evalStride *= impl.discretizations[coord].nsteps;
            }

            // equivalent to z += b*v, but handles rounding errors when the rhs
            // is very small
            z = std::fma(beta, impl.evals.at(evalIdx), z);

            // increment the offsets
            //
            // this is effectively the step that permutes [-1, -1, 2] --> [-1,  0, -1]
            {
                int pos = nCoords-1;
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
            int expectedIterations = 1 << (2*nCoords);
            if (cnt != expectedIterations) {
                std::stringstream msg;
                msg << "invalid number of permutations explored: expected = " << expectedIterations << ", got = " << cnt;
                OPENSIM_THROW(OpenSim::Exception, std::move(msg).str());
            }
        }

        return z;
    }

    // get the length of the path in the current state
    double Impl_GetPathLength(JorisFBP const& impl, SimTK::State const& s)
    {
        int nCoords = static_cast<int>(impl.coords.size());

        // get the input value of each coordinate in the current state
        std::array<double, g_MaxCoordsThatCanBeInterpolated> inputVals{};
        for (int coord = 0; coord < nCoords; ++coord) {
            inputVals[coord] = impl.coords[coord]->getValue(s);
        }

        return Impl_GetPathLength(impl, inputVals.data(), nCoords);
    }

    // get the *derivative* of the path length with respect to the given Coordinate index
    // (in impl.coords)
    //static double Impl_GetPathLengthDerivative(OpenSim::FunctionBasedPath::Impl const& impl,
    //                                           SimTK::State const& s,
    //                                           int coordIdx) {

    //    SimTK_ASSERT_ALWAYS(!impl.coords.empty(), "FBPs require at least one coordinate to affect the path");
    //    SimTK_ASSERT_ALWAYS(coordIdx != -1, "coord index must be valid");
    //    SimTK_ASSERT_ALWAYS(coordIdx < static_cast<int>(impl.coords.size()), "coord index must be valid");

    //    int nCoords = static_cast<int>(impl.coords.size());

    //    // get the input value of each coordinate in the current state
    //    std::array<double, g_MaxCoordsThatCanBeInterpolated> inputVals{};
    //    for (int coord = 0; coord < nCoords; ++coord) {
    //        inputVals[coord] = impl.coords[coord]->getValue(s);
    //    }

    //    // compute value at current point
    //    double v1 = Impl_GetPathLength(impl, inputVals.data(), nCoords);

    //    static constexpr double h = 0.00001;

    //    // alter the input value for the to-be-derived coordinate *slightly* and recompute
    //    inputVals[coordIdx] += h;
    //    double v2 = Impl_GetPathLength(impl, inputVals.data(), nCoords);

    //    // the derivative is how much the output changed when the input was altered
    //    // slightly (this is a poor-man's discrete derivative method)
    //    return (v2 - v1) / h;
    //}
    double Impl_GetPathLengthDerivative(JorisFBP const& impl,
                                        SimTK::State const& s,
                                        int coordIdx)
    {
        int nCoords = static_cast<int>(impl.coords.size());

        SimTK_ASSERT_ALWAYS(!impl.coords.empty(), "FBPs require at least one coordinate to affect the path");
        SimTK_ASSERT_ALWAYS(nCoords == static_cast<int>(impl.coords.size()), "You must call this function with the correct number of (precomputed) coordinate values");

        // get the input value of each coordinate in the current state
        std::array<double, g_MaxCoordsThatCanBeInterpolated> inputVals{};
        for (int coord = 0; coord < nCoords; ++coord) {
            inputVals[coord] = impl.coords[coord]->getValue(s);
        }

        // compute:
        //
        // - the index of the first discretization step *before* the input value
        //
        // - the polynomial of the curve at that step, given its fractional distance
        //   toward the next step
        using Polynomial = std::array<double, 4>;
        std::array<int, g_MaxCoordsThatCanBeInterpolated> closestDiscretizationSteps;
        std::array<Polynomial, g_MaxCoordsThatCanBeInterpolated> betas;
        for (int coord = 0; coord < nCoords; ++coord) {
            double inputVal = inputVals[coord];
            Discretization const& disc = impl.discretizations[coord];
            double step = (disc.end - disc.begin) / (disc.nsteps - 1);

            // compute index of first discretization step *before* the input value and
            // the fraction that the input value is towards the *next* discretization step
            int idx;
            double frac;
            if (inputVal < disc.begin+step) {
                idx = 1;
                frac = 0.0;
            } else if (inputVal > disc.end-2*step) {
                idx = disc.nsteps-3;
                frac = 0.0;
            } else {
                // solve for `n`: inputVal = begin + n*step
                double n = (inputVal - disc.begin) / step;
                double wholePart;
                double fractionalPart = std::modf(n, &wholePart);

                idx = static_cast<int>(wholePart);
                frac = fractionalPart;
            }

            // compute polynomial based on fraction the point is toward the next point
            double frac2 = frac*frac;
            double frac3 = frac2*frac;
            double frac4 = frac3*frac;
            double fracMinusOne = frac - 1;
            double fracMinusOne3 = fracMinusOne*fracMinusOne*fracMinusOne;

            Polynomial p;
            if (coord == coordIdx){
                // derivative
                p[0] = 5*frac4 - 10*frac3 + 4.5*frac2 + frac - 0.5;
                p[1] = -15*frac4 + 30*frac3 - 13.5*frac2 - 2*frac;
                p[2] = 15*frac4 - 30*frac3 + 13.5*frac2 + frac + 0.5;
                p[3] = frac2*(-5*frac2 + 10*frac - 4.5);
            } else {
                // 'normal' spline function
                p[0] =  0.5 * fracMinusOne3*frac*(2*frac + 1);
                p[1] = -0.5 * (frac - 1)*(6*frac4 - 9*frac3 + 2*frac + 2);
                p[2] =  0.5 * frac*(6*frac4 - 15*frac3 + 9*frac2 + frac + 1);
                p[3] = -0.5 * (frac - 1)*frac3*(2*frac - 3);
            }

            closestDiscretizationSteps[coord] = idx;
            betas[coord] = p;
        }

        std::array<int, g_MaxCoordsThatCanBeInterpolated> dimIdxOffsets;
        for (int coord = 0; coord < nCoords; ++coord) {
            dimIdxOffsets[coord] = -1;
        }

        double z = 0.0;
        int cnt = 0;
        while (dimIdxOffsets[0] < 3) {

            double beta = 1.0;
            int evalStride = 1;
            int evalIdx = 0;

            for (int coord = nCoords-1; coord >= 0; --coord) {
                int offset = dimIdxOffsets[coord];  // -1, 0, 1, or 2
                int closestStep = closestDiscretizationSteps[coord];
                int step = closestStep + offset;

                beta *= betas[coord][offset+1];
                evalIdx += evalStride * step;
                evalStride *= impl.discretizations[coord].nsteps;
            }

            double gridSize = (impl.discretizations[coordIdx].end-impl.discretizations[coordIdx].begin)/impl.discretizations[coordIdx].nsteps;
            z = std::fma(beta, impl.evals.at(evalIdx)/gridSize, z);
            {
                int pos = nCoords-1;
                ++dimIdxOffsets[pos];  // perform least-significant increment (may overflow)
                while (pos > 0 && dimIdxOffsets[pos] > 2) {  // handle overflows + carry propagation
                    dimIdxOffsets[pos] = -1;  // overflow
                    ++dimIdxOffsets[pos-1];  // carry propagation
                    --pos;
                }
            }

            ++cnt;
        }
        {
            int expectedIterations = 1 << (2*nCoords);
            if (cnt != expectedIterations) {
                std::stringstream msg;
                msg << "invalid number of permutations explored: expected = " << expectedIterations << ", got = " << cnt;
                OPENSIM_THROW(OpenSim::Exception, std::move(msg).str());
            }
        }

        return z;
    }

    // get the *derivative* of the path length with respect to the given Coordinate
    double Impl_GetPathLengthDerivative(JorisFBP const& impl,
                                        SimTK::State const& s,
                                        OpenSim::Coordinate const& c)
    {
        // figure out the index of the coordinate being referred to
        int coordIdx = -1;
        for (int i = 0; i < static_cast<int>(impl.coords.size()); ++i) {
            if (impl.coords[i] == &c) {
                coordIdx = i;
                break;
            }
        }

        // ensure the coordinate was actually found, or this alg will break
        if (coordIdx == -1) {
            std::stringstream msg;
            msg << "could not find coordiante '" << c.getName() << "' in the set of coordinates the FunctionBasedPath handles. Coordinates handled by this path are: ";
            char const* delim = "";
            for (OpenSim::Coordinate const* c : impl.coords) {
                msg << delim << c->getName();
                delim = ", ";
            }
            OPENSIM_THROW(OpenSim::Exception, std::move(msg).str());
        }

        // use the "raw" (non-lookup) version of this function with the index
        return Impl_GetPathLengthDerivative(impl, s, coordIdx);
    }

    // get the lengthening speed of the path in the current state
    double Impl_GetLengtheningSpeed(const JorisFBP& impl,
                                    const SimTK::State& state)
    {
        double lengtheningSpeed = 0.0;
        for (int coordIdx = 0; coordIdx < static_cast<int>(impl.coords.size()); ++coordIdx) {
            double deriv = Impl_GetPathLengthDerivative(impl, state, coordIdx);
            double coordSpeedVal = impl.coords[coordIdx]->getSpeedValue(state);

            lengtheningSpeed = std::fma(deriv, coordSpeedVal, lengtheningSpeed);
        }
        return lengtheningSpeed;
    }

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

    JorisFBP::JorisFBP()
    {
        constructProperty_FunctionBasedPathDiscretizationSet(FunctionBasedPathDiscretizationSet{});
        constructProperty_Evaluations();
    }

    double JorisFBP::getLength(const SimTK::State& s) const
    {
        return Impl_GetPathLength(*this, s);
    }

    double JorisFBP::getLengtheningSpeed(const SimTK::State& s) const
    {
        return Impl_GetLengtheningSpeed(*this, s);
    }

    double JorisFBP::computeMomentArm(const SimTK::State& s, const OpenSim::Coordinate& aCoord) const
    {
        return Impl_GetPathLengthDerivative(*this, s, aCoord);
    }

    void JorisFBP::addInEquivalentForces(const SimTK::State& state, double tension, SimTK::Vector_<SimTK::SpatialVec>&, SimTK::Vector& mobilityForces) const
    {
        const SimTK::SimbodyMatterSubsystem& matter = getModel().getMatterSubsystem();

        for (const OpenSim::Coordinate* coord :  coords) {
            double ma = computeMomentArm(state, *coord);
            double torqueOverCoord = -tension*ma;

            matter.addInMobilityForce(state,
                                      SimTK::MobilizedBodyIndex(coord->getBodyIndex()),
                                      SimTK::MobilizerUIndex(coord->getMobilizerQIndex()),
                                      torqueOverCoord,
                                      mobilityForces);
        }
    }

    void JorisFBP::extendFinalizeFromProperties()
    {
        Impl_InitFromFBPProperties(*this);
    }

    void JorisFBP::extendFinalizeConnections(OpenSim::Component& root)
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

        Impl_SetCoordinatePointersFromCoordinatePaths(*this, root);
    }

    OSIM_TEST(JorisFBP, CanBeDefaultConstructed)
    {
        JorisFBP fbp;
    }

    // TODO: port `Applications/FunctionBasedPathConversion/test/testFunctionBasedPathConversion.cpp`
    // see: `https://github.com/joris997/opensim-core/tree/interpolation/Applications/FunctionBasedPathConversion`
}

int main()
{
    return RunAllTests();
}
