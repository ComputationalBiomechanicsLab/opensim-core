#include "FunctionBasedPath.h"

#include <OpenSim/Simulation/Model/Model.h>
#include <OpenSim/Simulation/Model/PointBasedPath.h>
#include <OpenSim/Simulation/SimbodyEngine/Coordinate.h>
#include <SimTKcommon/internal/State.h>

#include <array>
#include <memory>
#include <vector>

static constexpr size_t g_MaxCoordsThatCanBeInterpolated = 8;  // important: upper limit that's used for stack allocations
static constexpr int g_MaxCoordsThatCanAffectPathDefault = static_cast<int>(g_MaxCoordsThatCanBeInterpolated);
static constexpr int g_NumProbingDiscretizationsDefault = 8;
static constexpr double g_MinProbingMomentArmChangeDefault = 0.001;
static constexpr int g_NumDiscretizationStepsPerDimensionDefault = 8;

// returns `true` if changing the supplied `Coordinate` changes the moment arm
// of the supplied `PointBasedPath` (PBP)
static bool coordAffectsPBP(
    OpenSim::PointBasedPath const& pbp,
    OpenSim::Coordinate const& c,
    SimTK::State& state,
    int numProbingSteps,
    double minMomentArmChangeRequired) {

    bool initialLockedState = c.getLocked(state);
    double initialValue = c.getValue(state);

    c.setLocked(state, false);

    double start = c.getRangeMin();
    double end = c.getRangeMax();
    double step = (end - start) / numProbingSteps;

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
static std::vector<OpenSim::Coordinate const*> coordsThatAffectPBP(OpenSim::Model const& model,
                                                                   OpenSim::PointBasedPath const& pbp,
                                                                   SimTK::State& st,
                                                                   int numProbingSteps,
                                                                   double minMomentArmChangeRequired) {
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
static Discretization discretizationForCoord(OpenSim::Coordinate const&, int numDiscretizationSteps) {
    // TODO: these hacky ranges were imported from the original code. Joris had this
    // commented out:
    //     dc.begin = std::max(c.getRangeMin(), -static_cast<double>(SimTK_PI));
    //     dc.end = std::min(c.getRangeMax(), static_cast<double>(SimTK_PI));
    SimTK_ASSERT_ALWAYS(numDiscretizationSteps >= 4, "need to supply more than 4 discretization steps");

    Discretization d;
    d.begin = -static_cast<double>(SimTK_PI)/2;
    d.end = static_cast<double>(SimTK_PI)/2;
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
static std::vector<double> computeEvaluationsFromPBP(OpenSim::PointBasedPath const& pbp,
                                                    SimTK::State& state,
                                                    OpenSim::Coordinate const** coords,
                                                    Discretization const* discs,
                                                    size_t ncoords) {
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

            double stepSz = (discr.end - discr.begin) / discr.nsteps;
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

struct OpenSim::FunctionBasedPath::Impl final {
    // direct pointers to each coordinate
    std::vector<OpenSim::Coordinate const*> coords;

    // absolute paths of each coordinate (1:1 with coords)
    std::vector<std::string> coordAbsPaths;

    // discretizations ranges for each coordinate (1:1 with coords)
    std::vector<Discretization> discretizations;

    // evaluations for each permutation of coordinates' discretizations
    std::vector<double> evals;

    // satisfies clone-ability requirement of SimTK::ClonePtr
    Impl* clone() const {
        auto p = std::unique_ptr<Impl>{new Impl{*this}};
        return p.release();
    }
};

// compute fresh implementation data from an existing PointBasedPath by
// evaluating it and fitting it to a function-based curve
//
// returns false if too many/too little coordinates affect the path
static bool Impl_ComputeFromPBP(OpenSim::FunctionBasedPath::Impl& impl,
                                const OpenSim::Model& model,
                                const OpenSim::PointBasedPath& pbp,
                                const OpenSim::FunctionBasedPath::FittingParams& params) {

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
static void Impl_InitFromFBPProperties(OpenSim::FunctionBasedPath::Impl& impl,
                                       OpenSim::FunctionBasedPath const& fbp) {

    OpenSim::FunctionBasedPathDiscretizationSet const& discSet = fbp.getProperty_FunctionBasedPathDiscretizationSet().getValue();

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
        OpenSim::FunctionBasedPathDiscretization const& disc = discSet[i];
        Discretization d;
        d.begin = disc.getProperty_x_begin().getValue();
        d.end = disc.getProperty_x_end().getValue();
        d.nsteps = disc.getProperty_num_points().getValue();
        impl.discretizations.push_back(d);
    }

    // set `evals` from evaluations property
    auto const& evalsProp = fbp.getProperty_Evaluations();
    impl.evals.clear();
    impl.evals.reserve(evalsProp.size());
    for (int i = 0; i < evalsProp.size(); ++i) {
        impl.evals.push_back(evalsProp.getValue(i));
    }
}

// ensure that the OpenSim::Coordinate* pointers held in Impl are up-to-date
//
// the pointers are there to reduce runtime path lookups
static void Impl_SetCoordinatePointersFromCoordinatePaths(OpenSim::FunctionBasedPath::Impl& impl,
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
static double Impl_GetPathLength(OpenSim::FunctionBasedPath::Impl const& impl,
                                 double const* inputVals,
                                 int nCoords) {

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
static double Impl_GetPathLength(OpenSim::FunctionBasedPath::Impl const& impl,
                                 SimTK::State const& s) {

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
static double Impl_GetPathLengthDerivative(OpenSim::FunctionBasedPath::Impl const& impl,
                                           SimTK::State const& s,
                                           int coordIdx) {

    SimTK_ASSERT_ALWAYS(!impl.coords.empty(), "FBPs require at least one coordinate to affect the path");
    SimTK_ASSERT_ALWAYS(coordIdx != -1, "coord index must be valid");
    SimTK_ASSERT_ALWAYS(coordIdx < static_cast<int>(impl.coords.size()), "coord index must be valid");

    int nCoords = static_cast<int>(impl.coords.size());

    // get the input value of each coordinate in the current state
    std::array<double, g_MaxCoordsThatCanBeInterpolated> inputVals{};
    for (int coord = 0; coord < nCoords; ++coord) {
        inputVals[coord] = impl.coords[coord]->getValue(s);
    }

    // compute value at current point
    double v1 = Impl_GetPathLength(impl, inputVals.data(), nCoords);

    static constexpr double h = 0.000001;

    // alter the input value for the to-be-derived coordinate *slightly* and recompute
    inputVals[coordIdx] += h;
    double v2 = Impl_GetPathLength(impl, inputVals.data(), nCoords);

    // the derivative is how much the output changed when the input was altered
    // slightly (this is a poor-man's discrete derivative method)
    return (v2 - v1) / h;
}

// get the *derivative* of the path length with respect to the given Coordinate
static double Impl_GetPathLengthDerivative(OpenSim::FunctionBasedPath::Impl const& impl,
                                           SimTK::State const& s,
                                           OpenSim::Coordinate const& c) {

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
static double Impl_GetLengtheningSpeed(const OpenSim::FunctionBasedPath::Impl& impl,
                                       const SimTK::State& state) {

    double lengtheningSpeed = 0.0;
    for (int coordIdx = 0; coordIdx < static_cast<int>(impl.coords.size()); ++coordIdx) {
        double deriv = Impl_GetPathLengthDerivative(impl, state, coordIdx);
        double coordSpeedVal = impl.coords[coordIdx]->getSpeedValue(state);

        lengtheningSpeed = std::fma(deriv, coordSpeedVal, lengtheningSpeed);
    }
    return lengtheningSpeed;
}

//----------------------------------------------------------------------------
//                               PUBLIC API
//----------------------------------------------------------------------------

OpenSim::FunctionBasedPath::FittingParams::FittingParams() :
    maxCoordsThatCanAffectPath{g_MaxCoordsThatCanAffectPathDefault},
    numProbingDiscretizations{g_NumProbingDiscretizationsDefault},
    minProbingMomentArmChange{g_MinProbingMomentArmChangeDefault},
    numDiscretizationStepsPerDimension{g_NumDiscretizationStepsPerDimensionDefault}
{
}

std::unique_ptr<OpenSim::FunctionBasedPath> OpenSim::FunctionBasedPath::fromPointBasedPath(
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

    std::unique_ptr<FunctionBasedPath> fbp{new FunctionBasedPath{}};

    // copy relevant data from source PBP into the output FBP
    fbp->upd_Appearance() = pbp.get_Appearance();
    fbp->setPathPointSet(pbp.getPathPointSet());
    fbp->setPathWrapSet(pbp.getWrapSet());

    OpenSim::FunctionBasedPath::Impl& impl = *fbp->_impl;

    // compute underlying impl data from the PBP
    if (!Impl_ComputeFromPBP(impl, model, pbp, params)) {
        return nullptr;
    }

    // write impl discretizations into the `Discretizations` property
    OpenSim::FunctionBasedPathDiscretizationSet& set = fbp->updProperty_FunctionBasedPathDiscretizationSet().updValue();
    for (size_t i = 0; i < impl.coords.size(); ++i) {
        auto disc = std::unique_ptr<OpenSim::FunctionBasedPathDiscretization>{new FunctionBasedPathDiscretization{}};
        disc->set_x_begin(impl.discretizations[i].begin);
        disc->set_x_end(impl.discretizations[i].end);
        disc->set_num_points(impl.discretizations[i].nsteps);
        disc->set_coordinate_abspath(impl.coordAbsPaths[i]);
        set.adoptAndAppend(disc.release());
    }

    // write evals into `Evaluations` property
    auto& evalsProp = fbp->updProperty_Evaluations();
    for (double eval : fbp->_impl->evals) {
        evalsProp.appendValue(eval);
    }

    return fbp;
}

OpenSim::FunctionBasedPath::FunctionBasedPath() : GeometryPath{}, _impl{new Impl{}}
{
    constructProperty_FunctionBasedPathDiscretizationSet(FunctionBasedPathDiscretizationSet{});
    constructProperty_Evaluations();
}
OpenSim::FunctionBasedPath::FunctionBasedPath(const FunctionBasedPath&) = default;
OpenSim::FunctionBasedPath::FunctionBasedPath(FunctionBasedPath&&) = default;
OpenSim::FunctionBasedPath& OpenSim::FunctionBasedPath::operator=(const FunctionBasedPath&) = default;
OpenSim::FunctionBasedPath& OpenSim::FunctionBasedPath::operator=(FunctionBasedPath&&) = default;
OpenSim::FunctionBasedPath::~FunctionBasedPath() noexcept = default;

void OpenSim::FunctionBasedPath::extendFinalizeFromProperties()
{
    Impl_InitFromFBPProperties(*_impl, *this);
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

    Impl_SetCoordinatePointersFromCoordinatePaths(*_impl, root);
}

double OpenSim::FunctionBasedPath::getLength(const SimTK::State& s) const
{
    if (isCacheVariableValid(s, _lengthCV)) {
        return getCacheVariableValue(s, _lengthCV);
    }

    double v = Impl_GetPathLength(*_impl, s);
    setCacheVariableValue(s, _lengthCV, v);
    return v;
}

void OpenSim::FunctionBasedPath::setLength(const SimTK::State& s, double length) const
{
    setCacheVariableValue(s, _lengthCV, length);
}

double OpenSim::FunctionBasedPath::getLengtheningSpeed(const SimTK::State& s) const
{
    if (isCacheVariableValid(s, _speedCV)) {
        return getCacheVariableValue(s, _speedCV);
    }

    double v = Impl_GetLengtheningSpeed(*_impl, s);
    setCacheVariableValue(s, _speedCV, v);
    return v;
}

void OpenSim::FunctionBasedPath::setLengtheningSpeed(const SimTK::State& s, double speed) const
{
    setCacheVariableValue(s, _speedCV, speed);
}

double OpenSim::FunctionBasedPath::computeMomentArm(const SimTK::State& s,
                                                    const Coordinate& aCoord) const
{
    return Impl_GetPathLengthDerivative(*_impl, s, aCoord);
}

/* add in the equivalent spatial forces on bodies for an applied tension
    along the GeometryPath to a set of bodyForces */
void OpenSim::FunctionBasedPath::addInEquivalentForces(const SimTK::State& s,
    const double& tension,
    SimTK::Vector_<SimTK::SpatialVec>& bodyForces,
    SimTK::Vector& mobilityForces) const
{
    const SimTK::SimbodyMatterSubsystem& matter = getModel().getMatterSubsystem();

    for (const OpenSim::Coordinate* coord : _impl->coords) {
        double ma = computeMomentArm(s, *coord);
        double torqueOverCoord = -tension*ma;

        matter.addInMobilityForce(s,
                                  SimTK::MobilizedBodyIndex(coord->getBodyIndex()),
                                  SimTK::MobilizerUIndex(coord->getMobilizerQIndex()),
                                  torqueOverCoord,
                                  mobilityForces);
    }
}

void OpenSim::FunctionBasedPath::generateDecorations(
        bool fixed,
        const ModelDisplayHints& hints,
        const SimTK::State& state,
        SimTK::Array_<SimTK::DecorativeGeometry>& appendToThis) const
{
    static bool shownOnce = []() {
        OpenSim::log_warn("Tried to call `generateDecorations` on a `FunctionBasedPath`. This function call will be ignored (there is no easy way to visualize function-based paths)");
        return true;
    }();
    (void)shownOnce;
}


void OpenSim::FunctionBasedPath::computePath(const SimTK::State& s) const
{
    OPENSIM_THROW(Exception, "Tried to call `computePath` on a `FunctionBasedPath`. You cannot call this method on a `GeometryPath` that is a `FunctionBasedPath` (it isn't path-based). Either remove this function call or replace the `FunctionBasedPath` with a `PointBasedPath` in the model");
}


void OpenSim::FunctionBasedPath::computeLengtheningSpeed(const SimTK::State& s) const
{
    double lengtheningspeed = getLengtheningSpeed(s);
    setLengtheningSpeed(s, lengtheningspeed);
}

double OpenSim::FunctionBasedPath::calcLengthAfterPathComputation(
        const SimTK::State& s,
        const Array<AbstractPathPoint*>& currentPath) const
{
    return getLength(s);
}

double OpenSim::FunctionBasedPath::calcPathLengthChange(const SimTK::State& s,
                                                        const WrapObject& wo,
                                                        const WrapResult& wr,
                                                        const Array<AbstractPathPoint*>& path) const
{
    OPENSIM_THROW(Exception, "Tried to call `calcPathLengthChange` on a `FunctionBasedPath`. You cannot call this method on a `GeometryPath` that is a `FunctionBasedPath` (it isn't path-based). Either remove this function call or replace the `FunctionBasedPath` with a `PointBasedPath` in the model");
}

void OpenSim::FunctionBasedPath::getPointForceDirections(const SimTK::State& s,
                                                         OpenSim::Array<PointForceDirection*> *rPFDs) const
{
    OPENSIM_THROW(Exception, "Tried to call `getPointForceDirections` on a `FunctionBasedPath`. You cannot call this method on a `GeometryPath` that is a `FunctionBasedPath` (it isn't path-based). Either remove this function call or replace the `FunctionBasedPath` with a `PointBasedPath` in the model");
}

const OpenSim::Array<OpenSim::AbstractPathPoint*>& OpenSim::FunctionBasedPath::getCurrentPath(const SimTK::State& s) const
{
    OPENSIM_THROW(Exception, "Tried to call `getCurrentPath` on a `FunctionBasedPath`. You cannot call this method on a `GeometryPath` that is a `FunctionBasedPath` (it does not contain a path). Either remove this function call or replace the `FunctionBasedPath` with a `PointBasedPath` in the model");
}

OpenSim::AbstractPathPoint* OpenSim::FunctionBasedPath::addPathPoint(
        const SimTK::State& s,
        int index,
        const PhysicalFrame& frame)
{
    OPENSIM_THROW(Exception, "Tried to call `addPathPoint` on a `FunctionBasedPath`. You cannot call this method on a `GeometryPath` that is a `FunctionBasedPath` (it does not contain a path). Either remove this function call or replace the `FunctionBasedPath` with a `PointBasedPath` in the model");
}

OpenSim::AbstractPathPoint* OpenSim::FunctionBasedPath::appendNewPathPoint(
        const std::string& proposedName,
        const PhysicalFrame& frame,
        const SimTK::Vec3& locationOnFrame)
{
    OPENSIM_THROW(Exception, "Tried to call `appendNewPathPoint` on a `FunctionBasedPath`. You cannot call this method on a `GeometryPath` that is a `FunctionBasedPath` (it does not contain a path). Either remove this function call or replace the `FunctionBasedPath` with a `PointBasedPath` in the model");
}

bool OpenSim::FunctionBasedPath::canDeletePathPoint(int index)
{
    OPENSIM_THROW(Exception, "Tried to call `canDeletePathPoint` on a `FunctionBasedPath`. You cannot call this method on a `GeometryPath` that is a `FunctionBasedPath` (it does not contain a path). Either remove this function call or replace the `FunctionBasedPath` with a `PointBasedPath` in the model");
}

bool OpenSim::FunctionBasedPath::deletePathPoint(const SimTK::State& s,
                                                 int index)
{
    OPENSIM_THROW(Exception, "Tried to call `deletePathPoint` on a `FunctionBasedPath`. You cannot call this method on a `GeometryPath` that is a `FunctionBasedPath` (it does not contain a path). Either remove this function call or replace the `FunctionBasedPath` with a `PointBasedPath` in the model");
}

bool OpenSim::FunctionBasedPath::replacePathPoint(const SimTK::State& s,
                                                  AbstractPathPoint* oldPathPoint,
                                                  AbstractPathPoint* newPathPoint)
{
    OPENSIM_THROW(Exception, "Tried to call `replacePathPoint` on a `FunctionBasedPath`. You cannot call this method on a `GeometryPath` that is a `FunctionBasedPath` (it does not contain a path). Either remove this function call or replace the `FunctionBasedPath` with a `PointBasedPath` in the model");
}

void OpenSim::FunctionBasedPath::addPathWrap(WrapObject& aWrapObject)
{
    OPENSIM_THROW(Exception, "Tried to call `addPathWrap` on a `FunctionBasedPath`. You cannot call this method on a `GeometryPath` that is a `FunctionBasedPath` (it does not contain wrapping objects). Either remove this function call or replace the `FunctionBasedPath` with a `PointBasedPath` in the model");
}

void OpenSim::FunctionBasedPath::deletePathWrap(const SimTK::State& s,
                                                int index)
{
    OPENSIM_THROW(Exception, "Tried to call `deletePathWrap` on a `FunctionBasedPath`. You cannot call this method on a `GeometryPath` that is a `FunctionBasedPath` (it does not contain wrapping objects). Either remove this function call or replace the `FunctionBasedPath` with a `PointBasedPath` in the model");
}

void OpenSim::FunctionBasedPath::moveUpPathWrap(const SimTK::State& s,
                                                int index)
{
    OPENSIM_THROW(Exception, "Tried to call `moveUpPathWrap` on a `FunctionBasedPath`. You cannot call this method on a `GeometryPath` that is a `FunctionBasedPath` (it does not contain wrapping objects). Either remove this function call or replace the `FunctionBasedPath` with a `PointBasedPath` in the model");
}

void OpenSim::FunctionBasedPath::moveDownPathWrap(const SimTK::State& s,
                                                  int index)
{
    OPENSIM_THROW(Exception, "Tried to call `moveDownPathWrap` on a `FunctionBasedPath`. You cannot call this method on a `GeometryPath` that is a `FunctionBasedPath` (it does not contain wrapping objects). Either remove this function call or replace the `FunctionBasedPath` with a `PointBasedPath` in the model");
}

void OpenSim::FunctionBasedPath::applyWrapObjects(const SimTK::State& s,
                                                  Array<AbstractPathPoint*>& path ) const
{
    OPENSIM_THROW(Exception, "Tried to call `applyWrapObjects` on a `FunctionBasedPath`. You cannot call this method on a `GeometryPath` that is a `FunctionBasedPath` (it doesn't know how to apply wrap objects to a sequence of path points). Either remove this function call or replace the `FunctionBasedPath` with a `PointBasedPath` in the model");
}

void OpenSim::FunctionBasedPath::namePathPoints(int aStartingIndex)
{
    static bool shownOnce = []() {
        OpenSim::log_warn("Tried to call `namePathPoints` on a `FunctionBasedPath`. This function call will be ignored (`FunctionBasedPath`s do not contain path points)");
        return true;
    }();
    (void)shownOnce;
}

void OpenSim::FunctionBasedPath::placeNewPathPoint(const SimTK::State& s,
                                                   SimTK::Vec3& aOffset,
                                                   int index,
                                                   const PhysicalFrame& frame)
{
    OPENSIM_THROW(Exception, "Tried to call `placeNwPathPoint` on a `FunctionBasedPath`. You cannot call this method on a `GeometryPath` that is a `FunctionBasedPath` (it has no path points). Either remove this function call or replace the `FunctionBasedPath` with a `PointBasedPath` in the model");
}
