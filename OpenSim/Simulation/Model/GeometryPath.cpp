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

//=============================================================================
// INCLUDES
//=============================================================================
#include "GeometryPath.h"
#include "ConditionalPathPoint.h"
#include "MovingPathPoint.h"
#include "PointForceDirection.h"
#include <OpenSim/Simulation/Wrap/PathWrap.h>
#include "Model.h"

//=============================================================================
// STATICS
//=============================================================================
using namespace std;
using namespace OpenSim;
using namespace SimTK;
using SimTK::Vec3;

//=============================================================================
// CONSTRUCTOR(S) AND DESTRUCTOR
//=============================================================================
//_____________________________________________________________________________
/*
 * Default constructor.
 */
GeometryPath::GeometryPath() :
    ModelComponent(),
    _preScaleLength(0.0)
{
    setAuthors("Peter Loan");
    constructProperties();
 }

//_____________________________________________________________________________
/*
* Perform set up functions after model has been deserialized or copied.
*
*/
void GeometryPath::extendFinalizeFromProperties()
{
    Super::extendFinalizeFromProperties();

    for (int i = 0; i < get_PathWrapSet().getSize(); ++i) {
        if (upd_PathWrapSet()[i].getName().empty()) {
            std::stringstream label;
            label << "pathwrap_" << i;
            upd_PathWrapSet()[i].setName(label.str());
        }
    }
}

void GeometryPath::extendConnectToModel(Model& aModel)
{
    Super::extendConnectToModel(aModel);

    OPENSIM_THROW_IF_FRMOBJ(get_PathPointSet().getSize() < 2,
        InvalidPropertyValue,
        getProperty_PathPointSet().getName(),
        "A valid path must be connected to a model by at least two PathPoints.")

    // Name the path points based on the current path
    // (i.e., the set of currently active points is numbered
    // 1, 2, 3, ...).
    namePathPoints(0);
}

//_____________________________________________________________________________
/*
 * Create the SimTK state, discrete and/or cache for this GeometryPath.
 */
 void GeometryPath::extendAddToSystem(SimTK::MultibodySystem& system) const 
{
    Super::extendAddToSystem(system);

    // Allocate cache entries to save the current length and speed(=d/dt length)
    // of the path in the cache. Length depends only on q's so will be valid
    // after Position stage, speed requires u's also so valid at Velocity stage.
    addCacheVariable<double>("length", 0.0, SimTK::Stage::Position);
    addCacheVariable<double>("speed", 0.0, SimTK::Stage::Velocity);
    // Cache the set of points currently defining this path.
    Array<AbstractPathPoint *> pathPrototype;
    addCacheVariable<Array<AbstractPathPoint *> >
        ("current_path", pathPrototype, SimTK::Stage::Position);

    // We consider this cache entry valid any time after it has been created
    // and first marked valid, and we won't ever invalidate it.
    addCacheVariable<SimTK::Vec3>("color", get_Appearance().get_color(), 
                                  SimTK::Stage::Topology);
}

 void GeometryPath::extendInitStateFromProperties(SimTK::State& s) const
{
    Super::extendInitStateFromProperties(s);
    markCacheVariableValid(s, "color"); // it is OK at its default value
}

//------------------------------------------------------------------------------
//                         GENERATE DECORATIONS
//------------------------------------------------------------------------------
// The GeometryPath takes care of drawing itself here, using information it
// can extract from the supplied state, including position information and
// color information that may have been calculated as late as Stage::Dynamics.
// For example, muscles may want the color to reflect activation level and 
// other path-using components might want to use forces (tension). We will
// ensure that the state has been realized to Stage::Dynamics before looking
// at it. (It is only guaranteed to be at Stage::Position here.)
void GeometryPath::
generateDecorations(bool fixed, const ModelDisplayHints& hints, 
                    const SimTK::State& state, 
                    SimTK::Array_<SimTK::DecorativeGeometry>& appendToThis) const
{        
    // There is no fixed geometry to generate here.
    if (fixed) { return; }

    const Array<AbstractPathPoint*>& pathPoints = getCurrentPath(state);

    assert(pathPoints.size() > 1);

    const AbstractPathPoint* lastPoint = pathPoints[0];
    MobilizedBodyIndex mbix(0);

    Vec3 lastPos = lastPoint->getLocationInGround(state);
    if (hints.get_show_path_points())
        DefaultGeometry::drawPathPoint(mbix, lastPos, getColor(state), appendToThis);

    Vec3 pos;

    for (int i = 1; i < pathPoints.getSize(); ++i) {
        AbstractPathPoint* point = pathPoints[i];
        PathWrapPoint* pwp = dynamic_cast<PathWrapPoint*>(point);

        if (pwp) {
            // A PathWrapPoint provides points on the wrapping surface as Vec3s
            Array<Vec3>& surfacePoints = pwp->getWrapPath();
            // The surface points are expressed w.r.t. the wrap surface's body frame.
            // Transform the surface points into the ground reference frame to draw
            // the surface point as the wrapping portion of the GeometryPath
            const Transform& X_BG = pwp->getParentFrame().getTransformInGround(state);
            // Cycle through each surface point and draw it the Ground frame
            for (int j = 0; j<surfacePoints.getSize(); ++j) {
                // transform the surface point into the Ground reference frame
                pos = X_BG*surfacePoints[j];
                if (hints.get_show_path_points())
                    DefaultGeometry::drawPathPoint(mbix, pos, getColor(state),
                        appendToThis);
                // Line segments will be in ground frame
                appendToThis.push_back(DecorativeLine(lastPos, pos)
                    .setLineThickness(4)
                    .setColor(getColor(state)).setBodyId(0).setIndexOnBody(j));
                lastPos = pos;
            }
        } 
        else { // otherwise a regular PathPoint so just draw its location
            pos = point->getLocationInGround(state);
            if (hints.get_show_path_points())
                DefaultGeometry::drawPathPoint(mbix, pos, getColor(state),
                    appendToThis);
            // Line segments will be in ground frame
            appendToThis.push_back(DecorativeLine(lastPos, pos)
                .setLineThickness(4)
                .setColor(getColor(state)).setBodyId(0).setIndexOnBody(i));
            lastPos = pos;
        }
    }
}

//_____________________________________________________________________________
/*
 * Connect properties to local pointers.
 */
void GeometryPath::constructProperties()
{
    constructProperty_PathPointSet(PathPointSet());

    constructProperty_PathWrapSet(PathWrapSet());
    
    Appearance appearance;
    appearance.set_color(SimTK::Gray);
    constructProperty_Appearance(appearance);
}

//_____________________________________________________________________________
/*
 * Name the path points based on their position in the set. To keep the
 * names up to date, this method should be called every time the path changes.
 *
 * @param aStartingIndex The index of the first path point to name.
 */
void GeometryPath::namePathPoints(int aStartingIndex)
{
    char indx[5];
    for (int i = aStartingIndex; i < get_PathPointSet().getSize(); i++)
    {
        sprintf(indx,"%d",i+1);
        AbstractPathPoint& point = get_PathPointSet().get(i);
        if (point.getName()=="" && hasOwner()) {
            point.setName(getOwner().getName() + "-P" + indx);
        }
    }
}

//_____________________________________________________________________________
/*
 * get the current path of the path
 *
 * @return The array of currently active path points.
 * 
 */
const OpenSim::Array <AbstractPathPoint*> & GeometryPath::
getCurrentPath(const SimTK::State& s)  const
{
    computePath(s);   // compute checks if path needs to be recomputed
    return getCacheVariableValue< Array<AbstractPathPoint*> >(s, "current_path");
}

// get the path as PointForceDirections directions 
// CAUTION: the return points are heap allocated; you must delete them yourself! 
// (TODO: that is really lame)
void GeometryPath::
getPointForceDirections(const SimTK::State& s, 
                        OpenSim::Array<PointForceDirection*> *rPFDs) const
{
    int i;
    AbstractPathPoint* start;
    AbstractPathPoint* end;
    const OpenSim::PhysicalFrame* startBody;
    const OpenSim::PhysicalFrame* endBody;
    const Array<AbstractPathPoint*>& currentPath = getCurrentPath(s);

    int np = currentPath.getSize();
    rPFDs->ensureCapacity(np);
    
    for (i = 0; i < np; i++) {
        PointForceDirection *pfd = 
            new PointForceDirection(currentPath[i]->getLocation(s), 
                                    currentPath[i]->getParentFrame(), Vec3(0));
        rPFDs->append(pfd);
    }

    for (i = 0; i < np-1; i++) {
        start = currentPath[i];
        end = currentPath[i+1];
        startBody = &start->getParentFrame();
        endBody = &end->getParentFrame();

        if (startBody != endBody)
        {
            Vec3 posStart, posEnd;
            Vec3 direction(0);

            // Find the positions of start and end in the inertial frame.
            //engine.getPosition(s, start->getParentFrame(), start->getLocation(), posStart);
            posStart = start->getLocationInGround(s);
            
            //engine.getPosition(s, end->getParentFrame(), end->getLocation(), posEnd);
            posEnd = end->getLocationInGround(s);

            // Form a vector from start to end, in the inertial frame.
            direction = (posEnd - posStart);

            // Check that the two points are not coincident.
            // This can happen due to infeasible wrapping of the path,
            // when the origin or insertion enters the wrapping surface.
            // This is a temporary fix, since the wrap algorithm should
            // return NaN for the points and/or throw an Exception- aseth
            if (direction.norm() < SimTK::SignificantReal){
                direction = direction*SimTK::NaN;
            }
            else{
                direction = direction.normalize();
            }

            // Get resultant direction at each point 
            rPFDs->get(i)->addToDirection(direction);
            rPFDs->get(i+1)->addToDirection(-direction);
        }
    }
}

/* add in the equivalent spatial forces on bodies for an applied tension 
    along the GeometryPath to a set of bodyForces */
void GeometryPath::addInEquivalentForces(const SimTK::State& s,
    const double& tension, 
    SimTK::Vector_<SimTK::SpatialVec>& bodyForces,
    SimTK::Vector& mobilityForces) const
{
    AbstractPathPoint* start = NULL;
    AbstractPathPoint* end = NULL;
    const SimTK::MobilizedBody* bo = NULL;
    const SimTK::MobilizedBody* bf = NULL;
    const Array<AbstractPathPoint*>& currentPath = getCurrentPath(s);
    int np = currentPath.getSize();

    const SimTK::SimbodyMatterSubsystem& matter = 
                                        getModel().getMatterSubsystem();

    // start point, end point,  direction, and force vectors in ground
    Vec3 po(0), pf(0), dir(0), force(0);
    // partial velocity of point in body expressed in ground 
    Vec3 dPodq_G(0), dPfdq_G(0);

    // gen force (torque) due to moving point under tension
    double fo, ff;

    for (int i = 0; i < np-1; ++i) {
        start = currentPath[i];
        end = currentPath[i+1];

        bo = &start->getParentFrame().getMobilizedBody();
        bf = &end->getParentFrame().getMobilizedBody();

        if (bo != bf) {
            // Find the positions of start and end in the inertial frame.
            po = start->getLocationInGround(s);
            pf = end->getLocationInGround(s);

            // Form a vector from start to end, in the inertial frame.
            dir = (pf - po);

            // Check that the two points are not coincident.
            // This can happen due to infeasible wrapping of the path,
            // when the origin or insertion enters the wrapping surface.
            // This is a temporary fix, since the wrap algorithm should
            // return NaN for the points and/or throw an Exception- aseth
            if (dir.norm() < SimTK::SignificantReal){
                dir = dir*SimTK::NaN;
            }
            else{
                dir = dir.normalize();
            }

            force = tension*dir;

            const MovingPathPoint* mppo =
                dynamic_cast<MovingPathPoint *>(start);

            // do the same for the end point of this segment of the path
            const MovingPathPoint* mppf =
                dynamic_cast<MovingPathPoint *>(end);

            // add in the tension point forces to body forces
            if (mppo) {// moving path point location is a function of the state
                // transform of the frame of the point to the base mobilized body
                auto X_BF = mppo->getParentFrame().findTransformInBaseFrame();
                bo->applyForceToBodyPoint(s, X_BF*mppo->getLocation(s), force,
                    bodyForces);
            }
            else {
                // transform of the frame of the point to the base mobilized body
                auto X_BF = start->getParentFrame().findTransformInBaseFrame();
                bo->applyForceToBodyPoint(s, X_BF*start->getLocation(s), force,
                    bodyForces);
            }

            if (mppf) {// moving path point location is a function of the state
                // transform of the frame of the point to the base mobilized body
                auto X_BF = mppf->getParentFrame().findTransformInBaseFrame();
                bf->applyForceToBodyPoint(s, X_BF*mppf->getLocation(s), -force,
                    bodyForces);
            }
            else {
                // transform of the frame of the point to the base mobilized body
                auto X_BF = end->getParentFrame().findTransformInBaseFrame();
                bf->applyForceToBodyPoint(s, X_BF*end->getLocation(s), -force,
                    bodyForces);
            }

            // Now account for the work being done by virtue of the moving
            // path point motion relative to the body it is on
            if(mppo){
                // torque (genforce) contribution due to relative movement 
                // of a via point w.r.t. the body it is connected to.
                dPodq_G = bo->expressVectorInGroundFrame(s, start->getdPointdQ(s));
                fo = ~dPodq_G*force;            

                // get the mobilized body the coordinate is couple to.
                const SimTK::MobilizedBody& mpbod =
                    matter.getMobilizedBody(mppo->getXCoordinate().getBodyIndex());

                // apply the generalized (mobility) force to the coordinate's body
                mpbod.applyOneMobilityForce(s, 
                    mppo->getXCoordinate().getMobilizerQIndex(), 
                    fo, mobilityForces);
            }

            if(mppf){
                dPfdq_G = bf->expressVectorInGroundFrame(s, end->getdPointdQ(s));
                ff = ~dPfdq_G*(-force);

                // get the mobilized body the coordinate is couple to.
                const SimTK::MobilizedBody& mpbod =
                    matter.getMobilizedBody(mppf->getXCoordinate().getBodyIndex());

                mpbod.applyOneMobilityForce(s, 
                    mppf->getXCoordinate().getMobilizerQIndex(), 
                    ff, mobilityForces);
            }
        }       
    }
}

//_____________________________________________________________________________
/*
 * Update the geometric representation of the path.
 * The resulting geometry is maintained at the VisibleObject layer.
 * This function should not be made public. It is called internally
 * by compute() only when the path has changed.
 * 
 */
void GeometryPath::updateGeometry(const SimTK::State& s) const
{
    // Check if the current path needs to recomputed.
    computePath(s);
}

//=============================================================================
// GET
//=============================================================================
//-----------------------------------------------------------------------------
// LENGTH
//-----------------------------------------------------------------------------
//_____________________________________________________________________________
/*
 * Compute the total length of the path.
 *
 * @return Total length of the path.
 */
double GeometryPath::getLength( const SimTK::State& s) const
{
    computePath(s);  // compute checks if path needs to be recomputed
    return( getCacheVariableValue<double>(s, "length") );
}

void GeometryPath::setLength( const SimTK::State& s, double length ) const
{
    setCacheVariableValue<double>(s, "length", length); 
}

void GeometryPath::setColor(const SimTK::State& s, const SimTK::Vec3& color) const
{
    setCacheVariableValue<SimTK::Vec3>(s, "color", color);
}

Vec3 GeometryPath::getColor(const SimTK::State& s) const
{
    return getCacheVariableValue<SimTK::Vec3>(s, "color");
}

//_____________________________________________________________________________
/*
 * Compute the lengthening speed of the path.
 *
 * @return lengthening speed of the path.
 */
double GeometryPath::getLengtheningSpeed( const SimTK::State& s) const
{
    computeLengtheningSpeed(s);
    return getCacheVariableValue<double>(s, "speed");
}
void GeometryPath::setLengtheningSpeed( const SimTK::State& s, double speed ) const
{
    setCacheVariableValue<double>(s, "speed", speed);    
}

void GeometryPath::setPreScaleLength( const SimTK::State& s, double length ) {
    _preScaleLength = length;
}
double GeometryPath::getPreScaleLength( const SimTK::State& s) const {
    return _preScaleLength;
}

//=============================================================================
// UTILITY
//=============================================================================
//_____________________________________________________________________________
/*
 * Add a new path point, with default location, to the path.
 *
 * @param aIndex The position in the pathPointSet to put the new point in.
 * @param frame The frame to attach the point to.
 * @return Pointer to the newly created path point.
 */
AbstractPathPoint* GeometryPath::
addPathPoint(const SimTK::State& s, int aIndex, const PhysicalFrame& frame)
{
    PathPoint* newPoint = new PathPoint();
    newPoint->setParentFrame(frame);
    Vec3 newLocation(0.0); 
    // Note: placeNewPathPoint() returns a location by reference.
    // It computes a new location according to the index where the new path point 
    // will be inserted along the path(among the other path points).
    placeNewPathPoint(s, newLocation, aIndex, frame);
    // Now set computed new location into the newPoint
    newPoint->setLocation(newLocation);
    upd_PathPointSet().insert(aIndex, newPoint);

    // Rename the path points starting at this new one.
    namePathPoints(aIndex);

    // Update start point and end point in the wrap instances so that they
    // refer to the same path points they did before the new point
    // was added. These indices are 1-based.
    aIndex++;
    for (int i=0; i<get_PathWrapSet().getSize(); i++) {
        int startPoint = get_PathWrapSet().get(i).getStartPoint();
        int endPoint = get_PathWrapSet().get(i).getEndPoint();
        if (startPoint != -1 && aIndex <= startPoint)
            get_PathWrapSet().get(i).setStartPoint(s,startPoint + 1);
        if (endPoint != -1 && aIndex <= endPoint)
            get_PathWrapSet().get(i).setEndPoint(s,endPoint + 1);
    }

    return newPoint;
}

AbstractPathPoint* GeometryPath::
appendNewPathPoint(const std::string& proposedName, 
                   const PhysicalFrame& frame,
                   const SimTK::Vec3& aPositionOnBody)
{
    PathPoint* newPoint = new PathPoint();
    newPoint->setParentFrame(frame);
    newPoint->setName(proposedName);
    newPoint->setLocation(aPositionOnBody);
    upd_PathPointSet().adoptAndAppend(newPoint);

    return newPoint;
}

//_____________________________________________________________________________
/*
 * Determine an appropriate default XYZ location for a new path point.
 * Note that this method is internal and should not be called directly on a new 
 * point as the point is not really added to the path (done in addPathPoint() 
 * instead)
 * @param aOffset The XYZ location to be determined.
 * @param aIndex The position in the pathPointSet to put the new point in.
 * @param frame The body to attach the point to.
 */
void GeometryPath::
placeNewPathPoint(const SimTK::State& s, SimTK::Vec3& aOffset, int aIndex, 
                  const PhysicalFrame& frame)
{
    // The location of the point is determined by moving a 'distance' from 'base' 
    // along a vector from 'start' to 'end.' 'base' is the existing path point 
    // that is in or closest to the index aIndex. 'start' and 'end' are existing
    // path points--which ones depends on where the new point is being added. 
    // 'distance' is 0.5 for points added to the middle of a path (so the point
    // appears halfway between the two adjacent points), and 0.2 for points that
    // are added to either end of the path. If there is only one point in the 
    // path, the new point is put 0.01 units away in all three dimensions.
    if (get_PathPointSet().getSize() > 1) {
        int start, end, base;
        double distance;
        if (aIndex == 0) {
            start = 1;
            end = 0;
            base = end;
            distance = 0.2;
        } else if (aIndex >= get_PathPointSet().getSize()) {
            start = aIndex - 2;
            end = aIndex - 1;
            base = end;
            distance = 0.2;
        } else {
            start = aIndex;
            end = aIndex - 1;
            base = start;
            distance = 0.5;
        }

        const Vec3 startPt = get_PathPointSet()[start].getLocation(s);
        const Vec3 endPt = get_PathPointSet()[end].getLocation(s);
        const Vec3 basePt = get_PathPointSet()[base].getLocation(s);

        Vec3 startPt2 = get_PathPointSet()[start].getParentFrame()
            .findStationLocationInAnotherFrame(s, startPt, frame);

        Vec3 endPt2 = get_PathPointSet()[end].getParentFrame()
            .findStationLocationInAnotherFrame(s, endPt, frame);

        aOffset = basePt + distance * (endPt2 - startPt2);
    } else if (get_PathPointSet().getSize() == 1) {
        aOffset= get_PathPointSet()[0].getLocation(s) + 0.01;
    }
    else {  // first point, do nothing?
    }
}

//_____________________________________________________________________________
/*
 * See if a path point can be deleted. All paths must have at least two
 * active path points to define the path.
 *
 * @param aIndex The index of the point to delete.
 * @return Whether or not the point can be deleted.
 */
bool GeometryPath::canDeletePathPoint( int aIndex)
{
    // A path point can be deleted only if there would remain
    // at least two other fixed points.
    int numOtherFixedPoints = 0;
    for (int i = 0; i < get_PathPointSet().getSize(); i++) {
        if (i != aIndex) {
            if (!(  get_PathPointSet().get(i).getConcreteClassName()
                  ==("ConditionalPathPoint")))
                numOtherFixedPoints++;
        }
    }

    if (numOtherFixedPoints >= 2)
        return true;

    return false;
}

//_____________________________________________________________________________
/*
 * Delete a path point.
 *
 * @param aIndex The index of the point to delete.
 * @return Whether or not the point was deleted.
 */
bool GeometryPath::deletePathPoint(const SimTK::State& s, int aIndex)
{
    if (canDeletePathPoint(aIndex) == false)
        return false;

    upd_PathPointSet().remove(aIndex);

    // rename the path points starting at the deleted position
    namePathPoints(aIndex);

    // Update start point and end point in the wrap instances so that they
    // refer to the same path points they did before the point was
    // deleted. These indices are 1-based. If the point deleted is start
    // point or end point, the path wrap range is made smaller by one point.
    aIndex++;
    for (int i=0; i<get_PathWrapSet().getSize(); i++) {
        int startPoint = get_PathWrapSet().get(i).getStartPoint();
        int endPoint   = get_PathWrapSet().get(i).getEndPoint();

        if (   (startPoint != -1 && aIndex < startPoint) 
            || (startPoint > get_PathPointSet().getSize()))
            get_PathWrapSet().get(i).setStartPoint(s, startPoint - 1);

        if (   endPoint > 1 
            && aIndex <= endPoint 
            && (   (endPoint > startPoint) 
                || (endPoint > get_PathPointSet().getSize())))
            get_PathWrapSet().get(i).setEndPoint(s, endPoint - 1);
    }

    return true;
}

//_____________________________________________________________________________
/*
 * Replace a path point in the set with another point. The new one is made a
 * member of all the same groups as the old one, and is inserted in the same
 * place the old one occupied.
 *
 *  @param aOldPathPoint Path point to remove.
 *  @param aNewPathPoint Path point to add.
 */
bool GeometryPath::
replacePathPoint(const SimTK::State& s, AbstractPathPoint* aOldPathPoint, 
                 AbstractPathPoint* aNewPathPoint) 
{
    if (aOldPathPoint != NULL && aNewPathPoint != NULL) {
        int count = 0;
        int index = get_PathPointSet().getIndex(aOldPathPoint);
        // If you're switching from non-via to via, check to make sure that the
        // path will be left with at least 2 non-via points.
        ConditionalPathPoint* oldVia = 
            dynamic_cast<ConditionalPathPoint*>(aOldPathPoint);
        ConditionalPathPoint* newVia = 
            dynamic_cast<ConditionalPathPoint*>(aNewPathPoint);
        if (oldVia == NULL && newVia != NULL) {
            for (int i=0; i<get_PathPointSet().getSize(); i++) {
                if (i != index) {
                    if (dynamic_cast<ConditionalPathPoint*>
                                        (&get_PathPointSet().get(i)) == NULL)
                        count++;
                }
            }
        } else {
            count = 2;
        }
        if (count >= 2 && index >= 0) {
            upd_PathPointSet().set(index, aNewPathPoint, true);
            //computePath(s);
            return true;
        }
    }
    return false;
}

//_____________________________________________________________________________
/*
 * Create a new wrap instance and add it to the set.
 *
 * @param aWrapObject The wrap object to use in the new wrap instance.
 */
void GeometryPath::addPathWrap(WrapObject& aWrapObject)
{
    PathWrap* newWrap = new PathWrap();
    newWrap->setWrapObject(aWrapObject);
    newWrap->setMethod(PathWrap::hybrid);
    upd_PathWrapSet().adoptAndAppend(newWrap);
    finalizeFromProperties();
}

//_____________________________________________________________________________
/*
 * Move a wrap instance up in the list. Changing the order of wrap instances for
 * a path may affect how the path wraps over the wrap objects.
 *
 * @param aIndex The index of the wrap instance to move up.
 */
void GeometryPath::moveUpPathWrap(const SimTK::State& s, int aIndex)
{
    if (aIndex > 0) {
        // Make sure wrap object is not deleted by remove().
        upd_PathWrapSet().setMemoryOwner(false); 

        PathWrap& wrap = get_PathWrapSet().get(aIndex);
        upd_PathWrapSet().remove(aIndex);
        upd_PathWrapSet().insert(aIndex - 1, &wrap);
        upd_PathWrapSet().setMemoryOwner(true);
    }
}

//_____________________________________________________________________________
/*
 * Move a wrap instance down in the list. Changing the order of wrap instances
 * for a path may affect how the path wraps over the wrap objects.
 *
 * @param aIndex The index of the wrap instance to move down.
 */
void GeometryPath::moveDownPathWrap(const SimTK::State& s, int aIndex)
{
    if (aIndex < get_PathWrapSet().getSize() - 1) {
        // Make sure wrap object is not deleted by remove().
        upd_PathWrapSet().setMemoryOwner(false);

        PathWrap& wrap = get_PathWrapSet().get(aIndex);
        upd_PathWrapSet().remove(aIndex);
        upd_PathWrapSet().insert(aIndex + 1, &wrap);
        upd_PathWrapSet().setMemoryOwner(true);
    }
}

//_____________________________________________________________________________
/*
 * Delete a wrap instance.
 *
 * @param aIndex The index of the wrap instance to delete.
 */
void GeometryPath::deletePathWrap(const SimTK::State& s, int aIndex)
{
    upd_PathWrapSet().remove(aIndex);

}

//==============================================================================
// SCALING
//==============================================================================
void GeometryPath::
extendPreScale(const SimTK::State& s, const ScaleSet& scaleSet)
{
    Super::extendPreScale(s, scaleSet);
    setPreScaleLength(s, getLength(s));
}

void GeometryPath::
extendPostScale(const SimTK::State& s, const ScaleSet& scaleSet)
{
    Super::extendPostScale(s, scaleSet);
    computePath(s);
}

//--------------------------------------------------------------------------
// COMPUTATIONS
//--------------------------------------------------------------------------
//=============================================================================
// PATH, WRAPPING, AND MOMENT ARM
//=============================================================================
//_____________________________________________________________________________
/*
 * Calculate the current path.
 */
void GeometryPath::computePath(const SimTK::State& s) const
{
    if (isCacheVariableValid(s, "current_path"))  {
        return;
    }

    Array<AbstractPathPoint*>& currentPath =
        updCacheVariableValue<Array<AbstractPathPoint*> >(s, "current_path");

    // Clear the current path.
    currentPath.setSize(0);

    // Add the active fixed and moving via points to the path.
    const PathPointSet& points = get_PathPointSet();
    for (int i = 0; i < points.getSize(); i++) {
        AbstractPathPoint* p = &points[i];
        if (p->isActive(s)) {
            currentPath.append(p); // <--- !!!!BAD
        }
    }
  
    // Use the current path so far to check for intersection with wrap objects, 
    // which may add additional points to the path.
    applyWrapObjects(s, currentPath);
    calcLengthAfterPathComputation(s, currentPath);

    markCacheVariableValid(s, "current_path");
}

//_____________________________________________________________________________
/*
 * Compute lengthening speed of the path.
 */
void GeometryPath::computeLengtheningSpeed(const SimTK::State& s) const
{
    if (isCacheVariableValid(s, "speed"))
        return;

    const Array<AbstractPathPoint*>& currentPath = getCurrentPath(s);

    double speed = 0.0;
    
    for (int i = 0; i < currentPath.getSize() - 1; i++) {
        speed += currentPath[i]->calcSpeedBetween(s, *currentPath[i+1]);
    }

    setLengtheningSpeed(s, speed);
}

//_____________________________________________________________________________
/*
 * Apply the wrap objects to the current path.
 */
void GeometryPath::
applyWrapObjects(const SimTK::State& s, Array<AbstractPathPoint*>& path) const 
{
    const PathWrapSet& pathWraps = get_PathWrapSet();
    const int numPathWraps = pathWraps.getSize();

    if (numPathWraps < 1) {
        return;
    }

    // Set the initial order to be the order they are listed in the path.
    std::vector<int> order(numPathWraps);
    for (int i = 0; i < numPathWraps; i++) {
        order[i] = i;
    }

    // If there is only one wrap object, calculate the wrapping only once.
    // If there are two or more objects, perform up to 8 iterations where
    // the result from one wrap object is used as the starting point for
    // the next wrap.
    const int maxIterations = numPathWraps == 1 ? 1 : 8;
    double computedLength = std::numeric_limits<double>::infinity();

    for (int iteration = 0; iteration < maxIterations; iteration++) {
        std::vector<WrapObject::WrapAction> result(numPathWraps);

        for (int pathWrapIdx = 0; pathWrapIdx < numPathWraps; pathWrapIdx++) {
            PathWrap& pathWrap = pathWraps.get(order[pathWrapIdx]);
            const WrapObject& wrapObject = *pathWrap.getWrapObject();

            // First remove this object's wrapping points from the current path.
            for (int j = 0; j < path.getSize(); j++) {
                const PathWrapPoint* pwp = &pathWrap.getWrapPoint1();
                if (path.get(j) == pwp) {
                    path.remove(j); // remove the first wrap point
                    path.remove(j); // remove the second wrap point
                    break;
                }
            }

            // skip computing this wrap if it is not active
            if (!wrapObject.get_active()) {
                continue;
            }

            // Find the start + end (inclusive) range of points in the caller-provided `path` that:
            // - start on an active point
            // - end on an active point
            // - account for via points
            int firstActivePointIdx = -1;
            int lastActivePointIdx = -1;
            {
                // 1. startPoint and endPoint are 1-based, where any number <1 is
                //    equivalent to the first (0) or last (PathPointSet.size() - 1)
                //    element.
                const int wrapStart = pathWrap.getStartPoint() > 0 ?
                                      pathWrap.getStartPoint() - 1 :
                                      0;
                const int wrapEnd = pathWrap.getEndPoint() > 0 ?
                                    pathWrap.getEndPoint() - 1 :
                                    get_PathPointSet().getSize() - 1;

                // 2. find the first point that is active in PathPointSet
                const AbstractPathPoint* firstActivePoint = nullptr;
                for (int i = wrapStart; i <= wrapEnd; i++) {
                    const AbstractPathPoint* p = &get_PathPointSet().get(i);
                    if (p->isActive(s)) {
                        firstActivePoint = p;
                        break;
                    }
                }
                if (firstActivePoint == nullptr) {
                    return;  // there are no active points in PathPointSet
                }

                // 3. find the last point that is active in PathPointSet
                const AbstractPathPoint* lastActivePoint = nullptr;
                for (int i = wrapEnd; i >= wrapStart; i--) {
                    const AbstractPathPoint* p = &get_PathPointSet().get(i);
                    if (p->isActive(s)) {
                        lastActivePoint = p;
                        break;
                    }
                }
                if (lastActivePoint == nullptr) {
                    return; // there are no active points in PathPointSet
                }

                // 4. Now find the indices of the active points, taken from PathPointSet,
                //    in the caller-provided `path` Array.
                for (int i = 0; i < path.getSize(); i++) {
                    const AbstractPathPoint* p = path.get(i);
                    if (p == firstActivePoint) {
                        firstActivePointIdx = i;
                    }
                    if (p == lastActivePoint) {
                        lastActivePointIdx = i;
                    }
                }
                if (firstActivePointIdx == -1 || lastActivePointIdx == -1) {
                    return;  // this should never happen.
                }
            }

            WrapResult bestWrap;
            bestWrap.wrap_pts.setSize(0);

            // Compute best wrapping for the PathWrap.
            //
            // Iterate over adjacent points (path segments) in the range and consider each path segment
            // for wrapping over the wrap object.
            //
            // The "best" wrap is the one that changes the path segment length the least.
            double smallestLengthChange = std::numeric_limits<double>::infinity();
            for (int pt1 = firstActivePointIdx; pt1 < lastActivePointIdx; pt1++) {
                AbstractPathPoint& p1 = *path.get(pt1);
                AbstractPathPoint& p2 = *path.get(pt1 + 1);

                // If the segments are auto-wrap points on same wrap object, skip
                // checking them for wrapping.
                if (p1.getWrapObject() != nullptr &&
                    p2.getWrapObject() != nullptr &&
                    p1.getWrapObject() != p2.getWrapObject()) {
                    continue;
                }

                WrapResult wr;
                wr.startPoint = pt1;
                wr.endPoint   = pt1 + 1;

                result[pathWrapIdx] = wrapObject.wrapPathSegment(s, p1, p2, pathWrap, wr);

                // mandatoryWrap:
                //
                // The segment actually intersected the wrap object. In the mandatory case,
                // the alg *must* choose this segment as the "best".
                //
                // note: this implementation means that only the first-intersecting segment
                //       is taken as the mandatory wrap (considered an ill-conditioned case).
                if (result[pathWrapIdx] == WrapObject::WrapAction::mandatoryWrap) {
                    bestWrap = wr;
                    pathWrap.setPreviousWrap(bestWrap);  // store best wrap for possible use next time

                    break; // stop trying other segments
                }

                // wrapped:
                //
                // The segment wrapped over the object. Unlike the mandatory case, the alg should
                // only choose this result if the length change is minimized by choosing it. Other
                // segments should also be attempted.
                if (result[pathWrapIdx] == WrapObject::WrapAction::wrapped) {
                    double lengthChange = calcPathLengthChange(s, wrapObject, wr, path);
                    if (lengthChange < smallestLengthChange) {
                        smallestLengthChange = lengthChange;

                        bestWrap = wr;
                        pathWrap.setPreviousWrap(wr);  // store best wrap for possible use next time
                    }
                    continue;  // keep trying other segments
                }
            }

            // Deallocate previous wrapping points if necessary.
            pathWrap.updWrapPoint2().getWrapPath().setSize(0);

            bool bestWrapFound = bestWrap.wrap_pts.getSize() > 0;
            if (bestWrapFound) {
                // If wrapping did occur, copy wrap info into the PathStruct.
                PathWrapPoint& p1 = pathWrap.updWrapPoint1();
                PathWrapPoint& p2 = pathWrap.updWrapPoint2();

                p1.getWrapPath().setSize(0);
                p2.getWrapPath() = bestWrap.wrap_pts;

                // In OpenSim, all conversion to/from the wrap object's
                // reference frame will be performed inside
                // wrapPathSegment(). Thus, all points in this function will
                // be in their respective body reference frames.
                // for (j = 0; j < wrapPath.getSize(); j++){
                //    convert_from_wrap_object_frame(wo, wrapPath.get(j));
                //    convert(ms->modelnum, wrapPath.get(j), wo->segment,
                //            ms->ground_segment);
                // }

                p1.setWrapLength(0.0);
                p2.setWrapLength(bestWrap.wrap_path_length);

                p1.setLocation(bestWrap.r1);
                p2.setLocation(bestWrap.r2);

                // Now insert the two new wrapping points into mp[] array.
                path.insert(bestWrap.endPoint, &p1);
                path.insert(bestWrap.endPoint + 1, &p2);
            } else {
                pathWrap.resetPreviousWrap();
                pathWrap.updWrapPoint2().getWrapPath().setSize(0);
            }
        }

        const double previousLength = computedLength;
        computedLength = calcLengthAfterPathComputation(s, path);

        if (std::abs(computedLength - previousLength) < 0.0005) {
            // the path's length did not change above some tolerance limit,
            // stop iterating (it's converged).
            break;
        }

        // edge-case:
        //
        // if the first wrap did not wrap, but the second one did not wrap
        // because a point was inside the object, switch the order of the first
        // two objects for the next iteration.
        if (iteration == 0 &&
            numPathWraps > 1 &&
            result[0] == WrapObject::WrapAction::noWrap &&
            result[1] == WrapObject::WrapAction::insideRadius) {

            std::swap(order[0], order[1]);

            // remove wrap object 0 from the list of path points
            for (int j = 0; j < path.getSize(); j++) {
                const PathWrapPoint* pwp = &pathWraps.get(0).getWrapPoint1();
                if (path.get(j) == pwp) {
                    path.remove(j); // remove the first wrap point
                    path.remove(j); // remove the second wrap point
                    break;
                }
            }
        }
    }
}

//_____________________________________________________________________________
/*
 * _calc_path_length_change - given the output of a successful path wrap
 * over a wrap object, determine the percent change in length of the
 * path segment incurred by wrapping.
 */
double GeometryPath::
calcPathLengthChange(const SimTK::State& s, const WrapObject& wo, 
                     const WrapResult& wr, const Array<AbstractPathPoint*>& path)  const
{
    const AbstractPathPoint* pt1 = path.get(wr.startPoint);
    const AbstractPathPoint* pt2 = path.get(wr.endPoint);

    double straight_length = pt1->calcDistanceBetween(s, *pt2);

    double wrap_length = pt1->calcDistanceBetween(s, wo.getFrame(), wr.r1);
    wrap_length += wr.wrap_path_length;
    wrap_length += pt2->calcDistanceBetween(s, wo.getFrame(), wr.r2);

    return wrap_length - straight_length; // return absolute diff, not relative
}

//_____________________________________________________________________________
/*
 * Compute the total length of the path. This function
 * assumes that the path has already been updated.
 */
double GeometryPath::
calcLengthAfterPathComputation(const SimTK::State& s, 
                               const Array<AbstractPathPoint*>& currentPath) const
{
    if (currentPath.getSize() < 2) {
        return 0.0;
    }

    double length = 0.0;
    for (int i = 0; i < currentPath.getSize() - 1; i++) {
        const AbstractPathPoint* p1 = currentPath[i];
        const AbstractPathPoint* p2 = currentPath[i+1];

        // If both points are wrap points on the same wrap object, then this
        // path segment wraps over the surface of a wrap object, so just add in 
        // the pre-calculated length.
        if (p1->getWrapObject() != nullptr &&
            p2->getWrapObject() != nullptr &&
            p1->getWrapObject() == p2->getWrapObject())
        {
            const PathWrapPoint* smwp = dynamic_cast<const PathWrapPoint*>(p2);
            if (smwp != nullptr) {
                length += smwp->getWrapLength();
            }
        } else {
            length += p1->calcDistanceBetween(s, *p2);
        }
    }

    setLength(s, length);
    return length;
}

//_____________________________________________________________________________
/*
 * Compute the path's moment arms for  specified coordinate.
 *
 * @param aCoord, the coordinate
 */   
double GeometryPath::
computeMomentArm(const SimTK::State& s, const Coordinate& aCoord) const
{
    if (!_maSolver)
        const_cast<Self*>(this)->_maSolver.reset(new MomentArmSolver(*_model));

    return _maSolver->solve(s, aCoord,  *this);
}

//_____________________________________________________________________________
// Override default implementation by object to intercept and fix the XML node
// underneath the model to match current version.
void GeometryPath::updateFromXMLNode(SimTK::Xml::Element& aNode, int versionNumber)
{
    if (versionNumber < XMLDocument::getLatestVersion()) {
        if (versionNumber < 30516) {
            // Create Appearance node under GeometryPath
            SimTK::Xml::Element appearanceElement("Appearance");
            aNode.appendNode(appearanceElement);
            SimTK::Xml::element_iterator visObjectIter = aNode.element_begin("VisibleObject");
            if (visObjectIter != aNode.element_end()) {
                SimTK::Xml::element_iterator oldPrefIter = visObjectIter->element_begin("display_preference");
                // old display_preference was used only for hide/show other options unsupported
                if (oldPrefIter != visObjectIter->element_end()) {
                    int oldPrefAsInt = 4;
                    oldPrefIter->getValueAs<int>(oldPrefAsInt);
                    if (oldPrefAsInt == 0) { // Hidden => Appearance/Visible
                        SimTK::Xml::Element visibleElement("visible");
                        visibleElement.setValueAs<bool>(false);
                        appearanceElement.insertNodeAfter(appearanceElement.element_end(), visibleElement);
                    }
                }
            }
            // If default_color existed, copy it to Appearance/color
            SimTK::Xml::element_iterator defaultColorIter = aNode.element_begin("default_color");
            if (defaultColorIter != aNode.element_end()) {
                // Move default_color to Appearance/color
                SimTK::Xml::Element colorElement("color");
                const SimTK::String& colorAsString = defaultColorIter->getValue();
                colorElement.setValue(colorAsString);
                appearanceElement.appendNode(colorElement);
            }
        }
    }
    // Call base class now assuming aNode has been corrected for current version
    Super::updateFromXMLNode(aNode, versionNumber);
}
