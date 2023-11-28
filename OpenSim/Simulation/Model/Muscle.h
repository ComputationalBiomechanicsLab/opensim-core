#ifndef OPENSIM_MUSCLE_H_
#define OPENSIM_MUSCLE_H_
/* -------------------------------------------------------------------------- *
 *                             OpenSim:  Muscle.h                             *
 * -------------------------------------------------------------------------- *
 * The OpenSim API is a toolkit for musculoskeletal modeling and simulation.  *
 * See http://opensim.stanford.edu and the NOTICE file for more information.  *
 * OpenSim is developed at Stanford University and supported by the US        *
 * National Institutes of Health (U54 GM072970, R24 HD065690) and by DARPA    *
 * through the Warrior Web program.                                           *
 *                                                                            *
 * Copyright (c) 2005-2017 Stanford University and the Authors                *
 * Author(s): Ajay Seth, Matthew Millard                                      *
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


// INCLUDE
#include "PathActuator.h"

#ifdef SWIG
    #ifdef OSIMSIMULATION_API
        #undef OSIMSIMULATION_API
        #define OSIMSIMULATION_API
    #endif
#endif

namespace OpenSim {

//==============================================================================
//                              Muscle Exceptions
//==============================================================================
class MuscleCannotEquilibrate : public Exception {
public:
    MuscleCannotEquilibrate(const std::string& file,
                            size_t line,
                            const std::string& func,
                            const Object& obj,
                            const std::string& detail) :
        Exception(file, line, func, obj, detail) {
        std::string msg = "Unable to compute equilibrium for this muscle.\n";
        msg += "Please verify that the initial activation is valid and that ";
        msg += "the length of the musculotendon actuator doesn't produce a ";
        msg += "pennation angle of 90 degrees or a negative fiber length.";
        addMessage(msg);
    }
};

//==============================================================================
//                                    Muscle
//==============================================================================
/**
 * A base class for modeling a muscle-tendon actuator. It defines muscle parameters
 * and methods to PathActuator, but does not implement all of the necessary methods,
 * and remains an abstract class. The path information for a muscle is contained
 * in PathActuator, and the force-generating behavior should be defined in
 * the derived classes.
 *
 * This class defines a subset of muscle models that include an active fiber
 * (contractile element) in series with a tendon. This class defines common 
 * data members and handles the geometry of a unipennate fiber in connection
 * with a tendon. No states are assumed, but concrete classes are free to
 * add whatever states are necessary to describe the specific behavior of a
 * muscle.
 *
 * @author Ajay Seth, Matt Millard
 *
 * (Based on earlier work by Peter Loan and Frank C. Anderson.)
 */
class OSIMSIMULATION_API Muscle : public PathActuator {
OpenSim_DECLARE_ABSTRACT_OBJECT(Muscle, PathActuator);
public:
//=============================================================================
// PROPERTIES
//=============================================================================
    OpenSim_DECLARE_PROPERTY(max_isometric_force, double,
        "Maximum isometric force that the fibers can generate");
    OpenSim_DECLARE_PROPERTY(optimal_fiber_length, double,
        "Optimal length of the muscle fibers");
    OpenSim_DECLARE_PROPERTY(tendon_slack_length, double,
        "Resting length of the tendon");
    OpenSim_DECLARE_PROPERTY(pennation_angle_at_optimal, double,
        "Angle between tendon and fibers at optimal fiber length expressed in radians");
    OpenSim_DECLARE_PROPERTY(max_contraction_velocity, double,
        "Maximum contraction velocity of the fibers, in optimal fiberlengths/second");
    OpenSim_DECLARE_PROPERTY(ignore_tendon_compliance, bool,
        "Compute muscle dynamics ignoring tendon compliance. Tendon is assumed to be rigid.");
    OpenSim_DECLARE_PROPERTY(ignore_activation_dynamics, bool,
        "Compute muscle dynamics ignoring activation dynamics. Activation is equivalent to excitation.");

//=============================================================================
// OUTPUTS
//=============================================================================
    OpenSim_DECLARE_OUTPUT(excitation, double, getExcitation,
            SimTK::Stage::Dynamics);
    OpenSim_DECLARE_OUTPUT(activation, double, getActivation,
            SimTK::Stage::Dynamics);
    OpenSim_DECLARE_OUTPUT(fiber_length, double, getFiberLength,
            SimTK::Stage::Position);
    OpenSim_DECLARE_OUTPUT(pennation_angle, double, getPennationAngle,
            SimTK::Stage::Position);
    OpenSim_DECLARE_OUTPUT(cos_pennation_angle, double, getCosPennationAngle,
            SimTK::Stage::Position);
    OpenSim_DECLARE_OUTPUT(tendon_length, double, getTendonLength,
            SimTK::Stage::Position);
    OpenSim_DECLARE_OUTPUT(normalized_fiber_length, double,
            getNormalizedFiberLength, SimTK::Stage::Position);
    OpenSim_DECLARE_OUTPUT(fiber_length_along_tendon, double,
            getFiberLengthAlongTendon, SimTK::Stage::Position);
    OpenSim_DECLARE_OUTPUT(tendon_strain, double, getTendonStrain,
            SimTK::Stage::Position);
    OpenSim_DECLARE_OUTPUT(passive_force_multiplier, double,
            getPassiveForceMultiplier, SimTK::Stage::Position);
    OpenSim_DECLARE_OUTPUT(active_force_length_multiplier, double,
            getActiveForceLengthMultiplier, SimTK::Stage::Position);
    OpenSim_DECLARE_OUTPUT(fiber_velocity, double, getFiberVelocity,
            SimTK::Stage::Velocity);
    OpenSim_DECLARE_OUTPUT(normalized_fiber_velocity, double,
            getNormalizedFiberVelocity, SimTK::Stage::Velocity);
    OpenSim_DECLARE_OUTPUT(fiber_velocity_along_tendon, double,
            getFiberVelocityAlongTendon, SimTK::Stage::Velocity);
    OpenSim_DECLARE_OUTPUT(tendon_velocity, double, getTendonVelocity,
            SimTK::Stage::Velocity);
    OpenSim_DECLARE_OUTPUT(force_velocity_multiplier, double,
            getForceVelocityMultiplier, SimTK::Stage::Velocity);
    OpenSim_DECLARE_OUTPUT(pennation_angular_velocity, double,
            getPennationAngularVelocity, SimTK::Stage::Velocity);
    OpenSim_DECLARE_OUTPUT(fiber_force, double, getFiberForce,
            SimTK::Stage::Dynamics);
    OpenSim_DECLARE_OUTPUT(fiber_force_along_tendon, double,
            getFiberForceAlongTendon, SimTK::Stage::Dynamics);
    OpenSim_DECLARE_OUTPUT(active_fiber_force, double, getActiveFiberForce,
            SimTK::Stage::Dynamics);
    OpenSim_DECLARE_OUTPUT(passive_fiber_force, double, getPassiveFiberForce,
            SimTK::Stage::Dynamics);
    OpenSim_DECLARE_OUTPUT(active_fiber_force_along_tendon, double,
            getActiveFiberForceAlongTendon, SimTK::Stage::Dynamics);
    OpenSim_DECLARE_OUTPUT(passive_fiber_force_along_tendon, double,
            getPassiveFiberForceAlongTendon, SimTK::Stage::Dynamics);
    OpenSim_DECLARE_OUTPUT(tendon_force, double, getTendonForce,
            SimTK::Stage::Dynamics);
    OpenSim_DECLARE_OUTPUT(fiber_stiffness, double, getFiberStiffness,
            SimTK::Stage::Dynamics);
    OpenSim_DECLARE_OUTPUT(fiber_stiffness_along_tendon, double,
            getFiberStiffnessAlongTendon, SimTK::Stage::Dynamics);
    OpenSim_DECLARE_OUTPUT(tendon_stiffness, double, getTendonStiffness,
            SimTK::Stage::Dynamics);
    OpenSim_DECLARE_OUTPUT(muscle_stiffness, double, getMuscleStiffness,
            SimTK::Stage::Dynamics);
    OpenSim_DECLARE_OUTPUT(fiber_active_power, double, getFiberActivePower,
            SimTK::Stage::Dynamics);
    OpenSim_DECLARE_OUTPUT(fiber_passive_power, double, getFiberPassivePower,
            SimTK::Stage::Dynamics);
    OpenSim_DECLARE_OUTPUT(tendon_power, double, getTendonPower,
            SimTK::Stage::Dynamics);
    OpenSim_DECLARE_OUTPUT(muscle_power, double, getMusclePower,
            SimTK::Stage::Dynamics);

//=============================================================================
// PUBLIC METHODS
//=============================================================================
    /** @name Constructors and Destructor
     */ 
    //@{
    /** Default constructor. */
    Muscle();

    // default destructor, copy constructor and copy assignment
    
    //--------------------------------------------------------------------------
    // MUSCLE PARAMETER ACCESSORS
    //--------------------------------------------------------------------------
    /** @name Muscle Parameters Access Methods
     */ 
    //@{
    /** get/set the maximum isometric force (in N) that the fibers can generate */
    double getMaxIsometricForce() const; 
    void setMaxIsometricForce(double maxIsometricForce);

    /** get/set the optimal length (in m) of the muscle fibers (lumped as a single fiber) */
    double getOptimalFiberLength() const;
    void setOptimalFiberLength(double optimalFiberLength);

    /** get/set the resting (slack) length (in m) of the tendon that is in series with the muscle fiber */
    double getTendonSlackLength() const;
    void setTendonSlackLength(double tendonSlackLength);

    /** get/set the angle (in radians) between fibers at their optimal fiber length and the tendon */
    double getPennationAngleAtOptimalFiberLength() const;
    void setPennationAngleAtOptimalFiberLength(double pennationAngle);
    
    /** get/set the maximum contraction velocity of the fibers, in optimal fiber-lengths per second */
    double getMaxContractionVelocity() const;
    void setMaxContractionVelocity(double maxContractionVelocity);

    // End of Muscle Parameter Accessors.
    //@} 

    //--------------------------------------------------------------------------
    // MUSCLE STATE DEPENDENT ACCESSORS
    //--------------------------------------------------------------------------
    /** @name Muscle State Dependent Access Methods
     *  Get quantities of interest common to all muscles
     */ 
    //@{

    /** Get/set Modeling (runtime) option to ignore tendon compliance when 
    computing muscle dynamics. This does not directly modify the persistent
    property value. **/
    bool getIgnoreTendonCompliance(const SimTK::State& s) const;
    void setIgnoreTendonCompliance(SimTK::State& s, bool ignore) const;

    /** Get/set Modeling (runtime) option to ignore activation dynamics when 
    computing muscle dynamics. This does not directly modify the persistent
    property value. **/
    bool getIgnoreActivationDynamics(const SimTK::State& s) const;
    void setIgnoreActivationDynamics(SimTK::State& s, bool ignore) const;

    /** get the activation level of the muscle, which modulates the active force
        of the muscle and has a normalized (0 to 1) value 
        Note: method remains virtual to permit override by deprecated muscles. */
    virtual double getActivation(const SimTK::State& s) const;

    /** get the current working fiber length (m) for the muscle */
    double getFiberLength(const SimTK::State& s) const;
    /** get the current pennation angle (radians) between the fiber and tendon at the current fiber length  */
    double getPennationAngle(const SimTK::State& s) const;
    /** get the cosine of the current pennation angle (radians) between the fiber and tendon at the current fiber length  */
    double getCosPennationAngle(const SimTK::State& s) const;
    /** get the current tendon length (m)  given the current joint angles and fiber length */
    double getTendonLength(const SimTK::State& s) const;
    /** get the current normalized fiber length (fiber_length/optimal_fiber_length) */
    double getNormalizedFiberLength(const SimTK::State& s) const;
    /** get the current fiber length (m) projected (*cos(pennationAngle)) onto the tendon direction */
    double getFiberLengthAlongTendon(const SimTK::State& s) const;
    /** get the current tendon strain (delta_l/tendon_slack_length is dimensionless)  */
    double getTendonStrain(const SimTK::State& s) const;

    /** the potential energy (J) stored in the fiber due to its parallel elastic element */
    double getFiberPotentialEnergy(const SimTK::State& s) const;
    /** the potential energy (J) stored in the tendon */    
    double getTendonPotentialEnergy(const SimTK::State& s) const;
    /** the total potential energy (J) stored in the muscle */  
    double getMusclePotentialEnergy(const SimTK::State& s) const;
    
    /** get the passive fiber (parallel elastic element) force multiplier */
    double getPassiveForceMultiplier(const SimTK::State& s) const;
    /** get the active fiber (contractile element) force multiplier due to current fiber length */
    double getActiveForceLengthMultiplier(const SimTK::State& s) const;

    /** get current fiber velocity (m/s) positive is lengthening */
    double getFiberVelocity(const SimTK::State& s) const;
    /** get normalized fiber velocity. This is the fiber velocity in m/s divided by
    the maximum contraction velocity expressed in m/s; therefore, this quantity is
    dimensionless and generally lies in the range [-1, 1]. */
    double getNormalizedFiberVelocity(const SimTK::State& s) const;
    /** get the current fiber velocity (m/s) projected onto the tendon direction */
    double getFiberVelocityAlongTendon(const SimTK::State& s) const;
    /** get pennation angular velocity (radians/s) */
    double getPennationAngularVelocity(const SimTK::State& s) const;
    /** get the tendon velocity (m/s) positive is lengthening */
    double getTendonVelocity(const SimTK::State& s) const;
    /** get the dimensionless multiplier resulting from the fiber's force-velocity curve */
    double getForceVelocityMultiplier(const SimTK::State& s) const;

    /** get the current fiber force (N) applied to the tendon */
    double getFiberForce(const SimTK::State& s) const;
    /**get the force of the fiber (N/m) along the direction of the tendon*/
    double getFiberForceAlongTendon(const SimTK::State& s) const;
    /** get the current active fiber force (N) due to activation*force_length*force_velocity relationships */
    double getActiveFiberForce(const SimTK::State& s) const;
    /** get the total force applied by all passive elements in the fiber (N) */
    double getPassiveFiberForce(const SimTK::State& s) const;
    /** get the current active fiber force (N) projected onto the tendon direction */
    double getActiveFiberForceAlongTendon(const SimTK::State& s) const;
    /** get the total force applied by all passive elements in the fiber (N)
        projected onto the tendon direction */
    double getPassiveFiberForceAlongTendon(const SimTK::State& s) const;
    /** get the current tendon force (N) applied to bones */
    double getTendonForce(const SimTK::State& s) const;

    /** get the current fiber stiffness (N/m) defined as the partial derivative
        of fiber force with respect to fiber length */
    double getFiberStiffness(const SimTK::State& s) const;
    /**get the stiffness of the fiber (N/m) along the direction of the tendon,
    that is the partial derivative of the fiber force along the tendon with
    respect to small changes in fiber length along the tendon*/
    double getFiberStiffnessAlongTendon(const SimTK::State& s) const;
    /** get the current tendon stiffness (N/m) defined as the partial derivative
        of tendon force with respect to tendon length */
    double getTendonStiffness(const SimTK::State& s) const;
    /** get the current muscle stiffness (N/m) defined as the partial derivative
        of muscle force with respect to muscle length */
    double getMuscleStiffness(const SimTK::State& s) const;

    /** get the current active fiber power (W) */
    double getFiberActivePower(const SimTK::State& s) const;
    /** get the current passive fiber power (W) */
    double getFiberPassivePower(const SimTK::State& s) const;

    /** get the current tendon power (W) */
    double getTendonPower(const SimTK::State& s) const;
    /** get the current muscle power (W) */
    double getMusclePower(const SimTK::State& s) const;
    
    /** get the stress in the muscle (part of the Actuator interface as well) */
    double getStress(const SimTK::State& s) const override;
    
    /** set the excitation (control) for this muscle. NOTE if controllers are connected to the
        muscle and are adding in their controls, and setExcitation is called after the model's
        computeControls(), then setExcitation will override the controller values. If called 
        before computeControls, then controller value(s) are added to the excitation set here. */
    void setExcitation(SimTK::State& s, double excitation) const;
    double getExcitation(const SimTK::State& s) const;


    /** DEPRECATED: only for backward compatibility */
    virtual void setActivation(SimTK::State& s, double activation) const = 0;

    // End of Muscle's State Dependent Accessors.
    //@} 

    /** Actuator interface for a muscle computes the tension in the muscle
        and applied by the tendon to bones (i.e. not the fiber force) */
    double computeActuation(const SimTK::State& s) const override = 0;

    /** @name Muscle initialization 
     */ 
    //@{
    /** Find and set the equilibrium state of the muscle (if any) */
    void computeEquilibrium(SimTK::State& s) const override final {
        return computeInitialFiberEquilibrium(s);
    }
    // End of Muscle's State Dependent Accessors.
    //@} 

    ///@cond
    //--------------------------------------------------------------------------
    // Estimate the muscle force for a given activation based on a rigid tendon 
    // assumption and neglecting passive fiber force. This provides a linear 
    // relationship between activation and force. This is used by CMC and 
    // StaticOptimization to solve the muscle force redundancy problem.
    //--------------------------------------------------------------------------    
    virtual double calcInextensibleTendonActiveFiberForce(SimTK::State& s, 
                                                  double aActivation) const;
    ///@endcond
//=============================================================================
// PROTECTED METHODS
//=============================================================================
protected:
    struct MuscleLengthInfo;
    struct FiberVelocityInfo;
    struct MuscleDynamicsInfo;
    struct MusclePotentialEnergyInfo;

    /** Developer Access to intermediate values calculate by the muscle model */
    const MuscleLengthInfo& getMuscleLengthInfo(const SimTK::State& s) const;
    MuscleLengthInfo& updMuscleLengthInfo(const SimTK::State& s) const;

    const FiberVelocityInfo& getFiberVelocityInfo(const SimTK::State& s) const;
    FiberVelocityInfo& updFiberVelocityInfo(const SimTK::State& s) const;

    const MuscleDynamicsInfo& getMuscleDynamicsInfo(const SimTK::State& s) const;
    MuscleDynamicsInfo& updMuscleDynamicsInfo(const SimTK::State& s) const;

    const MusclePotentialEnergyInfo& getMusclePotentialEnergyInfo(const SimTK::State& s) const;
    MusclePotentialEnergyInfo& updMusclePotentialEnergyInfo(const SimTK::State& s) const;

    //--------------------------------------------------------------------------
    // CALCULATIONS
    //--------------------------------------------------------------------------
    /** @name Muscle State Dependent Calculations
     *  Developers must override these methods to implement the desired behavior
     *  of their muscle models. Unless you are augmenting the behavior
     *  of an existing muscle class or writing a new derived class, you do not
     *  have access to these methods. 
     */ 
    //@{
    /** calculate muscle's position related values such fiber and tendon lengths,
        normalized lengths, pennation angle, etc... */
    virtual void calcMuscleLengthInfo(const SimTK::State& s, 
        MuscleLengthInfo& mli) const;

    /** calculate muscle's fiber velocity and pennation angular velocity, etc... */
    virtual void calcFiberVelocityInfo(const SimTK::State& s, 
        FiberVelocityInfo& fvi) const;

    /** calculate muscle's active and passive force-length, force-velocity, 
        tendon force, relationships and their related values */
    virtual void  calcMuscleDynamicsInfo(const SimTK::State& s, 
        MuscleDynamicsInfo& mdi) const;

    /** calculate muscle's fiber and tendon potential energy */
    virtual void calcMusclePotentialEnergyInfo(const SimTK::State& s,
        MusclePotentialEnergyInfo& mpei) const;

    /** This function modifies the fiber length in the supplied state such that  
    the fiber and tendon are developing the same force, taking activation and 
    velocity into account. This routine can assume that the state contains a
    meaningful estimate of muscle activation, joint positions, and joint 
    velocities. For example, this can produce fiber lengths suited to 
    beginning a forward dynamics simulation. 
    computeFiberEquilibriumAtZeroVelocity(). */
    virtual void computeInitialFiberEquilibrium(SimTK::State& s) const = 0;

    // End of Muscle's State Related Calculations.
    //@} 

    //--------------------------------------------------------------------------
    // PARENT INTERFACES
    //--------------------------------------------------------------------------
    /** @name Interfaces imposed by parent classes
     */ 
    //@{

    /** Force interface applies tension to bodies, and Muscle also checks that 
        applied muscle tension is not negative. */
    void computeForce(const SimTK::State& state, 
                      SimTK::Vector_<SimTK::SpatialVec>& bodyForces, 
                      SimTK::Vector& generalizedForce) const override;

    /** Potential energy stored by the muscle */
    double computePotentialEnergy(const SimTK::State& state) const override;

    /** Override PathActuator virtual to calculate a preferred color for the 
    muscle path based on activation. **/
    SimTK::Vec3 computePathColor(const SimTK::State& state) const override;
    
    /** Model Component creation interface */
    void extendConnectToModel(Model& aModel) override;
    void extendAddToSystem(SimTK::MultibodySystem& system) const override;
    void extendSetPropertiesFromState(const SimTK::State &s) override;
    void extendInitStateFromProperties(SimTK::State& state) const override;
    
    // Update the display geometry attached to the muscle
    virtual void updateGeometry(const SimTK::State& s);
    // End of Interfaces imposed by parent classes.
    //@} 


private:
    void setNull();
    void constructProperties();
    void copyData(const Muscle &aMuscle);

    //--------------------------------------------------------------------------
    // Implement Object interface.
    //--------------------------------------------------------------------------
    /** Override of the default implementation to account for versioning. */
    void updateFromXMLNode(SimTK::Xml::Element& aNode, int versionNumber=-1) override;


//=============================================================================
// DATA
//=============================================================================
protected:

    /** The assumed fixed muscle-width from which the fiber pennation angle can
        be calculated. */
    double _muscleWidth;


    struct MuscleOutputLength {
        double fiberLength            = SimTK::NaN;
        double fiberLengthAlongTendon = SimTK::NaN;
        double normFiberLength        = SimTK::NaN;

        double tendonStrain     = SimTK::NaN;
        double tendonLength     = SimTK::NaN;
        double normTendonLength = SimTK::NaN;

        double pennationAngle    = SimTK::NaN;
        double cosPennationAngle = SimTK::NaN;
        double sinPennationAngle = SimTK::NaN;

        double fiberPassiveForceLengthMultiplier = SimTK::NaN;
        double fiberActiveForceLengthMultiplier  = SimTK::NaN;
    };

    struct MuscleOutputVelocity {
        double fiberVelocity            = SimTK::NaN;
        double fiberVelocityAlongTendon = SimTK::NaN;
        double normFiberVelocity        = SimTK::NaN;

        double pennationAngularVelocity = SimTK::NaN;
        double tendonVelocity           = SimTK::NaN;
        double normTendonVelocity       = SimTK::NaN;

        double fiberForceVelocityMultiplier = SimTK::NaN;
    };

    struct MuscleOutputForce {
        double activation            = SimTK::NaN;
        double fiberForce            = SimTK::NaN;
        double fiberForceAlongTendon = SimTK::NaN;
        double normFiberForce        = SimTK::NaN;
        double activeFiberForce      = SimTK::NaN;
        double passiveFiberForce     = SimTK::NaN;

        double tendonForce = SimTK::NaN;
        double normTendonForce = SimTK::NaN;

        double fiberStiffness = SimTK::NaN;
        double fiberStiffnessAlongTendon = SimTK::NaN;
        double tendonStiffness = SimTK::NaN;
        double muscleStiffness = SimTK::NaN;

        double fiberActivePower = SimTK::NaN;
        double fiberPassivePower = SimTK::NaN;
        double tendonPower = SimTK::NaN;
        double musclePower = SimTK::NaN;
    };

    virtual const MuscleOutputLength& calcOutputLength(const SimTK::State& s) const
    {
        return getMuscleLengthInfo(s);
    }

    struct MuscleLengthInfo : MuscleOutputLength {
        SimTK::Vector userDefinedLengthExtras;//NA        NA
        MuscleLengthInfo():
            MuscleOutputLength(),
            userDefinedLengthExtras(0, SimTK::NaN){}
        friend std::ostream& operator<<(std::ostream& o, 
            const MuscleLengthInfo& mli) {
            o << "Muscle::MuscleLengthInfo should not be serialized!" 
              << std::endl;
            return o;
        }
    };

    struct FiberVelocityInfo: MuscleOutputVelocity {
        SimTK::Vector userDefinedVelocityExtras;//NA                  NA

        FiberVelocityInfo(): 
            MuscleOutputVelocity(),
            userDefinedVelocityExtras(0,SimTK::NaN){};
        friend std::ostream& operator<<(std::ostream& o, 
            const FiberVelocityInfo& fvi) {
            o << "Muscle::FiberVelocityInfo should not be serialized!" 
              << std::endl;
            return o;
        }
    };

    virtual const MuscleOutputVelocity& calcOutputVelocity(const SimTK::State& s) const
    {
        return getFiberVelocityInfo(s);
    }

    struct MuscleDynamicsInfo: MuscleOutputForce {
        SimTK::Vector userDefinedDynamicsExtras; //NA          NA

        MuscleDynamicsInfo(): 
            MuscleOutputForce(),
            userDefinedDynamicsExtras(0, SimTK::NaN){};
        friend std::ostream& operator<<(std::ostream& o, 
            const MuscleDynamicsInfo& mdi) {
            o << "Muscle::MuscleDynamicsInfo should not be serialized!" 
              << std::endl;
            return o;
        }
    };

    virtual const MuscleOutputForce& calcMuscleOutputForce(const SimTK::State& s) const
    {
        return getMuscleDynamicsInfo(s);
    }

    /**
        MusclePotentialEnergyInfo contains quantities related to the potential
        energy of the muscle (fiber + tendon) complex.
        
        The function that populates this struct, calcMusclePotentialEnergyInfo, can
        be called when position information is known. This function is
        dependent on calcMuscleLengthInfo.

        NAME                     DIMENSION              UNITS
        fiberPotentialEnergy      force*distance         J (Nm)   [1]
        tendonPotentialEnergy     force*distance         J (Nm)   [2]
        musclePotentialEnergy     force*distance         J (Nm)   [3]

        userDefinedPotentialEnergyExtras                         [4]

        [4] This vector is left for the muscle modeler to populate with any
            computationally expensive quantities that are computed in 
            calcMusclePotentialEnergyInfo, that might be useful for others to
            access.

    */
    struct MusclePotentialEnergyInfo {              //DIMENSION             UNITS
        double fiberPotentialEnergy;     //force*distance    J (Nm)     
        double tendonPotentialEnergy;    //force*distance    J (Nm)     
        double musclePotentialEnergy;    //force*distance    J (Nm)

        SimTK::Vector userDefinedPotentialEnergyExtras;//NA                  NA

        MusclePotentialEnergyInfo(): 
            fiberPotentialEnergy(SimTK::NaN),
            tendonPotentialEnergy(SimTK::NaN),
            musclePotentialEnergy(SimTK::NaN), 
            userDefinedPotentialEnergyExtras(0,SimTK::NaN){};
        friend std::ostream& operator<<(std::ostream& o, 
            const MusclePotentialEnergyInfo& fvi) {
            o << "Muscle::MusclePotentialEnergyInfo should not be serialized!" 
              << std::endl;
            return o;
        }
    };


    /** to support deprecated muscles */
    double _maxIsometricForce;
    double _optimalFiberLength;
    double _pennationAngleAtOptimal;
    double _tendonSlackLength;

    mutable CacheVariable<Muscle::MuscleLengthInfo> _lengthInfoCV;
    mutable CacheVariable<Muscle::FiberVelocityInfo> _velInfoCV;
    mutable CacheVariable<Muscle::MuscleDynamicsInfo> _dynamicsInfoCV;
    mutable CacheVariable<Muscle::MusclePotentialEnergyInfo> _potentialEnergyInfoCV;

//=============================================================================
};  // END of class Muscle
//=============================================================================
//=============================================================================

} // end of namespace OpenSim

#endif // OPENSIM_MUSCLE_H_
