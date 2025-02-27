#ifndef OPENSIM_MY_CABLE_H_
#define OPENSIM_MY_CABLE_H_

#include "OpenSim/Common/Component.h"
#include "OpenSim/Common/ComponentSocket.h"
#include "OpenSim/Common/Object.h"
#include "OpenSim/Common/Property.h"
#include "OpenSim/Common/Set.h"
#include "OpenSim/Simulation/Model/AbstractGeometryPath.h"
#include "OpenSim/Simulation/Model/Model.h"
#include "OpenSim/Simulation/Model/ModelComponent.h"
#include "OpenSim/Simulation/Model/PhysicalFrame.h"
#include "OpenSim/Simulation/Model/PhysicalOffsetFrame.h"
#include "OpenSim/Simulation/Model/Station.h"
#include "OpenSim/Simulation/Wrap/WrapCylinder.h"
#include "OpenSim/Simulation/Wrap/WrapObject.h"
#include "OpenSim/Simulation/osimSimulationDLL.h"
#include "PathPointSet.h"
#include <memory>
#include <simmath/internal/ContactGeometry.h>
#include <stdexcept>

#include <SimTKcommon/Scalar.h>
#include <SimTKcommon/internal/ExceptionMacros.h>

#include <simbody/internal/CableSpan.h>
#include <simbody/internal/MultibodySystem.h>
#include <simbody/internal/SimbodyMatterSubsystem.h>
#include <simbody/internal/common.h>

#include <OpenSim/Simulation/MomentArmSolver.h>
#include <OpenSim/Simulation/Wrap/PathWrapSet.h>

namespace OpenSim {

enum class NodeKind
{
    AttachmentPoint,
    Obstacle,
};

class OSIMSIMULATION_API ANode : public Component {
    OpenSim_DECLARE_ABSTRACT_OBJECT(ANode, Component);
    public:
    virtual NodeKind getKind() const = 0;
};

class OSIMSIMULATION_API NodeSet : public Set<ANode> {
    OpenSim_DECLARE_CONCRETE_OBJECT(NodeSet, Set<ANode>);

public:
    /** Use Super's constructors. @see Set */
    using Super::Super;
};

// An attachment point node.
class OSIMSIMULATION_API CableAttachmentPoint : public ANode {
    OpenSim_DECLARE_CONCRETE_OBJECT(CableAttachmentPoint, ANode);

public:
    OpenSim_DECLARE_UNNAMED_PROPERTY(Station, "TODO");

    CableAttachmentPoint() {
        constructProperties();
    }

    CableAttachmentPoint(
            const std::string& name,
            const PhysicalFrame& frame,
            const SimTK::Vec3& station
            ) {
        constructProperties();

        setName(name);
        upd_Station().setParentFrame(frame);
        upd_Station().set_location(station);
        std::cout << getName() << " --> socket " << get_Station().getParentFrame().getName() << "\n";
    }

    void constructProperties() {
        constructProperty_Station(Station());
    }

    NodeKind getKind() const override {
        return NodeKind::AttachmentPoint;
    }
};

// A surface obstacle node.
class OSIMSIMULATION_API CableObstacle : public ANode {
    OpenSim_DECLARE_CONCRETE_OBJECT(CableObstacle, ANode);

public:
    OpenSim_DECLARE_PROPERTY(contactPointHintInS, SimTK::Vec3,
            "A hint for the initial contact point.");

    OpenSim_DECLARE_SOCKET(wrapObject, WrapObject, "TODO");

    CableObstacle() { constructProperties(); }

    CableObstacle(const std::string& name, const WrapObject& wrapObject,
            const SimTK::Vec3& contactPointHintInS) {
        setName(name);
        constructProperties();
        set_contactPointHintInS(contactPointHintInS);
        connectSocket_wrapObject(wrapObject);
    }

    explicit CableObstacle(const WrapObject& wrapObject) {
        constructProperties();
        connectSocket_wrapObject(wrapObject);
    }

    void constructProperties() {
        constructProperty_contactPointHintInS(SimTK::Vec3{SimTK::NaN});
    }

    const WrapObject& getWrapObject() const {
        return getConnectee<WrapObject>("wrapObject");
    }

    NodeKind getKind() const override {
        return NodeKind::Obstacle;
    }
};

class OSIMSIMULATION_API MyCable : public AbstractGeometryPath {
    OpenSim_DECLARE_CONCRETE_OBJECT(MyCable, AbstractGeometryPath);

public:
    OpenSim_DECLARE_PROPERTY(accuracy, SimTK::Real, "TODO");
    OpenSim_DECLARE_PROPERTY(smoothness, SimTK::Real, "TODO");

    OpenSim_DECLARE_UNNAMED_PROPERTY(NodeSet, "TODO");

    MyCable() {
        constructProperties();
    }

    explicit MyCable(const std::string& name) {
        setName(name);
        constructProperties();
    }

    int getNumNodes() const
    {
        return get_NodeSet().getSize();
    }

    static Model convert(const std::string& modelPath);

    // See PhysicalOffsetFrame.h
    void extendAddToSystem(SimTK::MultibodySystem& system) const override;

    double getLength(const SimTK::State& s) const override {
        SimTK::Real length = 0.;
        for (const SimTK::CableSpan& cable : m_cableSpans) {
            length += cable.calcLength(s);
        }
        return length;
    }

    double getLengtheningSpeed(const SimTK::State& s) const override {
        SimTK::Real lengthDot = 0.;
        for (const SimTK::CableSpan& cable : m_cableSpans) {
            lengthDot += cable.calcLengthDot(s);
        }
        return lengthDot;
    }

    void addInEquivalentForces(const SimTK::State& state, const double& tension,
            SimTK::Vector_<SimTK::SpatialVec>& bodyForces,
            SimTK::Vector& mobilityForces) const override {
        for (const SimTK::CableSpan& cable : m_cableSpans) {
            cable.applyBodyForces(state, tension, bodyForces);
        }
    }

    double computeMomentArm(
            const SimTK::State& s, const Coordinate& aCoord) const override {
        throw std::runtime_error("Not yet implemented");
    }

    bool isVisualPath() const override {
        throw std::runtime_error("Not yet implemented");
    }

    void constructProperties() {
        constructProperty_accuracy(1e-6);
        constructProperty_smoothness(0.1 / 180. * SimTK::Pi);
        constructProperty_NodeSet(NodeSet());
    }

private:
    std::vector<SimTK::CableSpan> m_cableSpans;
};

} // namespace OpenSim

#endif
