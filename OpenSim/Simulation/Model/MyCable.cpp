#include "MyCable.h"

#include "OpenSim/Simulation/Model/AbstractGeometryPath.h"
#include "OpenSim/Simulation/Model/GeometryPath.h"
#include "OpenSim/Simulation/Model/MovingPathPoint.h"
#include "OpenSim/Simulation/Model/PathActuator.h"
#include "OpenSim/Simulation/Model/PathPoint.h"
#include "OpenSim/Simulation/Model/PathSpring.h"
#include "OpenSim/Simulation/Model/PhysicalFrame.h"
#include "OpenSim/Simulation/Model/PhysicalOffsetFrame.h"
#include "OpenSim/Simulation/Wrap/PathWrapPoint.h"
#include "OpenSim/Simulation/Wrap/WrapCylinder.h"
#include "OpenSim/Simulation/Wrap/WrapEllipsoid.h"
#include "OpenSim/Simulation/Wrap/WrapObject.h"
#include "OpenSim/Simulation/Wrap/WrapSphere.h"
#include "OpenSim/Simulation/Wrap/WrapTorus.h"
#include <algorithm>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>

#include <SimTKcommon/internal/ExceptionMacros.h>
#include <SimTKcommon/internal/Transform.h>

#include <simbody/internal/CableSpan.h>
#include <simbody/internal/MultibodySystem.h>

using namespace OpenSim;

namespace {

const Frame& FindBaseFrame(const Station& station) {
    return station.getParentFrame().findBaseFrame();
}

// helper: returns the location of the `Station` w.r.t. its base frame
SimTK::Vec3 GetLocationInBaseFrame(const CableAttachmentPoint& point) {
    return point.get_Station().getParentFrame().findTransformInBaseFrame() *
           point.get_Station().get_location();
}

SimTK::MobilizedBodyIndex GetMobilizedBodyIndex(
        const CableAttachmentPoint& point) {
    return point.get_Station().getParentFrame().getMobilizedBodyIndex();
}

SimTK::MobilizedBodyIndex GetMobilizedBodyIndex(const WrapObject& wrapObject) {
    return wrapObject.getFrame().getMobilizedBodyIndex();
}

SimTK::Transform GetTransformInBaseFrame(const WrapObject& wrapObject) {
    return wrapObject.getFrame().findTransformInBaseFrame().compose(
            wrapObject.getTransform());
}

template <typename T> bool printType(const Model& model) {
    std::cout << "LIST TYPE" << T::getClassName() << std::endl;
    bool found = false;
    for (const T& c : model.getComponentList<T>()) {
        std::cout << c.getClassName() << " : " << c.getName() << " @ "
                  << c.getAbsolutePathString() << std::endl;
        found = true;
    }
    return found;
}

template <typename T>
T& findMyComponent(Model& model, const std::string& name) {
    for (const T& c : model.getComponentList<T>()) {
        if (c.getName() == name) {
            return model.updComponent<T>(c.getAbsolutePathString());
        }
    }
    printType<T>(model);
    std::cout << "looking for: " << name << std::endl;
    throw std::runtime_error("component not found");
}

void printPathPointsInfo(const SimTK::State& s, const PathActuator& a) {
    std::cout << "Conversion of AbstractGeometryPath in " << a.getName()
              << std::endl;
    const OpenSim::Array<AbstractPathPoint*>& pts =
            a.getGeometryPath().getCurrentPath(s);
    for (int i = 0; i < pts.getSize(); ++i) {
        std::cout << "p[" << i << "] = ";
        if (dynamic_cast<const PathPoint*>(pts[i])) {
            std::cout << "PathPoint(" << pts[i]->getName() << ")" << std::endl;
            continue;
        }
        if (dynamic_cast<const PathWrapPoint*>(pts[i])) {
            std::cout << "PathWrap(" << pts[i]->getName() << ")" << std::endl;
            continue;
        }
        if (dynamic_cast<const MovingPathPoint*>(pts[i])) {
            throw std::runtime_error("contains moving path point");
        }
        throw std::runtime_error("unknown path point kind");
    }
}

} // namespace

namespace {
struct NodeData {
    std::string name;
    NodeKind kind;

    // AttachmentPoint data
    SimTK::Vec3 station{SimTK::NaN};
    std::string pathToPhysicalFrame;

    // Surface obstacle
    SimTK::Vec3 contactPointHintInS{SimTK::NaN};
    std::string pathToWrapObject;
};

bool isPathPoint(const AbstractPathPoint* pt) {
    return dynamic_cast<const PathPoint*>(pt);
};

bool isPathWrapPoint(const AbstractPathPoint* pt) {
    return dynamic_cast<const PathWrapPoint*>(pt);
};

std::vector<NodeData> readNodeData(
        const GeometryPath& g, const SimTK::State& s) {

    std::vector<NodeData> data;

    const OpenSim::Array<AbstractPathPoint*>& pts = g.getCurrentPath(s);
    const int numPathPoints = g.getPathPointSet().getSize();
    const int numObjects = g.getWrapSet().getSize();

    SimTK_ASSERT1_ALWAYS(pts.getSize() >= 2,
            "Need at least two attachment points, but got %i", pts.getSize());
    SimTK_ASSERT_ALWAYS(isPathPoint(pts[0]),
            "First path point should be an attachment point");
    SimTK_ASSERT_ALWAYS(isPathPoint(pts[pts.getSize() - 1]),
            "Last path point should be an attachment point");

    // Compute the contact point hints.
    std::vector<SimTK::Vec3> contactPointHintInS(
            pts.getSize(), SimTK::Vec3{SimTK::NaN});
    for (int i = 0; i < pts.getSize(); ++i) {
        if (isPathWrapPoint(pts[i]) && isPathWrapPoint(pts[i + 1])) {
            const PathWrapPoint& pathWrap =
                    *dynamic_cast<const PathWrapPoint*>(pts[i]);
            const WrapObject* wrapObject = pathWrap.getWrapObject();

            const SimTK::Transform& X_FS = wrapObject->getTransform();
            const SimTK::Transform& X_GF =
                    wrapObject->getFrame().getTransformInGround(s);
            const SimTK::Transform& X_GS = X_GF.compose(X_FS);
            contactPointHintInS[i] = X_GS.shiftBaseStationToFrame(
                    (pts[i]->getLocationInGround(s) +
                            pts[i + 1]->getLocationInGround(s)) /
                    2.);
            ++i;
        }
    }

    // Keep track of number of attachment points and surfaces.
    int count = 0;
    int countAttachmentPoints = 0;
    for (int i = 0; i < pts.getSize(); ++i) {
        NodeData node;

        // Add an attachment point.
        if (isPathPoint(pts[i])) {
            std::cout << "Add attachment point" << std::endl;

            const PathPoint& pathPoint =
                    *dynamic_cast<const PathPoint*>(pts[i]);

            node.kind = NodeKind::AttachmentPoint;
            node.name = pathPoint.getName();
            node.pathToPhysicalFrame =
                    pathPoint.getParentFrame().getAbsolutePathString();
            node.station = pathPoint.get_location();

            data.push_back(node);
            ++countAttachmentPoints;
            ++count;
            continue;
        }

        // Add an surface obstacle.
        if (isPathWrapPoint(pts[i])) {
            std::cout << "Add surface obstacle" << std::endl;

            const PathWrapPoint& pathWrap =
                    *dynamic_cast<const PathWrapPoint*>(pts[i]);
            const PathWrapPoint& nextPathWrap =
                    *dynamic_cast<const PathWrapPoint*>(pts[i + 1]);

            const WrapObject* wrapObject = pathWrap.getWrapObject();
            const WrapObject* nextWrapObject = nextPathWrap.getWrapObject();

            SimTK_ASSERT2_ALWAYS(
                    wrapObject->getAbsolutePathString() ==
                            nextWrapObject->getAbsolutePathString(),
                    "stop: %s == %s",
                    wrapObject->getAbsolutePathString().c_str(),
                    nextWrapObject->getAbsolutePathString().c_str());

            node.pathToWrapObject = wrapObject->getAbsolutePathString();
            node.contactPointHintInS = contactPointHintInS[i];
            node.kind = NodeKind::Obstacle;
            node.name = wrapObject->getName();

            data.push_back(node);
            ++count;
            ++i;
            continue;
        }

        if (dynamic_cast<const MovingPathPoint*>(pts[i])) {
            std::cout << "MovingPathPoint " << pts[i]->getName() << " @ " << pts[i]->getAbsolutePathString() << std::endl;
            throw std::runtime_error("contains moving path point");
        }

        throw std::runtime_error("unknown path point kind");
    }

    SimTK_ASSERT2_ALWAYS(countAttachmentPoints == numPathPoints,
            "num attachment points (=%i) does not match num path "
            "points "
            "(=%i)",
            countAttachmentPoints, numPathPoints);

    if (count != numPathPoints + numObjects) {
        std::cout << "WARNING: " << g.getName()
                  << " --> Not all wrap obstacles appear to have been "
                     "added!\n";
    }

    return data;
};

void appendNode(
        Model& model, const std::string& cablePath, const NodeData& node) {
    MyCable& cable = model.updComponent<MyCable>(cablePath);

    int numNodesPrev = cable.getNumNodes();

    if (node.kind == NodeKind::AttachmentPoint) {
        std::cout << "node[" << numNodesPrev << "]: Adding attachment point "
                  << node.name << "\n";
        cable.upd_NodeSet().adoptAndAppend(new CableAttachmentPoint());
        CableAttachmentPoint& point = dynamic_cast<CableAttachmentPoint&>(
                cable.upd_NodeSet()[numNodesPrev]);

        point.setName(node.name);
        point.upd_Station().set_location(node.station);
        point.upd_Station().setParentFrame(
                model.getComponent<PhysicalFrame>(node.pathToPhysicalFrame));
    }

    if (node.kind == NodeKind::Obstacle) {
        std::cout << "node[" << numNodesPrev << "]: Adding obstacle "
                  << node.name << "\n";
        cable.upd_NodeSet().adoptAndAppend(new CableObstacle());
        CableObstacle& obstacle =
                dynamic_cast<CableObstacle&>(cable.upd_NodeSet()[numNodesPrev]);

        obstacle.setName(node.name);
        obstacle.set_contactPointHintInS(node.contactPointHintInS);
        obstacle.connectSocket_wrapObject(
                model.getComponent<WrapObject>(node.pathToWrapObject));
    }

    SimTK_ASSERT_ALWAYS(
            cable.getNumNodes() == numNodesPrev + 1, "failed to add node");

    model.finalizeFromProperties();
    model.finalizeConnections();
}

void appendNodes(Model& model, const std::string& cablePath,
        const std::vector<NodeData>& nodes) {
    for (const NodeData& node : nodes) { appendNode(model, cablePath, node); }
}

void replaceGeometryPath(Model& model, const std::string& actuatorPath) {
    std::cout << "Init state" << std::endl;
    SimTK::State& s = model.initSystem();
    model.realizeReport(s);

    std::cout << "Read node data" << std::endl;

    auto isPathActuator = [&]() -> bool {
        return dynamic_cast<const PathActuator*>(
                &model.getComponent(actuatorPath));
    };
    auto isPathSpring = [&]() -> bool {
        return dynamic_cast<const PathSpring*>(
                &model.getComponent(actuatorPath));
    };
    auto getGeometryPath = [&]() -> std::string {
        if (isPathActuator()) {
            return model.getComponent<PathActuator>(actuatorPath)
                    .getPath()
                    .getAbsolutePathString();
        }
        if (isPathSpring()) {
            return model.getComponent<PathSpring>(actuatorPath)
                    .getPath()
                    .getAbsolutePathString();
        }
        throw std::runtime_error("unknown actuator");
    };

    const std::vector<NodeData> nodes = readNodeData(
            model.getComponent<GeometryPath>(getGeometryPath()), s);

    SimTK_ASSERT_ALWAYS(!nodes.empty(), "failed to read nodes");

    for (const NodeData& node : nodes) {
        if (node.kind == NodeKind::AttachmentPoint) {
            std::cout << "AttachmentPoint " << node.name << ", "
                      << node.pathToPhysicalFrame << "\n";
        }
        if (node.kind == NodeKind::Obstacle) {
            std::cout << "Obstacle " << node.name << ", "
                      << node.pathToWrapObject << "\n";
        }
    }

    std::cout << "Replace geometry path" << std::endl;
    if (isPathActuator()) {
        model.updComponent<PathActuator>(actuatorPath)
                .updProperty_path()
                .setValueAsObject(MyCable());
    }
    if (isPathSpring()) {
        model.updComponent<PathSpring>(actuatorPath)
                .updProperty_path()
                .setValueAsObject(MyCable());
    }

    model.finalizeFromProperties();
    model.finalizeConnections();

    std::cout << "Append nodes to geometry path" << std::endl;
    appendNodes(model, getGeometryPath(), nodes);
}

void replaceAllGeometryPath(Model& model) {
    model.finalizeFromProperties();
    model.finalizeConnections();

    auto ContainsMovingPathPoints = [&](const GeometryPath& g) ->bool {
        SimTK::State& s = model.initSystem();
        model.realizeReport(s);
        const OpenSim::Array<AbstractPathPoint*>& pts = g.getCurrentPath(s);
        for (int i = 0; i < pts.getSize(); ++i) {
            if (dynamic_cast<const MovingPathPoint*>(pts.get(i))) {
                return true;
            }
        }
        return false;
    };

    auto replaceFirst1 = [&]() -> bool {
        for (const PathActuator& c : model.getComponentList<PathActuator>()) {
            if (dynamic_cast<const MyCable*>(&c.getPath())) { continue; }
            if (!dynamic_cast<const GeometryPath*>(&c.getPath())) {
                std::cout << c.getName()
                          << " does not have a GeometryPath but a "
                          << c.getClassName() << "\n";
                continue;
            }
            if (ContainsMovingPathPoints(c.getGeometryPath())) {
                std::cout << "WARNING: Detected MovingPathPoint, skip conversion of " << c.getName() << std::endl;
                continue;
            }
            // Make sure the path is computed.

            replaceGeometryPath(model, c.getAbsolutePathString());
            return true;
        }
        return false;
    };

    auto replaceFirst2 = [&]() -> bool {
        for (const PathSpring& c : model.getComponentList<PathSpring>()) {
            if (dynamic_cast<const MyCable*>(&c.getPath())) { continue; }
            if (!dynamic_cast<const GeometryPath*>(&c.getPath())) {
                std::cout << c.getName()
                          << " does not have a GeometryPath but a "
                          << c.getClassName() << "\n";
                continue;
            }
            if (ContainsMovingPathPoints(c.getGeometryPath())) {
                std::cout << "WARNING: Detected MovingPathPoint, skip conversion of " << c.getName() << std::endl;
                continue;
            }
            // Make sure the path is computed.

            replaceGeometryPath(model, c.getAbsolutePathString());
            return true;
        }
        return false;
    };

    while (replaceFirst1() || replaceFirst2()) {}
}

} // namespace

void MyCable::extendAddToSystem(SimTK::MultibodySystem& system) const {
    Super::extendAddToSystem(system);

    std::vector<SimTK::CableSpan>& cable =
            const_cast<std::vector<SimTK::CableSpan>&>(m_cableSpans);
    cable.clear();

    SimTK::CableSubsystem& subsystem = _model->updCableSubsystem();

    const SimTK::Vec3 invalidStation{SimTK::NaN};
    const SimTK::MobilizedBodyIndex invalidBody =
            SimTK::MobilizedBodyIndex::Invalid();

    const NodeSet& nodes = get_NodeSet();
    SimTK_ASSERT_ALWAYS(nodes.getSize() >= 2, "must have atleast two nodes");
    SimTK_ASSERT_ALWAYS(nodes.get(0).getKind() == NodeKind::AttachmentPoint,
            "first node must be attachmentPoint");
    SimTK_ASSERT_ALWAYS(nodes.get(nodes.getSize() - 1).getKind() ==
                                NodeKind::AttachmentPoint,
            "last node must be attachmentPoint");

    auto asPoint = [&](int i) -> const CableAttachmentPoint& {
        return dynamic_cast<const CableAttachmentPoint&>(nodes.get(i));
    };
    auto asObstacle = [&](int i) -> const CableObstacle& {
        return dynamic_cast<const CableObstacle&>(nodes.get(i));
    };

    // Process the first node.
    {
        cable.emplace_back(subsystem, GetMobilizedBodyIndex(asPoint(0)),
                GetLocationInBaseFrame(asPoint(0)), invalidBody,
                invalidStation);
    }

    for (int i = 1; i < nodes.getSize(); ++i) {
        // Process an attachment point.
        if (nodes.get(i).getKind() == NodeKind::AttachmentPoint) {

            // Finalize the cable span.
            cable.back().setTerminationBodyIndex(
                    GetMobilizedBodyIndex(asPoint(i)));
            cable.back().setTerminationStation(
                    GetLocationInBaseFrame(asPoint(i)));

            // Stop at last node.
            if (i == nodes.getSize() - 1) { break; }

            // Start a new cable span.
            cable.emplace_back(subsystem, GetMobilizedBodyIndex(asPoint(i)),
                    GetLocationInBaseFrame(asPoint(i)), invalidBody,
                    invalidStation);
        }

        // Process an obstacle.
        if (nodes.get(i).getKind() == NodeKind::Obstacle) {
            // Add the obstacle.
            const WrapObject& wrapObject = asObstacle(i).getWrapObject();
            cable.back().addObstacle(GetMobilizedBodyIndex(wrapObject),
                    GetTransformInBaseFrame(wrapObject),
                    wrapObject.getContactGeometry(),
                    asObstacle(i).get_contactPointHintInS());
        }
    }

    for (SimTK::CableSpan& cableSpan : cable) {
        cableSpan.setCurveSegmentAccuracy(get_accuracy());
        cableSpan.setSmoothnessTolerance(get_smoothness());
    }
}

Model MyCable::convert(const std::string& modelPath) {
    Model model(modelPath);

    /* { */
    /*     SimTK::State& s = model.initSystem(); */
    /*     model.realizeReport(s); */
    /*     for (const GeometryPath& c: model.getComponentList<GeometryPath>()) { */
    /*         const OpenSim::Array<AbstractPathPoint*>& pts = c.getCurrentPath(s); */
    /*         for (int i = 0; i < pts.getSize(); ++i) { */
    /*         const AbstractPathPoint* p = pts.get(i); */
    /*             std::cout << "p[" << i << "] : " << p->getName() << " @ " << p->getAbsolutePathString() << p->getClassName() << "\n"; */
    /*         } */
    /*     } */
    /* } */

    /* const bool foundMoving = printType<MovingPathPoint>(model); */
    /* SimTK_ASSERT_ALWAYS(!foundMoving, "Found moving path points!"); */

    replaceAllGeometryPath(model);

    printType<CableAttachmentPoint>(model);
    printType<CableObstacle>(model);
    printType<MyCable>(model);

    model.finalizeFromProperties();
    model.finalizeConnections();
    return model;
}
