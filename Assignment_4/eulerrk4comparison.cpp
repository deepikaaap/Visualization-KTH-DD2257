/*********************************************************************
 *  Author  : Himangshu Saikia, Wiebke Koepp
 *  Init    : Tuesday, September 19, 2017 - 15:08:24
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <inviwo/core/datastructures/geometry/basicmesh.h>
#include <inviwo/core/interaction/events/mouseevent.h>
#include <labstreamlines/eulerrk4comparison.h>
#include <labstreamlines/integrator.h>

namespace inviwo {

// The Class Identifier has to be globally unique. Use a reverse DNS naming scheme
const ProcessorInfo EulerRK4Comparison::processorInfo_{
    "org.inviwo.EulerRK4Comparison",  // Class identifier
    "Euler RK4 Comparison",           // Display name
    "KTH Lab",                        // Category
    CodeState::Experimental,          // Code state
    Tags::None,                       // Tags
};

const ProcessorInfo EulerRK4Comparison::getProcessorInfo() const { return processorInfo_; }

EulerRK4Comparison::EulerRK4Comparison()
    : Processor()
    , inData("inData")
    , outMesh("meshOut")
    , propStartPoint("startPoint", "Start Point", vec2(0.5, 0.5), vec2(0), vec2(1024), vec2(0.5))
    , mouseMoveStart("mouseMoveStart", "Move Start", [this](Event* e) { eventMoveStart(e); },
                     MouseButton::Left, MouseState::Press | MouseState::Move)
    , stepSize("stepSize", "Step size", 1, 0, 2, 1)
    , integrationSteps("integrationSteps", "Number of integrations", 1, 1, 100, 1)
// TODO: Initialize additional properties
// propertyName("propertyIdentifier", "Display Name of the Propery",
// default value (optional), minimum value (optional), maximum value (optional), increment
// (optional)); propertyIdentifier cannot have spaces
{
    // Register Ports
    addPort(outMesh);
    addPort(inData);

    // Register Properties
    addProperty(propStartPoint);
    addProperty(mouseMoveStart);

    // TODO: Register additional properties
    // addProperty(propertyName);
    addProperty(stepSize);
    addProperty(integrationSteps);
}

void EulerRK4Comparison::eventMoveStart(Event* event) {
    if (!inData.hasData()) return;
    auto mouseEvent = static_cast<MouseEvent*>(event);
    vec2 mousePos = mouseEvent->posNormalized();

    // Map to range [0,1]^2
    mousePos = mousePos * 2 - vec2(1, 1);

    // Update starting point
    propStartPoint.set(mousePos);
    event->markAsUsed();
}

void EulerRK4Comparison::process() {
    // Get input
    if (!inData.hasData()) {
        return;
    }
    auto vol = inData.getData();

    // Retreive data in a form that we can access it
    const VectorField2 vectorField = VectorField2::createFieldFromVolume(vol);
    // The start point should be inside the volume (set maximum to the upper right corner)
    auto bboxMin = vectorField.getBBoxMin();
    propStartPoint.setMinValue(bboxMin);
    propStartPoint.setMaxValue(vectorField.getBBoxMax());

    // Initialize mesh, vertices and index buffers for the two streamlines and the points
    auto mesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> vertices;
	std::vector<BasicMesh::Vertex> vertices_rk4;
    auto indexBufferEuler = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
    auto indexBufferRK = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
    auto indexBufferPoints = mesh->addIndexBuffer(DrawType::Points, ConnectivityType::None);

    // Draw start point
    dvec2 startPoint = propStartPoint.get();
    vertices.push_back({vec3(startPoint.x, startPoint.y, 0), vec3(1), vec3(1), vec4(0, 0, 0, 1)});
    indexBufferPoints->add(static_cast<std::uint32_t>(0));
    indexBufferEuler->add(static_cast<std::uint32_t>(0));
    indexBufferRK->add(static_cast<std::uint32_t>(0));

    // TODO: Implement the Euler and Runge-Kutta of 4th order integration schemes
    // and then integrate forward for a specified number of integration steps and a given stepsize
    // (these should be additional properties of the processor)
    //Integrator integrator = Integrator();

    dvec2 startPoint_rk4 = propStartPoint.get();
    vertices_rk4.push_back({vec3(startPoint_rk4.x, startPoint_rk4.y, 0), vec3(1), vec3(1), vec4(1, 1, 1, 0)});
    for (int i = 0; i < integrationSteps; i++) {

        dvec2 new_pos_euler = Integrator::Euler(vectorField, startPoint, stepSize, 0, 0);
        
        indexBufferEuler->add(static_cast<std::uint32_t>(i));
        indexBufferPoints->add(static_cast<std::uint32_t>(i));
        vertices.push_back({vec3(new_pos_euler.x, new_pos_euler.y, 0), vec3(1), vec3(1), vec4(0, 0, 1, 1)});
        startPoint = new_pos_euler;
    }
    indexBufferPoints->add(static_cast<std::uint32_t>(integrationSteps));
    vertices.push_back({vec3(startPoint.x, startPoint.y, 0), vec3(1), vec3(1), vec4(0, 0, 0, 0)});
    mesh->addVertices(vertices);

    for (int i = 0; i < integrationSteps; i++) {

        dvec2 new_pos_rk4 = Integrator::RK4(vectorField, startPoint_rk4, stepSize, 0, 0);
        indexBufferRK->add(static_cast<std::uint32_t>(i+integrationSteps+1));
        indexBufferPoints->add(static_cast<std::uint32_t>(i+integrationSteps));
        vertices_rk4.push_back(
            {vec3(new_pos_rk4.x, new_pos_rk4.y, 0), vec3(1), vec3(1), vec4(0, 0, 0, 1)});
        startPoint_rk4 = new_pos_rk4;
    }

	
    mesh->addVertices(vertices_rk4);
    outMesh.setData(mesh);
}

}  // namespace inviwo
