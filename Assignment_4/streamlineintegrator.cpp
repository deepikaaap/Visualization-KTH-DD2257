/*********************************************************************
 *  Author  : Himangshu Saikia
 *  Init    : Tuesday, September 19, 2017 - 15:08:33
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <inviwo/core/interaction/events/mouseevent.h>
#include <inviwo/core/util/utilities.h>
#include <labstreamlines/integrator.h>
#include <labstreamlines/streamlineintegrator.h>
#include <labutils/scalarvectorfield.h>

namespace inviwo {

// The Class Identifier has to be globally unique. Use a reverse DNS naming scheme
const ProcessorInfo StreamlineIntegrator::processorInfo_{
    "org.inviwo.StreamlineIntegrator",  // Class identifier
    "Streamline Integrator",            // Display name
    "KTH Lab",                          // Category
    CodeState::Experimental,            // Code state
    Tags::None,                         // Tags
};

const ProcessorInfo StreamlineIntegrator::getProcessorInfo() const { return processorInfo_; }

StreamlineIntegrator::StreamlineIntegrator()
    : Processor()
    , inData("volIn")
    , outMesh("meshOut")
    , propStartPoint("startPoint", "Start Point", vec2(0.5, 0.5), vec2(-1), vec2(1), vec2(0.1))
    , propSeedMode("seedMode", "Seeds")
    , mouseMoveStart("mouseMoveStart", "Move Start", [this](Event *e) { eventMoveStart(e); },
                     MouseButton::Left, MouseState::Press | MouseState::Move)
    , stepSize("stepSize", "Step size", 1, 0, 2, 1)
    , integrationSteps("integrationSteps", "Number of integrations", 1, 1, 100, 1)
    , backwardDirection("backwardDirection", "Backward Direction")
    , integrateDirectionField("integrateDirectionField", "Integrate Direction Field")
    , maxArcLength("maxArcLength", "Max Arc Length", 1, 0, sqrt(8), 1)
    , maxIntegrationSteps("maxIntegrationSteps", "Max Integration Steps")
    , speedThreshold("speedThreshold", "Speed Threshold")
// TODO: Initialize additional properties
// propertyName("propertyIdentifier", "Display Name of the Propery",
// default value (optional), minimum value (optional), maximum value (optional),
// increment (optional)); propertyIdentifier cannot have spaces
{
    // Register Ports
    addPort(inData);
    addPort(outMesh);

    // Register Properties
    propSeedMode.addOption("one", "Single Start Point", 0);
    propSeedMode.addOption("multiple", "Multiple Seeds", 1);
    addProperty(propSeedMode);
    addProperty(propStartPoint);
    addProperty(mouseMoveStart);
    addProperty(stepSize);
    addProperty(integrationSteps);
    addProperty(backwardDirection);
    addProperty(integrateDirectionField);
    addProperty(maxArcLength);
    addProperty(speedThreshold);
    addProperty(maxIntegrationSteps);

    // TODO: Register additional properties
    // addProperty(propertyName);

    // Show properties for a single seed and hide properties for multiple seeds
    // (TODO)
    propSeedMode.onChange([this]() {
        if (propSeedMode.get() == 0) {
            util::show(propStartPoint, mouseMoveStart);
            // util::hide(...)
        } else {
            util::hide(propStartPoint, mouseMoveStart);
            // util::show(...)
        }
    });
}

void StreamlineIntegrator::eventMoveStart(Event *event) {
    if (!inData.hasData()) return;
    auto mouseEvent = static_cast<MouseEvent *>(event);
    vec2 mousePos = mouseEvent->posNormalized();

    // Map to range [0,1]^2
    mousePos = mousePos * 2 - vec2(1, 1);

    // Update starting point
    propStartPoint.set(mousePos);
    event->markAsUsed();
}

void StreamlineIntegrator::process() {
    // Get input
    if (!inData.hasData()) {
        return;
    }
    auto vol = inData.getData();

    // Retreive data in a form that we can access it
    auto vectorField = VectorField2::createFieldFromVolume(vol);

    // The start point should be inside the volume (set maximum to the upper
    // right corner)
    auto mesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> vertices;

    if (propSeedMode.get() == 0) {
        auto indexBufferPoints = mesh->addIndexBuffer(DrawType::Points, ConnectivityType::None);
        // Draw start point
        // Was originally vec2, changing to dvec2 to match the other part of the code.
        dvec2 startPoint = propStartPoint.get();
        vertices.push_back(
            {vec3(startPoint.x, startPoint.y, 0), vec3(0), vec3(0), vec4(0, 0, 0, 1)});
        indexBufferPoints->add(static_cast<std::uint32_t>(0));
        // TODO: Create one stream line from the given start point
        DrawStreamLine(vectorField, startPoint, indexBufferPoints.get(), vertices);
    } else {
        // TODO: Seed multiple stream lines either randomly or using a uniform grid
        // (TODO: Bonus, sample randomly according to magnitude of the vector field)
    }

    mesh->addVertices(vertices);
    outMesh.setData(mesh);
}

void StreamlineIntegrator::DrawStreamLine(const VectorField2 &vectorField, dvec2 &position,
                                          IndexBufferRAM *indexBufferPoints,
                                          std::vector<BasicMesh::Vertex> &vertices) {
    // Starting index may have to be changed while drawing multiple lines
    for (int i = 0; i < integrationSteps; i++) {
        dvec2 new_pos = Integrator::RK4(vectorField, position, stepSize, backwardDirection,
                                        integrateDirectionField);
		// length of the arc from the start point to the current position
        double arclength = sqrt(pow(new_pos.y - propStartPoint.get().y, 2) +
                                pow(new_pos.x - propStartPoint.get().x, 2));
		// check for arc length being within the limits
        if (arclength > maxArcLength) {
            break;
        }
		// check for the new position being within the boundaries of the vector field
        if ((vectorField.getBBoxMin()[0] < new_pos.x < vectorField.getBBoxMax()[0]) ||
            (vectorField.getBBoxMin()[1] < new_pos.x < vectorField.getBBoxMax()[1])) {
            break;
        }
		// check to stop integration for zero vector field
        if (new_pos == position) {
            break;
		}
		// check to stop integration for slow velocity(being below the speedThreshold)
        if (double(stepSize) * (sqrt(pow(new_pos.x - position.x, 2) + pow(new_pos.y - position.y, 2)) < speedThreshold)) {
            break;
		}
       
		indexBufferPoints->add(static_cast<std::uint32_t>(i));
        vertices.push_back({vec3(new_pos.x, new_pos.y, 0), vec3(1), vec3(1), vec4(0, 0, 1, 1)});
        position = new_pos;
    }
}
}  // namespace inviwo
