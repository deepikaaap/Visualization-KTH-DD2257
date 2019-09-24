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
    , integrationSteps("integrationSteps", "Number of integrations", 1, 0, 100, 1)
    , backwardDirection("backwardDirection", "Backward Direction")
    , integrateDirectionField("integrateDirectionField", "Integrate Direction Field")
    , maxArcLength("maxArcLength", "Max Arc Length", 1, 0, sqrt(8), 1)
    , maxIntegrationSteps("maxIntegrationSteps", "Max Integration Steps")
    , speedThreshold("speedThreshold", "Speed Threshold")
    , nSeedLines("nSeedLines", "Number of stream lines", 1, 1, 100, 1)
    , uniform("uniform", "Use uniform")
    , xSeed("xSeed", "# x-grids", 1, 0, 50, 1)
    , ySeed("ySeed", "# y-grids", 1, 0, 50, 1)

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
    addProperty(uniform);
    addProperty(nSeedLines);
    addProperty(xSeed);
    addProperty(ySeed);

    // TODO: Register additional properties
    // addProperty(propertyName);

    // Show properties for a single seed and hide properties for multiple seeds
    // (TODO)
    if (propSeedMode.get() == 0) {
        util::hide(nSeedLines, xSeed, ySeed, uniform);
        util::show(propStartPoint, mouseMoveStart);
    }
    uniform.onChange([this]() {
        if (uniform) {
            util::show(xSeed, ySeed);
            util::hide(nSeedLines);
        } else {
            util::hide(xSeed, ySeed);
            util::show(nSeedLines);
        }
    });
    propSeedMode.onChange([this]() {
        if (propSeedMode.get() == 0) {
            util::hide(nSeedLines, uniform, xSeed, ySeed);
            util::show(propStartPoint, mouseMoveStart);

        } else {
            util::hide(propStartPoint, mouseMoveStart);
            uniform.set(false);
            util::show(nSeedLines, uniform);
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
        auto indexBufferLines = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
        // Draw start point
        // Was originally vec2, changing to dvec2 to match the other part of the code.
        dvec2 startPoint = propStartPoint.get();
        vertices.push_back(
            {vec3(startPoint.x, startPoint.y, 0), vec3(0), vec3(0), vec4(0, 0, 0, 1)});
        indexBufferPoints->add(static_cast<std::uint32_t>(0));
        indexBufferLines->add(static_cast<std::uint32_t>(0));
        // TODO: Create one stream line from the given start point
        DrawStreamLine(vectorField, startPoint, indexBufferPoints.get(), indexBufferLines.get(),
                       vertices);
    }

    else if (!uniform) {
        srand(1);
        for (int i = 0; i < nSeedLines; i++) {
            auto indexBufferPoints = mesh->addIndexBuffer(DrawType::Points, ConnectivityType::None);
            auto indexBufferLines = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
            double start_x, start_y;

            double range_x = (vectorField.getBBoxMax()[0] - vectorField.getBBoxMin()[0]);
            double range_y = (vectorField.getBBoxMax()[1] - vectorField.getBBoxMin()[1]);
            start_x = vectorField.getBBoxMin()[0] + (double)rand() / (double)RAND_MAX * range_x;
            start_y = vectorField.getBBoxMin()[1] + (double)rand() / (double)RAND_MAX * range_y;
            dvec2 startPoint = {start_x, start_y};
            indexBufferPoints->add(static_cast<std::uint32_t>(vertices.size()));
            indexBufferLines->add(static_cast<std::uint32_t>(vertices.size()));
            vertices.push_back(
                {vec3(startPoint.x, startPoint.y, 0), vec3(0), vec3(0), vec4(0, 0, 0, 1)});

            DrawStreamLine(vectorField, startPoint, indexBufferPoints.get(), indexBufferLines.get(),
                           vertices);
        }
        // TODO: Seed multiple stream lines either randomly or using a uniform grid
        // (TODO: Bonus, sample randomly according to magnitude of the vector field)
    } else {
        double x_range = vectorField.getBBoxMax()[0] - vectorField.getMinValue()[0];
        double y_range = vectorField.getBBoxMax()[1] - vectorField.getMinValue()[1];

        for (int i = 0; i < xSeed; i++) {
            double start_x = x_range / (xSeed + 1) * (i + 1) + vectorField.getMinValue()[0];
            for (int j = 0; j < ySeed; j++) {
                auto indexBufferPoints =
                    mesh->addIndexBuffer(DrawType::Points, ConnectivityType::None);
                auto indexBufferLines =
                    mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);

                double start_y = y_range / (ySeed + 1) * (j + 1) + vectorField.getMinValue()[1];
                dvec2 startPoint = {start_x, start_y};
                indexBufferPoints->add(static_cast<std::uint32_t>(vertices.size()));
                indexBufferLines->add(static_cast<std::uint32_t>(vertices.size()));
                vertices.push_back(
                    {vec3(startPoint.x, startPoint.y, 0), vec3(0), vec3(0), vec4(0, 0, 0, 1)});
                DrawStreamLine(vectorField, startPoint, indexBufferPoints.get(),
                               indexBufferLines.get(), vertices);
            }
        }
    }

    mesh->addVertices(vertices);
    outMesh.setData(mesh);
}

void StreamlineIntegrator::DrawStreamLine(const VectorField2 &vectorField, const dvec2 &position,
                                          IndexBufferRAM *indexBufferPoints,
                                          IndexBufferRAM *indexBufferLines,
                                          std::vector<BasicMesh::Vertex> &vertices) {
    double arclength = 0.0;
    dvec2 pos = position;
    for (int i = 0; i < integrationSteps; i++) {
        dvec2 new_pos =
            Integrator::RK4(vectorField, pos, stepSize, backwardDirection, integrateDirectionField);
        LogProcessorInfo(i);
        LogProcessorInfo(new_pos);

        arclength += sqrt(pow(new_pos.y - pos.y, 2) + pow(new_pos.x - pos.x, 2));

        if (arclength > maxArcLength) {
            break;
        }
        if (!vectorField.isInside(new_pos)) {
            break;
        }
        if (new_pos == pos) {
            break;
        }
        if (sqrt(pow(new_pos.x - pos.x, 2) + pow(new_pos.y - pos.y, 2)) / double(stepSize) <
            speedThreshold) {
            break;
        }
        indexBufferPoints->add(static_cast<std::uint32_t>(vertices.size()));
        indexBufferLines->add(static_cast<std::uint32_t>(vertices.size()));
        vertices.push_back({vec3(new_pos.x, new_pos.y, 0), vec3(1), vec3(1), vec4(0, 0, 1, 1)});

        pos = new_pos;
    }
}
}  // namespace inviwo
