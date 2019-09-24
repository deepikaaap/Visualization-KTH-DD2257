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
    , maxArcLength("maxArcLength", "Max Arc Length",1,0,sqrt(8),1)
    , speedThreshold("speedThreshold", "Speed Threshold")
    , samplingType("samplingType", "Sampling")
    , nSeedLines("nSeedLines", "Number of stream lines", 1, 1, 1000, 1)
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
        samplingType.addOption("random", "Random sampling", 0);
        samplingType.addOption("uniform", "Uniform sampling", 1);
        samplingType.addOption("magnitude", "Magnitude based sampling", 2);
        addProperty(samplingType);
        addProperty(nSeedLines);
        addProperty(xSeed);
        addProperty(ySeed); 
        
        // TODO: Register additional properties
        // addProperty(propertyName);
        
        // Show properties for a single seed and hide properties for multiple seeds
        // (TODO)
        if (propSeedMode.get() == 0) {
            util::hide(nSeedLines, xSeed,ySeed, samplingType);
            util::show(propStartPoint, mouseMoveStart);
        }
        samplingType.onChange([this](){
            if(samplingType.get()==1){
                util::show(xSeed,ySeed);
                util::hide(nSeedLines);
            }
            else{
                util::hide(xSeed,ySeed);
                util::show(nSeedLines);
                
            }
            
        });
        propSeedMode.onChange([this]() {
            if (propSeedMode.get() == 0) {
                util::hide(nSeedLines, samplingType, xSeed, ySeed);
                util::show(propStartPoint, mouseMoveStart);
                
            } else {
                util::hide(propStartPoint, mouseMoveStart);
                util::show(nSeedLines, samplingType);
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
            vertices.push_back({vec3(startPoint.x, startPoint.y, 0), vec3(0), vec3(0), vec4(0, 0, 0, 1)});
            indexBufferPoints->add(static_cast<std::uint32_t>(0));
            indexBufferLines->add(static_cast<std::uint32_t>(0));
            // TODO: Create one stream line from the given start point
            DrawStreamLine(vectorField, startPoint, indexBufferPoints.get(), indexBufferLines.get(),vertices);
        }
        
        else if (samplingType.get()==0){
            srand(1);
            for(int i = 0; i <nSeedLines; i++) {
                auto indexBufferPoints = mesh->addIndexBuffer(DrawType::Points, ConnectivityType::None);
                auto indexBufferLines = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
                dvec2 startPoint = generateRandomSample(vectorField);
                indexBufferPoints->add(static_cast<std::uint32_t>(vertices.size()));
                indexBufferLines->add(static_cast<std::uint32_t>(vertices.size()));
                vertices.push_back({vec3(startPoint.x, startPoint.y, 0), vec3(0), vec3(0), vec4(0, 0, 0, 1)});
                
                DrawStreamLine(vectorField, startPoint, indexBufferPoints.get(),indexBufferLines.get(), vertices);
            }
            // TODO: Seed multiple stream lines either randomly or using a uniform grid
            // (TODO: Bonus, sample randomly according to magnitude of the vector field)
        }
        else if(samplingType.get()==1) {
            double x_range = vectorField.getBBoxMax()[0]-vectorField.getMinValue()[0];
            double y_range = vectorField.getBBoxMax()[1]-vectorField.getMinValue()[1];
            
            for(int i=0; i< xSeed; i++){
                double start_x = x_range / (xSeed+1) * (i+1) + vectorField.getMinValue()[0];
                for(int j=0; j< ySeed; j++){
                    auto indexBufferPoints = mesh->addIndexBuffer(DrawType::Points, ConnectivityType::None);
                    auto indexBufferLines = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
                    
                    double start_y = y_range / (ySeed+1) * (j+1) + vectorField.getMinValue()[1];
                    dvec2 startPoint = {start_x, start_y};
                    indexBufferPoints->add(static_cast<std::uint32_t>(vertices.size()));
                    indexBufferLines->add(static_cast<std::uint32_t>(vertices.size()));
                    vertices.push_back({vec3(startPoint.x, startPoint.y, 0), vec3(0), vec3(0), vec4(0, 0, 0, 1)});
                    DrawStreamLine(vectorField, startPoint, indexBufferPoints.get(),indexBufferLines.get(), vertices);
                }
            }
        }
        else {
            double totalMagnitude = getTotalMagnitude(vectorField);
            int counter = 0;
            while(counter < nSeedLines) {
                dvec2 samplePoint = generateRandomSample(vectorField);
                dvec2 value_samplePoint = vectorField.interpolate(samplePoint);
                double normalizedMagnitude = sqrt(pow(value_samplePoint[0],2) + pow(value_samplePoint[1],2)) / totalMagnitude;
                double randNum = (double)rand()/(double) RAND_MAX;
                auto indexBufferPoints = mesh->addIndexBuffer(DrawType::Points, ConnectivityType::None);
                auto indexBufferLines = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
                
                if(randNum < normalizedMagnitude) {
                    indexBufferPoints->add(static_cast<std::uint32_t>(vertices.size()));
                    indexBufferLines->add(static_cast<std::uint32_t>(vertices.size()));
                    vertices.push_back({vec3(samplePoint.x, samplePoint.y, 0), vec3(0), vec3(0), vec4(0, 0, 0, 1)});
                    
                    DrawStreamLine(vectorField, samplePoint, indexBufferPoints.get(),indexBufferLines.get(), vertices);
                    counter++;
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
            dvec2 new_pos = Integrator::RK4(vectorField, pos, stepSize,
                                            backwardDirection, integrateDirectionField);
            
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
            if (sqrt(pow(new_pos.x - pos.x,2) + pow(new_pos.y-pos.y,2)) / double(stepSize) < speedThreshold) {
                break;
            }
            indexBufferPoints->add(static_cast<std::uint32_t>(vertices.size()));
            indexBufferLines->add(static_cast<std::uint32_t>(vertices.size()));
            vertices.push_back({vec3(new_pos.x, new_pos.y, 0), vec3(1), vec3(1), vec4(0, 0, 1, 1)});
            
            pos = new_pos;
        }
    }
    
    double StreamlineIntegrator::getTotalMagnitude(const VectorField2 &vectorField){
        double maxMagnitude = 0.0;
        double magnitude;
        dvec2 maxPos;
        for(int i = 0; i < vectorField.getNumVerticesPerDim()[0]; i++) {
            for (int j = 0; j < vectorField.getNumVerticesPerDim()[1]; j++) {
                dvec2 pos =vectorField.getValueAtVertex({i,j});
                magnitude = sqrt(pow(pos[0],2)+ pow(pos[1],2));
                if(magnitude > maxMagnitude) {
                    maxMagnitude = magnitude;
                    maxPos = {i,j};
                }
            }
        }
        return maxMagnitude;
    }
    
    dvec2 StreamlineIntegrator::generateRandomSample(const VectorField2 &vectorField){
        double start_x, start_y;
        double range_x = (vectorField.getBBoxMax()[0] - vectorField.getBBoxMin()[0]);
        double range_y = (vectorField.getBBoxMax()[1] - vectorField.getBBoxMin()[1]);
        start_x = vectorField.getBBoxMin()[0]+ (double)rand()/(double) RAND_MAX * range_x;
        start_y = vectorField.getBBoxMin()[1]+ (double)rand()/(double) RAND_MAX * range_y;
        return {start_x, start_y};
    }
    
}  // namespace inviwo
