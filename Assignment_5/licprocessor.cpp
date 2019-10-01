/*********************************************************************
 *  Author  : Himangshu Saikia
 *  Init    : Monday, October 02, 2017 - 13:31:17
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <inviwo/core/datastructures/volume/volumeram.h>
#include <lablic/licprocessor.h>
#include <labstreamlines/integrator.h>
#include "cmath"

namespace inviwo {

// The Class Identifier has to be globally unique. Use a reverse DNS naming scheme
const ProcessorInfo LICProcessor::processorInfo_{
    "org.inviwo.LICProcessor",  // Class identifier
    "LICProcessor",             // Display name
    "KTH Labs",                 // Category
    CodeState::Experimental,    // Code state
    Tags::None,                 // Tags
};

const ProcessorInfo LICProcessor::getProcessorInfo() const { return processorInfo_; }

LICProcessor::LICProcessor()
    : Processor()
    , volumeIn_("volIn")
    , noiseTexIn_("noiseTexIn")
    , licOut_("licOut")
    , kernelSize("kernelSize", "Kernel Size", 5, 1, 200,1)
    , fastLIC("fastLIC", "Use FastLIC")
// TODO: Register additional properties
{
    // Register ports
    addPort(volumeIn_);
    addPort(noiseTexIn_);
    addPort(licOut_);

    // Register properties
    // TODO: Register additional properties
    addProperty(kernelSize);
    addProperty(fastLIC);
    fastLIC.set(true);
}

void LICProcessor::process() {
    // Get input
    if (!volumeIn_.hasData()) {
        return;
    }

    if (!noiseTexIn_.hasData()) {
        return;
    }
    
    auto vol = volumeIn_.getData();
    const VectorField2 vectorField = VectorField2::createFieldFromVolume(vol);
    vectorFieldDims_ = vol->getDimensions();

    auto tex = noiseTexIn_.getData();
    const RGBAImage texture = RGBAImage::createFromImage(tex);
    texDims_ = tex->getDimensions();

    // Prepare the output, it has the same dimensions as the texture and rgba values in [0,255]
    auto outImage = std::make_shared<Image>(texDims_, DataVec4UInt8::get());
    RGBAImage licImage(outImage);

    std::vector<std::vector<double>> licTexture(texDims_.x, std::vector<double>(texDims_.y));
    std::vector<std::vector<bool>> visited(texDims_.x, std::vector<bool>(texDims_.y, false));

    // TODO: Implement LIC and FastLIC
    // This code instead just creates a black image

    std::vector<dvec2> listForward;
    std::vector<dvec2> listBackward;
    double vf_x;
    double vf_y;
    double pxWidth = (double)(vectorField.getBBoxMax().x-vectorField.getBBoxMin().x) / texDims_.x, pxHeight = (double) (vectorField.getBBoxMax().y-vectorField.getBBoxMin().y) / texDims_.y;
    double step = pxWidth < pxHeight ? pxWidth: pxHeight;
    
    std::vector<std::vector<int>> interVarForward, interVarBackward;
    int pixelX, pixelY;
    for (auto i = 0; i < (int) texDims_.x; i++) {
        for (auto j = 0; j < (int) texDims_.y; j++) {
            //int val = int(licTexture[i][j]);
            // licImage.setPixel(size2_t(i, j), dvec4(val, val, val, 255));
            // or
            //licImage.setPixelGrayScale(size2_t(i, j), val);
            //double vf_x = i * (double)vectorFieldDims_.x/texDims_.x;
            //double vf_y = j * (double)vectorFieldDims_.y/texDims_.y;
            
            vf_x = (((double)i) / (texDims_.x-1)) * (vectorField.getBBoxMax().x - vectorField.getBBoxMin().x) + vectorField.getBBoxMin().x;
            vf_y = (((double)j) / (texDims_.y-1)) * (vectorField.getBBoxMax().y - vectorField.getBBoxMin().y) + vectorField.getBBoxMin().y;
            if(fastLIC == false) {
                listForward = Integrator::StreamLines(vectorField, {vf_x, vf_y}, step, 0, true, 0.0, kernelSize, false);
                listBackward = Integrator::StreamLines(vectorField, {vf_x, vf_y}, step, 1, true, 0.0, kernelSize, false);
                licTexture[i][j] = convolution(vectorField, listForward, listBackward, texture);
            }
            else {
                if(visited[i][j] == false) {
                    listForward = Integrator::StreamLines(vectorField, {vf_x, vf_y}, step, 0, true, 0.0, 2*kernelSize+1, true);
                    listBackward = Integrator::StreamLines(vectorField, {vf_x, vf_y}, step, 1, true, 0.0, 2*kernelSize+1, true);
                
                    /*
                    std::reverse(listBackward.begin() + 1, listBackward.end());
                    listBackward.insert(listBackward.end(), listForward.begin(), listForward.end());
                    listBackward = mapVectorToTextureField(listBackward,texDims_,vectorField);
                    */
                    
                    //LogProcessorInfo("Before mapVector");
                    
                    
                    interVarForward = mapVectorToTextureField(listForward,texDims_,vectorField);
                    interVarBackward = mapVectorToTextureField(listBackward,texDims_,vectorField);
                    
                    for(int k=0; k < kernelSize && k < (int)interVarForward.size(); k++){
                        pixelX = interVarForward[k][0];
                        pixelY = interVarForward[k][1];
                        if (visited[pixelX][pixelY]) {
                            continue;
                        }
                        licTexture[pixelX][pixelY] = LinearIntegralConvolution(texture, interVarForward, interVarBackward, k);
                        
                        visited[pixelX][pixelY] = true;
                    }
                    for(int k=1; k < kernelSize && k < (int)interVarBackward.size(); k++){
                        pixelX = interVarBackward[k][0];
                        pixelY = interVarBackward[k][1];
                        if (visited[pixelX][pixelY]) {
                            continue;
                        }
                        licTexture[pixelX][pixelY] = LinearIntegralConvolution(texture, interVarForward, interVarBackward, k);
                        visited[pixelX][pixelY] = true;
                    }
                    
                }
                
            }
                     
            int val = int(licTexture[i][j]);
            licImage.setPixel(size2_t(i, j), dvec4(val, val, val, 255));
        }
    }
                     

    licOut_.setData(outImage);
}
    
    std::vector<std::vector<int>> LICProcessor::mapVectorToTextureField(const std::vector<dvec2>& vectorList, const dvec2 &dims,const VectorField2 &vectorField ){
        std::vector<std::vector<int>> pixelList((int)vectorList.size(), std::vector<int>(2));
        double mapX, mapY;
        for(int i=0; i < (int)vectorList.size();i++){
            mapX = (vectorList[i].x - vectorField.getBBoxMin().x) * (double) (texDims_.x - 1) / (vectorField.getBBoxMax().x - vectorField.getBBoxMin().x);
            mapY = (vectorList[i].y - vectorField.getBBoxMin().y) * (double) (texDims_.y - 1) / (vectorField.getBBoxMax().y - vectorField.getBBoxMin().y);
            
            pixelList[i][0] = (int)round(mapX);
            pixelList[i][1] = (int)round(mapY);
             
        }
        return pixelList;
    }
   
    double LICProcessor::LinearIntegralConvolution(const RGBAImage& inputImage, const std::vector<std::vector<int>>& forwardList, const std::vector<std::vector<int>>& backwardList,const int& filterStartIndex){
        double conv=0.0;
        double kernelValue;
        for (int i = filterStartIndex; i < filterStartIndex+kernelSize && i<(int)forwardList.size(); i++) {
            
            kernelValue = boxKernel(i - filterStartIndex);
            conv += kernelValue * inputImage.readPixelGrayScale({(int) forwardList[i][0], (int) forwardList[i][1]});
        }
        for (int i = filterStartIndex; i < filterStartIndex+kernelSize && i<(int)backwardList.size(); i++) {
            kernelValue = boxKernel(i - filterStartIndex);
            conv += kernelValue * inputImage.readPixelGrayScale({(int) backwardList[i][0], (int) backwardList[i][1]});
        }
        return conv;
    }
    
    
    double LICProcessor::convolution(const VectorField2& vectorField, const std::vector<dvec2>& forwardList, const std::vector<dvec2>& backwardList, const RGBAImage& inputImage) {
        double conv = 0.0;
        double mapX;
        double mapY;
        for (int i = 0; i < (int) forwardList.size(); i++) {
            double kernelValue = boxKernel(i);
            mapX = (forwardList[i].x - vectorField.getBBoxMin().x) * (double) (texDims_.x - 1) / (vectorField.getBBoxMax().x - vectorField.getBBoxMin().x);
            mapY = (forwardList[i].y - vectorField.getBBoxMin().y) * (double) (texDims_.y - 1) / (vectorField.getBBoxMax().y - vectorField.getBBoxMin().y);
            conv += kernelValue * inputImage.readPixelGrayScale({(int) mapX, (int) mapY});
        }
        
        for (int i = 1; i < (int) backwardList.size(); i++) {
            double kernelValue = boxKernel(i);
            mapX = (backwardList[i].x - vectorField.getBBoxMin().x) * (double) (texDims_.x -1) / (vectorField.getBBoxMax().x - vectorField.getBBoxMin().x);
            mapY = (backwardList[i].y - vectorField.getBBoxMin().y) * (double) (texDims_.y -1 ) / (vectorField.getBBoxMax().y - vectorField.getBBoxMin().y);
            conv += kernelValue * inputImage.readPixelGrayScale({(int) mapX, (int) mapY});
        }
        
        return conv;
    }
    
    
    double LICProcessor::boxKernel(const int index) {
        if (index > kernelSize)
            return 0.0;
        else
            return 1.0 / (2*kernelSize);   //Maybe 2*kernelSize + 1.
    }

}  // namespace inviwo
