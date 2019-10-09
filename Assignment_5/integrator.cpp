/*********************************************************************
 *  Author  : Himangshu Saikia
 *  Init    : Wednesday, September 20, 2017 - 12:04:15
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <labstreamlines/integrator.h>

namespace inviwo {
    
    // TODO: Implement a single integration step here
    
    std::vector<dvec2> Integrator::StreamLines(const VectorField2& vectorField, const dvec2& position,
            const double stepSize, bool backwardDirection,
            bool integrateDirectionField, double speedThreshold, int kernelSize) {
        
        std::vector<dvec2> streamLinePoints;
        double speed = DBL_MAX;
        dvec2 pos = position;
        dvec2 new_pos;
        int count = 0;
        
        while (vectorField.isInside(pos) && speed > speedThreshold && count <= kernelSize) {
            streamLinePoints.push_back(pos);
            
             new_pos = Integrator::RK4(vectorField, pos, stepSize, backwardDirection,
                                            integrateDirectionField);
            speed = sqrt(pow(new_pos.x - pos.x, 2) + pow(new_pos.y - pos.y, 2)) / double(stepSize);
            pos = new_pos;
            
            count++;
        }
        return streamLinePoints;
     }
    
    dvec2 Integrator::Euler(const VectorField2& vectorField, const dvec2& position,
                            const double stepSize, bool backwardDirection,
                            bool integrateDirectionField) {
        // Access the vector field with vectorField.interpolate(...)
        dvec2 euler_vec;
        if (backwardDirection) {
            euler_vec = -1 * vectorField.interpolate(position);
        } else {
            euler_vec = vectorField.interpolate(position);
        }
        if (integrateDirectionField) {
            euler_vec = euler_vec / sqrt(pow(euler_vec.x, 2) + pow(euler_vec.y, 2));
        }
        return position + stepSize * euler_vec;
    }
    
    dvec2 Integrator::RK4(const VectorField2& vectorField, const dvec2& position, const double stepSize,
                          bool backwardDirection, bool integrateDirectionField) {
        dvec2 v1, v2, v3, v4;
        if (backwardDirection) {
            v1 = -1 * vectorField.interpolate(position);
            v2 = -1 * vectorField.interpolate(position + stepSize / 2 * v1);
            v3 = -1 * vectorField.interpolate(position + stepSize / 2 * v2);
            v4 = -1 * vectorField.interpolate(position + stepSize * v3);
            
        } else {
            v1 = vectorField.interpolate(position);
            v2 = vectorField.interpolate(position + stepSize / 2 * v1);
            v3 = vectorField.interpolate(position + stepSize / 2 * v2);
            v4 = vectorField.interpolate(position + stepSize * v3);
        }
        if (integrateDirectionField) {
            if (sqrt(pow(v1.x, 2) + pow(v1.y, 2)) != 0) {
                v1 = v1 / sqrt(pow(v1.x, 2) + pow(v1.y, 2));
            }
            if (sqrt(pow(v2.x, 2) + pow(v2.y, 2)) != 0) {
                v2 = v2 / sqrt(pow(v2.x, 2) + pow(v2.y, 2));
            }
            if (sqrt(pow(v3.x, 2) + pow(v3.y, 2)) != 0) {
                v3 = v3 / sqrt(pow(v3.x, 2) + pow(v3.y, 2));
            }
            if (sqrt(pow(v4.x, 2) + pow(v4.y, 2)) != 0) {
                v4 = v4 / sqrt(pow(v4.x, 2) + pow(v4.y, 2));
            }
        }
        return position + stepSize * (1.0 / 6 * v1 + 1.0 / 3 * v2 + 1.0 / 3 * v3 + 1.0 / 6 * v4);
    }
}  // namespace inviwo
