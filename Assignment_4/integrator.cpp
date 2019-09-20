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

dvec2 Integrator::Euler(const VectorField2& vectorField, const dvec2& position, const double stepSize) {
     //Access the vector field with vectorField.interpolate(...)
    
    return position + stepSize * vectorField.interpolate(position);
    
 }

 dvec2 Integrator::RK4(const VectorField2& vectorField, const dvec2& position, const double stepSize) {
     dvec2 v1 = vectorField.interpolate(position);
     dvec2 v2 = vectorField.interpolate(position + stepSize / 2*v1);
     dvec2 v3 = vectorField.interpolate(position + stepSize / 2*v2);
     dvec2 v4 = vectorField.interpolate(position + stepSize*v3);
     
     return position + stepSize * (1.0/6 * v1 + 1.0/3 * v2 + 1.0/3 * v3 + 1.0/6 *v4);
 }

}  // namespace inviwo
