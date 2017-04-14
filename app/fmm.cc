#include <complex>

namespace fmm{


// Particle to Particle, compute normal potential 
void P2P(body * target, vector<body*>& sources){
  double gravpotential = 0.0;
  point_t position = {0.,0.,0.};
  for(auto srci : sources){
    point_t vecposition = target->getPosition()-srci->getPosition(); 
    double dist = distance(target->getPosition(),srci->getPosition());
    // Avoid the target itself 
    if(dist > 0.0){
      gravpotentiel += target->getMass()* 1./dist; 
      position -= vecposition; 
    }// if 
  } // for
}// P2P

// Particle to Multipole 
void P2M(point_t& Mposition, ){
  complex_t Ymn[P*P], YnmTheta[P*P];
  
}


} // namespace fmm
