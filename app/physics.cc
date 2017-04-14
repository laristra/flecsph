#include "physics.h"

namespace physics{

double kernel(double r, double h){
  double rh = r/h;
  double result = 1.0/(M_PI*pow(h,3));
  if (0.0 <= rh && rh < 1.0) {
    result *= 1.0 - (3.0/2.0) * pow(rh,2) + (3.0/4.0) * pow(rh,3); 
    return result; 
  }else if (1.0 <= rh && rh < 2.0) {
    result *= (1.0/4.0) * pow(2-rh, 3);
    return result;
  }
  return 0.0;
} // kernel

point_t gradKernel(point_t vecP, double h){
  double coeff = 1.0/(M_PI*pow(h,4));
  double r = sqrt(vecP[0]*vecP[0]+vecP[1]*vecP[1]+vecP[2]*vecP[2]);
  double rh = r/h;

  point_t result = {0,0,0};
  if(0.0 <= rh && rh < 1.0){
    result = coeff*vecP;
    result *= ((-3.0/h)+(9.0*r/(4.0*h*h)));
  }else if(1.0 <= rh && rh < 2.0){
    result = coeff*vecP;
    result *= ((-3.0/r)+(3.0/h)+(-3.0*r/(4.0*h*h)));
  }
  return result;

} // gradKernel 

void computeDensity(tree_policy::body* source, std::vector<tree_policy::body*>& neighbors){
  double density = 0;
  assert(neighbors.size()>0);
  for(auto nb : neighbors){
    double dist = flecsi::distance(source->getPosition(),nb->getPosition());
    assert(dist>=0.0);
    double kernelresult = (1./2.)*(
        kernel(dist,source->getSmoothinglength())+
        kernel(dist,nb->getSmoothinglength()));
    assert(kernelresult>=0.0);
    density += kernelresult*nb->getMass();
  } // for
  assert(density>0);
  source->setDensity(density);
} // computeDensity

void computePressure(body* source){
  source->setPressure(source->getEntropy()*pow(source->getDensity(),kHeatRatio));
} // computePressure

void computeSoundspeed(body* source){
  source->setSoundspeed(pow(kHeatRatio*source->getPressure()/source->getDensity(),1./2.));
} // computeSoundspeed

double mu(body* source, body* nb){
  double result = 0.0;
  double h_ij = (1./2.)*(source->getSmoothinglength()+nb->getSmoothinglength());
  
  space_vector_t vecVelocity = flecsi::point_to_vector(source->getVelocity() - nb->getVelocity());
  space_vector_t vecPosition = flecsi::point_to_vector(source->getPosition() - nb->getPosition());
  double dotproduct = flecsi::dot(vecVelocity,vecPosition);
  
  if(dotproduct >= 0.0)
    return result;

  // Should add norm to space_vector
  double dist = flecsi::distance(source->getPosition(),nb->getPosition());
  result = dotproduct / (h_ij*(((dist*dist)/(h_ij*h_ij))+kViscEta+kViscEta));
  assert(result < 0.0);
  return result; 
} // mu

void computeHydro(body* source, std::vector<body*>& neighbors){
  point_t hydro = {0,0,0};
  for(auto nb : neighbors){

    // Artificial viscosity
    double density_ij = (1./2.)*(source->getDensity()+nb->getDensity());
    double soundspeed_ij = (1./2.)*(source->getSoundspeed()+nb->getSoundspeed());
    double mu_ij = mu(source,nb);
    double viscosity = (-kViscAlpha*mu_ij*soundspeed_ij+kViscBeta*mu_ij*mu_ij)/density_ij;
    assert(viscosity>=0.0);

    // Hydro force
    point_t vecPosition = source->getPosition()-nb->getPosition();
    double pressureDensity = source->getPressure()/(source->getDensity()*source->getDensity())
      + nb->getPressure()/(nb->getDensity()*nb->getDensity());
    point_t sourcekernelgradient = gradKernel(vecPosition,source->getSmoothinglength());
    point_t nbkernelgradient = gradKernel(vecPosition,nb->getSmoothinglength());
    point_t resultkernelgradient = (1./2.)*(sourcekernelgradient+nbkernelgradient);

    hydro += source->getMass()*nb->getMass()*(pressureDensity+viscosity)
        *resultkernelgradient;
  }
  source->setHydroForce(hydro);
}

void computeAcceleration(body* source, double dt){
  point_t acceleration = {0.,0.,0.};
  acceleration += source->getHydroForce();
  acceleration += source->getGravForce();
  acceleration /= source->getMass();
  
  // Use it in the new velocity 
  point_t velocity = {0.,0.,0.};
  velocity = source->getVelocityhalf() + dt/2.0*acceleration;

  source->setAcceleration(acceleration);
  source->setVelocity(velocity);
}

void moveBody(body* source, double dt){
  point_t position;
  point_t velocityhalf; 
  velocityhalf = source->getVelocity()+dt/2.0*source->getAcceleration();
  position = source->getPosition()+ dt*source->getVelocityhalf();
  source->setPosition(position);
  source->setVelocityhalf(velocityhalf);
}

void computeGrav(body* source, std::vector<body*>& neighbor){
  point_t grav = {0.,0.,0.};
  for(auto nb : neighbor){
    double dist = flecsi::distance(source->getPosition(),nb->getPosition());
    if(dist > 0.0){
      point_t vecPosition = nb->getPosition() - source->getPosition();
      grav += kGravConstant*source->getMass()*nb->getMass()/(dist*dist*dist) * vecPosition;
    }
  }
  source->setGravForce(grav);
}

} // namespace 
