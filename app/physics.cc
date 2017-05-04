#include "physics.h"

namespace physics{

double dt = 0.0;


// Basic spline kernel for SPH 
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

// Gradient of spline kernel 
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



// Compute density based on neighbors
void computeDensity(body_holder* srch, std::vector<body_holder*>& nbsh)
{
  body* source = srch->getBody();
  double density = 0;
  assert(nbsh.size()>0);
  for(auto nbh : nbsh){
    body* nb = nbh->getBody();
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



// Computre pressure on a body 
// This function does not need the neighbors
void computePressure(body_holder* srch)
{
  body* source = srch->getBody();
  source->setPressure(source->getEntropy()*
      pow(source->getDensity(),kHeatRatio));
} // computePressure

// Compute the sound speed on a body 
// This function does not need the neighbors
void computeSoundspeed(body_holder* srch)
{
  body* source = srch->getBody();
  source->setSoundspeed(pow(kHeatRatio*
        source->getPressure()/source->getDensity(),1./2.));
} // computeSoundspeed

// Compute MU using a source and a neighbor
// This function is used in artificial viscosity (Hydro Force)
// and to computeDt
double mu(body* source, body* nb){
  double result = 0.0;
  double h_ij = (1./2.)*
    (source->getSmoothinglength()+nb->getSmoothinglength()); 
  space_vector_t vecVelocity = flecsi::point_to_vector(
      source->getVelocity() - nb->getVelocity());
  space_vector_t vecPosition = flecsi::point_to_vector(
      source->getPosition() - nb->getPosition());
  double dotproduct = flecsi::dot(vecVelocity,vecPosition);
  if(dotproduct >= 0.0)
    return result;
  // Should add norm to space_vector
  double dist = flecsi::distance(source->getPosition(),nb->getPosition());
  result = dotproduct / (h_ij*(((dist*dist)/(h_ij*h_ij))+kViscEta+kViscEta));
  assert(result < 0.0);
  return result; 
} // mu

// Conpute the Hydrostatic force 
void computeHydro(body_holder* srch, std::vector<body_holder*>& ngbsh){
  body* source = srch->getBody();

  point_t hydro = {0,0,0};
  for(auto nbh : ngbsh){
    body* nb = nbh->getBody();

    // Artificial viscosity
    double density_ij = (1./2.)*(source->getDensity()+nb->getDensity());
    double soundspeed_ij = (1./2.)*
      (source->getSoundspeed()+nb->getSoundspeed());
    double mu_ij = mu(source,nb);
    double viscosity = (-kViscAlpha*mu_ij*soundspeed_ij+kViscBeta*mu_ij*mu_ij)
      /density_ij;
    assert(viscosity>=0.0);

    // Hydro force
    point_t vecPosition = source->getPosition()-nb->getPosition();
    double pressureDensity = source->getPressure()/(source->getDensity()*
        source->getDensity())
      + nb->getPressure()/(nb->getDensity()*nb->getDensity());
    point_t sourcekernelgradient = gradKernel(
        vecPosition,source->getSmoothinglength());
    point_t nbkernelgradient = gradKernel(
        vecPosition,nb->getSmoothinglength());
    point_t resultkernelgradient = (1./2.)*
      (sourcekernelgradient+nbkernelgradient);

    hydro += source->getMass()*nb->getMass()*(pressureDensity+viscosity)
        *resultkernelgradient;
  }
  source->setHydroForce(hydro);
} // computeHydro

// Compute acceleration for the FGrav and FHydro 
// If we consider the relaxation step, we need FRoche or FRot 
void computeAcceleration(body_holder* srch, double dt){
  body* source = srch->getBody();
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

// moveBody, apply the acceleration to compute new velocity 
// and position. Using a leapfrog scheme 
void moveBody(body_holder* srch, double dt){
  body* source = srch->getBody();
  point_t position;
  point_t velocityhalf; 
  velocityhalf = source->getVelocity()+dt/2.0*source->getAcceleration();
  position = source->getPosition()+ dt*source->getVelocityhalf();
  source->setPosition(position);
  source->setVelocityhalf(velocityhalf);
}// moveBody

// Compute the gravitation for using all the others particles
// This is the N^2 version of the algorithm
void computeGrav(body_holder* srch, std::vector<body_holder*>& ngbh){
  body* source = srch->getBody();
  point_t grav = {0.,0.,0.};
  for(auto nb : ngbh){
    double dist = flecsi::distance(source->getPosition(),nb->getPosition());
    if(dist > 0.0){
      point_t vecPosition = nb->getPosition() - source->getPosition();
      grav += kGravConstant*source->getMass()*nb->getMass()/
        (dist*dist*dist) * vecPosition;
    }
  }
  source->setGravForce(grav);
}// computeGrav

double computeDt(body_holder* srch, std::vector<body_holder*>& ngbhs)
{
  body* source = srch->getBody();
  // First compute the dt based on acceleration norm
  double accelNorm = 0.0;
  for(size_t i=0;i<gdimension;++i)
    accelNorm += source->getAcceleration()[i]*source->getAcceleration()[i];
  accelNorm = sqrt(accelNorm);
  double dt1 = sqrt(source->getSmoothinglength()/accelNorm); 
 
  // Second based on max mu 
  double max_mu_ij = -9999;
  for(auto nbh: ngbhs)
  {
    body* nb = nbh->getBody();
    double local_mu = fabs(mu(source,nb));
    max_mu_ij = std::max(max_mu_ij,local_mu);  
  }
  double dt2 = source->getSmoothinglength()/
    (source->getSoundspeed()+
     1.2*kViscAlpha*source->getSoundspeed()+
     1.2+kViscBeta*max_mu_ij);
  dt2 *= kCoeffDt;

  return std::min(dt1,dt2);
}// computeDt

} // namespace 
