/*~--------------------------------------------------------------------------~*
 * 
 * /@@@@@@@@  @@           @@@@@@   @@@@@@@@ @@@@@@@  @@      @@
 * /@@/////  /@@          @@////@@ @@////// /@@////@@/@@     /@@
 * /@@       /@@  @@@@@  @@    // /@@       /@@   /@@/@@     /@@
 * /@@@@@@@  /@@ @@///@@/@@       /@@@@@@@@@/@@@@@@@ /@@@@@@@@@@
 * /@@////   /@@/@@@@@@@/@@       ////////@@/@@////  /@@//////@@
 * /@@       /@@/@@//// //@@    @@       /@@/@@      /@@     /@@
 * /@@       @@@//@@@@@@ //@@@@@@  @@@@@@@@ /@@      /@@     /@@
 * //       ///  //////   //////  ////////  //       //      //  
 *
 *~--------------------------------------------------------------------------~*/

/**
 * @file sodtube.h
 * @author Julien Loiseau
 * @date April 2017
 * @brief Physics functions for the 1D Sod Tube implementation
 */

#include <vector>

#include "sodtube.h"

inline
point_t
operator*(
    const point_t p,const double val)
{
  point_t p1 = p;
  for(size_t i=0;i<gdimension;++i)
    p1[i]*=val;
  return p1;
}

inline
point_t
operator*(
    const double val, const point_t p)
{
  point_t p1 = p;
  for(size_t i=0;i<gdimension;++i)
    p1[i]*=val;
  return p1;
}

inline
point_t
operator*(
    const point_t p1, const point_t p2)
{
  point_t p3 = p1;
  for(size_t i=0;i<gdimension;++i)
    p3[i]*=p2[i];
  return p3;
}

inline 
bool 
operator==(
    const point_t& p1, const point_t& p2
    )
{
  for(size_t i=0;i<gdimension;++i)
    if(p1[i] != p2[i])
      return false;
  return true;
}


namespace sodtube{

const double kDt = 0.0025;

void randomDataSodTube1D(
    std::vector<std::pair<entity_key_t,body>>& bodies,
    int& nbodies,
    int& totalnbodies,
    int rank,
    int size)
{
  srand(time(NULL)*rank);
  // Default 2000 bodies 
  if(totalnbodies == 0)
    totalnbodies = 400;
  // Each process generates its bodies 
  nbodies = totalnbodies/size;
  if(rank == size-1)
    nbodies = totalnbodies-((size-1)*(totalnbodies/size));

  //double procposition = 1./(double)size;
  //double particlespace = 1/(double)totalnbodies;
  double distance = 0.0026566;
  double middle = totalnbodies*distance/2.;
  double procstart = distance*(totalnbodies/size)*rank;
  double procstop = procstart+(distance*totalnbodies/size);
  std::cout<<rank<<": ["<<procstart<<";"<<procstop<<"["<<std::endl;
  double position = procstart;
  double velocity = 0.0;
  double internalenergy = 2.5;
  double density = 1;
  point_t acceleration = {0};
  double mass = 2.65e-3;
  double smoothinglength = 1.0e-2;
  for(int i=0;i<nbodies;++i)
  {
    // Create empty body 
    body bi;
    if(position > middle)
    {
      internalenergy = 2;
      density = 0.125;
      mass = 3.3125e-4;
    }
    bi.setPosition(point_t{position});
    bi.setDensity(density);
    bi.setInternalenergy(internalenergy);
    bi.setMass(mass);
    bi.setVelocity(point_t{velocity});
    bi.setSmoothinglength(smoothinglength);
    bi.setAcceleration(acceleration);
    // Positions between [procstart;procstop[
    bodies.push_back(std::make_pair(entity_key_t::null(),bi));
    position += distance;
  } 
}

double kernel(double dist, double h)
{
  double sigma = 2./(3.*h);
  double q = dist/h;
  assert(q>=0);
  if(q<1)
    return sigma*(0.25*(2.0-q)*(2.0-q)*(2.0-q)-(1.0-q)*(1.0-q)*(1.0-q));
  if(q>=1 && q<2)
    return sigma*(0.25*(2.0-q)*(2.0-q)*(2.0-q));
  return 0.0;
}

double gradkernel(double dist, double h)
{
  double sigma = 2./(3.*h);
  double q = dist/h;
  assert(q>=0);
  if(q<1)
    return -sigma*(-0.75*(2.0-q)*(2.0-q)+3.*(1.0-q)*(1.0-q))/h;
  if(q>=1 && q<2)
    return -sigma*(-0.75*(2.0-q)*(2.0-q))/h;
  return 0.0;
}

void computeDensityApply(
    body_holder * nb, 
    body_holder * src)
{
  body * srcb = src->getBody();
  assert(srcb!=nullptr);
  double density = srcb->getDensity();
  // I am at least only in the vector 
  body * nbb = nb->getBody();
  assert(nbb!=nullptr);
    
  double dist = flecsi::distance(srcb->getPosition(),nbb->getPosition());
  double kern = kernel(dist,srcb->getSmoothinglength());
  density += nbb->getMass()*kern;
  
  srcb->setDensity(density);
}

void computeDensity(body_holder * src, std::vector<body_holder*>& neighb )
{
  body * srcb = src->getBody();
  assert(srcb!=nullptr);
  double density = 0.0;
  // I am at least only in the vector 
  assert(neighb.size() > 0);
  for(auto nb: neighb)
  {
    body * nbb = nb->getBody();
    assert(nbb!=nullptr);
    
    double dist = flecsi::distance(srcb->getPosition(),nbb->getPosition());
    double kern = kernel(dist,srcb->getSmoothinglength());
    density += nbb->getMass()*kern;
  }
  srcb->setDensity(density);
}

void computePressureSoundSpeed(body_holder * src)
{
  body * srcb = src->getBody();
  assert(srcb!=nullptr);
  double pressure = (1.4-1.0)*srcb->getDensity()*srcb->getInternalenergy();
  srcb->setPressure(pressure);
  double cs = sqrt((1.4-1.0)*srcb->getInternalenergy());
  srcb->setSoundspeed(cs);
}

void computeViscosity(body_holder * src, std::vector<body_holder*>& neighb )
{
  body * srcb = src->getBody();
  assert(srcb!=nullptr);
  double alpha = 2.0;
  double beta = 1.0;
  point_t acc = srcb->getAcceleration();
  double dudt = srcb->getDudt();
  assert(neighb.size()!=0);
  for(auto nb: neighb)
  {
    body * nbb = nb->getBody();
    assert(nbb!=nullptr);
    if(srcb->getPosition() == nbb->getPosition())
      continue;
    
    double dist = flecsi::distance(srcb->getPosition(),nbb->getPosition());
    point_t xv = (srcb->getPosition()-nbb->getPosition())
      *(srcb->getVelocity()-nbb->getVelocity());
    point_t rhat = (srcb->getPosition()-nbb->getPosition())/dist;
    point_t Pi = {}; 
    if(xv[0] < 0.0)
    {
      double average_c = 0.5 * (srcb->getSoundspeed() + nbb->getSoundspeed());
      double rhobar = 0.5*(srcb->getDensity()+nbb->getDensity());
      double hbar = 0.5*(srcb->getSmoothinglength()+nbb->getSmoothinglength());
      point_t mu = hbar*xv/(dist*dist + 0.01*(hbar*hbar));
      Pi = (-alpha*average_c*mu+beta*mu*mu)/rhobar;
    }

    double gradkern = gradkernel(dist,srcb->getSmoothinglength());
    acc += nbb->getMass()*Pi*gradkern*rhat;
    point_t vij = srcb->getVelocity() - nbb->getVelocity();
    dudt += (-0.5*Pi*nbb->getMass()*vij*gradkern*rhat)[0];

  }
  srcb->setDudt(dudt);
  srcb->setAcceleration(acc);
}

void computeAcceleration(body_holder * src, std::vector<body_holder*>& neighb)
{
  body * srcb = src->getBody();
  assert(srcb!=nullptr);
  point_t acc = {0};
  double dudt = 0;
  assert(neighb.size()!=0);
  for(auto nb: neighb)
  {
    body * nbb = nb->getBody();
    assert(nbb!=nullptr);
    
    if(nbb->getPosition() == srcb->getPosition())
      continue;
    
    double dist = flecsi::distance(srcb->getPosition(),nbb->getPosition());
    point_t rhat = (srcb->getPosition()-nbb->getPosition())/dist;
    double gradkern = gradkernel(dist,srcb->getSmoothinglength());
    acc += nbb->getMass()*
      (srcb->getPressure()/pow(srcb->getDensity(),2)+
       nbb->getPressure()/pow(nbb->getDensity(),2))*gradkern*rhat;
    point_t vij = srcb->getVelocity()-nbb->getVelocity();
    dudt += (-srcb->getPressure()/pow(srcb->getDensity(),2)*nbb->getMass()*vij
        *gradkern*rhat)[0];
  }
  srcb->setDudt(dudt);
  srcb->setAcceleration(acc);
}

void moveParticle(body_holder * src,std::array<point_t,2>& range)
{
  body * srcb = src->getBody();
  assert(srcb!=nullptr);
  if(srcb->getPosition()[0]>0.1
      && srcb->getPosition()[0]<0.9){
    srcb->setPosition(srcb->getPosition()+srcb->getVelocity()*kDt);
    srcb->setVelocity(srcb->getVelocity()+srcb->getAcceleration()*kDt);
    srcb->setInternalenergy(srcb->getInternalenergy()+srcb->getDudt()*kDt);
  }

}

}// namespace sodtube


