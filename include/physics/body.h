/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Los Alamos National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/

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
 * @file body.h
 * @author Julien Loiseau
 * @date April 2017
 * @brief Representation of a body, a particle for our SPH implementation 
 */

#ifndef body_h
#define body_h

#define OUTPUT

#include "tree_topology.h" // CHANGE TO FLECSI ONE!!!
#include "flecsi/geometry/point.h"
#include "flecsi/geometry/space_vector.h"

#include "user.h"

class body{

  static const size_t dimension = gdimension;
  using element_t = type_t; 
  using point_t = flecsi::point<element_t, dimension>;
  
public:
  
  body(const point_t& position, 
      const point_t& velocity, 
      const point_t& velocityhalf, 
      const point_t& acceleration,
      const double density, 
      const double pressure, 
      const double entropy, 
      const double electronFraction,
      const double mass,
      const double smoothinglength
  ):  position_(position), 
      velocity_(velocity),
      velocityhalf_(velocityhalf),
      acceleration_(acceleration),
      density_(density),
      pressure_(pressure),
      //entropy_(entropy),
      mass_(mass),
      //electronFraction_(electronFraction),
      smoothinglength_(smoothinglength),
      soundspeed_(0.0)
      ,internalenergy_(0.0)
      ,totalenergy_(0.0)
      ,dudt_(0.0)
      ,dedt_(0.0)
      ,adiabatic_(0.0)
      ,dadt_(0.0)
      //gravforce_(point_t{}),
      //hydroforce_(point_t{})
  {};

  body()
  {};

  const point_t& coordinates() const{return position_;}
  const point_t& getPosition() const{return position_;}
  double getMass() const{return mass_;}
  double getSmoothinglength() const{return smoothinglength_;}
  double getPressure() const{return pressure_;}
  double getSoundspeed() const{return soundspeed_;}
  //double getEntropy() const{return entropy_;}
  //double getElectronFraction() const{return electronFraction_;}
  double getDensity() const{return density_;}
  point_t getVelocity() const{return velocity_;}
  //point_t getHydroForce() const{return hydroforce_;}
  //point_t getGravForce() const{return gravforce_;}
  point_t getVelocityhalf() const{return velocityhalf_;}
  point_t getAcceleration() const{return acceleration_;}

  point_t getLinMomentum() const { 
    point_t res = {};
    for(size_t i = 0 ; i < dimension; ++i){
      res[i] = velocity_[i] * mass_;
    }
    return res;
  };
  //double getDudt(){return dudt_;};
  int64_t getId(){return id_;};
  double getDt(){return dt_;};
  int getType(){return type_;}; 

  bool is_wall(){return type_ == 1;};

  void setPosition(point_t position){position_ = position;}
  void setAcceleration(point_t acceleration){acceleration_ = acceleration;}
  void setVelocity(point_t velocity){velocity_ = velocity;}
  void setVelocityhalf(point_t velocityhalf){velocityhalf_ = velocityhalf;}
  //void setGravForce(point_t gravforce){gravforce_ = gravforce;}
  //void setHydroForce(point_t hydroforce){hydroforce_ = hydroforce;}
  void setSoundspeed(double soundspeed){soundspeed_ = soundspeed;}
  void setPressure(double pressure){pressure_ = pressure;}
  void setDensity(double density){density_ = density;}
  void setMass(double mass){mass_ = mass;};
  //void setLinMomentum(point_t lin_momentum){lin_momentum_ = lin_momentum;}
  void setSmoothinglength(double smoothinglength)
    {smoothinglength_=smoothinglength;};
 // void setDudt(double dudt){dudt_ = dudt;};
  void setDt(double dt){dt_ = dt;};
  void setId(int64_t id){id_ = id;};
  void setType(int type){type_ = type;}; 

  // Dependent of the problem 
    double getInternalenergy() const{return internalenergy_;}
    void setInternalenergy(double internalenergy)
        {internalenergy_=internalenergy;}
    double getTotalenergy() const{return totalenergy_;}
    void setTotalenergy(double totalenergy) {totalenergy_=totalenergy;}
    void setDudt(double dudt){dudt_ = dudt;};
    void setDedt(double dedt){dedt_ = dedt;};
    double getDudt(){return dudt_;};
    double getDedt(){return dudt_;};
    double getAdiabatic() const{return adiabatic_;}
    void setAdiabatic(double adiabatic){adiabatic_ = adiabatic;};
    double getDadt() const{return dadt_;};
    void setDadt(double dadt){dadt_ = dadt;};
  #ifdef VERLET
    double getDensityNM1(){return densityNM1_;};
    point_t getVelocityCor(){return velocityCor_;};
    double getDensityDt(){return densityDt_;};
    void setDensityNM1(double densityNM1){densityNM1_ = densityNM1;};
    void setVelocityCor(point_t velocityCor){velocityCor_ = velocityCor;};
    void setDensityDt(double densityDt){densityDt_ = densityDt;};
    point_t getVelocityTmp(){return velocityTmp_;};
    double getDensityTmp(){return densityTmp_;};
    point_t getVelocityNM1(){return velocityNM1_;};
    void setVelocityTmp(point_t velocityTmp){velocityTmp_ = velocityTmp;};
    void setDensityTmp(double densityTmp){densityTmp_ = densityTmp;};
    void setVelocityNM1(point_t velocityNM1){velocityNM1_=velocityNM1;};
    double getMaxVisc(){return maxVisc_;};
    void setMaxVisc(double maxVisc){maxVisc_ = maxVisc;};
  #endif


  friend std::ostream& operator<<(std::ostream& os, const body& b){
    // TODO change regarding to dimension 
    os << std::setprecision(10); 
    os << "Particle: Pos: " <<b.position_ << " rho: " << b.density_; 
    os << " h: " << b.smoothinglength_;
    os << " P: " << b.pressure_;
    os << " v: " << b.velocity_ ;//<< " VelH: " << b.velocityhalf_;
    os << " m: " << b.mass_;
    #ifdef INTERNAL_ENERGY
      os << " u: " << b.internalenergy_;
    #endif 
    os << " cs: " << b.soundspeed_;
    //os << " Force: hyd: " << b.hydroforce_;
    //os << " grav: " << b.gravforce_;
    os << " a: " << b.acceleration_;
    os << " id: " << b.id_; 
    return os;
  }      

private:
  point_t position_;
  point_t velocity_;
  point_t velocityhalf_;
  point_t acceleration_;
  double density_;
  double pressure_; 
  double entropy_;
  double mass_;
  double smoothinglength_; 
  double soundspeed_;
  double internalenergy_;
  double totalenergy_;
  double dudt_;
  double dedt_;
  double adiabatic_; 
  double dadt_;
  #ifdef VERLET
  double densityNM1_;
  point_t velocityCor_;
  double densityDt_;
  point_t velocityTmp_;
  double densityTmp_;
  point_t velocityNM1_;
  double maxVisc_; 
  #endif 
  //point_t lin_momentum_; //TODO : Need to check
  //double dudt_;
  //point_t gravforce_;
  //point_t hydroforce_;
  double dt_;
  int64_t id_;
  int type_; 
}; // class body 
  
#endif // body_h

