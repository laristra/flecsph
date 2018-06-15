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

#warning "CHANGE TO FLECSI ONE"
#include "tree_topology.h"
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
      const point_t& velocityTmp,
      const point_t& velocityCor,
      const point_t& velocityNM1, 
      const point_t& acceleration,
      const double density, 
      const double densityDt,
      const double densityTmp,
      const double densityNM1,
      const double pressure, 
      const double mass,
      const double smoothinglength,
      const double soundspeed, 
      const double maxvisc,
      const double dt, 
      const int64_t id, 
      const int type
  ):  position_(position), 
      velocity_(velocity),
      velocityTmp_(velocityTmp),
      velocityCor_(velocityCor),
      velocityNM1_(velocityNM1),
      acceleration_(acceleration),
      density_(density),
      densityDt_(densityDt),
      densityTmp_(densityTmp),
      densityNM1_(densityNM1),
      pressure_(pressure),
      mass_(mass),
      smoothinglength_(smoothinglength),
      soundspeed_(soundspeed),
      maxvisc_(maxvisc),
      dt_(dt),
      id_(id),
      type_(type)
  {};

  body()
  {};

  const point_t& coordinates() const{return position_;}
  const point_t& getPosition() const{return position_;}
  double getMass() const{return mass_;}
  double getSmoothinglength() const{return smoothinglength_;}
  double getPressure() const{return pressure_;}
  double getSoundspeed() const{return soundspeed_;}
  
  double getDensity() const{return density_;}
  double getDensityNM1() const{return densityNM1_;}
  double getDensityTmp() const{return densityTmp_;}
  double getDensityDt() const{return densityDt_;}

  point_t getVelocity() const{return velocity_;}
  point_t getVelocityCor() const{return velocityCor_;}
  point_t getVelocityTmp() const{return velocityTmp_;}
  point_t getVelocityNM1() const{return velocityNM1_;}

  point_t getAcceleration() const{return acceleration_;}
  int64_t getId(){return id_;};
  double getDt(){return dt_;};
  int getType(){return type_;}; 
  double getMaxVisc(){return maxvisc_;}

  bool is_wall(){return type_ == 1;};

  void setPosition(point_t position){position_ = position;}
  void setAcceleration(point_t acceleration){acceleration_ = acceleration;}
  
  void setVelocity(point_t velocity){velocity_ = velocity;}
  void setVelocityTmp(point_t velocityTmp){velocityTmp_ = velocityTmp;}
  void setVelocityNM1(point_t velocityNM1){velocityNM1_ = velocityNM1;}
  void setVelocityCor(point_t velocityCor){velocityCor_ = velocityCor;}
  
  void setSoundspeed(double soundspeed){soundspeed_ = soundspeed;}
  void setPressure(double pressure){pressure_ = pressure;}

  void setDensity(double density){density_ = density;}
  void setDensityDt(double densityDt){densityDt_ = densityDt;}
  void setDensityTmp(double densityTmp){densityTmp_ = densityTmp;}
  void setDensityNM1(double densityNM1){densityNM1_ = densityNM1;}
  
  void setMass(double mass){mass_ = mass;};
  void setSmoothinglength(double smoothinglength)
    {smoothinglength_=smoothinglength;};
  void setDt(double dt){dt_ = dt;};
  void setId(int64_t id){id_ = id;};
  void setType(int type){type_ = type;}; 
  void setMaxVisc(double maxvisc){maxvisc_ = maxvisc;}

  friend std::ostream& operator<<(std::ostream& os, const body& b){
    // TODO change regarding to dimension 
    os << std::setprecision(10); 
    os << "Particle: Pos: " <<b.position_ << " rho: " << b.density_; 
    os << " h: " << b.smoothinglength_;
    os << " P: " << b.pressure_;
    os << " v: " << b.velocity_ ;//<< " VelH: " << b.velocityhalf_;
    os << " m: " << b.mass_;
    os << " cs: " << b.soundspeed_;
    os << " a: " << b.acceleration_;
    os << " id: " << b.id_; 
    return os;
  }      

private:
  point_t position_;
  point_t velocity_;
  point_t velocityTmp_;
  point_t velocityCor_;
  point_t velocityNM1_;
  point_t acceleration_;
  double density_;
  double densityDt_; 
  double densityTmp_;
  double densityNM1_;
  double pressure_; 
  double mass_;
  double smoothinglength_; 
  double soundspeed_;
  double maxvisc_;
  double dt_;
  int64_t id_;
  int type_; 
}; // class body 
  
#endif // body_h

