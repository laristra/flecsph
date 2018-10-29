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

#include "tree_topology/tree_topology.h" // CHANGE TO FLECSI ONE!!!
#include "flecsi/geometry/point.h"
#include "flecsi/geometry/space_vector.h"

#include "user.h"

enum particle_type_t {NORMAL = 0 ,WALL = 1};

class body{

  static const size_t dimension = gdimension;
  using element_t = type_t; 
  using point_t = flecsi::point__<element_t, dimension>;
 

public:
  
  body(const point_t& position, 
      const point_t& velocity, 
      const point_t& velocityhalf, 
      const point_t& acceleration,
      const double density, 
      const double pressure, 
      const double entropy, 
      const double electronfraction,
      const double mass,
      const double smoothinglength
  ):  position_(position), 
      velocity_(velocity),
      velocityhalf_(velocityhalf),
      acceleration_(acceleration),
      density_(density),
      pressure_(pressure),
      entropy_(entropy),
      electronfraction_(electronfraction),
      mass_(mass),
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

   body(): type_(NORMAL)
   {};

  const point_t& coordinates() const{return position_;}   
  const point_t& getPosition() const{return position_;}
  double getMass() const{return mass_;}
  double getSmoothinglength() const{return smoothinglength_;}
  double getPressure() const{return pressure_;}
  double getSoundspeed() const{return soundspeed_;}
  double getEntropy() const{return entropy_;}
  double getElectronfraction() const{return electronfraction_;}
  double getDensity() const{return density_;}
  point_t getVelocity() const{return velocity_;}
  point_t getVelocityhalf() const{return velocityhalf_;}
  point_t getAcceleration() const{return acceleration_;}
  uint64_t neighbors(){return neighbors_;}
  particle_type_t type(){return type_;};

  point_t getLinMomentum() const { 
    point_t res = {};
    for(size_t i = 0 ; i < dimension; ++i){
      res[i] = velocity_[i] * mass_;
    }
    return res;
  };
  flecsi::topology::entity_id_t getId(){return id_;};
  flecsi::topology::entity_id_t id(){return id_;};
  double getDt(){return dt_;};
  double getMumax(){return mumax_;}
  int getType(){return type_;}; 

  bool is_wall(){return type_ == 1;};

  void set_neighbors(uint64_t neighbors){neighbors_ = neighbors;}
  void setPosition(point_t position){position_ = position;}
  void setAcceleration(point_t acceleration){acceleration_ = acceleration;}
  void setVelocity(point_t velocity){velocity_ = velocity;}
  void setVelocityhalf(point_t velocityhalf){velocityhalf_ = velocityhalf;}
  void setSoundspeed(double soundspeed){soundspeed_ = soundspeed;}
  void setPressure(double pressure){pressure_ = pressure;}
  void setEntropy(double entropy){entropy_ = entropy;}
  void setElectronfraction(double electronfraction){electronfraction_ = electronfraction;}
  void setDensity(double density){density_ = density;}
  void setMass(double mass){mass_ = mass;};
  void setSmoothinglength(double smoothinglength)
    {smoothinglength_=smoothinglength;};
  void setDt(double dt){dt_ = dt;};
  void setMumax(double mumax){mumax_ = mumax;};
  void setId(flecsi::topology::entity_id_t id){id_ = id;};
  void setType(particle_type_t type){type_ = type;};
  void setType(int type){type_= static_cast<particle_type_t>(type);};

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
  double electronfraction_;
  double mass_;
  double smoothinglength_; 
  double soundspeed_;
  double internalenergy_;
  double totalenergy_;
  double dudt_;
  double dedt_;
  double adiabatic_; 
  double dadt_;
  double dt_;
  double mumax_;
  flecsi::topology::entity_id_t id_;
  particle_type_t type_;
  int64_t neighbors_;  
}; // class body 
  
#endif // body_h

