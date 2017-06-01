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

#include "flecsi/topology/tree_topology.h"
#include "flecsi/geometry/point.h"
#include "flecsi/geometry/space_vector.h"

#define OUTPUT

static const size_t gdimension = 3;
using type_t = double;

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
      const double mass,
      const double smoothinglength
  ):  position_(position), 
      velocity_(velocity),
      velocityhalf_(velocityhalf),
      acceleration_(acceleration),
      density_(density),
      pressure_(pressure),
      entropy_(entropy),
      mass_(mass),
      smoothinglength_(smoothinglength),
      soundspeed_(0),
      gravforce_(point_t{}),
      hydroforce_(point_t{})
  {};

  body()
  {};

  const point_t& coordinates() const{return position_;}
  const point_t& getPosition() const{return position_;}
  double getMass() const{return mass_;}
  double getSmoothinglength() const{return smoothinglength_;}
  double getPressure() const{return pressure_;}
  double getSoundspeed() const{return soundspeed_;}
  double getEntropy() const{return entropy_;}
  double getDensity() const{return density_;}
  point_t getVelocity() const{return velocity_;}
  point_t getHydroForce() const{return hydroforce_;}
  point_t getGravForce() const{return gravforce_;}
  point_t getVelocityhalf() const{return velocityhalf_;}
  point_t getAcceleration() const{return acceleration_;}
  double getInternalenergy() const{return internalenergy_;}
  double getDudt(){return dudt_;};

  void setPosition(point_t position){position_ = position;}
  void setAcceleration(point_t acceleration){acceleration_ = acceleration;}
  void setVelocity(point_t velocity){velocity_ = velocity;}
  void setVelocityhalf(point_t velocityhalf){velocityhalf_ = velocityhalf;}
  void setGravForce(point_t gravforce){gravforce_ = gravforce;}
  void setHydroForce(point_t hydroforce){hydroforce_ = hydroforce;}
  void setSoundspeed(double soundspeed){soundspeed_ = soundspeed;}
  void setPressure(double pressure){assert(pressure > 0);pressure_ = pressure;}
  void setDensity(double density){assert(density >= 0);density_ = density;}
  void setMass(double mass){mass_ = mass;};
  void setInternalenergy(double internalenergy)
    {internalenergy_=internalenergy;};
  void setSmoothinglength(double smoothinglength)
    {smoothinglength_=smoothinglength;};
  void setDudt(double dudt){dudt_ = dudt;};
  void setDt(double dt){dt_ = dt;};
  void setId(int64_t id){id_ = id;};

  friend std::ostream& operator<<(std::ostream& os, const body& b){
    // TODO change regarding to dimension 
    os << "Particle: Pos: " <<b.position_ << " Dens: " << b.density_; 
    os << " H: " << b.smoothinglength_;
    os << " Pres: " << b.pressure_;
    os << " Vel: " << b.velocity_ << " VelH: " << b.velocityhalf_;
    os << " Mass: " << b.mass_;
    os << " Force: hyd: " << b.hydroforce_;
    os << " grav: " << b.gravforce_;
    os << "          Acc: " << b.acceleration_;
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
  double dudt_;
  point_t gravforce_;
  point_t hydroforce_;
  double dt_;
  int64_t id_;

}; // class body 
  
#endif // body_h

