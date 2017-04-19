#ifndef tree_h
#define tree_h

#include <vector>

#include "flecsi/topology/tree_topology.h"
#include "flecsi/geometry/point.h"
#include "flecsi/geometry/space_vector.h"

using namespace flecsi;

namespace flecsi{
namespace execution{
void specialization_driver(int argc, char * argv[]);
void driver(int argc, char*argv[]); 
} // namespace execution
} // namespace flecsi 

class tree_policy{
public:
  using tree_t = flecsi::topology::tree_topology<tree_policy>;
  using branch_int_t = uint64_t;
  static const size_t dimension = 3;
  using element_t = double; 
  using point_t = flecsi::point<element_t, dimension>;
  using space_vector_t = flecsi::space_vector<element_t,dimension>;

  class body  : public flecsi::topology::tree_entity<branch_int_t, dimension>{
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
        gravforce_(point_t(0.,0.,0.)),
        hydroforce_(point_t(0.,0.,0.))
    {};

    const point_t& coordinates() const{
      return position_; 
    }

    const point_t& getPosition() const{
      return position_;
    }
  
    double getMass() const{
      return mass_;
    }

    double getSmoothinglength() const{
      return smoothinglength_;
    }

    double getPressure() const{
      return pressure_; 
    }

    double getSoundspeed() const{
      return soundspeed_; 
    }

    double getEntropy() const{
      return entropy_; 
    }

    double getDensity() const{
      return density_;
    }

    point_t getVelocity() const{
      return velocity_;
    }
  
    point_t getHydroForce() const{
      return hydroforce_;
    }

    point_t getGravForce() const{
      return gravforce_;
    }

    point_t getVelocityhalf() const{
      return velocityhalf_;
    }

    point_t getAcceleration() const{
      return acceleration_;
    }

    void setPosition(point_t position){
      position_ = position;
    }

    void setAcceleration(point_t acceleration){
      acceleration_ = acceleration;
    }

    void setVelocity(point_t velocity){
      velocity_ = velocity; 
    }

    void setVelocityhalf(point_t velocityhalf){
      velocityhalf_ = velocityhalf;
    }

    void setGravForce(point_t gravforce){
      gravforce_ = gravforce;
    }

    void setHydroForce(point_t hydroforce){
      hydroforce_ = hydroforce; 
    }

    void setSoundspeed(double soundspeed){
      soundspeed_ = soundspeed;
    }

    void setPressure(double pressure){
      assert(pressure > 0);
      pressure_ = pressure;
    }

    void setDensity(double density){
      assert(density >= 0);
      density_ = density;
    }

    friend std::ostream& operator<<(std::ostream& os, const body& b){
      os << "Particle: Pos: " <<b.position_ << " Dens: " << b.density_; 
      os << " H: " << b.smoothinglength_;
      os << " Pres: " << b.pressure_ << std::endl;
      os << "          Vel: " << b.velocity_ << " VelH: " << b.velocityhalf_;
      os << " Mass: " << b.mass_ << std::endl;
      os << "   Force: hyd: " << b.hydroforce_;
      os << " grav: " << b.gravforce_ << std::endl;
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
    point_t gravforce_;
    point_t hydroforce_;
  }; // class body 
  
  using entity_t = body;

  class branch : public flecsi::topology::tree_branch<branch_int_t,dimension>{
  public:
    branch(){}

    void insert(body* ent){
      ents_.push_back(ent);
      if(ents_.size() > 8){
        refine();
      }
    } // insert
    
    auto begin(){
      return ents_.begin();
    }

    auto end(){
      return ents_.end();
    }

    auto clear(){
      ents_.clear();
    }

    void remove(body* ent){
      auto itr = find(ents_.begin(), ents_.end(), ent); 
      assert(itr!=ents_.end());
      ents_.erase(itr);
      if(ents_.empty()){
        coarsen();
      } 
    }

    point_t 
    coordinates(
        const std::array<flecsi::point<element_t, dimension>,2>& range) const{
      point_t p;
      branch_id_t bid = id(); 
      bid.coordinates(range,p);
      return p;
    }

   private:
    std::vector<body*> ents_;
  }; // class branch 

  bool should_coarsen(branch* parent){
    return true;
  }

  using branch_t = branch;

}; // class tree_policy

using tree_topology_t = flecsi::topology::tree_topology<tree_policy>;
using body = tree_topology_t::body;
using point_t = tree_topology_t::point_t;
using branch_t = tree_topology_t::branch_t;
using branch_id_t = tree_topology_t::branch_id_t;
using space_vector_t = tree_topology_t::space_vector_t;



#endif // tree_h
