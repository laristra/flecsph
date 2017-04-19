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

using entity_id_t = flecsi::topology::entity_id_t;


class entity_key_t
{
  public:
    using int_t = uint64_t;
    static const size_t dimension = 3;
    static constexpr size_t bits = sizeof(int_t) * 8;
    static constexpr size_t max_depth = (bits - 1)/dimension;

    entity_key_t()
    : id_(0)
    {}

    entity_key_t(
      const std::array<point_t, 2>& range,
      const point_t& p)
    : id_(int_t(1) << max_depth * dimension + (bits - 1) % dimension)
    {
      std::array<int_t, dimension> coords;
      for(size_t i = 0; i < dimension; ++i)
      {
        double min = range[0][i];
        double scale = range[1][i] - min;
        coords[i] = (p[i] - min)/scale * (int_t(1) << (bits - 1)/dimension);
      }
      size_t k = 0;
      for(size_t i = 0; i < max_depth; ++i)
      {
        for(size_t j = 0; j < dimension; ++j)
        {
          int_t bit = (coords[j] & int_t(1) << i) >> i;
          id_ |= bit << (k * dimension + j);
        }
        ++k;
      }
    }

  constexpr entity_key_t(const entity_key_t& bid)
  : id_(bid.id_)
  {}

  static 
  constexpr
  entity_key_t
  null()
  {
    return entity_key_t(0);
  }

  constexpr
  bool 
  is_null() const
  {
    return id_ == int_t(0);
  }

  entity_key_t&
  operator=(
      const entity_key_t& ek
  ) 
  {
    id_ = ek.id_; 
    return *this;
  }

  constexpr
  bool 
  operator==(
      const entity_key_t& ek
  ) const 
  {
    return id_ == ek.id_; 
  }

  constexpr 
  bool 
  operator!=(
    const entity_key_t& ek
  ) const
  {
    return id_ != ek.id_; 
  }

  void
  output_(
    std::ostream& ostr
  ) const
  {
    ostr << std::oct << id_ << std::dec;
    //constexpr int_t mask = ((int_t(1) << dimension) - 1) << bits - dimension;
    //size_t d = max_depth;
    //int_t id = id_;
    //
    //while((id & mask) == int_t(0))
    //{
    //  --d;
    //  id <<= dimension;
    //}
    //if(d == 0)
    //{
    //  ostr << "<root>";
    //  return;
    //}
    //id <<= 1 + (bits - 1) % dimension;
    //for(size_t i = 1; i <= d; ++i)
    //{
    //  int_t val = (id & mask) >> (bits - dimension);
    //  ostr << std::oct << val << std::dec; 
    //  //ostr << i << ":" << std::bitset<dimension>(val) << " ";
    //  id <<= dimension;
    //}
  }

  bool
  operator<(
    const entity_key_t& bid
  ) const
  {
    return id_ < bid.id_;
  }

private:
  int_t id_;
  constexpr
  entity_key_t(
    int_t id
  )
  :id_(id)
  {}
};


#endif // tree_h
