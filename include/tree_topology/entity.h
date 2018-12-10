#ifndef _mpisph_entity_h_
#define _mpisph_entity_h_

#include "tree_entity_id.h"

namespace flecsi{
namespace topology{

template<
  typename TT,
  typename T,
  size_t D>
class entity
{
  using element_t = TT;
  static const size_t dimension = D;
  using point_t = point__<element_t,dimension>;
  using key_t = key_id__<T,D>;

public:
  entity(){};

  entity
  (
    const point_t& coordinates,
    const element_t& mass,
    const entity_id_t& id
  ): coordinates_(coordinates), mass_(mass), id_(id)
  {};

  inline bool operator==(const entity& a)
  {
    return a.id_ == this->id_;
  }

  ~entity(){};

  point_t coordinates() const {return coordinates_;};
  element_t mass() const {return mass_;};
  key_t key() const {return key_;};
  element_t radius() const {radius_;};
  entity_id_t id() const {return id_;};

  void set_coordinates(const point_t& coordinates){coordinates_ = coordinates;};
  void set_mass(const element_t& mass){mass_ = mass;};
  void set_radius(const element_t& radius){radius_ = radius;};
  void set_key(const key_t& key){key_ = key;};
  void set_key(key_t& key){key_ = key;};
  void set_id(const entity_id_t& id){id_ = id;};

  friend std::ostream& operator<<(std::ostream& os, const entity& b){
    // TODO change regarding to dimension
    os << std::setprecision(10);
    os << "Particle: coord: " <<b.coordinates_;
    os << " mass: " <<b.mass_;
    os << " h: " << b.radius_;
    os << " id: " << b.id_;
    os << "key: "<< b.key_;
    return os;
  }

protected:
  point_t coordinates_;
  element_t mass_;
  entity_id_t id_;
  element_t radius_;
  key_t key_;
}; // class entity

} // namespace topology
} // namespace flecsi

#endif // _mpisph_entity_h_
