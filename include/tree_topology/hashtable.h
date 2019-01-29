#ifndef _hashtable_
#define _hashtable_

#include <vector>

/**
* @brief Class for hashtable
*/
template<
  typename KEY,
  typename TYPE>
class hashtable
{

public:

  /**
  * @brief Iterator on hashtable
  */
  typedef typename std::vector<std::pair<KEY,TYPE>>::iterator
    iterator;

  /**
  * @brief Creation of the hashtable: no collisions
  */
  hashtable(){
    collision_ = 0;
    ht_.resize(hash_size_);
    assert(ht_[0].first == 0);
  }

  ~hashtable(){
    ht_.clear();
  }

  /**
  * @brief Find a key in the hash table
  * Keys can be listed if conflict
  */
  typename std::vector<std::pair<KEY,TYPE>>::iterator
  find(const KEY& k){
    unsigned int index = hash_(k);
    auto it = ht_.begin() + index;
    while(it->first != KEY{} && it->first != k){
      ++it;
    }
    if(it->first != k) return ht_.end();
    return it;
  }

  /**
  * @brief Emplace an object in the hashtable, if conflict: advance
  */
  template<
    typename... ARGS>
  void
  emplace(const KEY& k, ARGS&&... args){
    unsigned int index = hash_(k);
    // Find an empty spot
    while(ht_[index].first != KEY{}){
      ++index;
      ++collision_;
    }
    ht_[index] = std::make_pair(k,std::forward<ARGS>(args)...);
    ++nelement_;
  }

  size_t collision(){return collision_;}

  void
  clear(){
    ht_.clear();
    ht_.resize(hash_size_);
    collision_ = 0;
    nelement_ = 0;
  }

  size_t
  size(){
    return nelement_;
  }

  iterator
  begin()
  {
    return ht_.begin();
  }

  iterator
  end()
  {
    return ht_.end();
  }

private:

  unsigned int
  hash_(const KEY& k){
    return k & hash_mask_;
  }

  const unsigned int hash_bit_ = 22;
  const size_t hash_size_ = 1<<hash_bit_;
  const unsigned int hash_mask_ = (1<<hash_bit_)-1;
  std::vector<std::pair<KEY,TYPE>> ht_;

  size_t collision_;
  size_t nelement_;
};

#endif // _hashtable_