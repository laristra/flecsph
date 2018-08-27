/*~-------------------------------------------------------------------~*
 * Copyright (c) 2018 Los Alamos National Security, LLC
 * All rights reserved.
 *~-------------------------------------------------------------------~*/

/**
 * @file nvector.h
 * @author Oleg Korobkin
 * @date 2018-01-05
 * @brief Handy vector class adapted for storing EoS tables
 */
#ifndef XEOS_NVECTOR_H_
#define XEOS_NVECTOR_H_

#include <vector>
#include <assert.h>

namespace xeos {

/*
 * Vector with multidimensional view: e.g., A(i,j,k) = 3.14;
 * Example:
 *   NVector<float> A(150);
 *   A.Reshape(10,15);
 *   A(8,12) = 230;
 *
 */
template<class S>
class NVector {
  typedef typename std::vector<S>   V;

 public:
  // provide nested types to expose structure
  typedef typename std::vector<S>::size_type  size_type;
  typedef S  value_type;

  // default constructor: 1D
  NVector<S> () {}

  // constructors of 1D, 2D and 3D zero-filled NVectors
  explicit // conversion from integer type
  NVector<S> (size_type Nx) : my_vector(Nx,0) { shape[0] = Nx; }

  NVector<S> (size_type Nx, size_type Ny) : my_vector(Nx*Ny,0) {
    num_dims = 2;
    shape[0] = Nx;
    shape[1] = Ny;
  }

  NVector<S> (size_type Nx, size_type Ny, size_type Nz)
      : my_vector(Nx*Ny*Nz,0) {
    num_dims = 3; shape[0] = Nx; shape[1] = Ny; shape[2] = Nz;
  }

  // conversion constructors in 1D, 2D and 3D (from vector)
  explicit // conversion from vector
  NVector<S> (V &_v) : my_vector(_v) {
    shape[0] = _v.size();
  }

  NVector<S> (size_type Nx, V &_v) : my_vector(_v) {
    assert(_v.size() == Nx);
    shape[0] = Nx;
  }

  NVector<S> (size_type Nx, size_type Ny, V &_v) : my_vector(_v) {
    assert(_v.size() == Nx*Ny);
    num_dims = 2;
    shape[0] = Nx;
    shape[1] = Ny;
  }

  NVector<S> (size_type Nx, size_type Ny, size_type Nz, V &_v)
      : my_vector(_v) {
    assert(_v.size() == Nx*Ny*Nz);
    num_dims = 3;
    shape[0] = Nx;
    shape[1] = Ny;
    shape[2] = Nz;
  }

  // copy constructor
  NVector<S> (const NVector<S> &_nv) : my_vector(_nv.my_vector) {
    num_dims = _nv.num_dims;
    for(int i=0;i<num_dims_max;i++)
      shape[i] = _nv.shape[i];
  }

  // move constructor
  NVector<S> (NVector<S> &&_nv) : my_vector(_nv.my_vector) {
    num_dims = _nv.num_dims;
    for(int i=0;i<num_dims_max;i++)
      shape[i] = _nv.shape[i];
  }

  // initializer-list constructor
  NVector<S> (std::initializer_list<S> il) : my_vector(il) {
    shape[0] = my_vector.size();
  }

  // destructor
  virtual ~NVector<S> () {}

  // useful stuff from the vector class:
  // begin/end, size(), push_back, operator[]
  auto Begin() { return begin(my_vector); }
  auto End()   { return end(my_vector); }

  friend auto begin(const NVector<S>& nv) { return begin(nv.my_vector); }
  friend auto begin(NVector<S>& nv)       { return begin(nv.my_vector); }

  friend auto end(const NVector<S>& nv) { return end(nv.my_vector); }
  friend auto end(NVector<S>& nv)       { return end(nv.my_vector); }

  size_type
  Size() const { return my_vector.size(); }

  void
  PushBack (const S elem) {
    assert(num_dims == 1); // forbid pushing if dimensions > 1
    my_vector.push_back(elem);
    ++shape[0];
  }

  S& operator[](const size_type i) { return my_vector.at(i); }

  // return shape
  const size_type*
  GetShape() const { return shape; }

  // return dimensions
  unsigned short
  GetNumDimensions() const { return num_dims; }

  // array elements access operator
  S& operator()(const size_type i) {
    assert (num_dims == 1);
    assert (i>=0 && i<shape[0]);
    return my_vector.at(i);
  }

  S& operator()(const size_type i, const size_type j) {
    assert (num_dims == 2);
    assert (i>=0 && i<shape[0]);
    assert (j>=0 && j<shape[1]);
    return my_vector.at(shape[1]*i+j);
  }

  S& operator()(const size_type i, const size_type j,
    const size_type k) {
    assert (num_dims == 3);
    assert (i>=0 && i<shape[0]);
    assert (j>=0 && j<shape[1]);
    assert (k>=0 && k<shape[2]);
    return my_vector.at(shape[2]*(shape[1]*i+j)+k);
  }

  // array elements access operator
  S operator()(const size_type i) const {
    assert (num_dims == 1);
    assert (i>=0 && i<shape[0]);
    return my_vector.at(i);
  }

  S operator()(const size_type i, const size_type j) const {
    assert (num_dims == 2);
    assert (i>=0 && i<shape[0]);
    assert (j>=0 && j<shape[1]);
    return my_vector.at(shape[1]*i+j);
  }

  S operator()(const size_type i, const size_type j,
    const size_type k) const {
    assert (num_dims == 3);
    assert (i>=0 && i<shape[0]);
    assert (j>=0 && j<shape[1]);
    assert (k>=0 && k<shape[2]);
    return my_vector.at(shape[2]*(shape[1]*i+j)+k);
  }

  // reshape operation
  void Reshape(const size_type Nx) {
    assert ( shape[0]*shape[1]*shape[2] == Nx);
    num_dims = 1;
    shape[0] = Nx;
    shape[1] = 1;
    shape[2] = 1;
  }

  void Reshape(const size_type Nx, const size_type Ny) {
    assert ( shape[0]*shape[1]*shape[2] == Nx*Ny);
    num_dims = 2;
    shape[0] = Nx;
    shape[1] = Ny;
    shape[2] = 1;
  }

  void Reshape(const size_type Nx, const size_type Ny,
      const size_type Nz) {
    assert ( shape[0]*shape[1]*shape[2] == Nx*Ny*Nz);
    num_dims = 3;
    shape[0] = Nx;
    shape[1] = Ny;
    shape[2] = Nz;
  }

  // resize
  void Resize(const size_type Nx) {
    num_dims = 1;
    my_vector.resize(Nx);
    shape[0] = Nx;
    shape[1] = 1;
    shape[2] = 1;
  }

  void Resize(const size_type Nx, const size_type Ny) {
    num_dims = 2;
    my_vector.resize(Nx*Ny);
    shape[0] = Nx;
    shape[1] = Ny;
    shape[2] = 1;
  }

  void Resize(const size_type Nx, const size_type Ny,
      const size_type Nz) {
    num_dims = 3;
    my_vector.resize(Nx*Ny*Nz);
    shape[0] = Nx;
    shape[1] = Ny;
    shape[2] = Nz;
  }

 private:
  unsigned short num_dims = 1; // one-dimensional by default
  static const unsigned short num_dims_max = 3;
  typename std::vector<S>::size_type shape[num_dims_max] = {1,1,1};

  V my_vector;

};

} //namespace xeos

#endif // XEOS_NVECTOR_H_

#if 0
#include <iostream>
int main() {
  using namespace xeos;
  using namespace std;
  NVector<float> A{1,2,3,4,5,6};
  cout << A.Size() << endl;
  A.Reshape(2,3);
  A(1,2) = 230;
  cout << A(0,2) << endl;
}
#endif


