/*~-------------------------------------------------------------------~*
 * Copyright (c) 2018 Los Alamos National Security, LLC
 * All rights reserved.
 *~-------------------------------------------------------------------~*/

/**
 * @file physquan.h
 * @author Oleg Korobkin
 * @date 2018-01-05
 * @brief Template classes to declare scalar and vector physical
 *        quantities in the form of distinct C++ types
 *
 * PhysicalQuantity: base class with units and dimensions
 * A   A
 * |   +---- PhysicalScalar<PhysKind,Type>
 * |         =
 * |         = PqDensity, PqPressure,
 * |           PqSpecificInternalEnergy etc.
 * |
 * +-------- PhysicalNVector<PhysKind,Type>
 *           =
 *           = NvDensity, NvEnthalph ...
 *
 * File includes two *.inc files:
 *  - pq_enum.inc: contains declaration of PhUnits and Pq enums;
 *  - quantities.inc: typedefs and dimensions of standard Pq's
 *
 */
#ifndef XEOS_PHYSUNITS_H_
#define XEOS_PHYSUNITS_H_

#include <string>
#include <iomanip>
#include <math.h>
#include "nvector.h"
#include <iostream>

namespace xeos {

//---------------------//
#include "pq_enum.inc" // defines enums PhUnits and Pq
//---------------------//

// for a given Pq kind, fills array with physical dimensions
// e.g. Pq::Energy -> {2, 1, -2, 0} [Length, Mass, Time, Temperature]
void
FillPhysicalDimensions (const Pq q, float retval[]);

/**
 * Base class to smoothly handle units conversions
 * TODO: how is it better than what is in Boost??
 */
class PhysicalQuantity {
 public:
  virtual ~PhysicalQuantity() {}

  // -- static members

  static PhUnits GetGlobalUnits() { return global_units; }
  static void    SetGlobalUnits(const PhUnits _u) { global_units = _u; }

  // -- non-static members
  virtual Pq     PqKind() const { return Pq::Undefined; }
  PhUnits        GetUnits() const { return val_units; }
  const float*   GetPhysicalDimensionsArray() const { return dim; }
  virtual void   ConvertTo(const PhUnits _u) { val_units=_u; }
  virtual std::string ToString() const { return "<PhysicalQuantity>"; }

  // friendly ostream operator
  friend std::ostream&
  operator<< (std::ostream &os, const PhysicalQuantity &q) {
    return os << q.ToString();
  }

 protected:
  // dimensions: length, mass, time, temperature (cgsK)
  static const int num_basic_units = 4;
  float dim[num_basic_units] = {0,0,0,0}; // default: dimensionless

  // default units: none
  PhUnits val_units = PhUnits::NONE;

 private:
  // static members
  static PhUnits global_units;

}; // class PhysicalQuantity


/**
 * Template class for scalar physical quantity types
 * E.g.:
 *   typedef PhysicalScalar<Pq::Temperature,double> PqTemperature;
 *   PqTemperature temp(0.3, PhUnits::NUCLEAR);
 *   temp.ConvertTo(PhUnits::SI);
 *   std::cout << (double)temp << std::endl;
 */
template<Pq K, typename S>
class PhysicalScalar : public PhysicalQuantity {

 public:
  // define nested type 'value_type'
  typedef S value_type;

  // -- constructors ---------------------------------------------------
  PhysicalScalar<K,S> ()
      : PhysicalScalar<K,S> ((S)0, global_units) {}

  PhysicalScalar<K,S> (const S _v)
      : PhysicalScalar<K,S> (_v, global_units) {}

  // --= workhorse constructor =--
  PhysicalScalar<K,S> (const S _v, const PhUnits _u) : val(_v) {
    FillPhysicalDimensions (K, dim);
    PhysicalQuantity::val_units = _u;

    // IMPORTANT!
    // if global units are not set, sets them to '_u'
    if (global_units == PhUnits::NONE)
      global_units = _u;
  }

  // maximal constructor for 'other' type, with dimensions;
  // 'other' type is useful for defining various physical constants
  // which often have wierd dimensions, e.g. [G] = L^3 M^-1 T^-2 etc.
  PhysicalScalar<K,S> ( const S _v,
      const std::initializer_list<float> _dim, const PhUnits _u)
      : PhysicalScalar<K,S> (_v, _u) {

    // allow this constructor only for quantities of type 'other'
    assert (K == Pq::Other);

    int i=0;
    for (auto it=begin(_dim);it<end(_dim);++it,++i)
      dim[i]= (*it);
  }

  // minimal constructor for 'other' type, w/o units
  PhysicalScalar<K,S> (const S _v,
      const std::initializer_list<float> _dim)
      : PhysicalScalar<K,S> (_v, _dim, global_units) {}

  // copy constructor
  PhysicalScalar<K,S> (const PhysicalScalar<K,S> &q)
      : PhysicalScalar<K,S> (q.val, q.GetUnits()) {}

  // destructor
  virtual ~PhysicalScalar<K,S>() override {}

  // -- operators ------------------------------------------------------

  // assignment operator: <- same type (PhysicalScalar<K,S>)
  // ABORT if units are different
  PhysicalScalar<K,S>&
  operator= (const PhysicalScalar<K,S> &q) {
    assert(GetUnits()==q.GetUnits());
    val = q.val;
  }

  // assignment operator: <- double: assigns only the value
  PhysicalScalar<K,S>&
  operator= (const S _v) { val = _v; }

  // explicit conversion to scalar type
  explicit
  operator S() const { return val; }

  // kindly return own kind
  virtual Pq
  PqKind() const override { return K; }

  // ---------------------------------------------------
  // other non-static members

  virtual void
  ConvertTo (const PhUnits new_units) override;

  // string representation
  virtual std::string
  ToString() const {
    std::stringstream buffer;
    buffer << std::scientific << std::setprecision(8);
    buffer << val << " [" << PqUnitsString(K,GetUnits()) << "]";
    return buffer.str();
  }

 private:
  S val;

}; // class PhysicalScalar


/**
 * Template class for vector physical quantities
 *  - K: physical quantity kind
 *  - S: value type
 * Example:
 *   typedef PhysicalScalar<Pq::Pressure,double> NvPressure;
 *   PhysicalQuantity::SetGlobalUnits(PhUnits::CGS);
 *   NvPressure pvec = {0.3, 1.4, 5.3, 6.12}; // [Ba]
 *   pvec.Reshape(2,2);
 *   pvec.ConvertTo(PhUnits::SI);
 *   std::cout << (double)pvec(1,1) << std::endl;
 *
 */
template<Pq K, class S>
class PhysicalNVector : public PhysicalQuantity, public NVector<S> {

  typedef NVector<S> V;

 public:
  // provide nested types to expose structure
  typedef typename V::size_type   size_type;
  typedef                     S   value_type;

  // default constructor
  PhysicalNVector<K,S> ()
      : PhysicalNVector<K,S> (PhysicalQuantity::GetGlobalUnits()) {}

  // constructor by units
  PhysicalNVector<K,S> (const PhUnits _u) : V() {
    PhysicalQuantity::val_units = _u;
    FillPhysicalDimensions(K,dim);
  }

  PhysicalNVector<K,S> (const V &_v)
      : PhysicalNVector<K,S> (_v, PhysicalQuantity::GetGlobalUnits()) {
  }

  // --= workhorse vector-lvalue copy constructor =--
  PhysicalNVector<K,S> (const V &_v,
      const PhUnits _u = PhysicalQuantity::GetGlobalUnits())
      : V(_v) {
    PhysicalQuantity::val_units = _u;
    FillPhysicalDimensions(K,dim);
  }

  // --= workhorse vector-rvalue copy constructor =--
  PhysicalNVector<K,S> (const V &&_v,
      const PhUnits _u = PhysicalQuantity::GetGlobalUnits())
      : V(_v) {
    PhysicalQuantity::val_units = _u;
    FillPhysicalDimensions(K,dim);
  }

  // initializer-list constructor
  PhysicalNVector<K,S> (std::initializer_list<S> il) : V(il) {
    PhysicalQuantity::val_units = PhysicalQuantity::GetGlobalUnits();
    FillPhysicalDimensions(K,dim);
  }

  // copy constructor / assignment operator
  PhysicalNVector<K,S>&
  operator= (const PhysicalNVector<K,S> &q) { /* TODO */ }

  // virtual destructor
  virtual ~PhysicalNVector<K,S> () {}

  // make sure to use all these functions from NVector
  using V::PushBack;
  using V::Begin;
  using V::End;
  using V::Size;
  using V::Resize;
  using V::operator[];
  using V::operator();

  friend auto
  begin(PhysicalNVector<K,S>& p) { return p.Begin(); }

  friend auto
  end(PhysicalNVector<K,S>& p) { return p.End(); }

  // return the kind of the physical quantity (Pq)
  virtual Pq
  PqKind() const override { return K; }

  // convert to a different system of units
  virtual void
  ConvertTo (const PhUnits new_units) override;

protected:
  using PhysicalQuantity::num_basic_units;
  using PhysicalQuantity::dim;

}; // class PhysicalNVector


//
// implementations (former physunits.cpp)
//
PhUnits PhysicalQuantity::global_units = PhUnits::NONE;

// convert to a different system of units
template <Pq K, class S>
void PhysicalScalar<K,S>::ConvertTo (const PhUnits new_units) {

  PhUnits old_units = GetUnits();
  if (new_units == old_units) // sanity check
    return;

  // conversion
  double cfactor_old = CgsConversionFactor(old_units, dim);
  double cfactor_new = CgsConversionFactor(new_units, dim);

  // compute the new value
  val *= cfactor_old / cfactor_new;

  // set new units
  PhysicalQuantity::val_units = new_units;
}

// PhysicalNVector: implementation of the units conversion method
template <Pq K, class V>
void PhysicalNVector<K,V>::ConvertTo (const PhUnits new_units) {

  PhUnits old_units = PhysicalQuantity::GetUnits();
  if (new_units == old_units) // sanity check
    return;

  // conversion
  double cfactor_old = CgsConversionFactor (old_units,
      PhysicalQuantity::dim);
  double cfactor_new = CgsConversionFactor (new_units,
      PhysicalQuantity::dim);

  // compute the new value
  for (int i=0; i<this-> Size(); i++)
    (*this)[i] *= cfactor_old / cfactor_new;

  // set new units
  PhysicalQuantity::val_units = new_units;
}

// ---------------------- //
#include "quantities.inc" // physical quantities database
// ---------------------- //

} // namespace xeos
#endif // XEOS_PHYSUNITS_H_

#if 0
#include <iostream>
int main() {
  using namespace xeos;
  typedef PhysicalNVector<Pq::Pressure,double> xPressure;
  PhysicalQuantity::SetGlobalUnits(PhUnits::CGS);
  xPressure pvec = {0.3, 1.4, 5.3, 6.12}; // [Ba]
  pvec.Reshape(2,2);
  pvec.ConvertTo(PhUnits::SI);
  std::cout << (double)pvec(1,1) << " "
            << UnitsToString(pvec.GetUnits()) << std::endl;
}
#endif
