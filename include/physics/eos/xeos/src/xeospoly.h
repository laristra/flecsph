/*~-------------------------------------------------------------------~*
 * Copyright (c) 2018 Los Alamos National Security, LLC
 * All rights reserved.
 *~-------------------------------------------------------------------~*/

/**
 * @file xeospoly.h
 * @author Oleg Korobkin
 * @date 2018-01-05
 * @brief Polytropic equation of state (analytic)
 */
#ifndef XEOS_XEOSPOLY_H_
#define XEOS_XEOSPOLY_H_

#include <math.h>
#include <iostream>
#include "xeosbase.h"

namespace xeos {

/**
 * 1D polytropic equation of state, P(rho) = K rho^Gamma
 */
class XeosPolytropic : public XeosBase {

 public:
  // import constructors from the base class
  using XeosBase::XeosBase;

  // default constructor
  XeosPolytropic()
      : XeosPolytropic(1.0,5./3.,PhysicalQuantity::GetGlobalUnits()) {}

  // workhorse constructor
  // Note: the K constant can be specified in any units;
  //       it is automatically converted to global units
  //       by the constructor.
  XeosPolytropic (const double _K, const double _G,
      const PhUnits _u = PhysicalQuantity::GetGlobalUnits())
      : K(_K, {(float)(3*_G-1), (float)(1-_G),-2,0}, _u), Gamma(_G) {
    K.ConvertTo(PhysicalQuantity::GetGlobalUnits());
  }

  // see explanation in the XeosBase class
  virtual std::vector<Pq> GetPrimaryQuantities() const override;
  virtual std::vector<Pq> GetAdditionalQuantities() const override;

  virtual void
  operator() (const int num_out,
      eos_in in_array[],
      eos_out out_array[]) override;

  virtual void
  DfDx (const int num_out, const Pq dX,
     eos_in  in_array[],
     eos_out dF_array[]) override { /* TODO */ }

  // scalar eos functions
  void operator() (const PqDensity& rho, PqPressure *pres) const;
  void operator() (const PqDensity& rho,
      PqSpecificInternalEnergy *eps) const;
  void operator() (const PqPressure& pres, PqDensity *rho) const;
  void operator() (const PqPressure& pres,
      PqSpecificInternalEnergy *eps) const;
  void operator() (const PqSpecificInternalEnergy& eps,
      PqDensity *rho) const;
  void operator() (const PqSpecificInternalEnergy& eps,
      PqPressure *pres) const;

  // vector eos functions
  void operator() (const NvDensity& rho, NvPressure *pres) const;
  void operator() (const NvPressure& pres, NvDensity *rho) const;

  // getters
  inline double
  GetGamma() const { return Gamma; }

  inline PqOther
  GetK() const { return K; }

  // consistency check declaration
  virtual bool
  ConsistencyCheck() const override;

 private:
  // eos parameters: K and Gamma in P = K*rho**Gamma
  double Gamma = 5./3.;
  PqOther K = {1.0,{4,-2./3.,-2,0},
      PhysicalQuantity::GetGlobalUnits()};

}; // class XeosPolytropic


std::vector<Pq>
XeosPolytropic::GetPrimaryQuantities() const {
  // note how there is neither temperature
  // nor entropy in this eos class
  return std::vector<Pq> {
    Pq::Density,
    Pq::Pressure,
    Pq::SpecificInternalEnergy
  };
}


std::vector<Pq>
XeosPolytropic::GetAdditionalQuantities() const {
  // no additional quantities, sorry
  return std::vector<Pq> {};
}

// compute equation of state at a single point
// TODO: work in progress...
void
XeosPolytropic::operator() (
    const int num_out,  // # of output arguments
    eos_in in_array[],
    eos_out out_array[]) {

  PqOther w(3.5, {1,0,0,0});
  out_array[0] = (&w);

  using namespace std;
  cout << "in[0].kind  = "
       << PqString((*in_array[0]).PqKind()) << endl;
  cout << "out_array[0].kind = "
       << PqString((*out_array[0]).PqKind()) << endl;
}


// overloaded () operators for different types of arguments
// rho --> pressure
void XeosPolytropic::operator() (
    const PqDensity& rho, PqPressure *pres) const {
  *pres = (double)K*pow((double)rho,Gamma);
}

// rho --> eps
void XeosPolytropic::operator() (
    const PqDensity& rho, PqSpecificInternalEnergy *eps) const {
  *eps = (double)K*pow((double)rho,Gamma-1.)/(Gamma-1.);
}

// pressure --> rho
void XeosPolytropic::operator() (
    const PqPressure& pres, PqDensity *rho) const {
  *rho = pow((double)pres/(double)K,1./Gamma);
}

// pressure --> eps (internal energy)
void XeosPolytropic::operator() (
    const PqPressure& pres,
    PqSpecificInternalEnergy *eps) const {
  *eps = pow((double)pres,1.-1./Gamma)
       * pow((double)K,1./Gamma)/(Gamma-1);
}

// eps --> density
void XeosPolytropic::operator() (
    const PqSpecificInternalEnergy& eps,
    PqDensity *rho) const {
  *rho = pow((Gamma-1)*(double)eps/(double)K, 1./(Gamma-1.));
}

// eps --> pressure
void XeosPolytropic::operator() (
    const PqSpecificInternalEnergy& eps,
    PqPressure *pres) const {
  *pres = pow((double)K,-1./(Gamma-1))
        * pow((Gamma-1)*(double)eps, Gamma/(Gamma-1));
}

// density[] -> pressure[]
void XeosPolytropic::operator() (
    const NvDensity& rho,
    NvPressure *pres) const {
  assert (rho.Size() == pres->Size());
  auto it = begin(rho);
  auto jt = begin(*pres);
  for (; it<end(rho); ++it, ++jt)
    (*jt) = (double)K*pow(*it,Gamma);
}

// pressure[] -> density[]
void XeosPolytropic::operator() (
    const NvPressure& pres,
    NvDensity *rho) const {
  assert (pres.Size() == rho->Size());
  auto it = begin(pres);
  auto jt = begin(*rho);
  for (; it<end(pres); ++it, ++jt)
    (*jt) = pow(*it/(double)K,1./Gamma);
}

// define consistency check
bool XeosPolytropic::ConsistencyCheck() const {

  bool retval;
  const double tol = 1e-15;
  PqDensity rho1(1.0), rho2;
  PqPressure P1, P2;
  PqSpecificInternalEnergy eps1, eps2;

  // disallow Gamma = 1 case (isothermal)
  assert (fabs(Gamma-1.) > tol);

  // rho <-> P
  (*this) (rho1, &P1);
  (*this) (P1, &rho2);
  (*this) (rho2, &P2);
  retval = abs((double)rho1 - (double)rho2) < tol;
  retval = retval && abs((double)P1 - (double)P2) < tol;

  // P <-> eps
  (*this) (P1, &eps1);
  (*this) (eps1, &P2);
  (*this) (P2, &eps2);
  retval = retval && (abs((double)P1 - (double)P2) < tol);
  retval = retval && (abs((double)eps1 - (double)eps2) < tol);

  // P <-> eps
  (*this) (rho1, &eps1);
  (*this) (eps1, &rho2);
  (*this) (rho2, &eps2);
  retval = retval && (abs((double)rho1 - (double)rho2) < tol);
  retval = retval && (abs((double)eps1 - (double)eps2) < tol);

  return retval;
}

} // namespace xeos
#endif // XEOS_XEOSPOLY_H_
