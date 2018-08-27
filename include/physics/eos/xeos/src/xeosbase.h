/*~-------------------------------------------------------------------~*
 * Copyright (c) 2018 Los Alamos National Security, LLC
 * All rights reserved.
 *~-------------------------------------------------------------------~*/

/**
 * @file xeosbase.h
 * @author Oleg Korobkin
 * @date 2018-01-05
 * @brief Abstract base classes for arbirtrary equations of state
 */
#ifndef XEOS_XEOSBASE_H_
#define XEOS_XEOSBASE_H_

#include "physquan.h"

namespace xeos {

// these types will be used in a universal eos function call (TODO)
typedef const PhysicalQuantity * const eos_in;
typedef PhysicalQuantity*             eos_out;

/** 
 * Abstract base class for all equations of state
 */
class XeosBase {

 public:
  // default constructor: does nothing tbh
  XeosBase() {}

  // get primary quantities which this EoS understands
  virtual std::vector<Pq>
  GetPrimaryQuantities() const = 0;

  // get additional quantities
  virtual std::vector<Pq>
  GetAdditionalQuantities() const = 0;

  // optional feature
  virtual bool
  ConsistencyCheck() const = 0;

  // universal multi-parameter eos function
  virtual void
  operator() (const int num_out, // number of output arguments
       eos_in  in_array[],       // input quantities
       eos_out out_array[]) = 0; // output quantities

  // universal multi-parameter eos first derivative(s):
  // computes one or more derivatives wrt X at a given point
  //
  // {dF/dX(a1,a2,...), dG/dX(a1,a2,...), ... }
  //
  virtual void
  DfDx (const int num_out,    // number of output arguments
     const Pq dX,             // differentiation variable X
     eos_in  in_array[],      // input parameters: {a1, a2, ...}
     eos_out dF_array[]) = 0; // array of {dFdX, dGdX, ...}

}; // class XeosBase


/**
 * Abstract base class for tabulated equations of state
 */
class XeosTabulated : public XeosBase {

 public:
  // Demand from every tabulated eos class an ability to
  // read eos table (or tables, if there are several).
  virtual void ReadEosTables() = 0;

  // This function returns a 1D grid along one of the primary
  // variables which parameterize the table. This can be either
  // a uniformly-spaced grid, or a custom grid like this:
  // T = {0.1, 0.2, 0.5, 1.0, 2.0, 10} [GK]
  virtual double
  GridPoint(const Pq var, const int i) const = 0;

  // Total number of points in the table
  // along the primary variable direction
  virtual const int
  GridSize(const Pq) const = 0;

  double EosTableMinval(const Pq var) const {
    return GridPoint(var, 0);
  }

  double EosTableMaxval(const Pq var) const {
    return GridPoint(var, GridSize(var)-1);
  }

 protected:
  std::string data_path; // where the stuff is

}; // class XeosTabulated

} // namespace xeos
#endif // XEOS_XEOSBASE_H_

