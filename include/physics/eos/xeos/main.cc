/*~-------------------------------------------------------------------~*
 * Copyright (c) 2018 Los Alamos National Security, LLC
 * All rights reserved.
 *~-------------------------------------------------------------------~*/

/**
 * @file main.cc
 * @author Oleg Korobkin
 * @date 2018-01-05
 * @brief XEOS library use cases and tests
 */
#include "main.h"
#include <iostream>
#include "xeospoly.h"

static void TestNVector();
static void TestUnits();
static void TestPhysicalNVector();
static void TestXeosPolytropic();

int main() {
  TestNVector();
  TestUnits();
  TestPhysicalNVector();
  TestXeosPolytropic();
  return -0;
}

//-------------------------------------------------------------
static void TestNVector() {

  using namespace std;
  using xeos::NVector;

  cout << "== Testing the class NVector<double> ==" << endl;
  cout << endl;
  cout << "1. creating a 3-dimensional NVector x: " << endl;
  cout << "   NVector<double>  x(12,100,10);" << endl;
  NVector<double>  x(12,100,10);
  const vector<double>::size_type *shape = x.GetShape();
  cout << "   x.num_dims = " << x.GetNumDimensions() << endl;
  cout << "   x.shape = {";
  for (int i=0;i<3;i++)
    cout << (int)shape[i] << " ";
  cout << "} " << endl;
  cout << "   x.Size() = " << x.Size() << endl;
  x[4] = 45;
  cout << "   x[4] = " << x[4] << endl;
  cout << endl;

  cout << "2. copy constructor:" << endl;
  cout << "   NVector<double> y = x;" << endl;
  NVector<double> y = x;
  cout << "   y.num_dims =  " << y.GetNumDimensions() << endl;
  shape = y.GetShape();
  cout << "   y.shape = {";
  for (int i=0;i<3;i++)
    cout << (int)shape[i] << " ";
  cout << "} " << endl;
  cout << "   y.Size() = " << y.Size() << endl;
  y[4] = 45;
  cout << "   y[4] = " << y[4] << endl;
  cout << endl;

  cout << "3. conversion constructors, vector -> NVector:" << endl;
  cout << " - vector<double> v(36, 314.);" << endl;
  cout << "   NVector<double> nv0(v);" << endl;
  vector<double> v(36, 314.);
  NVector<double> nv0(v);
  cout << "   nv0.num_dims =  " << nv0.GetNumDimensions() << endl;
  shape = nv0.GetShape();
  cout << "   nv0.shape = {";
  for (int i=0;i<3;i++)
    cout << (int)shape[i] << " ";
  cout << "} " << endl;
  cout << "   nv0.Size() = " << nv0.Size() << endl;
  cout << "   nv0[4] = " << nv0[4] << endl;
  cout << endl;

  cout << " - NVector<double> nv1(36, v);" << endl;
  NVector<double> nv1(36, v);
  cout << "   nv1.num_dims =  " << nv1.GetNumDimensions() << endl;
  shape = nv1.GetShape();
  cout << "   nv1.shape = {";
  for (int i=0;i<3;i++)
    cout << (int)shape[i] << " ";
  cout << "}" << endl;
  cout << "   nv1.Size() = " << nv1.Size() << endl;
  cout << "   nv1[4] = " << nv1[4] << endl;
  cout << endl;

  cout << " - NVector<double> nv2(9,4, v);" << endl;
  NVector<double> nv2(9,4, v);
  cout << "   nv2.num_dims =  " << nv2.GetNumDimensions() << endl;
  shape = nv2.GetShape();
  cout << "   nv2.shape = {";
  for (int i=0;i<3;i++)
    cout << (int)shape[i] << " ";
  cout << "}" << endl;
  cout << "   nv2.Size() = " << nv2.Size() << endl;
  cout << "   nv2[4] = " << nv2[4] << endl;
  cout << endl;

  cout << " - NVector<double> nv(6,2,3, v);" << endl;
  NVector<double> nv(6,2,3, v);
  cout << "   nv.num_dims =  " << nv.GetNumDimensions() << endl;
  shape = nv.GetShape();
  cout << "   nv.shape = {";
  for (int i=0;i<3;i++)
    cout << (int)shape[i] << " ";
  cout << "} " << endl;
  cout << "   nv.Size() = " << nv.Size() << endl;
  cout << "   nv[4] = " << nv[4] << endl;
  cout << endl;

  cout << "5. initializer-list constructor:" << endl;
  cout << "   NVector<double> z {1,2,3,45,6.7};" << endl;
  NVector<double> z {1,2,3,45,6.7};
  cout << "   z.num_dims =  " << z.GetNumDimensions() << endl;
  shape = z.GetShape();
  cout << "   z.shape = {";
  for (int i=0;i<3;i++)
    cout << (int)shape[i] << " ";
  cout << "} " << endl;
  cout << "   z.Size() = " << z.Size() << endl;
  cout << "   z(4) = " << z(4) << "  # note the () access" << endl;
  cout << endl;

  cout << "6. testing Reshape() function:" << endl;
  cout << " - nv.Reshape(36): " << endl;
  nv.Reshape(36);
  cout << "   dimensions = " << nv.GetNumDimensions() << ", {";
  shape = nv.GetShape();
  for (int i=0;i<3;i++)
    cout << (int)shape[i] << " ";
  cout << "} " << "nv(30) = " << nv(30) << endl;

  cout << " - nv.Reshape(6,6);" << endl;
  nv.Reshape(6,6);
  cout << "   dimensions = " << nv.GetNumDimensions() << ", {";
  shape = nv.GetShape();
  for (int i=0;i<3;i++)
    cout << (int)shape[i] << " ";
  cout << "} " << "nv(5,5) = " << nv(5,5) << endl;

  cout << " - nv.Reshape(1,9,4): " << endl;
  nv.Reshape(1,9,4);
  cout << "   dimensions = " << nv.GetNumDimensions() << ", {";
  shape = nv.GetShape();
  for (int i=0;i<3;i++)
    cout << (int)shape[i] << " ";
  cout << "} " << "nv(0,3,5) = " << nv(0,5,3) << endl;
  cout << endl;

  cout << "7. testing assignment through operator():" << endl;
  cout << "   nv.Reshape(2,3,6); " << endl;
  cout << "   nv = { ";
  nv.Reshape(2,3,6);
  for (int i=0;i<2;i++) {
    for (int j=0;j<3;j++) {
      if(i || j) cout << "          ";
      for (int k=0;k<6;k++) {
        nv(i,j,k) = 100*(i+1) + 10*(j+1) + k+1;
        cout << nv(i,j,k) << " ";
      }
      if (i==1 && j==2)
        cout << "}";
      cout << endl;
    }
    cout << endl;
  }
} // test_nvector

static void TestUnits() {

  using namespace std;
  using namespace xeos;

  PhysicalQuantity::SetGlobalUnits (PhUnits::CGS);
  vector<PhUnits> all_units = {
    PhUnits::CGS, PhUnits::SI, PhUnits::NUCLEAR, PhUnits::GR
  };

  cout << "== Testing physical units conversion ==" << endl;
  cout << endl;

  cout << "1. physical constants in various units:" << endl
       << endl;
  PqOther G(6.67384e-11, {3,-1,-2,0}, PhUnits::SI);
  PhysicalScalar<Pq::Other,float> m_n(939.565379,
    {0,1,0,0}, PhUnits::NUCLEAR);
  PqOther clight(2.99792458e10, {1,0,-1,0});
  PqOther MeV(1.0, {2,1,-2,0}, PhUnits::NUCLEAR);
  PqOther hbar(197.3268718, {2,1,-1,0}, PhUnits::NUCLEAR);
  PqOther Rsun(696392e3, {1,0,0,0}, PhUnits::SI);
  PqOther kB(1.0, {2,1,-2,-1}, PhUnits::NUCLEAR);
  PqDensity rho_sat(2.3e+17,PhUnits::SI);
  PqPressure Patm(1010000.0);

  PhysicalQuantity *pqs[9] =
    { &G, &m_n, &clight, &MeV, &hbar, &Rsun, &kB, &rho_sat, &Patm };

  vector<string> names = {
    "gravitational constant",
    "neutron mass          ",
    "speed of light        ",
    "MeV                   ",
    "hbar                  ",
    "Rsun                  ",
    "Boltzmann's constant  ",
    "saturation density    ",
    "atmospheric pressure  "};

  for (int i=0;i<names.size(); i++) {
    cout << setw(25);
    cout << names[i];
    for (int j=0;j<all_units.size(); j++) {
      pqs[i]->ConvertTo(all_units[j]);
      if (j!=2)
        cout << fixed << " :" << setw(24);
      else
        cout << fixed << " :" << setw(30);
      cout << *pqs[i];
    }
    cout << endl;
  }
  cout << endl;

  cout << "2. polytropic constant K in P = K*rho^gam:" << endl;
  cout << "   float gam = 4./3.;" << endl;
  cout << "   K = 1 [cgs];" << endl;

  float gam = 4./3.;
  PqOther K(1.0, {3*gam-1,1-gam,-2,0});
  cout << "   K units: L^{3*gam-1} M^{1-gam} T^{-2} --> ";
  cout << fixed << std::setprecision(3);
  const float *Kdm = K.GetPhysicalDimensionsArray();
  cout << "{" << Kdm[0] << "," << Kdm[1] << ","<< Kdm[2] << "} ";
  cout << endl;

  cout << "   K.ConvertTo(PhUnits::NUCLEAR) --> ";
  K.ConvertTo(PhUnits::NUCLEAR);
  cout << "K = " << K << endl;
  cout << endl;

} // test_units


static void TestPhysicalNVector() {

  using namespace std;
  using namespace xeos;

  cout << "== Testing PhysicalNVector class ==" << endl;
  cout << "1. creating NvDensity with 12 elements, "
       << " then reshaping to 3x4:" << endl;
  cout << "   NvDensity dens(PhUnits::NUCLEAR);" << endl;
  cout << "   nv.Resize(12); " << endl;
  cout << "   nv.Reshape(3,4); " << endl;
  NvDensity dens(PhUnits::NUCLEAR);
  dens.Resize(12);
  dens.Reshape(3,4);
  for (int i=0;i<3;i++)
    for (int j=0;j<4;j++)
      dens(i,j) = i + 0.25*j;

  cout << fixed << setprecision(3) << "   dens = { ";
  for (int j=0;j<3;j++) {
    if(j) cout << "            ";
    for (int k=0;k<4;k++) {
      cout << dens(j,k) << " ";
    }
    if (j==2)
      cout << "} [nuc] == [MeV/c^2/fm^3]";
    cout << endl;
  }
  cout << endl << endl;

  cout << "2. converting to cgs units: " << endl;
  dens.ConvertTo(PhUnits::CGS);
  cout << scientific << setprecision(3) << "   dens = { ";
  for (int j=0;j<3;j++) {
    if(j) cout << "            ";
    for (int k=0;k<4;k++) {
      cout << dens(j,k) << " ";
    }
    if (j==2)
      cout << "} [g/cm^3]";
    cout << endl;
  }
  cout << endl << endl;

} // test_PhysicalNVector


static void TestXeosPolytropic() {

  using namespace std;
  using namespace xeos;

  cout << "== Polytropic equation of state ==" << endl;
  cout << "1. specifying K in different units" << endl;
  cout << "   PhysicalQuantity::SetGlobalUnits(PhUnits::GR);"
       << endl;

  PhysicalQuantity::SetGlobalUnits(PhUnits::GR);
  cout << "   XeosPolytropic eosPoly(10.,4./3.,PhUnits::CGS);" << endl;
  XeosPolytropic eosPoly(2e+14,4./3.,PhUnits::CGS);
  cout << "   K = 2e+14 [cgs] --> K = ";
  cout << (double)eosPoly.GetK() << endl;
  cout << endl;

  cout << "2. testing eos function: rho -> P" << endl;
  cout << "   PqDensity rho(5.0); // [geom]" << endl;
  cout << "   PqPressure P;       // [geom]" << endl;
  cout << "   eosPoly(rho, &P); --> P = ";
  PqDensity rho(5.0);
  PqPressure P;
  eosPoly(rho, &P);
  cout << P << endl;
  cout << endl;

  cout << "3. convert to cgs, find P->eps->rho, ";
  cout << "convert back to geom:"<< endl;
  cout << "   P.ConvertTo(PhUnits::CGS); --> P = ";
  P.ConvertTo(PhUnits::CGS);
  cout << P << endl;
  cout << "   rho.ConvertTo(PhUnits::CGS);" << endl;
  cout << "   set global units to cgs (see above);" << endl;
  cout << "   XeosPolytropic eosPolyCgs(2e+14, 4./3.);" << endl;
  cout << "   K = 2e+14 [cgs] --> K = ";
  rho.ConvertTo(PhUnits::CGS);
  PhysicalQuantity::SetGlobalUnits(PhUnits::CGS);
  XeosPolytropic eosPolyCgs(2e+14, 4./3.);
  cout << (double)eosPolyCgs.GetK() << endl;
  cout << "   PqSpecificInternalEnergy eps; // [erg/g]" << endl;
  cout << "   eosPolyCgs(P, &eps);   --> eps = ";
  PqSpecificInternalEnergy eps;
  eosPolyCgs(P, &eps);
  cout << eps << endl;
  cout << "   eosPolyCgs(eps, &rho); --> rho = ";
  eosPolyCgs(eps, &rho);
  cout << rho << endl;
  cout << "   rho.ConvertTo(PhUnits::GR); --> rho = ";
  rho.ConvertTo(PhUnits::GR);
  cout << rho << endl;
  cout << "   expected: rho = 5.0 [1/M^2]" << endl;
  cout << endl;

  cout << "4. computing pressure NVector from density" << endl;
  cout << "   NvDensity dens();" << endl;
  cout << "   dens.Resize(3,4);" << endl;
  NvDensity dens;
  dens.Resize(3,4);
  for (int i=0;i<3;i++)
    for (int j=0;j<4;j++)
      dens(i,j) = i + 0.25*j;

  cout << fixed << setprecision(3) << "   dens = { ";
  for (int j=0;j<3;j++) {
    if(j) cout << "            ";
    for (int k=0;k<4;k++) {
      cout << dens(j,k) << " ";
    }
    if (j==2)
      cout << "} [g/cm^3]";
    cout << endl;
  }
  cout << endl;

  cout << "   NvPressure pres;" << endl;
  cout << "   pres.Resize(3,4);" << endl;
  cout << "   eosPolyCgs(dens, &pres);" << endl;
  NvPressure pres;
  pres.Resize(3,4);
  eosPolyCgs(dens, &pres);

  cout << scientific << setprecision(3) << "   pres = { ";
  for (int j=0;j<3;j++) {
    if(j) cout << "            ";
    for (int k=0;k<4;k++) {
      cout << pres(j,k) << " ";
    }
    if (j==2)
      cout << "} [Ba]";
    cout << endl;
  }
  cout << endl;

  cout << "5. computing density vector from pressure" << endl;
  cout << "   eosPolyCgs(pres, &dens);" << endl;
  eosPolyCgs(pres, &dens);
  cout << fixed << setprecision(3) << "   dens = { ";
  for (int j=0;j<3;j++) {
    if(j) cout << "            ";
    for (int k=0;k<4;k++) {
      cout << dens(j,k) << " ";
    }
    if (j==2)
      cout << "} [g/cm^3]";
    cout << endl;
  }
  cout << "   - must be the same as in (4) above." << endl;
  cout << endl;

  cout << "6. eos consistency check" << endl;
  bool passed = eosPolyCgs.ConsistencyCheck();
  cout << "   eosPolyCgs.ConsistencyCheck() --> "
       << (passed?"true":"false") << " // "
       << (passed?"passed":"failed")
       << endl;
  cout << endl;


/*
  cout << endl << "== Testing operator() in xeos ==" << endl;
  eos_out qarg[1], zarg[2];
  zarg[1] = &P;
  qarg[0] = &dens;
  eos_poly(1, qarg, zarg);
  cout << "(*zarg[0]).kind = "
       << pq_string((*zarg[0]).pq_kind()) << endl;
  cout << "(*zarg[1]).kind = "
       << pq_string((*zarg[1]).pq_kind()) << endl;

  cout << endl << "== XeosPolytropic consistency check ==" << endl;
  cout << "consistency check passed: "
       << (eos_poly.consistency_check()?"true":"false") << endl;
 */
}

