/*~-------------------------------------------------------------------~*
 * Copyright (c) 2018 Los Alamos National Security, LLC
 * All rights reserved.
 *~-------------------------------------------------------------------~*/

/**
 * @file pq_enum.inc
 * @author Oleg Korobkin
 * @date 2018-01-05
 * @brief Contains enum PhUnits:: representing different physical
 *        units, conversion factor calculator for each units system
 *        to CGS(K), and enum class Pq:: of all physical quantities
 */

// Enumerator of physical units
enum class PhUnits {
  NONE,
  CGS,
  SI,
  NUCLEAR,
  GR
};


// Enumerator of physical quantities
enum class Pq {
  Undefined,
  Other,
  Density,
  Pressure,
  SpecificInternalEnergy,
  Temperature,
  Enthalpy,
  RestMassDensity,
  BaryonNumberDensity,
  ElectronFraction,
  ChemPotentialDiff,
  SoundSpeedSq,
  SpecificEntropy,
  NeutronMassFraction,
  ProtonMassFraction,
  DeuteronMassFraction,
  TritonMassFraction,
  HelionMassFraction,
  AlphaMassFraction,
  HeavyNucMassFraction,
  AverageAtomicMass,
  AverageNuclearCharge,
  Velocity,
  Mass
};

// conversion factor with respect to cgs
double CgsConversionFactor (const PhUnits u, const float d[])
{
  double retval = 0.0;
  switch (u) {
    case (PhUnits::CGS):
      retval = 1.0;
      break;

    case (PhUnits::SI):
      // [L] = 100 cm, [M] = 1000 g
      retval = pow(10., 2*d[0] + 3*d[1]);
      break;

    case (PhUnits::NUCLEAR):
      retval = 1.0               // fm      [1e-13 cm]
          * pow(1.782662,d[1])   // MeV/c^2 [1e-33 g]
          * pow(3.33564095,d[2]) // fm/c    [1e-24 s]
          * pow(1.1604505,d[3])  // MeV     [1e10 K]
          * pow(10.,-13*d[0]-27*d[1]-24*d[2]+10*d[3]);
      break;

    case (PhUnits::GR):
      retval =
          pow(1.47703,d[0])    // GMsun/c^2 [1e5 cm]
        * pow(1.9891, d[1])    // Msun      [1e33 g]
        * pow(4.92686,d[2])    // GMsun/c^3 [1e-6 s]
        * pow(10.,5*d[0]+33*d[1]-6*d[2]);
      break;

    default:
      assert(false); // such units not implemented
  }
  return retval;
}


// converts PhUnits into std::string
std::string
UnitsToString(const PhUnits u) {
  std::string str;
  switch(u) {
    case PhUnits::CGS:
      str = "cgsK";
      break;

    case PhUnits::SI:
      str = "SI";
      break;

    case PhUnits::NUCLEAR:
      str = "nuc";
      break;

    case PhUnits::GR:
      str = "geometric";
      break;

    default:
      str = "<STREAM UNDEFINED FOR THESE UNITS>";
  }
  return str;
}

