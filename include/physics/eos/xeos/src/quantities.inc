/*~-------------------------------------------------------------------~*
 * Copyright (c) 2018 Los Alamos National Security, LLC
 * All rights reserved.
 *~-------------------------------------------------------------------~*/

/**
 * @file main.cc
 * @author Oleg Korobkin
 * @date 2018-01-05
 * @brief This file is included into physquan.h: it contains
 *        shortcut typedefs and database of dimensions / output
 *        formats for known physical quantities. This is the
 *        file you want to edit in order to add a new physical
 *        quantity.
 *        Use :set syntax=cpp for syntax highlighting in vim.
 */

// --------------------------------------------------------------
// scalar and vector typedefs macro:
//   PQ_TYPEDEF(density) -> PqDensity, NvDensity
#define PQ_TYPEDEF(Q) \
    typedef PhysicalScalar<Pq::Q,double> Pq ## Q; \
    typedef PhysicalNVector<Pq::Q,double> Nv ## Q;

PQ_TYPEDEF(Other)
PQ_TYPEDEF(Density)
PQ_TYPEDEF(Pressure)
PQ_TYPEDEF(Temperature)
PQ_TYPEDEF(Enthalpy)
PQ_TYPEDEF(SpecificInternalEnergy)
PQ_TYPEDEF(RestMassDensity)
PQ_TYPEDEF(BaryonNumberDensity)
PQ_TYPEDEF(ElectronFraction)
PQ_TYPEDEF(ChemPotentialDiff)
PQ_TYPEDEF(SoundSpeedSq)
PQ_TYPEDEF(SpecificEntropy)
PQ_TYPEDEF(NeutronMassFraction)
PQ_TYPEDEF(ProtonMassFraction)
PQ_TYPEDEF(DeuteronMassFraction)
PQ_TYPEDEF(TritonMassFraction)
PQ_TYPEDEF(HelionMassFraction)
PQ_TYPEDEF(AlphaMassFraction)
PQ_TYPEDEF(HeavyNucMassFraction)
PQ_TYPEDEF(AverageAtomicMass)
PQ_TYPEDEF(AverageNuclearCharge)
PQ_TYPEDEF(Velocity)
PQ_TYPEDEF(Mass)

// return string identifying physical quantity
const std::string PqString(const Pq kind) {
  std::string s;
  switch(kind) {
    case Pq::AlphaMassFraction:
      s = "mass fraction of alphas";
      break;
    case Pq::AverageAtomicMass:
      s = "average atomic mass";
      break;
    case Pq::AverageNuclearCharge:
      s = "average nuclear charge";
      break;
    case Pq::BaryonNumberDensity:
      s = "baryon number density";
      break;
    case Pq::ChemPotentialDiff:
      s = "difference in n-p chemical potentials";
      break;
    case Pq::Density:
      s = "density";
      break;
    case Pq::DeuteronMassFraction:
      s = "mass fraction of deuterons";
      break;
    case Pq::ElectronFraction:
      s = "electron fraction";
      break;
    case Pq::Enthalpy:
      s = "enthalpy";
      break;
    case Pq::HeavyNucMassFraction:
      s = "mass fraction of heavy nuclei";
      break;
    case Pq::HelionMassFraction:
      s = "mass fraction of helions";
      break;
    case Pq::Mass:
      s = "mass";
      break;
    case Pq::NeutronMassFraction:
      s = "mass fraction of neutrons";
      break;
    case Pq::Other:
      s = "other";
      break;
    case Pq::Pressure:
      s = "pressure";
      break;
    case Pq::ProtonMassFraction:
      s = "mass fraction of protons";
      break;
    case Pq::RestMassDensity:
      s = "rest mass density";
      break;
    case Pq::SoundSpeedSq:
      s = "sound speed squared";
      break;
    case Pq::SpecificEntropy:
      s = "specific entropy";
      break;
    case Pq::SpecificInternalEnergy:
      s = "specific internal energy";
      break;
    case Pq::Temperature:
      s = "temperature";
      break;
    case Pq::TritonMassFraction:
      s = "mass fraction of tritons";
      break;
    case Pq::Velocity:
      s = "velocity";
      break;
    case Pq::Undefined:
      s = "undefined";
      break;
    default:
      s = "";
  }
  return s;
}

void
FillPhysicalDimensions (const Pq q, float retval[]) {

  const unsigned short ms = 1;
  const unsigned short ln = 0;
  const unsigned short tm = 2;
  const unsigned short Tp = 3;
  const unsigned short num_basic_units = 4;
  for(int i=0;i<num_basic_units;i++) retval[i] = 0;

  switch(q) {
    // dimensionless quantities and other
    case Pq::Other:
    case Pq::ElectronFraction:
    case Pq::NeutronMassFraction:
    case Pq::ProtonMassFraction:
    case Pq::DeuteronMassFraction:
    case Pq::TritonMassFraction:
    case Pq::HelionMassFraction:
    case Pq::AlphaMassFraction:
    case Pq::HeavyNucMassFraction:
    case Pq::AverageAtomicMass:
    case Pq::AverageNuclearCharge:
      break;

    case Pq::Density:
    case Pq::RestMassDensity:
      retval[ms] = 1;
      retval[ln] = -3;
      break;

    case Pq::BaryonNumberDensity:
      retval[ln] = -3;
      break;

    case Pq::Pressure:
      retval[ms] = 1;
      retval[ln] = -1;
      retval[tm] = -2;
      break;

    case Pq::Enthalpy:
    case Pq::SpecificInternalEnergy:
    case Pq::SoundSpeedSq:
    case Pq::ChemPotentialDiff:
      retval[ln] = 2;
      retval[tm] = -2;
      break;

    case Pq::SpecificEntropy:
      retval[ln] = 2;
      retval[tm] = -2;
      retval[Tp] = -1;
      break;

    case Pq::Mass:
      retval[ms] = 1;
      break;

    case Pq::Temperature:
      retval[Tp] = 1;
      break;

    case Pq::Velocity:
      retval[ln] = 1;
      retval[tm] = -1;
      break;

    default:
      assert(false); // unknown Pq
  }
}


// helper function for class 'Pq'
std::string
PqUnitsString(const Pq pq_kind, const PhUnits u) {
  std::string str;
  switch(pq_kind) {
    case Pq::Density:
    case Pq::RestMassDensity:
      switch(u) {
        case PhUnits::CGS:
          str = "g/cm^3";
          break;
        case PhUnits::SI:
          str = "kg/m^3";
          break;
        case PhUnits::NUCLEAR:
          str = "MeV/c^2/fm^3";
          break;
        case PhUnits::GR:
          str = "1/M^2";
          break;
        default:
          assert (false);
      }
      break;

    case Pq::Pressure:
      switch(u) {
        case PhUnits::CGS:
          str = "Ba";
          break;
        case PhUnits::SI:
          str = "Pa";
          break;
        case PhUnits::NUCLEAR:
          str = "MeV/fm^3";
          break;
        case PhUnits::GR:
          str = "1/M^2";
          break;
        default:
          assert (false);
      }
      break;

    case Pq::SpecificInternalEnergy:
    case Pq::Enthalpy:
    case Pq::ChemPotentialDiff:
      switch(u) {
        case PhUnits::CGS:
          str = "erg/g";
          break;
        case PhUnits::SI:
          str = "J/kg";
          break;
        case PhUnits::NUCLEAR:
          str = "MeV/baryon";
          break;
        case PhUnits::GR:
          str = "1";
          break;
        default:
          assert (false);
      }
      break;

    case Pq::SoundSpeedSq:
      switch(u) {
        case PhUnits::CGS:
          str = "cm^2/s^2";
          break;
        case PhUnits::SI:
          str = "m^2/s^2";
          break;
        case PhUnits::NUCLEAR:
          str = "c^2";
          break;
        case PhUnits::GR:
          str = "1";
          break;
        default:
          assert (false);
      }
      break;

    case Pq::SpecificEntropy:
      switch(u) {
        case PhUnits::CGS:
          str = "erg/g/K"; // TODO: check this!
          break;
        case PhUnits::SI:
          str = "J/kg/K";
          break;
        case PhUnits::NUCLEAR:
          str = "kB/baryon";
          break;
        case PhUnits::GR:
          str = "1/K";
          break;
        default:
          assert (false);
      }
      break;

    case Pq::Temperature:
      switch(u) {
        case PhUnits::CGS:
        case PhUnits::GR:
        case PhUnits::SI:
          str = "K";
          break;
        case PhUnits::NUCLEAR:
          str = "MeV";
        default:
          assert (false);
      }
      break;

    case Pq::ElectronFraction:
    case Pq::NeutronMassFraction:
    case Pq::ProtonMassFraction:
    case Pq::DeuteronMassFraction:
    case Pq::TritonMassFraction:
    case Pq::HelionMassFraction:
    case Pq::AlphaMassFraction:
    case Pq::HeavyNucMassFraction:
    case Pq::AverageAtomicMass:
    case Pq::AverageNuclearCharge:
      str = "1";
      break;

    case Pq::Undefined:
    case Pq::Other:
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
          str = "geom";
          break;
        default:
          assert (false);
      }
      break;
    
    case Pq::Mass:
      switch(u) {
        case PhUnits::CGS:
          str = "g";
          break;
        case PhUnits::SI:
          str = "kg";
          break;
        case PhUnits::NUCLEAR:
          str = "MeV/c^2";
          break;
        case PhUnits::GR:
          str = "M";
          break;
        default:
          assert (false);
      }
      break;

    case Pq::Velocity:
      switch(u) {
        case PhUnits::CGS:
          str = "cm/s";
          break;
        case PhUnits::SI:
          str = "m/s";
          break;
        case PhUnits::NUCLEAR:
          str = "1";
          break;
        case PhUnits::GR:
          str = "1";
          break;
        default:
          assert (false);
      }
      break;

    default:
      assert (false);
  } // switch pq_kind
  return str;
}

