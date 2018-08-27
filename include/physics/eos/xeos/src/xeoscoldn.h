/*~-------------------------------------------------------------------~*
 * Copyright (c) 2018 Los Alamos National Security, LLC
 * All rights reserved.
 *~-------------------------------------------------------------------~*/

/**
 * @file xeoscoldn.h
 * @author Oleg Korobkin
 * @date 2018-01-05
 * @brief XeosColdNuclear: tabulated 1D equation of state
 *        for cold nuclear matter
 */
#ifndef XEOS_XEOSCOLDN_H_
#define XEOS_XEOSCOLDN_H_

#include <math.h>
#include "physquan.h"
#include "xeosbase.h"
#include "filereader.h"

namespace xeos {

/*
 * enum for format types
 */
enum class XeosColdNuclearFormat {
  kFourColumnsCGS,
  kSixteenColumnsNuclear
};


/**
 * Single-parameter tabulated EoS for cold nuclear matter
 * TODO: implement details, give usage example here
 */
template <class FR>
class XeosColdNuclear : public XeosTabulated {
 public:
  typedef enum XeosColdNuclearFormat format_type;

  // default constructor
  XeosColdNuclear<FR>()
      : XeosColdNuclear<FR>( "",
        format_type::kFourColumnsCGS,
        PhysicalQuantity::GetGlobalUnits()) {}

  // constructor which omits units
  XeosColdNuclear<FR> (const std::string _path, const format_type _fmt)
      : XeosColdNuclear<FR>( _path, _fmt,
        PhysicalQuantity::GetGlobalUnits()) {}

  // workhorse constructor declaration
  XeosColdNuclear<FR> (const std::string _path, const format_type _fmt,
      const PhUnits _u);

  // destructor
  virtual ~XeosColdNuclear<FR> ();

  // report EoS units
  PhUnits GetUnits() const { return eos_units; }
  format_type GetFormat() const { return format; }

  virtual void ReadEosTables() override;

  virtual double
  GridPoint(const Pq var, const int i) const override;

  virtual const int
  GridSize(const Pq) const override { return num_records; }

  virtual std::vector<Pq>
  GetPrimaryQuantities() const override { return *primary_vars; }

  virtual std::vector<Pq>
  GetAdditionalQuantities() const override { return additional_vars; }

  virtual void
  operator() (const int num_out, eos_in in_array[],
      eos_out out_array[]) override;

  virtual void
  DfDx (const int num_out, const Pq dX, eos_in  in_array[],
      eos_out dF_array[]) override;

  virtual bool
  ConsistencyCheck() const override;
  double consistency_check_threshold = 2e-5;

 protected:
  // units in which to keep the tables
  PhUnits eos_units = PhysicalQuantity::GetGlobalUnits();

  using XeosTabulated::data_path;

 private:

  std::vector<Pq> *primary_vars;
  std::vector<Pq> additional_vars;
  FR freader;
  int num_records;
  int num_fields;
  format_type format;
  NVector<double> *table_nvs[16];
  PhysicalQuantity *table_pqs[16];

}; // class XeosColdNuclear


//
// XeosColdNuclear: workhorse constructor definition
//
template <class FR>
XeosColdNuclear<FR>::XeosColdNuclear (
    const std::string _path,
    const format_type _fmt,
    const PhUnits _u) {

  XeosTabulated::data_path = _path;
  format = _fmt;
  eos_units = _u;
  freader.SetFilename(_path);

  NvPressure               *pres_tab;
  NvDensity                 *rho_tab;
  NvEnthalpy               *enth_tab;
  NvBaryonNumberDensity    *nbar_tab;
  NvElectronFraction         *ye_tab;
  NvSpecificInternalEnergy  *eps_tab;
  NvSoundSpeedSq            *cs2_tab;
  NvChemPotentialDiff       *muh_tab;
  NvSpecificEntropy         *ent_tab;

  NvNeutronMassFraction      *Xn_tab;
  NvProtonMassFraction       *Xp_tab;
  NvDeuteronMassFraction     *Xd_tab;
  NvTritonMassFraction       *Xt_tab;
  NvHelionMassFraction       *Xh_tab;
  NvAlphaMassFraction        *Xa_tab;
  NvHeavyNucMassFraction    *Xhv_tab;
  NvAverageAtomicMass        *Ab_tab;
  NvAverageNuclearCharge     *Zb_tab;

  // setup primary and secondary variables
  switch(_fmt) {
    case format_type::kFourColumnsCGS:
      primary_vars = new std::vector<Pq> {
          Pq::Density,
          Pq::Pressure,
          Pq::Enthalpy,
          Pq::BaryonNumberDensity
      };

      num_fields = 4; // number of fields depends on format
      rho_tab = new NvDensity(PhUnits::CGS);
      pres_tab = new NvPressure(PhUnits::CGS);
      enth_tab = new NvEnthalpy(PhUnits::CGS);
      nbar_tab = new NvBaryonNumberDensity(PhUnits::CGS);

      table_nvs[0] = rho_tab;     table_pqs[0] = rho_tab;
      table_nvs[1] = pres_tab;    table_pqs[1] = pres_tab;
      table_nvs[2] = enth_tab;    table_pqs[2] = enth_tab;
      table_nvs[3] = nbar_tab;    table_pqs[3] = nbar_tab;
      break;

    case format_type::kSixteenColumnsNuclear:
      primary_vars = new std::vector<Pq> {
          Pq::Pressure,
          Pq::Density,
          Pq::ElectronFraction,
          Pq::SpecificInternalEnergy
      };

      additional_vars.push_back(Pq::ChemPotentialDiff);
      additional_vars.push_back(Pq::SoundSpeedSq);
      additional_vars.push_back(Pq::SpecificEntropy);
      additional_vars.push_back(Pq::NeutronMassFraction);

      additional_vars.push_back(Pq::ProtonMassFraction);
      additional_vars.push_back(Pq::DeuteronMassFraction);
      additional_vars.push_back(Pq::TritonMassFraction);
      additional_vars.push_back(Pq::HelionMassFraction);

      additional_vars.push_back(Pq::AlphaMassFraction);
      additional_vars.push_back(Pq::HeavyNucMassFraction);
      additional_vars.push_back(Pq::AverageAtomicMass);
      additional_vars.push_back(Pq::AverageNuclearCharge);

      num_fields = 16;
      pres_tab = new NvPressure(PhUnits::NUCLEAR);
      rho_tab  = new NvDensity(PhUnits::CGS);
      ye_tab   = new NvElectronFraction(PhUnits::NUCLEAR);
      eps_tab  = new NvSpecificInternalEnergy(PhUnits::NUCLEAR);

      muh_tab  = new NvChemPotentialDiff {PhUnits::NUCLEAR};
      cs2_tab  = new NvSoundSpeedSq {PhUnits::CGS};
      ent_tab  = new NvSpecificEntropy {PhUnits::NUCLEAR};
      Xn_tab   = new NvNeutronMassFraction {PhUnits::NUCLEAR};

      Xp_tab   = new NvProtonMassFraction {PhUnits::NUCLEAR};
      Xd_tab   = new NvDeuteronMassFraction {PhUnits::NUCLEAR};
      Xt_tab   = new NvTritonMassFraction {PhUnits::NUCLEAR};
      Xh_tab   = new NvHelionMassFraction {PhUnits::NUCLEAR};

      Xa_tab   = new NvAlphaMassFraction {PhUnits::NUCLEAR};
      Xhv_tab  = new NvHeavyNucMassFraction {PhUnits::NUCLEAR};
      Ab_tab   = new NvAverageAtomicMass {PhUnits::NUCLEAR};
      Zb_tab   = new NvAverageNuclearCharge {PhUnits::NUCLEAR};

      table_nvs[ 0] = pres_tab;     table_pqs[ 0] = pres_tab;
      table_nvs[ 1] = rho_tab;      table_pqs[ 1] = rho_tab;
      table_nvs[ 2] = ye_tab;       table_pqs[ 2] = ye_tab;
      table_nvs[ 3] = eps_tab;      table_pqs[ 3] = eps_tab;

      table_nvs[ 4] = muh_tab;      table_pqs[ 4] = muh_tab;
      table_nvs[ 5] = cs2_tab;      table_pqs[ 5] = cs2_tab;
      table_nvs[ 6] = ent_tab;      table_pqs[ 6] = ent_tab;
      table_nvs[ 7] = Xn_tab;       table_pqs[ 7] = Xn_tab;

      table_nvs[ 8] = Xp_tab;       table_pqs[ 8] = Xp_tab;
      table_nvs[ 9] = Xd_tab;       table_pqs[ 9] = Xd_tab;
      table_nvs[10] = Xt_tab;       table_pqs[10] = Xt_tab;
      table_nvs[11] = Xh_tab;       table_pqs[11] = Xh_tab;

      table_nvs[12] = Xa_tab;       table_pqs[12] = Xa_tab;
      table_nvs[13] = Xhv_tab;      table_pqs[13] = Xhv_tab;
      table_nvs[14] = Ab_tab;       table_pqs[14] = Ab_tab;
      table_nvs[15] = Zb_tab;       table_pqs[15] = Zb_tab;

      break;

    default:
      assert(false); // format not implemented
  }

  ReadEosTables();

  if (!ConsistencyCheck())
    std::cerr << "WARNING(XeosColdNuclear): equation of state"
              << "failed consistency check!" << std::endl;
}

//
// XeosColdNuclear: virtual destructor
//
template <class FR>
XeosColdNuclear<FR>::~XeosColdNuclear () {
  delete primary_vars;
}

//
// File reader (TODO: implement with templates)
//
template <class FR>
void XeosColdNuclear<FR>::ReadEosTables() {
  freader.Open();
  int file_len = freader.NumLines();
  int header_len;
  double fields[16];

  PqVelocity clight(1.0, PhUnits::GR); // light speed squared
  clight.ConvertTo(PhUnits::CGS);
  const double c2 = (double)clight * (double)clight;

  switch (format) {
    case format_type::kFourColumnsCGS:
      header_len = 1;
      break;

    case format_type::kSixteenColumnsNuclear:
      freader.Rewind();
      header_len = freader.SkipHashHeader();
      break;

    default:
      assert (false);
  }

  // establish number of records in 1D table and resize the table
  num_records = file_len - header_len;
  for (int i=0; i<num_fields; ++i)
    table_nvs[i]-> Resize(num_records);

  freader.Rewind();
  freader.SkipHeader(header_len);
  for (int i=0; i<num_records; ++i) {
    if (freader.ReadFields(num_fields,fields) < num_fields) {
      num_records = i;
      break; // skip empty / incomplete lines in the end of file
    }

    for (int j=0;j<num_fields;++j)
      (*table_nvs[j])(i)= fields[j];

    // format-specific quirks
    switch (format) {
      case format_type::kFourColumnsCGS:
        (*table_nvs[2])(i)= exp(fields[2]/c2) * c2;
        break;

      case format_type::kSixteenColumnsNuclear:
        (*table_nvs[1])(i)= pow(10.,fields[1]);
        break;

      default:
        assert (false);
    }
  }

  for (int i=0;i<num_fields;++i)
    table_pqs[i]-> ConvertTo(eos_units);

} // ReadEosTable()


//
// Return a point of the one-dimensional grid
//
template <class FR>
double XeosColdNuclear<FR>::GridPoint (
    const Pq var, const int i) const {
  double retval;
  int j;
  assert (i>=0 && i<num_records);

  for (j=0; j<num_fields; ++j) {
    if (var == table_pqs[j]-> PqKind())
      break;
  }
  if (j<num_fields)
    retval = (*table_nvs[j])(i);
  else
    assert (false); // physical quantity is not found

  return retval;
}

//
// General EoS function
//
template <class FR>
void XeosColdNuclear<FR>::operator() (const int num_out,
    eos_in in_array[],
    eos_out out_array[]) {
  /*      +
     TODO
   +      */
}

//
// General derivative(s), (dF/dX)_const
//
template <class FR>
void XeosColdNuclear<FR>::DfDx (const int num_out, const Pq dX,
    eos_in  in_array[],
    eos_out dF_array[]) {
  /*      +
     TODO
   +      */
}

template <class FR>
bool XeosColdNuclear<FR>::ConsistencyCheck() const {
  bool retval = true;

  double rho, pres, enth, nbar;
  PqVelocity clight(1.0, PhUnits::GR);
  PqMass m_baryon(1.66e-24, PhUnits::CGS);
  clight.ConvertTo(eos_units);
  m_baryon.ConvertTo(eos_units);
  double c2 = (double)clight * (double)clight;
  double mb = (double)m_baryon;

  switch(format) {
    case XeosColdNuclearFormat::kFourColumnsCGS:
      for (int i=0; i<num_records; ++i) {
        rho =  (*table_nvs[0])(i);
        pres = (*table_nvs[1])(i);
        enth = (*table_nvs[2])(i);
        nbar = (*table_nvs[3])(i);
        if (consistency_check_threshold <
            fabs(1.0 - (rho*c2 + pres)/(nbar*mb*enth))) {
          retval = false;
          std::cerr << std::scientific << std::setprecision(12);
          std::cerr << "XeosColdNuclear::ConsistencyCheck failed "
                    << "at grid point i=" << i << ": " << std::endl;
          std::cerr << "    (rho + pres/c2)/(nbar*mb)| = "
                    << (rho + pres/c2)/(nbar*mb) << std::endl;
          std::cerr << "    enth/c2 = " << enth/c2 << std::endl;
        }
      }
    break;

  case XeosColdNuclearFormat::kSixteenColumnsNuclear:
    /*      +
       TODO
     +      */
    retval = true;
    break;

  default:
    assert (false);
  } // switch (format)

  return retval;
}

} // namespace xeos
#endif // XEOS_XEOSCOLDN_H_

#if 0
#include <iostream>

int main() {
  using namespace std;
  using namespace xeos;

  PhysicalQuantity::SetGlobalUnits(PhUnits::NUCLEAR);

  XeosColdNuclear<AsciiFileReader>
      eosA("../data/rnsid/eosA",
      XeosColdNuclearFormat::kFourColumnsCGS);
  XeosColdNuclear<AsciiFileReader>
      eSFHo("../data/rnsid/sfho_0.1MeV_beta.txt",
      XeosColdNuclearFormat::kSixteenColumnsNuclear);

  cout << scientific << setprecision(12) << setw(18);
  for (int i=0;i<eSFHo.GridSize(Pq::Density);++i)
    cout << eSFHo.GridPoint(Pq::Density, i) << " "
         << eSFHo.GridPoint(Pq::Pressure,i) << " "
         << eSFHo.GridPoint(Pq::ElectronFraction,i) << " "
         << eSFHo.GridPoint(Pq::ChemPotentialDiff,i) << " "
         << eSFHo.GridPoint(Pq::SoundSpeedSq,i) << " "
         << eSFHo.GridPoint(Pq::SpecificEntropy,i) << " "
         << eSFHo.GridPoint(Pq::NeutronMassFraction,i) << " "
         << eSFHo.GridPoint(Pq::ProtonMassFraction,i) << " "
         << eSFHo.GridPoint(Pq::DeuteronMassFraction,i) << " "
         << eSFHo.GridPoint(Pq::TritonMassFraction,i) << " "
         << eSFHo.GridPoint(Pq::HelionMassFraction,i) << " "
         << eSFHo.GridPoint(Pq::AlphaMassFraction,i) << " "
         << eSFHo.GridPoint(Pq::HeavyNucMassFraction,i) << " "
         << eSFHo.GridPoint(Pq::AverageAtomicMass,i) << " "
         << eSFHo.GridPoint(Pq::AverageNuclearCharge,i) << " "
         << endl;
  cout << endl << endl;

  for (int i=0;i<eosA.GridSize(Pq::Density);++i)
    cout << eosA.GridPoint(Pq::Density, i) << " "
         << eosA.GridPoint(Pq::Pressure,i) << endl;

}
#endif
