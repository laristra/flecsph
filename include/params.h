/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Los Alamos National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/

 /*~--------------------------------------------------------------------------~*
 * 
 * /@@@@@@@@  @@           @@@@@@   @@@@@@@@ @@@@@@@  @@      @@
 * /@@/////  /@@          @@////@@ @@////// /@@////@@/@@     /@@
 * /@@       /@@  @@@@@  @@    // /@@       /@@   /@@/@@     /@@
 * /@@@@@@@  /@@ @@///@@/@@       /@@@@@@@@@/@@@@@@@ /@@@@@@@@@@
 * /@@////   /@@/@@@@@@@/@@       ////////@@/@@////  /@@//////@@
 * /@@       /@@/@@//// //@@    @@       /@@/@@      /@@     /@@
 * /@@       @@@//@@@@@@ //@@@@@@  @@@@@@@@ /@@      /@@     /@@
 * //       ///  //////   //////  ////////  //       //      //  
 *
 *~--------------------------------------------------------------------------~*/

/**
 * @file params.h
 * @author Oleg Korobkin
 * @date June 2018
 * @brief Handle global parameters: parsing parameter files, default vals etc. 
 */

#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <sstream>
#include <stdlib.h>
#include <assert.h>

#define MAXLEN_FNAME 255

namespace param {

  // --- list of global parameters ------------------------------------------
  int64_t nparticles;        //< global number of particles
  int    sodtest_num;        //< Sod test number (1..5)
  double poly_gamma;         //< polytropic index
  double initial_dt;         //< intial timestep
  char   initial_data_h5part[MAXLEN_FNAME]; //< h5part output filename
  // ---

/*!
   trim whitespace from the begging and end of line
 */
std::string trim(const std::string& str) {
  using namespace std;
  const size_t strBegin = str.find_first_not_of(" \t");
  if (strBegin == string::npos)
    return "";
  const size_t strEnd = str.find_last_not_of(" \t");
  return str.substr(strBegin, strEnd - strBegin + 1);
}


/*!
   setup global defaults
 */
void set_default_params(int rank, int size) {
  nparticles = 1000;
  poly_gamma = 1.4;
  sodtest_num = 1;
  initial_dt = 1.0;
  sprintf(initial_data_h5part,"%s","init.h5part");
}


/**
 * @brief Sets a global parameter param_name to the value, given by its
 *        string representation param_value
 */
void set_param(const std::string& param_name,
               const std::string& param_value) {
  using namespace std;
  bool run_job;

  // strip trailing comments
  size_t i2 = param_value.find("#");
  string str_value;
  if (i2 != string::npos) 
    str_value = trim(param_value.substr(0,i2));
  else
    str_value = param_value;

  // remove (single or double) quotes 
  if (str_value[0]=='\'' or str_value[0]=='\"') 
    str_value = str_value.substr(1,str_value.length()-2);
  istringstream iss(str_value);

  // boolean parameters ------------------------------
  bool lparam_value = (str_value == "yes" 
                    or str_value == "'yes'" 
                    or str_value == "\"yes\"");
  if (param_name == "run_job") {
    run_job = lparam_value;
  }
 
  // integer parameters ------------------------------
  else if (param_name == "nparticles") {
    iss >> nparticles;
  }
  else if (param_name == "sodtest_num") {
    iss >> sodtest_num;
  }

  // real parameters ---------------------------------
  else if (param_name == "poly_gamma") {
    iss >> poly_gamma;
  }
  else if (param_name == "intial_dt") {
    iss >> initial_dt;
  }

  // string parameters -------------------------------
  else if (param_name == "initial_data_h5part") {
    strcpy(initial_data_h5part, str_value.c_str());
  }
  // unknown parameter -------------------------------
  else { 
    clog_one(fatal) << "ERROR: unknown parameter " << param_name << endl;
    exit(2);
  }
 
  clog_one(trace) << param_name << ": " << param_value << endl;
}               

/**
 * @brief Parameter file parser
 *
 * Simulation parameter file contains strings of pairs:
 *  param = value
 *
 * assigining values to specific parameters. Parameters can have four 
 * different types: logical (boolean), integer, real and string. 
 * Parameter file can contain comments. Comments begin with '#'; 
 * everything after '#' is ignored. Strings can be embraced in single 
 * quotes (') or double quotes ("). Logical parameters are initialized 
 * with "yes" or "no".
 *
 * Example of a parameter file snippet:
 * --->>> sodtube_N100.par ---------------
  # Sodtube test #1 for 100 particles
  nparticles= 1000    # global number of particles
  poly_gamma = 1.4    # polytropic index
  sodtest_num  = 1    # which Sod test to generate
  output_singlefile = yes
  intial_data_h5part  = "hdf5_sodtube.h5part"
 * ---<<<  --------------------------------
 */
void read_params(std::string parfile) {
  using namespace std;
  ifstream infile;
  string line;
  infile.open(parfile.c_str());
  if (!infile) {
    cerr << "ERROR: Unable to open file " << parfile << endl;
    exit(1);
  }

  for (int ln=1; std::getline(infile,line); ++ln) {

    // skip comments (lines starting with '#' at any position)
    // and blank lines
    bool is_blank = true, is_comment = false;
    for (int i=0;i<line.length();i++) {
      char c = line[i];
      is_comment = (c == '#');
      is_blank = (c == ' ' || c == '\t');
      if (is_comment or not is_blank) break;
    }
    if (is_comment or is_blank) continue;

    // line must be in the form:
    // lhs = rhs
    size_t eqpos = line.find("=");
    if (eqpos == string::npos) {
      cerr << "ERROR in parameter file " << parfile
           << ", line #" << ln << endl;
      exit(1);
    }
    string lhs = line.substr(0,eqpos);
    string rhs = line.substr(eqpos+1);

    // strip whitespace
    lhs = trim(lhs);
    rhs = trim(rhs);
    if (lhs.length() == 0 or rhs.length() == 0) {
      cerr << "ERROR in parameter file " << parfile
           << ", line #" << ln << ": wrong format" << endl;
      exit(1);
    }
    set_param(lhs,rhs);
  }
  infile.close();
}

} // namespace params

