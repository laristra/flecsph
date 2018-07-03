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
 *
 * Parameters are essential for maintaining data in scientific numerical
 * experiments. They allow to launch and keep track of a battery of experiments
 * without recompiling the code. State-of-the-art simulation codes contain
 * hundreds, sometimes thousands of parameters. Clearly, it is impossible to
 * remember what are their optimal values, which also may vary case-tocase.
 * This is why it is important to store parameters in separate files.
 *
 * For FleCSPH, we picked a simple ASCII format, which is easy to read in any
 * system. The format contains pairs of "PARAM = VALUE". There are four types
 * of parameters: integer, real, string and boolean. Integer and real
 * parameters are intialized as usual, e.g.:
 *
 *   nparticles = 250
 *   poly_gamma = 1.4
 *
 * Boolean parameters are set with "yes" or "no", and string ones are set with
 * strings, which can be enclosed in single or double quotes:
 *
 *   output_conserved = yes                # produce scalar output
 *   initial_data_prefix = "sodtube_n10k"  # name of the input file
 *
 * Comments can be added to parameter files: everything which starts with '#',
 * is considered a comment and ignored. Blank spaces and empty lines are
 * ignored.
 *
 * All parameters in a parfile have exactly the same name as in the code, to
 * avoid confusion. All parameters are declared as const references in the
 * param:: namespace. It is also possible to #define a parameter instead to
 * ioptimze the performance. In this case, when a parameter has been defined via 
 * a macro, it needs to be commented out in a parametre file. 
 *
 * Example of #define-ing parameters:
 * 
 *  #define nparticles 1000
 *  #define sodtest_num  1
 *  #define poly_gamma  1.4
 *  #define initial_data_prefix "sodtube_n10k"
 * 
 * To introduce a new parameter:
 *  - add its declaration below using DECLARE_PARAM or DECLARE_STRING_PARAM 
 *    macro (see below);
 *  - add REAL_NUMERIC_PARAM, READ_BOOLEAN_PARAM or READ_STRING_PARAM
 *    to the set_param() function -- see examples below.
 *
 * TODO: introduce a set of macros to do both declaration and reading:
 *       REGISTER_BOOLEAN_PARAM, REGISTER_INTEGER_PARAM, 
 *       REGISTER_STRING_PARAM, REGISTER_REAL_PARAM
 *       instead of DECLARE/READ pair
 */

#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <sstream>
#include <stdlib.h>
#include <assert.h>
#include "cinchlog.h"

#ifndef PARAMS_H
#define PARAMS_H

//////////////////////////////////////////////////////////////////////
#define DECLARE_PARAM(PTYPE,PNAME,PDEF) \
  PTYPE _##PNAME = (PDEF); \
  const PTYPE & PNAME = _##PNAME;

#define STRING_MAXLEN 511
#define DECLARE_STRING_PARAM(PNAME,PDEF) \
  char _##PNAME[STRING_MAXLEN] = (PDEF); \
  const char * PNAME = _##PNAME;

#define DO_QUOTE(str) #str
#define QUOTE(str) DO_QUOTE(str)

#define READ_BOOLEAN_PARAM(PNAME) \
  if (param_name == QUOTE(PNAME)) { \
    (_##PNAME) = lparam_value; unknown_param=false;}

#define READ_NUMERIC_PARAM(PNAME) \
  if (param_name == QUOTE(PNAME)) { \
    iss >> (_##PNAME); unknown_param=false;}

#define READ_STRING_PARAM(PNAME) \
  if (param_name == QUOTE(PNAME)) { \
    strcpy(_##PNAME, str_value.c_str()); unknown_param = false;}

//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
///
/// Parameter #define-s: uncomment to fix parameters at compile time:
///
//#define nparticles 1000
//#define sodtest_num  1
//#define poly_gamma  1.4
//#define initial_data_prefix "initial_data"

namespace param {

/// final iteration (= total iterations + 1, if counting from 0)
#ifndef final_iteration
  DECLARE_PARAM(int64_t,final_iteration,10)
#endif

/// global number of particles
#ifndef nparticles
  DECLARE_PARAM(int64_t,nparticles,1000)
#endif

/// number of Sodtest to run (1..5)
#ifndef sodtest_num
  DECLARE_PARAM(unsigned short,sodtest_num,1)
#endif

/// polytropic index
#ifndef poly_gamma
  DECLARE_PARAM(double,poly_gamma,1.4)
#endif

/// inital timestep
#ifndef initial_dt
  DECLARE_PARAM(double,initial_dt,0.001)
#endif

/// file prefix for input and intiial data file[s]
#ifndef initial_data_prefix
  DECLARE_STRING_PARAM(initial_data_prefix,"initial_data")
#endif

/// a test boolean parameter
#ifndef run_job
  DECLARE_PARAM(bool,run_job,false)
#endif
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


/**
 * @brief Sets a global parameter param_name to the value, given by its
 *        string representation param_value
 */
void set_param(const std::string& param_name,
               const std::string& param_value) {
  using namespace std;
  bool run_job;
  bool unknown_param = true;

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
# ifndef run_job
    READ_BOOLEAN_PARAM(run_job)
# endif

  // integer parameters ------------------------------
# ifndef final_iteration
    READ_NUMERIC_PARAM(final_iteration)
# endif

# ifndef nparticles
    READ_NUMERIC_PARAM(nparticles)
# endif

# ifndef sodtest_num
    READ_NUMERIC_PARAM(sodtest_num)
# endif

  // real parameters ---------------------------------
# ifndef poly_gamma
    READ_NUMERIC_PARAM(poly_gamma)
# endif

# ifndef initial_dt
    READ_NUMERIC_PARAM(initial_dt)
# endif

  // string parameters -------------------------------
# ifndef initial_data_prefix
    READ_STRING_PARAM(initial_data_prefix)
# endif

  // unknown parameter -------------------------------
  if (unknown_param) {
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
 * Example parameter file:
 * --->>> sodtube_N100.par ---------------
  # Sodtube test #1 for 100 particles
  nparticles= 1000    # global number of particles
  poly_gamma = 1.4    # polytropic index
  sodtest_num  = 1    # which Sod test to generate
  output_singlefile = yes
  initial_data_prefix  = "sodtube_np10k"
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

#endif // PARAMS_H
