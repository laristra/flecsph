# STELLAR COLLAPSE EOS READER

Aug.28.2018

This is third party library that reads various tabulated EOS 
to simulate merging compact object binaries. EOSs can be found
[here](https://stellarcollapse.org/equationofstate)

It is orignated from EOS reader in `νbhlight` which is 
developed/maintained by Jonah Miller in LANL. Below is
basic description about code

* `decs.h` contains the public API for all implemented functions, 
along with the structs, includes, externs, and macros. 
* `constants.h` contains global constants. It is used for units.
* `eos_stellar_collapse.h` contains private variables and prototypes 
for internal structures for the table reader
* `eos.c` contains wrapper functions that incorporate tabulated, 
ideal gas, and polytrope EOS's. 
* `eos_poly.c` is polytrope EOS code
* `eos_gamma.c` is ideal gas EOS code
* `eos_stellar_collapse.c` is implementation of the polytrope reader
* `root_finding.c` is implementation of root finding for inverting table
* `util.c` implementation of some useful utility functions, such as 1D 
interpolation and setting units.

## Note for units
Unit system can vary with respect to each code. Here, we describe the
original unit system for this code. Every unit conversions can be done
easily.

This code basically adapts the unit system in `νbhligh`.
The idea is that each `X_unit` variable should be set to some value in CGS. 
Then dividing by that unit gives you code units and multiplying by it gives 
you CGS. The user specifies mass via `M_unit`. Then the length unit is either 
specified by `L_unit` or by multiplying a black hole mass by `G/c^2`. All other 
units are then derived from that so that code units are GR units with `G=c=1`.

One exception to the units scheme: In `νbhlight, it uses MeV as code units 
for temperature. The cgs units are Kelvin, but the code doesn't convert in 
and out of temperature units. Everything is in MeV.

# Usage

You just simply type:

```{engine=sh}
   make
```

Then it will generate `libTABEOS_SC.a` and `libTABEOS_RF.a` as libs

# Contact

If you have any questions to use this, please contact Jonah Miller (jonahm@lanl.gov), 
Hyun Lim (hylim1988@gmail.com), or Oleg Korobkin (korobkin@lanl.gov)
