# XEOS library

This project aims at creating a uniform user-friendly
C++ API to a variety of equations of state, from simple analytic 
ones such as a polytrope or an ideal gas, to more complex 
multi-parametric EoS tables like the ones provided by 
[CompOSE](https://compose.obspm.fr/).

Use case examples:

1. Set global physical unitrs to e.g. CGS
All physical quantities created hereafter will be in CGS
by default, unless specified otherwise.
```cpp
   PhysicalQuantity::SetGlobalUnits (PhUnits::CGS);
```


2. Create the simplest analytic EOS - polytropic.
Notice how you can specify polytropic constant in units which
are different from global. It is automatically converted to 
global units.
```cpp
  XeosPolytropic eos_poly = new XeosPolytropic(
   0.01,              // K in P = K*rho^Gamma 
   4./3.,             // Gamma
   PhUnits::NUCLEAR); // units of K (different from global)
```

3. Create an instance of more complex, tabulated CompOSE EOS with three parameters (TODO):
```cpp
  XeosCompose3D eos_compose = new XeosCompose3D(
   "data/compose/eos.thermo",  // path to the data file
   {Pq::Pressure,              // read only these variables
    Pq::ElectronFraction,
    Pq::Temperature,
    Pq::BaryonMassDensity},  
   PhUnits::NUCLEAR);  // CompOSE tables are in nuclear units
```

4. Equation of state uses physical quantities as input and output
parameters. It computes output using overloaded `operator()`, 
with the first argument being a const reference, and the second 
being a pointer. Consider the following example of using `eos_poly`
to compute pressure from density:
```cpp
  PqDensity  rho(12.0); // define density in [g/cm3]
  PqPressure P;         // declare pressure variable [dynes]
  eos_poly(rho, &P);    // arg1: input, arg2: output.
```

5. The following code compute pressure at {rho, eps, Ye} for the 
CompOSE eos (TODO):
```cpp
  PqSpecificInternalEnergy eps(0.01); // specify eps [erg/g]
  PqElectronFraction         Ye(0.5); // specify Ye  [mol/g]
  eos_compose(rho,eps,Ye, &P);        // args: 1-3: in, 4: out
```

6. All equation-of-state objects are derived from base class
`xeos_base`, which can be used to access both equations in a unified
manner, e.g. through a pointer:
```cpp  
  XeosBase *eosptr = { &eos_poly, &eos_compose };
  (*eosptr[0])(rho, &P);
  (*eosptr[1])(rho,eps,Ye, &P);
```

7. Equations of state classes also take vector quantities as 
arguments. Standard class for vector operations is 
`PhysicalNVector<Kind,ScalarType>`, which is a templated class
for representing N-dimensional arrays of physical quantities. 
Example of computing an array of temperatures:
```cpp
  typedef 
    PhysicalNVector<Pq::Temperature,double> NvTemperature;

  const int Np = 1000;
  NvTemperature Temps(Np);
  NvDensity rhos(Np);
  NvSpecificInternalEnergy epss(Np);
  NvElectronFraction Yess(Np);
  NvPressure Ps(Np);

  for (i=0; i<Np; ++i) {
    rhos[i] = <...> // fill in arrays
    epss[i] = <...>
    Yess[i] = <...>
  }

  // now everything is ready for an EoS function call, 
  // which is just like before for the scalars:
  eos_poly(rhos, &Ps);
  eos_compose(rhos, epss, Yess, &Temps);
```

8. TODO
Compute multiple quantities: some complex equations of
state have interface which allows computation of several
physical quantities at once. 
```cpp
  NvEntropy Ents(Np);
  eos_compose (3, (eos_in) {&rhos, &epss, &Yess},  // input arrays
                  (eos_out){&Temps, &Ps, &Ents);  // output arrays
```


