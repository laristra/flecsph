## Index of equations of state files

(collection of files as supplied with RNS code by N. Stergioulas)

| file name	| description                                                    |
| --------- |:-------------------------------------------------------------- |
| eosA		| Pandharipande neutron: A&B EOS A                               |
| eosAU		| WFF1 (denoted AV14+UVII in WFF) matched to Negele-Vautherin    |
| eosB		| Pandharipande hyperon: A&B EOS B                               |
| eosC		| Bethe-Johnson model 1: A&B EOS C                               |
| eosF		| Pandharipande (hyperon) and Arponen: A&B EOS F                 |
| eosFP		| Friedman & Pandharipande EOS.(Negele-Vautherin for             |
| 			| N<.1 fm^-3,then * BPS for N<.001 fm^-3                         |
| eosFPS	| Lorenz, Ravenhall and Pethick, 1993, PRL 70,379                |
| eosG		| Canuto-Chitre (solid): A&B EOS G                               |
| eosL		| Pandharipande and Smith (mean field): A&B EOS L                |
| eosN		| Walecka-Serot: A&B EOS N                                       |
| eosNV		| Negele & Vautherin, 1973, Nucl. Phys. A207, 298                |
| eosO		| Bowers, Gleeson, and Pedigo. A&B EOS O                         |
| eosUU		| WFF2 (denoted UV14+UVII in WFF) matched to Negele-Vautherin    |
| eosWNV	| 	WFF3 (denoted UV14+TNI in WFF) matched to Negele-Vautherin   |
| eosWS		| WFF3 (UV14+TNI) mathed to EOS FPS                              |

A&B = Arnett and Bowers (1977) APJS 33, 415
WFF = Wiringa, Fiks and Fabrocini (1988), Phys. Rev. C 38, 1010


The equation of state files are in the format required for
rotating star codes by Nikolaos Stergioulas:
`rns.c`	(rapidly rotating neutron stars in equilibrium)

The format required is as follows:
line 1) number of tabulated points in file
remaining lines) four columns:

1. column 1) energy density/c^2 (g/cm^3)
2. column 2) pressure (dynes/cm^2)
3. column 3) enthalpy (cm^2/s^2)  !!! Actually, c^2 log(h/c^2) !!!
4. column 4) baryon number density (1/cm^3)

Baryon mass is set to mbaryon = 1.66e-24 g.

Also, the database contains files in more detailed 16-column
format, representing beta-equilibrium configuration for low-temperature
1D slices of the following equations:
- SFHo, SFHx: Steiner, Hempel & Fischer (2013), ApJ 774:
  `sfho_0.1MeV_beta.dat`, 
  `sfhx_0.1MeV_beta.dat`
- DD2: Hempel et al. (2012), ApJ 748
  `dd2_0.1MeV_beta.dat`
- LS220: Lattimer & Swesty (1991), Nuc. Phys. A, 535
  `ls220_0.01MeV_beta.dat`,
  `ls220_0.1MeV_beta.dat`

