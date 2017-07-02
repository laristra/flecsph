# Usage guide 

For the input file we are using the H5hut format. 
Headers containts general informations like: 

- Number of particles: "nparticles"
- Dimension: "dimension"
- Timestep: "timestep"
- Is fixed timestep used ? "used_fixed_timestep"

Not implemented yet: 
- Physics constants ? 
- Different files for output ? See Oleg for that
- Different types of EOS

Then for each Step we save:

Header:
- Timestep "timestep"
- ???

Particles:
- Position X: "x"
- Position Y: "y" 
- Position Z: "z" 
- Velocity X: "vx"
- Velocity Y: "vy" 
- Velocity Z: "vz"
- Acceleration X: "ax"
- Acceleration Y: "ay"
- Acceleration Z: "az"
- Smoothing Length: "h"
- Density: "rho"
- Internal Energy: "u"
- Pressure: "P"
- Mass: "m"
- Id: "id" 
- Time step: "dt"
 
Not implemented yet:
- Particle type: "type"
- Electron fraction: "Ye"

Types are all double except for:

- id = int64_t
- nparticles = int64_t 
