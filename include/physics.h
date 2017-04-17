#ifndef physics_h 
#define physics_h 

#include <vector>

#include "tree.h"

namespace physics{

const double kHeatRatio = 2.0;
const double kViscAlpha = 2.0;
const double kViscBeta  = 2.0;
const double kViscEta = 0.01; 
const double kGravConstant = 1.0;
const double kCoeffDt = 0.1;


// Utils functions
double kernel(double, double);
double mu(body*, body*);

// Main functions
void computeDensity(body*, std::vector<body*>&);
void computeSoundspeed(body*);
void computePressure(body*);
void computeHydro(body*, std::vector<body*>&);
void computeAcceleration(body*,double);
void computeGrav(body*,std::vector<body*>&);
void moveBody(body*,double);

} // namespace

#endif // physics_h  
