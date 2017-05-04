#ifndef physics_h 
#define physics_h 

#include <vector>

#include "tree.h"

namespace physics{

extern double dt;
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
void computeDensity(body_holder*, std::vector<body_holder*>&);
void computeSoundspeed(body_holder*);
void computePressure(body_holder*);
void computeHydro(body_holder*, std::vector<body_holder*>&);
void computeAcceleration(body_holder*,double);
void computeGrav(body_holder*,std::vector<body_holder*>&);
void moveBody(body_holder*,double);
double computeDt(body_holder*,std::vector<body_holder*>&);

} // namespace

#endif // physics_h  
