#ifndef sodtube_h
#define sodtube_h

#include <vector>

#include "tree.h"


namespace sodtube{
  
  void randomDataSodTube1D(
      std::vector<std::pair<entity_key_t,body>>&,
      int&, int&, int, int);

  void computeDensity(body_holder*,std::vector<body_holder*>&);
  void computeDensityApply(body_holder*,body_holder*);

  void computePressureSoundSpeed(body_holder*);
  void computeViscosity(body_holder*,std::vector<body_holder*>&);
  void moveParticle(body_holder*,std::array<point_t,2>&);
  void computeAcceleration(body_holder*,std::vector<body_holder*>&);


} // namespace sod_tube

#endif // sodtube_h
