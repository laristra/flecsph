
#include <cinchdevel.h>
#include <cinchtest.h>

#include <iostream>
#include <cmath>
#include <mpi.h>

#include "tree_fmm.h"

using namespace ::testing;

namespace flecsi{
  namespace execution{
    void driver(int argc, char* argv[]){
    }
  }
}

TEST(tree_colorer, mpi_qsort){
  int rank;
  int size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  srand(time(NULL)*rank);
  tree_fmm<double,gdimension> tfmm;

  point_t sinkPosition(0.,0.,0.);
  point_t sourcePosition(1.,1.,1.);
  double sourceMass = 1.0;
  point_t fc = {};
  double jacobi[9]={};
  double hessian[27]={};

  rank|| clog(info) << "test pSink=" << sinkPosition
                 << " pSource="   << sourcePosition << std::endl;

  point_t fc_res(1./(3*sqrt(3)),1./(3*sqrt(3)),1./(3*sqrt(3)));

  double jacobi_res[9] = {
                        0.,1./(3*sqrt(3)),1./(3*sqrt(3)),
                        1./(3*sqrt(3)),0.,1./(3*sqrt(3)),
                        1./(3*sqrt(3)),1./(3*sqrt(3)),0.
                      };
  double hessian_res[27] = {
                       -0.2566,0.1283,0.1283,
                       0.1283,0.1283,0.32075,
                       0.1283,0.32075,0.1283,

                       0.1283,0.1283,0.32075,
                       0.1283,-0.2566,0.1283,
                       0.32075,0.1283,0.1283,

                       0.1283,0.32075,0.1283,
                       0.32075,0.1283,0.1283,
                       0.1283,0.1283,-0.2566
                      };

  tfmm.computeAcceleration(sinkPosition,sourcePosition,sourceMass,fc,
    jacobi,hessian);

  for(int i = 0; i<27; ++i){
    if(i<3){
      rank|| clog(info) << "fc i=" << i << ": code=" << fc[i] 
                     << " res=" << fc_res[i]<<std::endl;
      ASSERT_TRUE(fabs(fc[i]-fc_res[i]) < 1.0e-10);
    }
    if(i<9){
      rank|| clog(info) << "dfcdr i=" << i << ": code=" << jacobi[i] 
                     << " res=" << jacobi_res[i] << std::endl;
      ASSERT_TRUE(fabs(jacobi[i]-jacobi_res[i]) < 1.0e-10);
    }
    rank|| clog(info) << "dfcdrdr i=" << i << ": code=" << hessian[i]
                   << " res=" << hessian_res[i] << std::endl;
    ASSERT_TRUE(fabs(hessian[i]-hessian_res[i]) < 1.0e-4);

  }

  sinkPosition = point_t(1.,1.,1.);
  sourcePosition = point_t(0.,0.,0.);
  sourceMass = 1.0;
  fc = point_t{};
  memset(jacobi,0.,9*sizeof(double));
  memset(hessian,0.,sizeof(double)*27);


  rank|| clog(info) << std::endl;
  rank|| clog(info) << "test pSink=" << sinkPosition 
                 << " pSource=" << sourcePosition << std::endl;


  point_t fc_res_2 = point_t(-1./(3*sqrt(3)),-1./(3*sqrt(3)),-1./(3*sqrt(3)));

  double jacobi_res_2[9] = {
                  0.,1./(3*sqrt(3)),1./(3*sqrt(3)),
                  1./(3*sqrt(3)),0.,1./(3*sqrt(3)),
                  1./(3*sqrt(3)),1./(3*sqrt(3)),0.
                };

  double hessian_res_2[27] = {
                  0.2566,-0.1283,-0.1283,
                  -0.1283,-0.1283,-0.32075,
                  -0.1283,-0.32075,-0.1283,

                  -0.1283,-0.1283,-0.32075,
                  -0.1283,0.2566,-0.1283,
                  -0.32075,-0.1283,-0.1283,

                  -0.1283,-0.32075,-0.1283,
                  -0.32075,-0.1283,-0.1283,
                  -0.1283,-0.1283,0.2566
                };

  tfmm.computeAcceleration(sinkPosition,sourcePosition,sourceMass,fc,
    jacobi,hessian);

  for(int i = 0; i<27; ++i){
    if(i<3){
      rank|| clog(info) << "fc i=" << i << ": code=" << fc[i]
                      << " res=" << fc_res_2[i] << std::endl;
      ASSERT_TRUE(fabs(fc[i]-fc_res_2[i]) < 1.0e-4);
    }
    if(i<9){
      rank|| clog(info) << "dfcdr i=" << i << ": code=" << jacobi[i] 
                     << " res=" << jacobi_res_2[i] << std::endl;
      ASSERT_TRUE(fabs(jacobi[i]-jacobi_res_2[i]) < 1.0e-4);
    }
    rank|| clog(info) << "dfcdrdr i=" << i << ": code=" << hessian[i]
                   << " res=" << hessian_res_2[i] << std::endl;
    ASSERT_TRUE(fabs(hessian[i]-hessian_res_2[i]) < 1.0e-4);
  }


  sinkPosition = point_t(-0.1,-0.1,-0.1);
  sourcePosition = point_t(-0.5,-0.5,-0.5);
  sourceMass = 1.0;
  fc = point_t{};
  memset(jacobi,0.,9*sizeof(double));
  memset(hessian,0.,sizeof(double)*27);


  rank|| clog(info)<< std::endl
                << "test pSink=" << sinkPosition
                << " pSource="   << sourcePosition<<std::endl;


  point_t fc_res_3 = point_t(-1.20281,-1.20281,-1.20281);

  double jacobi_res_3[9] = {
                  0.,3.00703,3.00703,
                  3.00703,0.,3.00703,
                  3.00703,3.00703,0.
                };

  double hessian_res_3[27] = {
                  10.0234,-5.01172,-5.01172,
                  -5.01172,-5.01172,-12.5293,
                  -5.01172,-12.5293,-5.01172,

                  -5.01172,-5.01172,-12.5293,
                  -5.01172,10.0234,-5.01172,
                  -12.5293,-5.01172,-5.01172,

                  -5.01172,-12.5293,-5.01172,
                  -12.5293,-5.01172,-5.01172,
                  -5.01172,-5.01172,10.0234
                };

  tfmm.computeAcceleration(sinkPosition,sourcePosition,sourceMass,fc,
    jacobi,hessian);

    for(int i = 0; i<27; ++i){
    if(i<3){
      rank|| clog(info)<<"fc i="<<i<<": code="<<fc[i]
                    <<" res="<<fc_res_3[i]<<std::endl;
      ASSERT_TRUE(fabs(fc[i]-fc_res_3[i]) < 1.0e-4);
    }
    if(i<9){
      rank|| clog(info)<<"dfcdr i="<<i<<": code="<<jacobi[i]
                    <<" res="<<jacobi_res_3[i]<<std::endl;
      ASSERT_TRUE(fabs(jacobi[i]-jacobi_res_3[i]) < 1.0e-4);
    }
    rank|| clog(info)<<"dfcdrdr i="<<i<<": code="<<hessian[i]
                  <<" res="<<hessian_res_3[i]<<std::endl;
    ASSERT_TRUE(fabs(hessian[i]-hessian_res_3[i]) < 1.0e-4);
  }
}
