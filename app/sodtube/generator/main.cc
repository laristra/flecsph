#include "hdf5SimIO.h"

int main(int argc, char * argv[]){

  if(argc!=2){
    printf("./sodtube_generator [nParticles]\n");
    exit(-1);
  }

  int rank, size; 
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  int nparticles = atoi(argv[1]);
  int nparticlesproc = nparticles/size;
  if(rank==size-1){
    nparticlesproc = nparticles - nparticlesproc*(size-1);
  }

  if(rank==0){
    printf("Generating %d particles\n",nparticles);
  }

  double* x = new double[nparticlesproc];
  double* y = new double[nparticlesproc];
  double* z = new double[nparticlesproc];
  double* pressure = new double[nparticlesproc];


 
  Flecsi_Sim_IO::HDF5SimIO testDataSet; 
  testDataSet.setFilename("test.h5");
  testDataSet.createDataset(Flecsi_Sim_IO::unStructuredGrid,4,3);

  Flecsi_Sim_IO::Variable _x,_y,_z,_pressure;

  _x.createVariable("x",Flecsi_Sim_IO::point,"double",nparticlesproc,x);
  _y.createVariable("y",Flecsi_Sim_IO::point,"double",nparticlesproc,y);
  _z.createVariable("z",Flecsi_Sim_IO::point,"double",nparticlesproc,z);
  _pressure.createVariable("pressure",Flecsi_Sim_IO::point,"double",nparticlesproc,pressure);

  testDataSet.vars.push_back(_x);
  testDataSet.vars.push_back(_y);
  testDataSet.vars.push_back(_z);
  testDataSet.vars.push_back(_pressure);

  testDataSet.writePointData(0,MPI_COMM_WORLD);

  
  MPI_Finalize();
  return 0;
}
