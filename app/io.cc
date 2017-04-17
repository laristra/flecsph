#include "io.h"

namespace io{

// Read data from txt file 
// This is not the most beautiful way, but enough for testing
void inputDataTxt(
    std::vector<body*>& bodies, 
    const char * filename,
    tree_topology_t& t
  ){
  // Load data from text file and setup bodies 
  FILE * inputfile = fopen(filename,"r");
  size_t nbodies = 0;
  double x, y, z, velX, velY, velZ, velHX, velHY, velHZ;
  double accX, accY, accZ, smoothinglength, pressure;
  double entropy, density, mass, tmp, angularMoment;
  int tmpD;
  while(fscanf(inputfile,
      "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf"
      " %lf %lf %lf %lf %lf %lf %d %lf %lf %lf %lf %lf \n",
      &x,&y,&z,                       // positions
      &velX,&velY,&velZ,              // velocity
      &velHX,&velHY,&velHZ,           // velHalf 
      &accX,&accY,&accZ,              // acceleration
      &density,&pressure,&entropy,    // density, pressure, entropy 
      &mass,&smoothinglength,&tmp,    // mass, smoothing length, dt 
      &tmpD,&tmp,&tmp,                // step, totalTime, tmp
      &tmp,&tmp,&angularMoment        // tmp, tmp, angularMoment
  )!=EOF){
    point_t position = {x,y,z}; 
    point_t velocity = {velX,velY,velZ};
    point_t velocityhalf = {velHX,velHY,velHZ};
    point_t acceleration = {accX, accY, accZ};
    auto bi = t.make_entity(position,velocity,velocityhalf,
        acceleration,density,pressure,entropy,mass,smoothinglength);   
    bodies.push_back(bi);
    ++nbodies;
    if(nbodies < 5)
      std::cout << *bi << std::endl;
  } // while 
} // inputDataTxt

// Input data fro HDF5 File
void inputDataHDF5(
    std::vector<body*> bodies,
    const char * filename,
    tree_topology_t& t){
  // Data double buffer to read 
  double* data;
  hid_t file; 
  hid_t dset; 
  herr_t status; 
  size_t nbodies; 

  // Open the file 
  file = H5Fopen(filename,H5F_ACC_RDONLY,H5P_DEFAULT);
  
  // First get the number of bodies (Here size of the point array in hdf5)
  dset = H5Dopen(file,"nParticles",H5P_DEFAULT);
  //status = 

}// inputDataHDF5

// Output data in txt format 
void outputDataTxt(){
  
}// ouputDataTxt

// Output data in HDF5 format 
// Generate the associate XDMF file 
void outputDataHDF5(
    std::vector<body*>& bodies,
    const char* filename,
    int step,
    double dt
    ){
  hsize_t size = bodies.size(); 
  int nbodies = bodies.size(); 
  char outputfilename[128];
  hid_t dset, file, space;
  herr_t status; 
  int cpt;
  double *data;  
  data = (double*)malloc(sizeof(double)*size);
  sprintf(outputfilename,"%s%d.h5",filename,step); 
  file = H5Fcreate(outputfilename,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
  space = H5Screate_simple(1,&size,NULL);  
  
  // Positions 
  // X
  cpt = 0;
  for(auto bi: bodies)
    data[cpt++] =  bi->getPosition()[0];
  dset = H5Dcreate(file,"/X",H5T_NATIVE_DOUBLE,space,
      H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT); 
  status = H5Dwrite(dset,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  status = H5Dclose(dset);  
  // Y
  cpt = 0;
  for(auto bi: bodies)
    data[cpt++] =  bi->getPosition()[1];
  dset = H5Dcreate(file,"/Y",H5T_NATIVE_DOUBLE,space,
      H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT); 
  status = H5Dwrite(dset,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  status = H5Dclose(dset); 
  // Z
  cpt = 0;
  for(auto bi: bodies)
    data[cpt++] =  bi->getPosition()[2];
  dset = H5Dcreate(file,"/Z",H5T_NATIVE_DOUBLE,space,
      H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT); 
  status = H5Dwrite(dset,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  status = H5Dclose(dset);

  // Density 
  cpt = 0;
  for(auto bi: bodies)
    data[cpt++] =  bi->getDensity();
  dset = H5Dcreate(file,"/density",H5T_NATIVE_DOUBLE,space,
      H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT); 
  status = H5Dwrite(dset,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  status = H5Dclose(dset);

  // Time 
  status = H5Sclose(space); 
  size = 1;
  space = H5Screate_simple(1,&size,NULL); 
  dset = H5Dcreate(file,"/time",H5T_NATIVE_DOUBLE,space,
      H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  status = H5Dwrite(dset,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  status = H5Dclose(dset);  

  status = H5Sclose(space);
  status = H5Fclose(file);  
  free(data);
}// outputDataHDF5

} // namespace io


