/*~--------------------------------------------------------------------------~*
 * 
 * /@@@@@@@@  @@           @@@@@@   @@@@@@@@ @@@@@@@  @@      @@
 * /@@/////  /@@          @@////@@ @@////// /@@////@@/@@     /@@
 * /@@       /@@  @@@@@  @@    // /@@       /@@   /@@/@@     /@@
 * /@@@@@@@  /@@ @@///@@/@@       /@@@@@@@@@/@@@@@@@ /@@@@@@@@@@
 * /@@////   /@@/@@@@@@@/@@       ////////@@/@@////  /@@//////@@
 * /@@       /@@/@@//// //@@    @@       /@@/@@      /@@     /@@
 * /@@       @@@//@@@@@@ //@@@@@@  @@@@@@@@ /@@      /@@     /@@
 * //       ///  //////   //////  ////////  //       //      //  
 *
 *~--------------------------------------------------------------------------~*/

/**
 * @file main_driver.cc
 * @author Julien Loiseau
 * @date April 2017
 * @brief Specialization and Main driver used in FleCSI. 
 * The Specialization Driver is normally used to register data and the main 
 * code is in the Driver.  
 */

#include <iostream>
#include <numeric> // For accumulate
#include <iostream>

#include <mpi.h>
#include <legion.h>
#include <omp.h>

#include "flecsi/execution/execution.h"
#include "flecsi/data/data_client.h"
#include "flecsi/data/data.h"

#include <tree_colorer.h>
#include "io.h"
#include "sodtube.h"

namespace flecsi{
namespace execution{

void
mpi_init_task(int inputparticles){
  // TODO find a way to use the file name from the specialiszation_driver
  
  int rank;
  int size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  
  int nbodies = 0;
  int totaliters = 80;
  int iteroutput = 1;
  int totalnbodies = 0;
  double totaltime = 0.0;
  double maxtime = 10.0;
  double macangle = 0.5;
  double mcell = 1.0e-6;
  std::vector<std::pair<entity_key_t,body>> rbodies;
  std::array<point_t,2> range;
  std::vector<std::pair<entity_key_t,entity_key_t>> rangeproc;
  std::vector<std::pair<point_t,point_t>> rangeposproc;
  mpi_ghosts_t ghosts_data;
  std::vector<body_holder> recv_COM;
  
  tree_colorer<double,1> tcolorer;

  //printf("%d/%d, file %s\n",rank,size,filename);   
  // Read data from file, each process read a part of it 
  // For HDF5, no problems because we know the number of particles 
  // For txt format, work on the number of lines yet ... 
  totalnbodies = inputparticles;
  sodtube::randomDataSodTube1D(rbodies,nbodies,totalnbodies,rank,size);

  double smoothinglength = 0.0;
  int iter = 0;

#ifdef OUTPUT
  tcolorer.mpi_output_txt(rbodies,iter,"output_sodtube"); 
#endif

  ++iter; 
  do
  { 
    // Get the worst smothing length 
    if(rank==0)
      smoothinglength = rbodies[0].second.getSmoothinglength();

    MPI_Bcast(&smoothinglength,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    if(rank==0)
      std::cout<<"H="<<smoothinglength<<std::endl;

    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==0)
      std::cout<<std::endl<<"#### Iteration "<<iter<<std::endl;
    MPI_Barrier(MPI_COMM_WORLD);

    tcolorer.mpi_compute_range(rbodies,range,smoothinglength);

    // The bodies are loaded
    // Compute the key    
    for(auto& bi: rbodies){
      bi.first = entity_key_t(range,bi.second.coordinates());
    }
    // Check for duplicates keys, particles
    //assert(rbodies.end() == 
    //  std::unique(rbodies.begin(),rbodies.end(),
    //    [](const auto& left, const auto& right){
    //      return left.first == right.first;
    //    }));

    // Apply a distributed sort algorithm 
    tcolorer.mpi_qsort(rbodies,totalnbodies);
    //assert(rbodies.size() == (size_t)targetnbodies[rank]); 
    assert(rbodies.end() == std::unique(rbodies.begin(),rbodies.end(),
        [] (const auto& left, const auto& right){
          return left.second.coordinates()==right.second.coordinates() &&
            left.first == right.first;
        }));

    // Generate the local tree 
    tree_topology_t tree(range[0],range[1]);

    // create the bodies holders for my local bodies
    // They will be registred as exclusive
    std::vector<body_holder*> bodies;
    double lowmass = 10.0;
    for(auto& bi: rbodies){
      auto nbi = tree.make_entity(bi.second.getPosition(),&(bi.second),rank,
          bi.second.getMass());
      // Add them in the tree
      tree.insert(nbi);
      bodies.push_back(nbi);
    } 
    // Search and share the branches 
    //mpi_branches_exchange_useful(tree,rbodies,range,rangeproc);
    tcolorer.mpi_branches_exchange(
        tree,
        rbodies,
        rangeposproc,
        smoothinglength);

    // For test, gather the neighbors and loop over here
    tcolorer.mpi_compute_ghosts(tree,smoothinglength,range);
    tcolorer.mpi_refresh_ghosts(tree,range);

    //if(size==1)
    //  assert(ghosts_data.recvbodies.size() == 0);

#ifdef OUTPUTGRAPH
    if(iter==1){
      // Display current tree, with GHOSTS,SHARED,EXCL,NONLOCAL
      MPI_Barrier(MPI_COMM_WORLD);
      tcolorer.mpi_tree_traversal_graphviz(tree,range);
    }
#endif

    // Do the Sod Tube physics
    if(rank==0)
      std::cout<<"Density"<<std::flush; 
    for(auto& bi: bodies)
    {
      bi->getBody()->setDensity(0.0);
      tree.apply_in_radius(bi->coordinates(),2*smoothinglength,
          sodtube::computeDensityApply,bi);
    }
    if(rank==0)
      std::cout<<".done"<<std::endl;

    if(rank==0)
      std::cout<<"PressureSoundSpeed"<<std::flush; 
    for(auto& bi: bodies)
    {
      sodtube::computePressureSoundSpeed(bi);
    }
    if(rank==0)
      std::cout<<".done"<<std::endl;

    MPI_Barrier(MPI_COMM_WORLD);
    // Share data again, need the density and pressure of the ghosts 
    tcolorer.mpi_refresh_ghosts(tree,range);    

    if(rank==0)
      std::cout<<"Acceleration"<<std::flush; 
    for(auto& bi: bodies)
    {
      auto ents = tree.find_in_radius(bi->coordinates(),2*smoothinglength);
      // Apply physics
      auto vecents = ents.to_vec(); 
      sodtube::computeAcceleration(bi,vecents);
    }
    if(rank==0)
      std::cout<<".done"<<std::endl;

    if(rank==0)
      std::cout<<"Viscosity"<<std::flush; 
    for(auto& bi: bodies)
    {
      auto ents = tree.find_in_radius(bi->coordinates(),2*smoothinglength);
      // Apply physics
      auto vecents = ents.to_vec(); 
      sodtube::computeViscosity(bi,vecents);
    }
    if(rank==0)
      std::cout<<".done"<<std::endl;

    if(rank==0)
      std::cout<<"MoveParticles"<<std::flush; 
    for(auto& bi: bodies)
    {
      sodtube::moveParticle(bi,range);
    }
    if(rank==0)
      std::cout<<".done"<<std::endl;

#ifdef OUTPUT
    if(iter % iteroutput == 0)
    { 
      tcolorer.mpi_output_txt(rbodies,iter/iteroutput,"output_sodtube");
    }
    // output to see repartition
    /*char fn ame[64];
    sprintf(fname,"part_%d_%d.txt",rank,iter); 
    FILE * tmpf = fopen(fname,"w");
    fprintf(tmpf,"X Y Z\n");
    for(auto bi: rbodies)
    {
      fprintf(tmpf,"%.4f %.4f %.4f\n",bi.second.getPosition()[0],
          bi.second.getPosition()[1],bi.second.getPosition()[2]); 
    }
    fclose(tmpf);*/
#endif
    ++iter;
    
  }while(iter<totaliters);
}

flecsi_register_task(mpi_init_task,mpi,index);

void 
specialization_driver(int argc, char * argv[]){
  if (argc!=2)  {
    std::cerr << "Error not enough arguments\n"
        "Usage: tree <datafile>\n";
    exit(-1); 
  }

  std::cout << "In user specialization_driver" << std::endl;
  /*const char * filename = argv[1];*/
  /*std::string  filename(argv[1]);
  std::cout<<filename<<std::endl;*/
  flecsi_execute_task(mpi_init_task,mpi,index,atoi(argv[1])/*,filename*/); 
} // specialization driver

void 
driver(int argc,  char * argv[]){
  std::cout << "In user driver" << std::endl;
} // driver


} // namespace
} // namespace


