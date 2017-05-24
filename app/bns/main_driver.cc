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

#include "mpi_partition.h"
#include "physics.h"
#include "io.h"
#include "sodtube.h"

inline bool
operator==(
  const point_t& p1,
  const point_t& p2
)
{
  for(size_t i=0;i<gdimension;++i)
    if(p1[i]!=p2[i])
      return false;
  return true;
}

inline 
point_t 
operator+(
  const point_t& p, 
  const double& val)
{
  point_t pr = p; 
  for(size_t i=0;i<gdimension;++i)
    pr[i]+=val;
  return pr;
}

inline 
point_t 
operator-(
  const point_t& p, 
  const double& val)
{
  point_t pr = p; 
  for(size_t i=0;i<gdimension;++i)
    pr[i]-=val;
  return pr;
}

namespace flecsi{
namespace execution{

void
mpi_init_task(/*std::string sfilename*/int inputparticles){
  // TODO find a way to use the file name from the specialiszation_driver
  //std::cout<<sfilename<<std::endl;
  const char * filename = "../data/data_bns_4169.txt";
  //const char * filename = "../data/data_binary_rdy_16288.txt";
  
  int rank;
  int size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  
  int nbodies = 0;
  int totaliters = 100000;
  int iteroutput = 10;
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
  
  //printf("%d/%d, file %s\n",rank,size,filename);   
  // Read data from file, each process read a part of it 
  // For HDF5, no problems because we know the number of particles 
  // For txt format, work on the number of lines yet ... 
  io::inputDataTxtRange(rbodies,nbodies,totalnbodies,rank,size,filename); 

  double smoothinglength = 0.0;
  int iter = 0;

#ifdef OUTPUT
  mpi_output_txt(rbodies,iter); 
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

    mpi_compute_range(rbodies,range,smoothinglength);

    // The bodies are loaded
    // Compute  the key    
    for(auto&  bi: rbodies){
      bi.first = entity_key_t(range,bi.second.coordinates());
    }
    // Check for duplicates keys, particles
    //assert(rbodies.end() == 
    //  std::unique(rbodies.begin(),rbodies.end(),
    //    [](const auto& left, const auto& right){
    //      return left.first == right.first;
    //    }));

    // Apply a distributed sort algorithm 
    mpi_sort_unbalanced(rbodies,totalnbodies);
    //assert(rbodies.size() == (size_t)targetnbodies[rank]); 
    assert(rbodies.end() == std::unique(rbodies.begin(),rbodies.end(),
        [](const auto& left, const auto& right){
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
    mpi_branches_exchange_useful_positions(
        tree,
        rbodies,
        rangeposproc,
        smoothinglength);

    // For test, gather the neighbors and loop over here
    mpi_compute_ghosts(tree,smoothinglength,ghosts_data,range);
    mpi_refresh_ghosts(tree,ghosts_data,range);

    if(size==1)
      assert(ghosts_data.recvbodies.size() == 0);

#ifdef OUTPUTGRAPH
    if(iter==1){
      // Display current tree, with GHOSTS,SHARED,EXCL,NONLOCAL
      MPI_Barrier(MPI_COMM_WORLD);
      mpi_tree_traversal_graphviz(tree,range);
    }
#endif

    // Do the BNS Physics 
    
    // Compute Density
    if(rank == 0) 
      std::cout<<"Density"<<std::flush;
    for(auto& bi: bodies)
    {
      auto ents = tree.find_in_radius(bi->coordinates(),2*smoothinglength);
      auto vecents = ents.to_vec();
      physics::computeDensity(bi,vecents);
    }
    if(rank == 0)
      std::cout<<".done"<<std::endl<<std::flush;
  
    // Compute Pressure and SoundSpeed
    if(rank == 0)
      std::cout<<"Pressure"<<std::flush;
    for(auto& bi: bodies)
    {
      physics::computePressure(bi);
      physics::computeSoundspeed(bi);
    }
    if(rank==0)
      std::cout<<".done"<<std::endl<<std::flush;

    // Exchange data
    mpi_refresh_ghosts(tree,ghosts_data,range);

    // Compute Hydro 
    if(rank==0)
      std::cout<<"Hydro"<<std::flush;
    for(auto bi: bodies)
    {
      auto ents = tree.find_in_radius(bi->coordinates(),2*smoothinglength);
      auto vecents = ents.to_vec();
      physics::computeHydro(bi,vecents);
    }
    if(rank==0)
      std::cout<<".done"<<std::endl<<std::flush;


    if(size==1)
      assert(tree.entities().size() == (size_t)totalnbodies);
    // Compute and exchange the aggregate  
    if(rank==0)
      std::cout<<"Grav"<<std::flush;
    
#if 1
    // compute the COM data
    tree_traversal_com(tree);
    // Gather all the COM
    std::vector<mpi_cell> recvcells;
    std::vector<int> nrecvcells;
    mpi_exchange_cells(tree,recvcells,nrecvcells,mcell);
    mpi_compute_fmm(tree,recvcells,macangle);
    mpi_gather_cells(tree,recvcells,nrecvcells);
    mpi_gather_ghosts_com(tree,recvcells,nrecvcells,range);
    MPI_Barrier(MPI_COMM_WORLD);
#else 
    int nbody = 0;
    tree_traversal_com(tree);
    tree_traversal_grav(tree,tree.root(),mcell,macangle,nbody);
    assert(nbody == totalnbodies);
#endif

    if(rank==0)
      std::cout<<".done"<<std::endl<<std::flush;
    
    // Compute Acceleration 
    if(rank==0)
      std::cout<<"Acceleration"<<std::flush;
    for(auto bi: bodies)
    {
      physics::computeAcceleration(bi,physics::dt);
    }
    if(rank==0)
      std::cout<<".done"<<std::endl<<std::flush;

    // Compute the new DT
    physics::dt = 1;
    if(rank==0)
      std::cout<<"Dt"<<std::flush;
    for(auto bi: bodies)
    {
      // Compute the local dt for each body 
      // Use only the neighbors in the smoothing length 
      auto ents = tree.find_in_radius(bi->coordinates(),2*smoothinglength);
      auto vecents = ents.to_vec();
      double localdt = physics::computeDt(bi,vecents);
      physics::dt = std::min(physics::dt,localdt);

    }

    // Do a reduction for dt
    MPI_Allreduce(MPI_IN_PLACE,&physics::dt,1,MPI_DOUBLE,
        MPI_MIN,MPI_COMM_WORLD);
    //physics::dt = 1.0e-4;
    if(rank==0)
      std::cout<<".done"<<std::endl<<std::flush;
    assert(physics::dt<1&&physics::dt>0);

    if(rank==0)
      std::cout<<"dt="<<physics::dt<<std::endl<<std::flush;

    // Move the particles 
    if(rank==0)
      std::cout<<"MoveParticles"<<std::flush;    
    for(auto bi: bodies)
    {
      physics::moveBody(bi,physics::dt);
    }
    if(rank==0)
      std::cout<<".done"<<std::endl<<std::flush;

    // Add rotation 


#ifdef OUTPUT
    if(iter % iteroutput == 0)
    {
      mpi_output_txt(rbodies,iter/iteroutput);
    }
    // output to see repartition
    /*char fname[64];
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
    totaltime += physics::dt; 
    if(rank==0)
      std::cout<<"Time: "<<totaltime<<std::endl;

  }while(totaltime<maxtime);
}

flecsi_register_task(mpi_init_task,mpi,index);

void 
specialization_driver(int argc, char * argv[]){
  if (argc!=2) {
    std::cerr << "Error not enough arguments\n"
        "Usage: tree <datafile>\n";
    exit(-1); 
  }

  std::cout << "In user specialization_driver" << std::endl;
  /*const char * filename = argv[1];*/
  /*std::string filename(argv[1]);
  std::cout<<filename<<std::endl;*/
  flecsi_execute_task(mpi_init_task,mpi,index,atoi(argv[1])/*,filename*/); 
} // specialization driver

void 
driver(int argc, char * argv[]){
  std::cout << "In user driver" << std::endl;
} // driver


} // namespace
} // namespace


