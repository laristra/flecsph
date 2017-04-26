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

inline bool
operator==(
  const point_t& p1,
  const point_t& p2
)
{
  return p1[0]==p2[0]&&p1[1]==p2[1]&&p1[2]==p2[2];
}

namespace flecsi{
namespace execution{

void
mpi_init_task(/*std::string sfilename*/){
  // TODO find a way to use the file name from the specialiszation_driver
  //std::cout<<sfilename<<std::endl;
  const char * filename = "../data/data_test_40.txt";

  int rank; 
  int size; 
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank); 
  
  int nbodies = 0;
  int totalnbodies = 0;
  std::vector<std::pair<entity_key_t,body>> rbodies;
  
  printf("%d/%d, file %s\n",rank,size,filename);   
  // Read data from file, each process read a part of it 
  // For HDF5, no problems because we know the number of particles 
  // For txt format, work on the number of lines yet ... 
  io::inputDataTxtRange(rbodies,nbodies,totalnbodies,rank,size,filename); 

  std::array<point_t,2> range;
  mpi_compute_range(rbodies,range);

  // The bodies are loaded
  // Compute the key and sort them 
  for(auto& bi: rbodies){
    bi.first = entity_key_t(range,bi.second.coordinates());
  }

  // Check for duplicates keys, particles
  assert(rbodies.end() == 
      std::unique(rbodies.begin(),rbodies.end(),
        [](const auto& left, const auto& right){
          return left.first == right.first;
        }));
 
  // Target number of bodies for every process
  // The last one will takes more 
  std::vector<int> targetnbodies;
  for(int i=0;i<size;++i){
    if(i!=size-1){
      targetnbodies.push_back(totalnbodies/size);
    }else{ 
      targetnbodies.push_back(totalnbodies-((size-1)*(totalnbodies/size)));
    }
  }

  // Apply a distributed sort algorithm 
  mpi_sort(rbodies,targetnbodies);
  assert(rbodies.size() == (size_t)targetnbodies[rank]); 
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
  for(auto bi: rbodies){
    auto nbi = tree.make_entity(bi.second.getPosition(),&bi.second,rank);
    bodies.push_back(nbi); 
    // Add them in the tree
    tree.insert(bi);
  }

  // Search and share the branches 
  mpi_branches_exchange(tree);

  // Do the research of ghost and shared 
  std::vector<body_holder*> ghosts;
  std::vector<body_holder*> shared;
  
  // In this version consider constant H 
  // TODO handle different H, problem the holder does not contain the H 
  // yet, add it or find a way 
  double smoothinglength = bodies[0]->getBody()->getSmoothinglength(); 

  // 1. Search for ghost
  // Do the tree traversal for my bodies
  // The entities that are non local are ghost 
  if(rank==0)
    std::cout<<"Ghost search";
  for(auto bi: bodies)
  {
    auto neighbors = tree.find_in_radius(bi->coordinates(),
        2*smoothinglength);
    for(auto nb: neighbors)
      if(!nb->is_local()){
        nb->setLocality(GHOST);
        ghosts.push_back(nb);
      }
  }
  if(rank==0)
    std::cout<<".done"<<std::endl;


  // 2. Search for shared 
  // Search for the non local entity 
  // The entities found that are LOCAL are shared 
  if(rank==0)
    std::cout<<"Shared search";
  // Get all the entities in the tree
  auto allents = tree.entities();
  for(auto bi: allents)
  {
    // If this entity if not mine, search neighbors
    if(!bi->is_local()){
      auto neighbors = tree.find_in_radius(bi->coordinates(),
        2*smoothinglength);
      // Mark them as shared
    for(auto nb: neighbors)
      if(nb->is_local())
        nb->setLocality(SHARED);
    }
  }
  if(rank==0)
    std::cout<<".done"<<std::endl;
  
  // Display current tree, with GHOSTS,SHARED,EXCL,NONLOCAL
  MPI_Barrier(MPI_COMM_WORLD);
  mpi_tree_traversal_graphviz(tree,range);

  // Index everything

  // Register data and create the final tree ?
  //
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
  flecsi_execute_task(mpi_init_task,mpi,index/*,filename*/); 
} // specialization driver

void 
driver(int argc, char * argv[]){
  std::cout << "In user driver" << std::endl;
} // driver


} // namespace
} // namespace

