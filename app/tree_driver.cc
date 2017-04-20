#ifndef tree_driver_h
#define tree_driver_h

#include <iostream>
#include <numeric> // For accumulate
#include <iostream>

#include <mpi.h>
#include <legion.h>
#include <omp.h>

#include "flecsi/execution/execution.h"
#include "flecsi/data/data_client.h"
#include "flecsi/data/data.h"

#include "physics.h"
#include "io.h"

std::ostream&
operator<<(
  std::ostream& ostr,
  const entity_key_t& id
)
{
  id.output_(ostr);
  return ostr;
}

template<class I>
auto is_contiguous(I first, I last)
{ 
  auto test = true;
  auto const n = std::distance(first, last);
  for (auto i = 0; i < n && test; ++i) {
    test &= (*(std::next(first, i))).first == 
     (*(std::next(std::addressof(*first), i))).first;
  }        
  return test;        
}

// First sort, here a silly bubble 
void mpi_sort(std::vector<std::pair<entity_key_t,body>>& rbodies,
    std::vector<int> targetnbodies)
{
  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  // Normally the pivot position should be selected randomly 
  // Here we hae the assumption that at the begining the particles 
  // are approximetly well distributed then we consider:
  // p0 = 10000-11750
  // p1 = 11750-13500
  // p2 = 13500-15250
  // p3 = 15250-17777
  // In this case we does not need to share the pivots
  // In the random case we need to add a communication step
  std::vector<entity_key_t> pivots; 
  entity_key_t rangeperproc = (entity_key_t::first_key()
      -entity_key_t::last_key())
    /(size);  
  entity_key_t currentkey = entity_key_t::first_key();

  if(rank==0)
    std::cout << "start: "<<currentkey<<std::endl;
  
  for(int i=0;i<size-1;++i)
  {
    currentkey = currentkey + rangeperproc; 
    pivots.push_back(currentkey); 
    if(rank==0)
      std::cout << i << ": "<<currentkey<<std::endl;
  }

  // Then we create buckets based on this pivot 
  std::vector<std::vector<std::pair<entity_key_t,body>>> buckets;
  std::vector<int> bucketssize;
  bucketssize.resize(size);
  std::fill(bucketssize.begin(),bucketssize.end(),0);
  buckets.resize(size);
  
  for(std::vector<std::pair<entity_key_t,body>>::iterator 
      it =  rbodies.begin();
      it != rbodies.end();/*++it*/)
  {
    for(int i=0;i<size;++i){
      if(i == 0)
      {
        if(it->first < pivots[0]){
          buckets[0].push_back(std::move(*it)); 
          bucketssize[0]++;
          rbodies.erase(it);
        }
      }else if(i == size-1){
        // All the others elements should fit here ... test anyway
        if(it->first >= pivots[size-2]){
          buckets[i].push_back(std::move(*it));
          bucketssize[i]++;
          rbodies.erase(it);
        }
      }else{
        if(it->first <= pivots[i] && it->first > pivots[i-1]){
          buckets[i].push_back(std::move(*it));
          bucketssize[i]++;
          rbodies.erase(it);
        }
      }
    }
  }
  assert(rbodies.empty());

  //std::cout<<rank<<":"<<bucketssize[0]<<";"<<bucketssize[1]<<std::endl;

  // The vectors of vector are not contigous in memory
  // Let's append the vectors 
  std::vector<std::pair<entity_key_t,body>> sendbuffer;
  for(auto bucket: buckets)
  {
    std::move(bucket.begin(),bucket.end(),std::back_inserter(sendbuffer));
    bucket.erase(bucket.begin(),bucket.end());
  }
  assert(is_contiguous(sendbuffer.begin(),sendbuffer.end()));

  // Share the bucket 
  // First send the sizes of each bucket 
  std::vector<int> receivbucket;
  receivbucket.resize(size);
  MPI_Alltoall(&bucketssize[0],1,MPI_INT,&receivbucket[0],1,MPI_INT,
      MPI_COMM_WORLD);

  // Create the offset for alltoallv
  //  First receiv offset
  std::vector<int> receivoffsets;
  receivoffsets.resize(size);
  std::partial_sum(receivbucket.begin(),receivbucket.end(),&receivoffsets[0]); 
  // As we need an exscan, add a zero and delete the last element 
  receivoffsets.insert(receivoffsets.begin(),0);
  //receivoffsets.pop_back(); 
  
  // Then send offsets
  std::vector<int> sendoffsets;
  sendoffsets.resize(size);
  std::partial_sum(bucketssize.begin(),bucketssize.end(),&sendoffsets[0]);
  // As we need an exscan, add a zero and delete the last element 
  sendoffsets.insert(sendoffsets.begin(),0);
  //sendoffsets.pop_back();
 
  // The receiv buffer, work in th rbodies
  //std::vector<std::pair<entity_key_t,body>> receivbuffer; 
  rbodies.resize(receivoffsets.back(),
      std::make_pair(entity_key_t::null(),body()));

  // We need to keep the last value to generate rbodies before erasing
  receivoffsets.pop_back(); 
  sendoffsets.pop_back();


  // Create the MPI datatype associate with the pair 
  // First the entity_key_t which is an uint64_t
  int blocks[2] = {1,26};
  MPI_Datatype types[2] = {MPI_UINT64_T,MPI_DOUBLE}; 
  MPI_Aint displacements[2];
  MPI_Aint lb[2];

  MPI_Datatype pair_entity_body;
  MPI_Aint uintex,doublex;

  MPI_Type_get_extent(MPI_UINT64_T,&lb[0],&uintex);
  MPI_Type_get_extent(MPI_DOUBLE,&lb[1],&doublex);
  displacements[0] = static_cast<MPI_Aint>(0);
  displacements[1] = uintex;
  MPI_Type_create_struct(2,blocks,displacements,types,&pair_entity_body);
  MPI_Type_commit(&pair_entity_body);

  int datasize;
  MPI_Type_size(pair_entity_body,&datasize);

  // Use this array for the global buckets communication
  MPI_Alltoallv(&sendbuffer[0],&bucketssize[0],&sendoffsets[0],pair_entity_body,
      &rbodies[0],&receivbucket[0],&receivoffsets[0],pair_entity_body,
      MPI_COMM_WORLD);

  // Sort the incoming buffer 
  sort(rbodies.begin(),rbodies.end(),[](auto& left, auto &right)
  {
    return left.first<right.first;
  }); 

  std::vector<int> totalnbodies(size);
  int mybodies = rbodies.size();
  // Share the final array size of everybody 
  MPI_Allgather(&mybodies,1,MPI_INT,&totalnbodies[0],1,MPI_INT,MPI_COMM_WORLD);

  if(rank == 0)
    std::cout<<"Repartition: "<<totalnbodies[0]<<";"<<totalnbodies[1]<<";"<<totalnbodies[2]<<";"<<totalnbodies[3]<<std::endl;

  if(rank == 0)
    std::cout<<"Target: "<<targetnbodies[0]<<";"<<targetnbodies[1]<<";"<<targetnbodies[2]<<";"<<targetnbodies[3]<<std::endl;

  // Distribution using full right and then full left
  // First full right, normally the worst if size iteration
  for(int i=0;i<size;++i)
  {   
    std::vector<int> needs(size);
    // Compute the current needs  
    for(int i=0;i<size;++i){
      needs[i] = targetnbodies[i]-totalnbodies[i];
    }  
    // Look at the right if someone have the bodies I need 
    if(rank!=0 && needs[rank]>0 && totalnbodies[rank-1]>=needs[rank]){
      int nrecv = needs[rank];
      if(totalnbodies[rank-1]<nrecv)
        nrecv = totalnbodies[rank-1];
      std::vector<std::pair<entity_key_t,body>> recvbuffer(nrecv);
      MPI_Recv(&recvbuffer[0],nrecv,pair_entity_body,rank-1,0,MPI_COMM_WORLD,
          MPI_STATUS_IGNORE);
      // Add and keep sorted
      rbodies.insert(rbodies.end(),recvbuffer.begin(),
          recvbuffer.end());
      std::sort(rbodies.begin(),rbodies.end(),[]
          (auto& left,auto& right)
          {
            return left.first<right.first;
          });
    }
    if(rank!=size-1 && needs[rank+1]>0 && totalnbodies[rank]>=needs[rank+1]){
      int nsend = needs[rank+1];
      if(nsend>totalnbodies[rank])
        nsend = totalnbodies[rank];
      int position = rbodies.size() - nsend;
      MPI_Send(&(rbodies[position]),nsend,pair_entity_body,rank+1,0,
          MPI_COMM_WORLD);
      // Suppress the bodies 
      rbodies.erase(rbodies.end()-nsend,rbodies.end()); 
    }
    // Gather new array size
    mybodies = rbodies.size();
    // Share the final array size of everybody 
    MPI_Allgather(&mybodies,1,MPI_INT,&totalnbodies[0],1,
        MPI_INT,MPI_COMM_WORLD);
  }

  // Now go to left  
  while(!std::equal(totalnbodies.begin(),totalnbodies.end(),
        targetnbodies.begin()))
  {
    std::vector<int> needs(size);
    // Compute the current needs  
    for(int i=0;i<size;++i){
      needs[i] = targetnbodies[i]-totalnbodies[i];
    }  
    // Look at the right if someone have the bodies I need 
    if(rank!=size-1 && needs[rank]>0 && totalnbodies[rank+1]>=needs[rank]){
      int nrecv = needs[rank];
      if(totalnbodies[rank+1]<nrecv)
        nrecv = totalnbodies[rank+1];
      std::vector<std::pair<entity_key_t,body>> recvbuffer(nrecv);
      MPI_Recv(&recvbuffer[0],nrecv,pair_entity_body,rank+1,0,MPI_COMM_WORLD,
          MPI_STATUS_IGNORE);
      // Add and keep sorted
      rbodies.insert(rbodies.end(),recvbuffer.begin(),
          recvbuffer.end());
      std::sort(rbodies.begin(),rbodies.end(),[]
          (auto& left,auto& right)
          {
            return left.first<right.first;
          });
    }
    if(rank!=0 && needs[rank-1]>0 && totalnbodies[rank]>=needs[rank-1]){
      int nsend = needs[rank-1];
      if(nsend>totalnbodies[rank])
        nsend = totalnbodies[rank];
      int position = rbodies.size()-nsend;
      MPI_Send(&rbodies[position],nsend,pair_entity_body,rank-1,0,
          MPI_COMM_WORLD);
      // Suppress the bodies 
      rbodies.erase(rbodies.end()-nsend,rbodies.end()); 
    }
    // Gather new array size
    mybodies = rbodies.size();
    // Share the final array size of everybody 
    MPI_Allgather(&mybodies,1,MPI_INT,&totalnbodies[0],1,
        MPI_INT,MPI_COMM_WORLD);
  }
  if(rank == 0){
    std::cout<<"Repartition: ";
    for(auto num: totalnbodies)
      std::cout<<num<<";";
    std::cout<<std::endl;
  }
}

namespace flecsi{
namespace execution{

void
mpi_task(/*const char * filename*/){
  const char * filename = "../data/data_binary_8338.txt";
  //std::vector<body*> rbodies; // Body read by the process

  int rank; 
  int size; 
  int nbodies = 0;
  int totalnbodies = 0;
  tree_topology_t t;
  std::vector<std::pair<entity_key_t,body>> rbodies;

  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank); 
  printf("%d/%d, file %s\n",rank,size,filename);   
  // Read data from file, each process read a part of it 
  // For HDF5, no problems because we know the number of particles 
  // For txt format, work on the number of lines yet ... 
  io::inputDataTxtRange(rbodies,nbodies,totalnbodies,rank,size,filename); 
  //std::cout << "Read done" << std::endl;

  // Compute the range to compute the keys 
  double max[3] = {-9999,-9999,-9999};
  double min[3] = {9999,9999,9999};
  for(auto bi: rbodies){
    for(int i=0;i<3;++i){
        if(bi.second.coordinates()[i]>max[i])
          max[i] = bi.second.coordinates()[i];
        if(bi.second.coordinates()[i]<min[i])
          min[i] = bi.second.coordinates()[i];
      }
  }

  // Do the MPI Reduction 
  MPI_Allreduce(MPI_IN_PLACE,max,3,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD); 
  MPI_Allreduce(MPI_IN_PLACE,min,3,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD); 
 
  point_t minposition = {min[0],min[1],min[2]};
  point_t maxposition = {max[0],max[1],max[2]};

  if(rank==0)
    std::cout <<"boundaries: "<< minposition << maxposition << std::endl;

  std::array<point_t,2> range = {minposition,maxposition};

  // The bodies are loaded
  // Compute the key and sort them 
  for(auto& bi: rbodies){
    bi.first = entity_key_t(range,bi.second.coordinates());
    //std::cout << bi.first << std::endl;
  }
  // Sort based on keys 
  std::sort(rbodies.begin(),rbodies.end(), [](auto &left,auto &right){
      return left.first < right.first; 
  });

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

  // Get the first and last key, send to everyone
  //std::array<entity_key_t,2> locboundaries; 
  //locboundaries[0] = rbodies.front().first;
  //locboundaries[1] = rbodies.back().first; 
  //std::cout <<rank<<"/"<<size<<" f: "<<locboundaries[0]<<" l: "
  //  <<locboundaries[1]<<std::endl;

  // Share the boundary keys, here we know it is a 64 unsigned int 
  // Make a vector for the right number of pairs 
  //vector<entity_key_t> globalboundaries;

  //MPI_Allgather(&locboundaries[0],2,);  

  // Create one vector per process with its keys 

  // Send the informations vio MPI_Alltoallv

  // Generate the local tree 

  // Do the research of ghost and shared 
  
  // Index everything

  // Register data and create the final tree ?
  MPI_Barrier(MPI_COMM_WORLD);

}

flecsi_register_task(mpi_task,mpi,index);

void 
specialization_driver(int argc, char * argv[]){
  if (argc!=2) {
    std::cerr << "Error not enough arguments\n"
        "Usage: tree <datafile>\n";
    exit(-1); 
  }

  std::cout << "In user specialization_driver" << std::endl;
  /*const char * filename = argv[1];*/

  flecsi_execute_task(mpi_task,mpi,index/*,filename*/); 
} // specialization driver

void 
driver(int argc, char * argv[]){
  std::cout << "In user driver" << std::endl;
} // driver


} // namespace
} // namespace


#endif // tree_driver_h
