#include "mpi_partition.h"

#include <fstream>
#include <iostream>
  
std::ostream&
operator<<(
  std::ostream& ostr,
  const entity_key_t& id
)
{
  id.output_(ostr);
  return ostr;
}

inline 
bool
operator==(
    const point_t& p1, 
    const point_t& p2)
{
  for(size_t i=0;i<gdimension;++i)
    if(p1[i]!=p2[i])
      return false;
  return true;
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

// This method sort the particles over all the processes
// Then proceed to balancing the particles over the processes
void mpi_sort(std::vector<std::pair<entity_key_t,body>>& rbodies,
    std::vector<int> targetnbodies)
{
  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  std::sort(rbodies.begin(),rbodies.end(),[](auto& left, auto& right){
    return left.first == right.first;    
  });

  // If there is just one process, skip the distribution
  if(size==1)
    return;
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
  //std::cout<<entity_key_t::first_key()<<std::endl; 
  //std::cout<<entity_key_t::last_key()<<std::endl;
  entity_key_t currentkey = entity_key_t::first_key();

  //if(rank==0)
  //  std::cout << "start: "<<currentkey<<std::endl;
  
  for(int i=0;i<size-1;++i)
  {
    currentkey = currentkey + rangeperproc; 
    pivots.push_back(currentkey); 
    //if(rank==0)
    //  std::cout << i << ": "<<currentkey<<std::endl;
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
  int blocks[2] = {1,14};
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

  int typesize=0;
  MPI_Type_size(pair_entity_body,&typesize);
  //std::cout<<sizeof(std::pair<entity_key_t,body>)<<":"<<typesize<<std::endl;
  assert(sizeof(std::pair<entity_key_t,body>)==typesize); ;

  // Trnaform the offsets for bytes 
 /* for(auto& bs: bucketssize)
    bs*=sizeof(std::pair<entity_key_t,body>);
  for(auto& bs: sendoffsets)
    bs*=sizeof(std::pair<entity_key_t,body>);
  for(auto& bs: receivbucket)
    bs*=sizeof(std::pair<entity_key_t,body>);
  for(auto& bs: receivoffsets)
    bs*=sizeof(std::pair<entity_key_t,body>);*/



  // Use this array for the global buckets communication
  MPI_Alltoallv(&sendbuffer[0],&bucketssize[0],&sendoffsets[0],pair_entity_body,
      &rbodies[0],&receivbucket[0],&receivoffsets[0],pair_entity_body,
      MPI_COMM_WORLD);

  // Sort the incoming buffer 
  sort(rbodies.begin(),rbodies.end(),[](auto& left, auto &right)
  {
    return left.first<right.first;
  }); 

  // Check for duplicates
  assert(rbodies.end() == std::unique(rbodies.begin(),rbodies.end(),
      [](const auto& left, const auto& right ){ 
        return left.first == right.first;
      })
  );

  std::vector<int> totalnbodies(size);
  int mybodies = rbodies.size();
  // Share the final array size of everybody 
  MPI_Allgather(&mybodies,1,MPI_INT,&totalnbodies[0],1,MPI_INT,MPI_COMM_WORLD);

  if(rank == 0){
    std::cout<<"Repartition: ";
    for(auto num: totalnbodies)
      std::cout<<num<<";";
    std::cout<<std::endl;
  }

  //if(rank == 0){
  //  std::cout<<"Target: ";
  //  for(auto num: targetnbodies)
  //    std::cout<<num<<";";
  //  std::cout<<std::endl;
  //}

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
      int position = 0;
      MPI_Send(&rbodies[position],nsend,pair_entity_body,rank-1,0,
          MPI_COMM_WORLD);
      // Suppress the bodies 
      rbodies.erase(rbodies.begin(),rbodies.begin()+nsend); 
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

  // Display final particles 
  //for(auto bi: rbodies)
  //  std::cout<<rank<<": "<< bi.first <<std::endl;
}

void mpi_tree_traversal_graphviz(tree_topology_t & tree,
   std::array<point_t,2>& range)
{
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  char fname[64];
  sprintf(fname,"output_graphviz_%d.gv",rank);
  std::ofstream output;
  output.open(fname);
  output<<"digraph G {"<<std::endl;

  std::stack<branch_t*> stk;
  // Get root
  auto rt = tree.root();
  stk.push(rt);

  while(!stk.empty()){
    branch_t* cur = stk.top();
    stk.pop();
    if(!cur->is_leaf()){
      // Add the child to the stack and add for display 
      for(size_t i=0;i<(1<<gdimension);++i){
        auto br = tree.child(cur,i);
        //if(br){
          stk.push(br);
          if(gdimension == 3){
            output<<std::oct<<cur->id().value_()
              <<"->"<<br->id().value_()<<std::dec<<std::endl;
          }else if(gdimension == 1){
            output<<std::bitset<6>(cur->id().value_())<<"->"<<
              std::bitset<6>(br->id().value_())<<std::endl;
          }
        //}
      }
    }else{
      for(auto ent: *cur){
        entity_key_t key(range,ent->coordinates());
        output<<std::bitset<6>(cur->id().value_())<<
          "->"<<key<<std::endl;
        //fprintf(output,"\"%lo\" -> \"%lo\"\n",cur->id().value_(),
        //    key.truncate_value(17));
        switch (ent->getLocality()){
          case SHARED:
            output<<key<<"[shape=box,color=blue]"<<std::endl;
            break;
          case EXCL:
            output<<key<<" [shape=box,color=red]"<<std::endl;
            //fprintf(output,"\"%lo\" [shape=box,color=red]\n",
            //  key.truncate_value(17));
            break;
          case GHOST:
            output<<key<<" [shape=box,color=green]"<<std::endl;
            //fprintf(output,"\"%lo\" [shape=box,color=green]\n",
            //  key.truncate_value(17));
            break;
          default:
            output<<key<<" [shape=circle,color=black]"<<std::endl;
            //fprintf(output,"\"%lo\" [shape=circle,color=black]\n",
            //  key.truncate_value(17));
            break;
        }
      }
    } 
  }
  output<<"}"<<std::endl;
  output.close();
}

// The aim of this method is to shared the global informations on the tree 
// and the local informations with the process neighbors 
// Send all the body holders to everyone 
// Exchage all the tree informations
// Here we need to gather all the nodes and bodies info (key)
// The communication will propagate info in all the processes 
// For 4 processes :
// 0 <-> 1 2 <-> 3
// 0 <-> 2 1 <-> 3
void mpi_branches_exchange(tree_topology_t& tree)
{
  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  if(rank==0)
    std::cout<<"Branches repartition";
  int comms = log2(size);
  // TODO handle for non power of 2
  for(int comm=0;comm<comms;++comm)
  {
    int neighb = rank ^ (1 << comm);
    auto allents_ncontiguous = tree.entities();
    
    // Create a new vector to have contiguous memory TODO better way 
    std::vector<body_holder_mpi_t> allents;
    for(auto ent: allents_ncontiguous)
      allents.push_back(body_holder_mpi_t{ent->coordinates(),ent->getOwner()});
    
    // Send this to the neighb
    // 1. exchange to sizes
    int sendsize = allents.size()*sizeof(body_holder_mpi_t);
    int recvsize = 0;
    MPI_Sendrecv(&sendsize,1,MPI_INT,neighb,0,
      &recvsize,1,MPI_INT,neighb,MPI_ANY_TAG,
      MPI_COMM_WORLD,MPI_STATUS_IGNORE);

    std::vector<body_holder_mpi_t> recvents;
    recvents.resize(recvsize/sizeof(body_holder_mpi_t));
 
    MPI_Sendrecv(&allents[0],sendsize,MPI_BYTE,neighb,0,
      &recvents[0],recvsize,MPI_BYTE,neighb,MPI_ANY_TAG,
      MPI_COMM_WORLD,MPI_STATUS_IGNORE);

    // Add the entities in the tree, changing the ptr to nullptr
    for(auto bi: recvents)
    {
      assert(bi.owner!=rank);
      auto nbi = tree.make_entity(bi.position,nullptr,bi.owner);
      tree.insert(nbi);
    }
  }
  MPI_Barrier(MPI_COMM_WORLD); 
  if(rank == 0)
    std::cout<<".done"<<std::endl;
}
  
// compute the range, minposition and maxposition from a group of bodies
void 
mpi_compute_range(
    std::vector<std::pair<entity_key_t,body>>& bodies,
    std::array<point_t,2>& range
    )
{
  int rank, size; 
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  // Compute the range to compute the keys 
  double max[gdimension];
  for(size_t i=0;i<gdimension;++i)
    // TODO replace by the max value 
    max[i] = -9999;
  double min[gdimension];
  for(size_t i=0;i<gdimension;++i)
    // TODO replace by the max value 
    min[i] = 9999;
  for(auto bi: bodies){
    for(size_t i=0;i<gdimension;++i){
        if(bi.second.coordinates()[i]>max[i])
          max[i] = bi.second.coordinates()[i];
        if(bi.second.coordinates()[i]<min[i])
          min[i] = bi.second.coordinates()[i];
      }
  }

  // Do the MPI Reduction 
  MPI_Allreduce(MPI_IN_PLACE,max,gdimension,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD); 
  MPI_Allreduce(MPI_IN_PLACE,min,gdimension,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD); 
 
  point_t minposition; 
  point_t maxposition; 

  for(size_t i=0;i<gdimension;++i){
    minposition[i] = min[i]-0.1;
    maxposition[i] = max[i]+0.1;
  }

  if(rank==0)
    std::cout <<"boundaries: "<< minposition << maxposition << std::endl;
  
  range[0] = minposition;
  range[1] = maxposition;
}


void mpi_tree_traversal(tree_topology_t& tree)
{
  // Do a post order traversal 
}

// Go through the tree to see the structure

void mpi_gather_ghosts(
    tree_topology_t& tree,
    double smoothinglength,
    std::vector<body>& recvbodies)
{
  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  // 1. for each processes, bucket of send, bucket of receiv 
  std::vector<int> nsendholders;
  nsendholders.resize(size);
  std::fill(nsendholders.begin(),nsendholders.end(),0);
  std::vector<int> nrecvholders;
  nrecvholders.resize(size);
  std::fill(nrecvholders.begin(),nrecvholders.end(),0);
  std::vector<std::set<body_holder*>> sendholders;
  sendholders.resize(size);
  std::vector<std::set<body_holder*>> recvholders;
  recvholders.resize(size);

  // Clear the body vector 
  recvbodies.clear();

  int totalrecvbodies = 0;
  // 1.a bodies to send  
  auto treeents = tree.entities();
  for(auto bi: treeents)
  {
    if(!bi->is_local())
    {  
      auto gstneighbs = tree.find_in_radius(bi->coordinates(),
          2*smoothinglength);
      for(auto nb: gstneighbs)
      {
        // I f I am the owner, I can send data to this ghost owner
        if(nb->is_local())
        {
          //nsendholders[bi->getOwner()]++;
          sendholders[bi->getOwner()].insert(nb);
        }
      }
    }
  }
  // Unique, some bodies are required for several neighbs
  for(int i=0;i<size;++i){
    nsendholders[i] = sendholders[i].size();
  }
  
  //1.b body I will receive 
  for(auto bi: treeents)
  { 
    if(bi->is_local())
    {
      assert(bi->getOwner() == rank);
      auto bodiesneighbs = tree.find_in_radius(bi->coordinates(),
           2*smoothinglength);
        for(auto nb: bodiesneighbs)
       {
         if(!nb->is_local())
         {
          //nrecvholders[nb->getOwner()]++;
          recvholders[nb->getOwner()].insert(nb);
          //totalrecvbodies++;
         }
       }
    }
  }
  // Unique, some bodies are required for several neighbs
  for(int i=0;i<size;++i){
    nrecvholders[i] = recvholders[i].size();
    totalrecvbodies += nrecvholders[i];
  }

  // Make a vector with the recvholsters to be able to connect the pointer
  // at the end of the communication
  std::vector<body_holder*> totalrecvholders;
  for(auto proc: recvholders)
  {
    for(auto bi: proc)
    {
      totalrecvholders.push_back(bi);
    }
  }

  // Now gather the bodies data to send in a vector 
  std::vector<body> sendbodies;
  // Take the holders in the order 0 to n
  for(auto proc: sendholders)
  {
    for(auto bi: proc)
    {
      auto bodyholder = tree.get(bi->id());
      assert(bodyholder->getBody()!=nullptr);
      sendbodies.push_back(*(bodyholder->getBody()));
    }
  }

  // Prepare offsets for alltoallv
  std::vector<int> nrecvoffsets;
  nrecvoffsets.resize(size);
  std::fill(nrecvoffsets.begin(),nrecvoffsets.end(),0);
  std::vector<int> nsendoffsets;
  nsendoffsets.resize(size);
  std::fill(nsendoffsets.begin(),nsendoffsets.end(),0);

  for(int i=1;i<size;++i)
  {
    nrecvoffsets[i] = nrecvholders[i-1]+nrecvoffsets[i-1];
    nsendoffsets[i] = nsendholders[i-1]+nsendoffsets[i-1]; 
  }

  /*std::cout<<rank<<": recv "
    <<nrecvholders[0]<<":"
    <<nrecvholders[1]<<":"
    <<nrecvholders[2]<<":"
    <<nrecvholders[3]<<":"
    <<std::endl;
  std::cout<<rank<<": send "
    <<nsendholders[0]<<":"
    <<nsendholders[1]<<":"
    <<nsendholders[2]<<":"
    <<nsendholders[3]<<":"
    <<std::endl;

  std::cout<<rank<<": offrecv "
    <<nrecvoffsets[0]<<":"
    <<nrecvoffsets[1]<<":"
    <<nrecvoffsets[2]<<":"
    <<nrecvoffsets[3]<<":"
    <<std::endl;
  std::cout<<rank<<": offsend "
    <<nsendoffsets[0]<<":"
    <<nsendoffsets[1]<<":"
    <<nsendoffsets[2]<<":"
    <<nsendoffsets[3]<<":"
    <<std::endl;*/


  recvbodies.resize(totalrecvbodies);

  // Convert the offsets to byte
  for(int i=0;i<size;++i)
  {
    nsendholders[i]*=sizeof(body);
    nrecvholders[i]*=sizeof(body);
    nsendoffsets[i]*=sizeof(body);
    nrecvoffsets[i]*=sizeof(body);
  }
  
  /*std::cout<<"sizebody"<<sizeof(body)<<std::endl;

  std::cout<<rank<<":d recv "
    <<nrecvholders[0]<<":"
    <<nrecvholders[1]<<":"
    <<nrecvholders[2]<<":"
    <<nrecvholders[3]<<":"
    <<std::endl;
  std::cout<<rank<<":d send "
    <<nsendholders[0]<<":"
    <<nsendholders[1]<<":"
    <<nsendholders[2]<<":"
    <<nsendholders[3]<<":"
    <<std::endl;

  std::cout<<rank<<":d offrecv "
    <<nrecvoffsets[0]<<":"
    <<nrecvoffsets[1]<<":"
    <<nrecvoffsets[2]<<":"
    <<nrecvoffsets[3]<<":"
    <<std::endl;
  std::cout<<rank<<":d offsend "
    <<nsendoffsets[0]<<":"
    <<nsendoffsets[1]<<":"
    <<nsendoffsets[2]<<":"
    <<nsendoffsets[3]<<":"
    <<std::endl;*/


  MPI_Alltoallv(&sendbodies[0],&nsendholders[0],&nsendoffsets[0],MPI_BYTE,
      &recvbodies[0],&nrecvholders[0],&nrecvoffsets[0],MPI_BYTE,
      MPI_COMM_WORLD);

  // Sort the bodies based on position
  std::sort(recvbodies.begin(),recvbodies.end(),
      [](auto& left, auto& right){
        return left.coordinates()[0]<right.coordinates()[0];
      });
  // Sort the holders too 
  std::sort(totalrecvholders.begin(),totalrecvholders.end(),
      [](auto& left, auto& right){
        return left->coordinates()[0]<right->coordinates()[0];
      });


  // Then link the holders with these bodies
  auto it = recvbodies.begin(); 
  for(auto& bi: totalrecvholders)
  {
    auto bh = tree.get(bi->id());
    assert(bh->getLocality()==NONLOCAL||bh->getLocality()==GHOST);
    bh->setBody(&(*it));
    assert(bh->coordinates()==bh->getBody()->coordinates());
    //std::cout<<*(bh->getBody())<<std::endl;
    //bh->setLocality(NONLOCAL);
    ++it;
  }  
  
}

void mpi_output_txt(std::vector<std::pair<entity_key_t,body>>&rbodies,
    std::vector<int>& targetnbodies,
    int iter)
{

  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  int totalnbodies = std::accumulate(targetnbodies.begin(),
      targetnbodies.end(),0);
  //std::cout<<totalnbodies<<std::endl;
  // Display initial data 
  // Ouput the data
  // Gather on process 0 then output everythings
  std::vector<std::pair<entity_key_t,body>> gatheroutput;
  std::vector<int> recvcount;
  std::vector<int> recvoffsets;
  if(rank==0){
    recvcount.resize(size);
    recvoffsets.resize(size);
    std::fill(recvoffsets.begin(),recvoffsets.end(),0);
    for(int i=0;i<size;++i)
    {
      recvcount[i] = targetnbodies[i]*sizeof(std::pair<entity_key_t,body>);
      if(i < size-1){
        recvoffsets[i+1] = recvoffsets[i] + targetnbodies[i];
        //recvoffsets[i+1] *=  sizeof(std::pair<entity_key_t,body>);
      }
    }
    gatheroutput.resize(totalnbodies);
    for(auto& val: recvoffsets)
    {
      val *= sizeof(std::pair<entity_key_t,body>);
    } 
  }
  
   MPI_Gatherv(&rbodies[0],
      targetnbodies[rank]*sizeof(std::pair<entity_key_t,body>),
      MPI_BYTE,
      &gatheroutput[0],&recvcount[0],&recvoffsets[0],MPI_BYTE,
      0,MPI_COMM_WORLD);
  
  if(rank == 0)
  {
    char name[64];
    sprintf(name,"output_sod_%03d.txt",iter);
    std::cout<<"Output in file "<<name<<std::endl;
    //std::cout<<"Received: "<<gatheroutput.size()<<std::endl;
    FILE * file;
    file = fopen(name,"w");
    fprintf(file,"# X d p u\n");
    // Write in an output file 
    for(auto bi: gatheroutput)
    {
      fprintf(file,"%.10f %.10f %.10f %.10f\n",
          bi.second.getPosition()[0],bi.second.getDensity(),
          bi.second.getPressure(),bi.second.getInternalenergy());
    }
    fclose(file);
  }
}
