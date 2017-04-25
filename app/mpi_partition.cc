#include "mpi_partition.h"
  
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
  return p1[0]==p2[0]&&
    p1[1]==p2[1]&&
    p1[2]==p2[2];
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
  int blocks[2] = {1,24};
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

  // Use this array for the global buckets communication
  MPI_Alltoallv(&sendbuffer[0],&bucketssize[0],&sendoffsets[0],pair_entity_body,
      &rbodies[0],&receivbucket[0],&receivoffsets[0],pair_entity_body,
      MPI_COMM_WORLD);

  // Sort the incoming buffer 
  sort(rbodies.begin(),rbodies.end(),[](auto& left, auto &right)
  {
    return left.first<right.first;
  }); 

  for(auto it=rbodies.begin();it<rbodies.end();++it){
    auto nx = std::next(it);
    if((*nx).first == (*it).first){
      std::cout<<(*nx).first<<";"<<(*it).first<<std::endl;
      std::cout<<(*nx).second<<std::endl;
      std::cout<<(*it).second<<std::endl;

    }
  }

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

  if(rank == 0){
    std::cout<<"Target: ";
    for(auto num: targetnbodies)
      std::cout<<num<<";";
    std::cout<<std::endl;
  }

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

  FILE * output;
  char fname[64];
  sprintf(fname,"output_graphviz_%d.gv",rank);
  output = fopen(fname,"w");
  fprintf(output,"digraph G {");

  std::stack<branch_t*> stk;
  // Get root
  auto rt = tree.root();
  stk.push(rt);

  while(!stk.empty()){
    branch_t* cur = stk.top();
    stk.pop();
    if(!cur->is_leaf()){
      // Add the child to the stack and add for display 
      for(int i=0;i<8;++i){
        auto br = tree.child(cur,i);
        if(br){
          stk.push(br);
          fprintf(output,"\"%llo\" -> \"%llo\"\n",
              cur->id().value_(),br->id().value_());
        }
      }
    }else{
      for(auto ent: *cur){
        entity_key_t key(range,ent->coordinates());
        fprintf(output,"\"%llo\" -> \"%llo\"\n",cur->id().value_(),
            key.truncate_value(17));
        switch (ent->getLocality()){
          case SHARED:
            fprintf(output,"\"%llo\" [shape=box,color=blue]\n",
              key.truncate_value(17));
            break;
          case EXCL:
            fprintf(output,"\"%llo\" [shape=box,color=red]\n",
              key.truncate_value(17));
            break;
          case GHOST:
            fprintf(output,"\"%llo\" [shape=box,color=green]\n",
              key.truncate_value(17));
            break;
          default:
            fprintf(output,"\"%llo\" [shape=circle,color=black]\n",
              key.truncate_value(17));
            break;
        }
      }
    } 
  }
  fprintf(output,"}");
  fclose(output);
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
    //std::cout<<"comm: "<<rank<<"<>"<<neighb<<std::endl;
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

    //for(auto ent: allents)
    //  std::cout<<rank<<": SEND"<<ent.position
    //    <<":"<<ent.owner<<" S: "<<sizeof(ent)<<std::endl;
    //body_holder_mpi_t * ptr = &allents[0];
    //for(int i = 0 ; i < 20 ; ++i)
    //  std::cout<<rank<<": SEND"<<(ptr+i)->position
    //    <<":"<<(ptr+i)->owner<<std::endl;
    //std::cout<<rank<<": SIZE:"<<sendsize<<":"<<recvsize<<std::endl;

    MPI_Sendrecv(&allents[0],sendsize,MPI_BYTE,neighb,0,
      &recvents[0],recvsize,MPI_BYTE,neighb,MPI_ANY_TAG,
      MPI_COMM_WORLD,MPI_STATUS_IGNORE);

    //for(auto ent: recvents)
    //  std::cout<<rank<<": RECV"<<ent.position<<ent.owner<<std::endl;
    // Add the entities in the tree, changing the ptr to nullptr
    for(auto bi: recvents)
    {
      //if(bi.getOwner()!=rank){
        auto nbi = tree.make_entity(bi.position,nullptr,bi.owner);
        //std::cout<<rank<<": inserting: "<<*nbi<<std::endl;
        tree.insert(nbi);
      //}
    }
  }
  MPI_Barrier(MPI_COMM_WORLD); 
  if(rank == 0)
    std::cout<<".done"<<std::endl;
}

void mpi_tree_traversal(tree_topology_t& tree)
{
  // Do a post order traversal 
}

// Go through the tree to see the structure 
