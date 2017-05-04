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



void mpi_sort_unbalanced(
    std::vector<std::pair<entity_key_t,body>>& rbodies,
    int totalnbodies
    )
{
  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  // Sort the keys 
  std::sort(rbodies.begin(),rbodies.end(),
      [](auto& left, auto& right)
      {
        return left.first < right.first;
      });

  // If on process, done 
  if(size==1)
    return;

  // Create a vector for the samplers 
  std::vector<entity_key_t> keys_sample;
  // Number of elements for sampling 
  // In this implementation we share up to 256KB to 
  // the master. 
  size_t noct = 256*1024;
  size_t nsample = noct / sizeof(entity_key_t);
  if(rbodies.size()<nsample){nsample = rbodies.size();}
  int chuncksize = rbodies.size()/nsample;
  //std::cout<<rank<<": nsample: "<<nsample<<" chuncksize: "<<chuncksize<<std::endl;

  for(size_t i=0;i<nsample;++i)
  {
    keys_sample.push_back(rbodies[chuncksize*i].first);
  }
  assert(keys_sample.size()==nsample);

  std::vector<entity_key_t> master_keys;
  std::vector<int> master_recvcounts;
  std::vector<int> master_offsets;
  int master_nkeys = 0; 

  if(rank==0)
    master_recvcounts.resize(size);
  
  // Echange the number of samples
  MPI_Gather(&nsample,1,MPI_INT,
      &master_recvcounts[0],1,MPI_INT,0,MPI_COMM_WORLD);

  // Master 
  if(rank == 0)
  {
    master_offsets.resize(size); 
    master_nkeys = noct/sizeof(entity_key_t)*size;
    if(totalnbodies<master_nkeys){master_nkeys=totalnbodies;}
    // Number to receiv from each process
    for(int i=0;i<size;++i)
      master_recvcounts[i]*=sizeof(entity_key_t);
    std::partial_sum(master_recvcounts.begin(),master_recvcounts.end(),
        &master_offsets[0]); 
    master_offsets.insert(master_offsets.begin(),0);
    master_keys.resize(master_nkeys);
  }

  MPI_Gatherv(&keys_sample[0],nsample*sizeof(entity_key_t),MPI_BYTE,
      &master_keys[0],&master_recvcounts[0],&master_offsets[0],MPI_BYTE,
      0,MPI_COMM_WORLD);

  // Generate the splitters
  std::vector<entity_key_t> splitters; 
  splitters.resize(size-1);
  if(rank==0)
  {

    std::sort(master_keys.begin(),master_keys.end());
    std::cout<<entity_key_t::first_key()<<std::endl;
    chuncksize = master_nkeys/size;
    for(int i=0;i<size-1;++i){
      splitters[i] = master_keys[(i+1)*chuncksize];
      std::cout<<splitters[i]<<std::endl;
    }
    std::cout<<entity_key_t::last_key()<<std::endl;
  }

  // Bradcast the splitters 
  MPI_Bcast(&splitters[0],(size-1)*sizeof(entity_key_t),MPI_BYTE,
      0,MPI_COMM_WORLD);

  // The rbodies array is already sorted. We just need to determine the 
  // limits for each process
  std::vector<int> sendcount(size);
  std::fill(sendcount.begin(),sendcount.end(),0);
  for(auto bi: rbodies)
  {
    for(int i = 0; i< size;++i)
    {
      if(i == 0 && bi.first < splitters[i]){
        sendcount[0]++;
      }else if(i == size-1 && bi.first >= splitters[size-2]){
        sendcount[i]++;
      }else if(bi.first < splitters[i] && bi.first >= splitters[i-1]){
        sendcount[i]++;
      }
    }
  }
 
  // Share the bucket 
  // First send the sizes of each bucket 
  std::vector<int> recvcount(size);
  MPI_Alltoall(&sendcount[0],1,MPI_INT,&recvcount[0],1,MPI_INT,
      MPI_COMM_WORLD);

  // Create the offset for alltoallv
  //  First receiv offset
  std::vector<int> recvoffsets(size);
  std::partial_sum(recvcount.begin(),recvcount.end(),&recvoffsets[0]); 
  // As we need an exscan, add a zero
  recvoffsets.insert(recvoffsets.begin(),0);
  
  // Then send offsets
  std::vector<int> sendoffsets(size);
  std::partial_sum(sendcount.begin(),sendcount.end(),&sendoffsets[0]);
  // As we need an exscan, add a zero
  sendoffsets.insert(sendoffsets.begin(),0);
 
  // The receiv buffer, work in th rbodies
  std::vector<std::pair<entity_key_t,body>> recvbuffer; 
  recvbuffer.resize(recvoffsets.back(),
      std::make_pair(entity_key_t::null(),body()));

  // We need to keep the last value to generate rbodies before erasing
  //recvoffsets.pop_back(); 
  //sendoffsets.pop_back();

  // Trnaform the offsets for bytes 
  for(auto& bs: sendcount)
    bs*=sizeof(std::pair<entity_key_t,body>);
  for(auto& bs: sendoffsets)
    bs*=sizeof(std::pair<entity_key_t,body>);
  for(auto& bs: recvcount)
    bs*=sizeof(std::pair<entity_key_t,body>);
  for(auto& bs: recvoffsets)
    bs*=sizeof(std::pair<entity_key_t,body>);

  // Use this array for the global buckets communication
  MPI_Alltoallv(&rbodies[0],&sendcount[0],&sendoffsets[0],MPI_BYTE,
      &recvbuffer[0],&recvcount[0],&recvoffsets[0],MPI_BYTE,
      MPI_COMM_WORLD);

  //rbodies.clear();
  rbodies = recvbuffer;

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

  std::vector<int> totalprocbodies;
  totalprocbodies.resize(size);
  int mybodies = rbodies.size();
  // Share the final array size of everybody 
  MPI_Allgather(&mybodies,1,MPI_INT,&totalprocbodies[0],1,MPI_INT,
      MPI_COMM_WORLD);

  if(rank == 0){
    std::cout<<"Repartition: ";
    for(auto num: totalprocbodies)
      std::cout<<num<<";";
    std::cout<<std::endl;
  }
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

  // Trnaform the offsets for bytes 
  for(auto& bs: bucketssize)
    bs*=sizeof(std::pair<entity_key_t,body>);
  for(auto& bs: sendoffsets)
    bs*=sizeof(std::pair<entity_key_t,body>);
  for(auto& bs: receivbucket)
    bs*=sizeof(std::pair<entity_key_t,body>);
  for(auto& bs: receivoffsets)
    bs*=sizeof(std::pair<entity_key_t,body>);

  // Use this array for the global buckets communication
  MPI_Alltoallv(&sendbuffer[0],&bucketssize[0],&sendoffsets[0],MPI_BYTE,
      &rbodies[0],&receivbucket[0],&receivoffsets[0],MPI_BYTE,
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

  std::vector<int> totalnbodies;
  totalnbodies.resize(size);
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
    if(rank!=0 && needs[rank]>0 && totalnbodies[rank-1]>0){
      int nrecv = needs[rank];
      if(totalnbodies[rank-1]<nrecv)
        nrecv = totalnbodies[rank-1];
      std::vector<std::pair<entity_key_t,body>> recvbuffer(nrecv);
      MPI_Recv(&recvbuffer[0],nrecv*sizeof(std::pair<entity_key_t,body>),
          MPI_BYTE,rank-1,0,MPI_COMM_WORLD,
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
    if(rank!=size-1 && needs[rank+1]>0 && totalnbodies[rank]>0){
      int nsend = needs[rank+1];
      if(nsend>totalnbodies[rank])
        nsend = totalnbodies[rank];
      int position = rbodies.size() - nsend;
      MPI_Send(&(rbodies[position]),nsend*sizeof(std::pair<entity_key_t,body>),
          MPI_BYTE,rank+1,0,
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
      MPI_Recv(&recvbuffer[0],nrecv*sizeof(std::pair<entity_key_t,body>),
          MPI_BYTE,rank+1,0,MPI_COMM_WORLD,
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
      MPI_Send(&rbodies[position],nsend*sizeof(std::pair<entity_key_t,body>),
          MPI_BYTE,rank-1,0,
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
            output<<std::bitset<64>(cur->id().value_())<<"->"<<
              std::bitset<64>(br->id().value_())<<std::endl;
          }
        //}
      }
    }else{
      for(auto ent: *cur){
        entity_key_t key(range,ent->coordinates());
        output<<std::bitset<64>(cur->id().value_())<<
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

void
mpi_branches_exchange_useful(tree_topology_t& tree,
    std::vector<std::pair<entity_key_t,body>>& rbodies,
    std::array<point_t,2>& rangekeys,
    std::vector<std::pair<entity_key_t,entity_key_t>>& ranges)
{
  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Barrier(MPI_COMM_WORLD);
  if(rank==0)
    std::cout<<"Branches repartition" << std::flush;
  
  // Search for my min and max key needed
  std::pair<entity_key_t,entity_key_t> range;
  range.first = entity_key_t::last_key();
  range.second = entity_key_t::first_key();

  // Lowest key:
  auto minkey = std::min_element(rbodies.begin(),rbodies.end(),
      [](auto& left, auto& right)
      {
        return left.first < right.first; 
      });
  auto maxkey = std::max_element(rbodies.begin(),rbodies.end(),
      [](auto& left, auto& right)
      {
        return left.first < right.first; 
      });
  
  range.first = (*minkey).first;
  range.second = (*maxkey).first;
  //std::cout<<rank<<": "<<(*minkey).first<<"-"<<(*maxkey).first<<std::endl;

  // Point out of bounds
  auto outofbounds = [](point_t& p, point_t& limit) ->bool
  {
    for(size_t i=0;i<gdimension;++i)
      if(p[i]<limit[i])
        return true;
    return false;
  };

  // Get all my entities 
  for(auto ent: rbodies)
  {
    // TODO check validity for 3 dimensions
    point_t pos = ent.second.getPosition(); 
    double h = ent.second.getSmoothinglength()*2;

    // not considering less than the total range ! 
    point_t p = pos-h;
    if(!outofbounds(p,rangekeys[0]))
      if(entity_key_t(rangekeys,p) < range.first)
        range.first = entity_key_t(rangekeys,p);
    p = pos+h;
    if(!outofbounds(rangekeys[1],p))
      if(entity_key_t(rangekeys,p) > range.second)
        range.second = entity_key_t(rangekeys,p);
  }
  // Handle on the extremities for 0 and size-1
  if(range.first < entity_key_t::first_key())
    range.first = entity_key_t::first_key();

  if(range.second > entity_key_t::last_key())
    range.second = entity_key_t::last_key();

  //std::cout<<rank<<": "<<range.first<<" "<<range.second<<std::endl<<std::flush;
  //MPI_Barrier(MPI_COMM_WORLD);
  //exit(-1);

  // Gather the keys of everyone 
  if(ranges.size() == 0)
    ranges.resize(size);
  //std::vector<std::pair<entity_key_t,entity_key_t>> ranges(size);
  MPI_Allgather(&range,sizeof(std::pair<entity_key_t,entity_key_t>),
      MPI_BYTE,&ranges[0],sizeof(std::pair<entity_key_t,entity_key_t>),
      MPI_BYTE,MPI_COMM_WORLD);

  // Now generate the snedbufer, ordered by processes
  //  for the holders 
  std::vector<body_holder_mpi_t> sendbuffer;
  std::vector<int> sendcount(size);
  // To be able to search in the rbodies array 
  std::vector<std::pair<entity_key_t,body>> s1(size);
  std::vector<std::pair<entity_key_t,body>> s2(size);

  for(int i=0;i<size;++i)
  {
    s1[i].first=ranges[i].first;
    s2[i].first=ranges[i].second;
  }

  for(int i=0;i<size;++i)
  {
    if(i==rank)
    {
      sendcount[i]=0;
      continue;
    }
    auto lb = std::lower_bound(rbodies.begin(),rbodies.end(),s1[i],
        [](auto& left, auto& right)
        {
          return left.first < right.first; 
        });
    auto hb = std::upper_bound(rbodies.begin(),rbodies.end(),s2[i],
        [](auto& left, auto& right)
        {
          return left.first < right.first; 
        });
    // TODO check here, if last element is valid 
    if(lb != rbodies.end())
    {
      sendcount[i]=std::distance(lb,hb); 
      for(;lb!=hb;++lb)
        sendbuffer.push_back(
            body_holder_mpi_t{lb->second.coordinates(),rank,
            lb->second.getMass()}); 
    }else{
      sendcount[i]=0;
    } 
  }

  //int tot = 0;
  //for(auto val: sendcount)
  //{
  //  std::cout<<rank<<": sendto "<<tot<<"="<<val<<" from: "<<
  //    s1[tot].first<<";"<<s2[tot].first<<std::endl;
  //  tot++;
  //}

  //Count the elements to send 
  std::vector<int> sendoffsets(size);
  std::vector<int> recvcount(size);
  std::vector<int> recvoffsets(size);

  MPI_Alltoall(&sendcount[0],1,MPI_INT,
      &recvcount[0],1,MPI_INT,MPI_COMM_WORLD);

  int totalrecv = 0;
  for(int i=0;i<size;++i)
  {
    totalrecv += recvcount[i];
    sendcount[i]*=sizeof(body_holder_mpi_t);
    recvcount[i]*=sizeof(body_holder_mpi_t);
    if(i<size-1)
    {
      sendoffsets[i+1]=sendoffsets[i]+sendcount[i];
      recvoffsets[i+1]=recvoffsets[i]+recvcount[i];
    }
  }

  std::vector<body_holder_mpi_t> recvbuffer(totalrecv);
  MPI_Alltoallv(&sendbuffer[0],&sendcount[0],&sendoffsets[0],MPI_BYTE,
      &recvbuffer[0],&recvcount[0],&recvoffsets[0],MPI_BYTE,
      MPI_COMM_WORLD);

  // Add in the tree 
  for(auto bi: recvbuffer)
  {
    assert(bi.owner!=rank);
    auto nbi = tree.make_entity(bi.position,nullptr,bi.owner,bi.mass);
    tree.insert(nbi);
  }

  if(rank==0)
    std::cout<<".done"<<std::endl;


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
  MPI_Barrier(MPI_COMM_WORLD);
  if(rank==0)
    std::cout<<"Branches repartition";

//#if 0
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
      auto nbi = tree.make_entity(bi.position,nullptr,bi.owner,
          bi.mass);
      tree.insert(nbi);
    }
  }
//#endif

#if 0
  // Do a MPI_Allgatherv instead 
  // Get all my tree entities
  auto allents_ncontiguous = tree.entities();
    
  // Create a new vector to have contiguous memory TODO better way ?
  std::vector<body_holder_mpi_t> allents;
  for(auto ent: allents_ncontiguous)
  {
    allents.push_back(body_holder_mpi_t{ent->coordinates(),ent->getOwner()});
  } 
  // Prepare MPI_Allgatherv
  // 1. exchange the sizes
  int sendcount;
  std::vector<int> recvcount(size);
  sendcount = allents.size();
  MPI_Allgather(&sendcount,1,MPI_INT,
      &recvcount[0],1,MPI_INT,MPI_COMM_WORLD);
  // Compute recv offset and do transfert 
  std::vector<int> recvoffsets(size);
  recvoffsets[0] = 0;
  int totalrecv = 0;
  for(int i=0;i<size;++i)
  {
    totalrecv += recvcount[i];
    recvcount[i] = recvcount[i]*sizeof(body_holder_mpi_t);
    if(i<size-1)
      recvoffsets[i+1] = recvoffsets[i]+recvcount[i];
  } 
  std::vector<body_holder_mpi_t> recvbuffer(totalrecv);

  MPI_Allgatherv(&allents[0],sendcount*sizeof(body_holder_mpi_t),MPI_BYTE,
      &recvbuffer[0],&recvcount[0],&recvoffsets[0],MPI_BYTE,
      MPI_COMM_WORLD);

  // Add everything in my tree
  // changing the ptr to nullptr 
  for(auto bi: recvbuffer)
  {
    if(bi.owner != rank)
    {
      auto nbi = tree.make_entity(bi.position,nullptr,bi.owner);
      tree.insert(nbi);
    }
  }

#endif
  MPI_Barrier(MPI_COMM_WORLD); 
  if(rank == 0)
    std::cout<<".done"<<std::endl;
}
  
// compute the range, minposition and maxposition from a group of bodies
void 
mpi_compute_range(
    std::vector<std::pair<entity_key_t,body>>& bodies,
    std::array<point_t,2>& range,
    double smoothinglength
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
    minposition[i] = min[i]-2*smoothinglength;
    maxposition[i] = max[i]+2*smoothinglength;
  }

  if(rank==0)
    std::cout <<"boundaries: "<< minposition << maxposition << std::endl;
  
  range[0] = minposition;
  range[1] = maxposition;
}


void mpi_tree_traversal_com(tree_topology_t& tree)
{
  std::function<void(branch_t*)>traverse;
  traverse = [&tree,&traverse](branch_t * b)
  {
    //std::cout<< b->id() <<std::endl<<std::flush;
    double mass = 0.0;
    point_t com = {0,0,0};
    if(b->is_leaf())
    {
      for(auto child: *b)
      {
        // Only for local particles 
        if(child->is_local()){
          com += child->getMass()*child->getPosition();
          mass += child->getMass();
        }
      }
      com = com / mass;
    }
    else
    {
      for(int i=0;i<(1<<gdimension);++i)
      {
        auto branch = tree.child(b,i);
        traverse(branch);
        if(branch->getMass()!=0)
        {
          // Add this children position and coordinates
          com += branch->getMass()*branch->getPosition();
          mass += branch->getMass();
        }
      }
      com = com / mass;
    }
    b->setMass(mass);
    b->setPosition(com);
    //std::cout<<com<<";"<<mass<<std::endl;

  };
  traverse(tree.root());
}

void mpi_refresh_ghosts(
    tree_topology_t& tree, 
    mpi_ghosts_t& refresh,
    std::array<point_t,2>& range    
    )
{
  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  
  MPI_Barrier(MPI_COMM_WORLD);
  if(rank==0)
    std::cout<<"Refresh Ghosts" << std::flush;

  // Refresh the sendbodies with new data
  auto itsb = refresh.sendbodies.begin();
  for(auto proc: refresh.sendholders)
  {
    for(auto bi: proc)
    {
      assert(bi->getBody()!=nullptr);
      *itsb = *(bi->getBody());
      itsb++;
    }
  } 

  MPI_Alltoallv(&refresh.sendbodies[0],&refresh.nsendholders[0],
      &refresh.nsendoffsets[0],MPI_BYTE,
      &refresh.recvbodies[0],&refresh.nrecvholders[0],
      &refresh.nrecvoffsets[0],MPI_BYTE,
      MPI_COMM_WORLD);

  // Sort the bodies based on key
  std::sort(refresh.recvbodies.begin(),refresh.recvbodies.end(),
      [range](auto& left, auto& right){
        return entity_key_t(range,left.coordinates())<
        entity_key_t(range,right.coordinates());
      });
 
  // Then link the holders with these bodies
  auto it = refresh.recvbodies.begin(); 
  for(auto& bi: refresh.totalrecvholders)
  {
    auto bh = tree.get(bi->id());
    assert(bh->getLocality()==NONLOCAL||bh->getLocality()==GHOST);
    bh->setBody(&(*it));
    assert(bh->coordinates()==bh->getBody()->coordinates());
    //std::cout<<*(bh->getBody())<<std::endl;
    //bh->setLocality(NONLOCAL);
    ++it;
  }  
  
  MPI_Barrier(MPI_COMM_WORLD);
  if(rank==0)
    std::cout<<".done" << std::endl << std::flush;
}

void 
mpi_compute_ghosts(
    tree_topology_t& tree,
    double smoothinglength,
    mpi_ghosts_t& ghosts_data,
    std::array<point_t,2>& range)
{
  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  
  MPI_Barrier(MPI_COMM_WORLD);
  if(rank==0)
    std::cout<<"Compute Ghosts" << std::flush;

  // Clean the structure 
  ghosts_data.sendbodies.clear();
  ghosts_data.recvbodies.clear();
  ghosts_data.totalrecvholders.clear();
  ghosts_data.sendholders.clear();
  // For the first iteration init the count and offset size in the structure
  if(ghosts_data.nsendholders.size()==0)
  {
    ghosts_data.nsendholders.resize(size);
    ghosts_data.nsendoffsets.resize(size);
    ghosts_data.nrecvholders.resize(size);
    ghosts_data.nrecvoffsets.resize(size);
  }
  assert(ghosts_data.nsendholders.size()==(size_t)size);

  // 1. for each processes, bucket of send, bucket of receiv 
  ghosts_data.sendholders.resize(size);
  std::vector<std::set<body_holder*>> recvholders(size);

  // Considering the biggest h
  // TODO add a reduction over h 
  int totalrecvbodies = 0;
  int totalsendbodies = 0;
  auto treeents = tree.entities();
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
           // THIS IS TRUE BECAUSE WE CONSIDER THE BIGGEST SMOOTHING LENGTH
           // The distant particle will need mine
           ghosts_data.sendholders[nb->getOwner()].insert(bi);;
          // I will also need the distant particle 
          recvholders[nb->getOwner()].insert(nb);
         }
       }
    }
  }

  for(int i=0;i<size;++i){
    ghosts_data.nsendholders[i] = ghosts_data.sendholders[i].size();
    assert(ghosts_data.nsendholders[i]>=0);
    totalsendbodies += ghosts_data.nsendholders[i];
  }
  
  for(int i=0;i<size;++i){
    ghosts_data.nrecvholders[i] = recvholders[i].size();
    assert(ghosts_data.nrecvholders[i]>=0);
    totalrecvbodies += ghosts_data.nrecvholders[i];
  }
  
  // Make a vector with the recvholsters to be able to connect the pointer
  // at the end of the communication
  for(auto proc: recvholders){
    ghosts_data.totalrecvholders.insert(
        ghosts_data.totalrecvholders.end(),proc.begin(),proc.end());
  }

  // Now gather the bodies data to send in a vector 
  // Take the holders in the order 0 to n processes
  //for(auto proc: ghosts_data.sendholders)
  //{
  //  for(auto bi: proc)
  //  {
  //    auto bodyholder = tree.get(bi->id());
  //    assert(bodyholder->getBody()!=nullptr);
  //    ghosts_data.sendbodies.push_back(*(bodyholder->getBody()));
  //  }
  //}
  ghosts_data.sendbodies.resize(totalsendbodies);

  // Prepare offsets for alltoallv
  ghosts_data.nrecvoffsets[0]=0;
  ghosts_data.nsendoffsets[0]=0;

  for(int i=1;i<size;++i)
  {
    ghosts_data.nrecvoffsets[i] = ghosts_data.nrecvholders[i-1]+
      ghosts_data.nrecvoffsets[i-1];
    ghosts_data.nsendoffsets[i] = ghosts_data.nsendholders[i-1]+
      ghosts_data.nsendoffsets[i-1]; 
  }

  ghosts_data.recvbodies.resize(totalrecvbodies);

  // Convert the offsets to byte
  for(int i=0;i<size;++i)
  {
    ghosts_data.nsendholders[i]*=sizeof(body);
    assert(ghosts_data.nsendholders[i]>=0);
    ghosts_data.nrecvholders[i]*=sizeof(body);
    assert(ghosts_data.nrecvholders[i]>=0);
    ghosts_data.nsendoffsets[i]*=sizeof(body);
    assert(ghosts_data.nsendoffsets[i]>=0);
    ghosts_data.nrecvoffsets[i]*=sizeof(body);
    assert(ghosts_data.nrecvoffsets[i]>=0);
  }
  
  //std::cout<<rank<<": "<<  <<std::endl;

  // Sort the aolders once
  std::sort(ghosts_data.totalrecvholders.begin(),
            ghosts_data.totalrecvholders.end(),
      [range](auto& left, auto& right){
        return entity_key_t(range,left->coordinates())<
        entity_key_t(range,right->coordinates());
      });


  MPI_Barrier(MPI_COMM_WORLD);
  if(rank==0)
    std::cout<<".done"<<std::endl;
  
}

void mpi_output_txt(
    std::vector<std::pair<entity_key_t,body>>&rbodies,
    int iter)
{

  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  std::vector<int> processnbodies(size);
  processnbodies[rank] = rbodies.size();
  // Gather the number of bodies per process on 0
  MPI_Gather(&processnbodies[rank],1,MPI_INT,
      &processnbodies[0],1,MPI_INT,0,MPI_COMM_WORLD);

  int totalnbodies = std::accumulate(processnbodies.begin(),
      processnbodies.end(),0);
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
      recvcount[i] = processnbodies[i]*sizeof(std::pair<entity_key_t,body>);
      if(i < size-1){
        recvoffsets[i+1] = recvoffsets[i] + processnbodies[i];
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
      processnbodies[rank]*sizeof(std::pair<entity_key_t,body>),
      MPI_BYTE,
      &gatheroutput[0],&recvcount[0],&recvoffsets[0],MPI_BYTE,
      0,MPI_COMM_WORLD);
  
  if(rank == 0)
  {
    char name[64];
    sprintf(name,"output_sod_%05d.txt",iter);
    std::cout<<"Output in file "<<name<<std::endl;
    //std::cout<<"Received: "<<gatheroutput.size()<<std::endl;
    FILE * file;
    file = fopen(name,"w");
    fprintf(file,"# pX pY pZ d p u vX vY vZ\n");
    // Write in an output file 
    for(auto bi: gatheroutput)
    {
      if(gdimension==1)
        fprintf(file,"%.10f %.10f %.10f %.10f %.10f\n",
          bi.second.getPosition()[0],bi.second.getDensity(),
          bi.second.getPressure(),bi.second.getInternalenergy(),
          bi.second.getVelocity()[0]);
      if(gdimension==3)
        fprintf(file,"%.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f\n",
          bi.second.getPosition()[0],bi.second.getPosition()[1],
          bi.second.getPosition()[2],bi.second.getDensity(),
          bi.second.getPressure(),bi.second.getInternalenergy(),
          bi.second.getVelocity()[0],bi.second.getVelocity()[1],
          bi.second.getVelocity()[2]);

    }
    fclose(file);
  }
}

void 
mpi_gather_com(
    tree_topology_t& tree,
    std::array<point_t,2>& range,
    std::vector<std::pair<entity_key_t,entity_key_t>>& rangeproc,
    std::vector<body_holder>& recv_COM
    )
{
  int rank,size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);


  if(rank==0)
    std::cout<<"Gather COM"<<std::flush;
  MPI_Barrier(MPI_COMM_WORLD);

  // Reset the recv vector 
  recv_COM.clear();
  // Compute COM data for the current tree
  mpi_tree_traversal_com(tree);

  // Display my root data
  double mymass = tree.root()->getMass();
  //std::cout<<rank<<":"<<tree.root()->getPosition()<<";"
  //  <<tree.root()->getMass()<<std::endl;
  double total = 0.0;
  MPI_Reduce(&(mymass),&total,1,MPI_DOUBLE,
      MPI_SUM,0,MPI_COMM_WORLD); 
  //if(rank==0)
  //  std::cout<<"Total mass = "<<total<<std::endl;
  
  // Get the interesting branches for everyones. 
  // Do a traversal and create a vector for every other processes

  // Set a unique vector for everyone, to be contiguous 
  std::vector<body_holder> shared_COM;
  std::vector<int> shared_COM_size(size);
  std::vector<int> shared_COM_offsets(size);
  std::fill(shared_COM_size.begin(),shared_COM_size.end(),0);
  std::fill(shared_COM_offsets.begin(),shared_COM_offsets.end(),0);

  std::cout<<rank<<": keyrank="<<rangeproc[rank].first<<" : "<<rangeproc[rank].second<<std::endl;


  std::function<void(
    branch_t *,
    std::vector<body_holder>&,
    std::pair<entity_key_t,entity_key_t>&,
    int&)>traverse;

  traverse = [rank,&tree,&traverse,&range]
    (branch_t * b,
     std::vector<body_holder>& vbh,
     std::pair<entity_key_t,entity_key_t>& rangekeys,
     int& nelements)
  {
    // If in the range, go further 
    if(rangekeys.first<entity_key_t(range,b->getPosition())
        &&entity_key_t(range,b->getPosition())<rangekeys.second)
    {
      if(!b->is_leaf())
      {
        for(int i=0;i<(1<<gdimension);++i)
          traverse(tree.child(b,i),vbh,rangekeys,nelements);
      }else{
        // If I am a leaf check if one of my children is out of the range
        for(auto child: *b)
        {
          // Only for local particles 
          if(child->is_local()){
            if(rangekeys.first<entity_key_t(range,child->getPosition())
        &&entity_key_t(range,child->getPosition())<rangekeys.second)
            {
              vbh.push_back(body_holder(child->getPosition(),nullptr,rank,
                  child->getMass()));
              nelements++;
            }
          }   
        }
      }
    }else
    {
      // Mass = 0 just for elements I did not owe
      if(b->getMass()!=0)
      {
        // Not in range, stop going down and add in vector 
        vbh.push_back(body_holder(b->getPosition(),nullptr,rank,
              b->getMass()));
        nelements++;
      }
    } 
  };

  for(int i=0;i<size;++i)
  {
    int nelements = 0;
    // Do a tree traversal knowing the limits of every processes limits
    // Only take the keys that are outside
    //std::cout<<rank<<" search for "<<i<<": "<<rangeproc[i].first
    //  <<":"<<rangeproc[i].second<<std::endl;
    if(i!=rank)
      traverse(tree.root(),shared_COM,rangeproc[i],nelements);
    shared_COM_size[i]=nelements*sizeof(body_holder);
    if(i<size-2)
      shared_COM_offsets[i+1] = shared_COM_offsets[i]+shared_COM_size[i];
    //std::cout<<rank<<" for "<<i<<": "<<nelements<<std::endl;
  }

  std::vector<int> recvcount(size);
  std::vector<int> recvoffsets(size);
  // Share data size
  MPI_Alltoall(&shared_COM_size[0],1,MPI_INT,
      &recvcount[0],1,MPI_INT,MPI_COMM_WORLD);
  std::partial_sum(recvcount.begin(),recvcount.end(),&recvoffsets[0]); 
  // As we need an exscan, add a zero and delete the last element 
  recvoffsets.insert(recvoffsets.begin(),0);
 
  recv_COM.resize(recvoffsets.back()/sizeof(body_holder)); 
  // Share the data
  MPI_Alltoallv(&shared_COM[0],&shared_COM_size[0],&shared_COM_offsets[0],
      MPI_BYTE,&recv_COM[0],&recvcount[0],&recvoffsets[0],
      MPI_BYTE,MPI_COMM_WORLD); 

  if(rank==0)
    std::cout<<".done"<<std::endl<<std::flush;
  if(rank==0)
    std::cout<<"Total mass = "<<total<<std::endl;
  

  MPI_Barrier(MPI_COMM_WORLD);
}

