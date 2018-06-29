/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Los Alamos National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/

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
 * @file utils.h
 * @author Julien Loiseau
 * @date June 2018
 * @brief Function needed for MPI distribution of the bodies 
 */

#ifndef _mpisph_utils_
#define _mpisph_utils_

// Local version of assert to handle MPI abord
static void mpi_assert_fct(
  bool expression, 
  const char *file,
  int line);

#define mpi_assert( err ) (mpi_assert_fct(err,__FILE__,__LINE__))

// Local version of assert to handle MPI abord
static void mpi_assert_fct(
  bool expression, 
  const char *file,
  int line)
{
  if (!(expression)) {
     fprintf(stderr, "Failed assertion in %s in %d\n",file, line);
     MPI_Abort(MPI_COMM_WORLD, 1);
  }
}

namespace mpi_utils{



  /**
 * @brief      Simple version of all to all 
 * Use to generate the offsets and do the pre exchange 
 * Then realise the MPI_Alltoall call 
 *
 * @param[in]  sendcount   The sendcount
 * @param      sendbuffer  The sendbuffer
 * @param      recvbuffer  The recvbuffer
 *
 * @tparam     M           The type of data sent 
 */
  template<
    typename M>
  void 
  mpi_alltoallv(
      std::vector<int> sendcount,
      std::vector<M>& sendbuffer,
      std::vector<M>& recvbuffer
    )
  {
    int rank, size; 
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size); 

    std::vector<int> recvcount(size); 
    std::vector<int> recvoffsets(size); 
    std::vector<int> sendoffsets(size); 

    // Exchange the send count 
    MPI_Alltoall(&sendcount[0],1,MPI_INT,&recvcount[0],1,MPI_INT,
      MPI_COMM_WORLD);
    
    // Generate the send and recv offsets 
    std::partial_sum(recvcount.begin(),recvcount.end(),&recvoffsets[0]); 
    // As we need an exscan, add a zero
    recvoffsets.insert(recvoffsets.begin(),0);
    
    // Then send offsets
    std::partial_sum(sendcount.begin(),sendcount.end(),&sendoffsets[0]);
    // As we need an exscan, add a zero
    sendoffsets.insert(sendoffsets.begin(),0);
    
    // Set the recvbuffer to the right size
    recvbuffer.resize(recvoffsets.back()); 
    
    // Trnaform the offsets for bytes 
    for(int i=0;i<size;++i){
      sendcount[i] *= sizeof(M);
      recvcount[i] *= sizeof(M);
      sendoffsets[i] *= sizeof(M);
      recvoffsets[i] *= sizeof(M);
    } // for
    
    // Use this array for the global buckets communication
    MPI_Alltoallv(&sendbuffer[0],&sendcount[0],&sendoffsets[0],MPI_BYTE,
      &recvbuffer[0],&recvcount[0],&recvoffsets[0],MPI_BYTE,MPI_COMM_WORLD);
  }

  void reduce_min(
    double& value)
  {
    MPI_Allreduce(MPI_IN_PLACE,&value,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
  }

  void reduce_min(
    float& value)
  {
    MPI_Allreduce(MPI_IN_PLACE,&value,1,MPI_FLOAT,MPI_MIN,MPI_COMM_WORLD);
  }

  void reduce_min(
    int& value)
  {
    MPI_Allreduce(MPI_IN_PLACE,&value,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
  }


  void reduce_sum(
    double& value)
  {
    MPI_Allreduce(MPI_IN_PLACE,&value,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  }

  void reduce_sum(
    float& value)
  {
    MPI_Allreduce(MPI_IN_PLACE,&value,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
  }

  void reduce_sum(
    int& value)
  {
    MPI_Allreduce(MPI_IN_PLACE,&value,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  }

  void reduce_sum(
    point_t& value)
  {
    MPI_Allreduce(MPI_IN_PLACE,&value[0],gdimension,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  }


  /**
   * @brief      Export to a file the current tree in memory 
   * This is useful for small number of particles to help representing the tree 
   *
   * @param      tree   The tree to output
   * @param      range  The range of the particles, use to construct entity_keys
   */
  void mpi_tree_traversal_graphviz(
    tree_topology_t & tree/*,*/
    /*std::array<point_t,2>& range*/)
  {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    char fname[64];
    sprintf(fname,"output_graphviz_%02d.gv",rank);
    std::ofstream output;
    output.open(fname);
    output<<"digraph G {"<<std::endl<<"forcelabels=true;"<<std::endl;

    std::stack<branch_t*> stk;
    // Get root
    auto rt = tree.root();
    stk.push(rt);

    while(!stk.empty()){
      branch_t* cur = stk.top();
      stk.pop();
      if(!cur->is_leaf()){
        if(gdimension == 3){
          output<<std::oct<<cur->id().value_()<<" [label=\""<< 
          cur->id().value_()<<std::dec<< "\", xlabel=\"" << cur->sub_entities() 
          <<"\"];"<<std::endl;
        }
        if(gdimension == 2){
          output<<cur->id().value_()<<" [label=\""<< cur->id().value_() 
          << "\", xlabel=\"" << cur->sub_entities() <<"\"];"<<std::endl;

        }

        // Add the child to the stack and add for display 
        for(size_t i=0;i<(1<<gdimension);++i)
        {
          auto br = tree.child(cur,i);
          stk.push(br);
          if(gdimension == 3){
            output<<std::oct<<cur->id().value_()
              <<"->"<<br->id().value_()<<std::dec<<std::endl;
          }
          if(gdimension == 2){
            output<<cur->id().value_()
              <<"->"<<br->id().value_()<<std::dec<<std::endl;  
          }
        }
      }else{
        if(gdimension == 3){
          output<<std::oct<<cur->id().value_()<<" [label=\""<< 
          cur->id().value_() <<std::dec<< "\", xlabel=\"" << cur->sub_entities() 
          <<"\"];"<<std::endl;
        }
        if(gdimension == 2){
          output<<cur->id().value_()<<" [label=\""<< cur->id().value_() 
          << "\", xlabel=\"" << cur->sub_entities() <<"\"];"<<std::endl; 
        }    
        for(auto ent: *cur)
        {
          entity_key_t key(ent->coordinates());
          int64_t key_int = key.truncate_value(tree.max_depth()+2);
          if(gdimension == 3){
            output<<std::oct<<cur->id().value_()<<
              "->"<<key_int<<std::endl;
          }
          if(gdimension == 2){
            output<<cur->id().value_()<<
              "->"<<key_int<<std::endl;
          }
          switch (ent->getLocality())
          {
            case 2:
              output<<key_int<<" [shape=box,color=blue]"<<std::endl;
              break;
            case 3:
              output<<key_int<<" [shape=box,color=red]"<<std::endl;
              //fprintf(output,"\"%lo\" [shape=box,color=red]\n",
              //  key.truncate_value(17));
              break;
            case 1:
              output<<key_int<<" [shape=box,color=green]"<<std::endl;
              //fprintf(output,"\"%lo\" [shape=box,color=green]\n",
              //  key.truncate_value(17));
              break;
            default:
              output<<key_int<<" [shape=circle,color=black]"<<std::endl;
              //fprintf(output,"\"%lo\" [shape=circle,color=black]\n",
              //  key.truncate_value(17));
              break;
          }
          output<<std::dec;
        }
      } 
    }
    output<<"}"<<std::endl;
    output.close();
  }

  /**
 * @brief      Compute the local range of particles
 * range 0 = min range 1 = max
 *
 * @param      bodies  The bodies
 * @param      range   The range
 */
  void 
  local_range(
    std::vector<std::pair<entity_key_t,body>>& bodies,
    std::array<point_t,2>& range)
  {
    
    range[1] = bodies.back().second.coordinates();
    range[0] = bodies.back().second.coordinates();
    
    for(auto bi: bodies){
      for(size_t i=0;i<gdimension;++i){
        if(bi.second.coordinates()[i]>range[1][i])
          range[1][i] = bi.second.coordinates()[i];
        if(bi.second.coordinates()[i]<range[0][i])
          range[0][i] = bi.second.coordinates()[i];
      }
    }
  }


}; // utils

#endif // _mpisph_utils_
