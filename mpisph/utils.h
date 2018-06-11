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
 * @date April 2017
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

}; // utils

#endif // _mpisph_utils_