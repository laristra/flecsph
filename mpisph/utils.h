/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Triad National Security, LLC
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

#include <numeric>

#include "tree.h"

// Local version of assert to handle MPI abort
#define mpi_assert(assertion)                                                  \
  ((assertion) ? true                                                          \
               : (fprintf(stderr, "Failed assertion in %s in %d\n", __FILE__,  \
                          __LINE__) &&                                         \
                  fflush(stderr) && MPI_Abort(MPI_COMM_WORLD, 1)));

namespace mpi_utils {

/**
 * @brief Simple version of all gather
 * Send the size of arrays and then the all gather operation
 */
template <typename M>
void mpi_allgatherv(const std::vector<M> &send, std::vector<M> &recv,
                    std::vector<int> &count = 0) {
  int size, rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  count.clear();
  count.resize(size);

  int my_count = send.size();
  // Gather the total size
  MPI_Allgather(&my_count, 1, MPI_INT, &count[0], 1, MPI_INT, MPI_COMM_WORLD);
  // Convert to byte
  std::vector<int> count_byte(size);
  std::vector<int> offset_byte(size);
  int64_t total = 0L;

#pragma omp parallel for reduction(+ : total)
  for (int i = 0; i < size; ++i) {
    total += count[i];
    count_byte[i] = count[i] * sizeof(M);
  }
  std::partial_sum(count_byte.begin(), count_byte.end(), &offset_byte[0]);
  offset_byte.insert(offset_byte.begin(), 0);
  recv.resize(total);

  MPI_Allgatherv(&send[0], count_byte[rank], MPI_BYTE, &recv[0], &count_byte[0],
                 &offset_byte[0], MPI_BYTE, MPI_COMM_WORLD);
}

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
template <typename M>
void mpi_alltoallv(std::vector<int> sendcount, std::vector<M> &sendbuffer,
                   std::vector<M> &recvbuffer) {
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  std::vector<int> recvcount(size);
  std::vector<int> recvoffsets(size);
  std::vector<int> sendoffsets(size);

  // Exchange the send count
  MPI_Alltoall(&sendcount[0], 1, MPI_INT, &recvcount[0], 1, MPI_INT,
               MPI_COMM_WORLD);

  // Generate the send and recv offsets
  std::partial_sum(recvcount.begin(), recvcount.end(), &recvoffsets[0]);
  // As we need an exscan, add a zero
  recvoffsets.insert(recvoffsets.begin(), 0);

  // Then send offsets
  std::partial_sum(sendcount.begin(), sendcount.end(), &sendoffsets[0]);
  // As we need an exscan, add a zero
  sendoffsets.insert(sendoffsets.begin(), 0);

  // Set the recvbuffer to the right size
  recvbuffer.resize(recvoffsets.back());

  // Trnaform the offsets for bytes
#pragma omp parallel for
  for (int i = 0; i < size; ++i) {
    sendcount[i] *= sizeof(M);
    assert(sendcount[i] >= 0);
    recvcount[i] *= sizeof(M);
    assert(recvcount[i] >= 0);
    sendoffsets[i] *= sizeof(M);
    assert(sendoffsets[i] >= 0);
    recvoffsets[i] *= sizeof(M);
    assert(recvoffsets[i] >= 0);
  } // for

  // Use this array for the global buckets communication
  MPI_Alltoallv(&sendbuffer[0], &sendcount[0], &sendoffsets[0], MPI_BYTE,
                &recvbuffer[0], &recvcount[0], &recvoffsets[0], MPI_BYTE,
                MPI_COMM_WORLD);
}

template <typename M>
void mpi_alltoallv_p2p(std::vector<int> &sendcount, std::vector<M> &sendbuffer,
                       std::vector<M> &recvbuffer) {
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  std::vector<int> recvcount(size), recvoffsets(size), sendoffsets(size);
  // Exchange the send count
  MPI_Alltoall(&sendcount[0], 1, MPI_INT, &recvcount[0], 1, MPI_INT,
               MPI_COMM_WORLD);
  std::partial_sum(recvcount.begin(), recvcount.end(), &recvoffsets[0]);
  recvoffsets.insert(recvoffsets.begin(), 0);
  std::partial_sum(sendcount.begin(), sendcount.end(), &sendoffsets[0]);
  sendoffsets.insert(sendoffsets.begin(), 0);
  // Set the recvbuffer to the right size
  recvbuffer.resize(recvoffsets.back());
  // Transform the offsets for bytes
#pragma omp parallel for
  for (int i = 0; i < size; ++i) {
    sendcount[i] *= sizeof(M);
    assert(sendcount[i] >= 0);
    recvcount[i] *= sizeof(M);
    assert(recvcount[i] >= 0);
    sendoffsets[i] *= sizeof(M);
    assert(sendoffsets[i] >= 0);
    recvoffsets[i] *= sizeof(M);
    assert(recvoffsets[i] >= 0);
  } // for
  std::vector<MPI_Status> status(size);
  std::vector<MPI_Request> request(size);
#pragma omp parallel for
  for (int i = 0; i < size; ++i) {
    if (sendcount[i] != 0) {
      char *start = (char *)&(sendbuffer[0]);
      MPI_Isend(start + sendoffsets[i], sendcount[i], MPI_BYTE, i, 0,
                MPI_COMM_WORLD, &request[i]);
    }
  }
#pragma omp parallel for
  for (int i = 0; i < size; ++i) {
    if (recvcount[i] != 0) {
      char *start = (char *)&(recvbuffer[0]);
      MPI_Recv(start + recvoffsets[i], recvcount[i], MPI_BYTE, i, MPI_ANY_TAG,
               MPI_COMM_WORLD, &status[i]);
    }
    if (sendcount[i] != 0) {
      MPI_Wait(&request[i], &status[i]);
    }
  }
} // mpi_alltoallv_p2p

template <typename M>
void mpi_alltoallv_p2p(std::vector<int> &sendcount,
                       std::vector<std::vector<M>> &sendbuffer,
                       std::vector<std::vector<M>> &recvbuffer) {
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  std::vector<int> recvcount(size), recvoffsets(size), sendoffsets(size);
  // Exchange the send count
  MPI_Alltoall(&sendcount[0], 1, MPI_INT, &recvcount[0], 1, MPI_INT,
               MPI_COMM_WORLD);
  // Set the recvbuffer to the right size
  // recvbuffer.resize(recvoffsets.back());
  // Transform the offsets for bytes
#pragma omp parallel for
  for (int i = 0; i < size; ++i) {
    recvbuffer[i].resize(recvcount[i]);
    sendcount[i] *= sizeof(M);
    assert(sendcount[i] >= 0);
    recvcount[i] *= sizeof(M);
    assert(recvcount[i] >= 0);
  } // for
  std::vector<MPI_Status> status(size);
  std::vector<MPI_Request> request(size);
#pragma omp parallel for
  for (int i = 0; i < size; ++i) {
    if (sendcount[i] != 0 && rank != i) {
      MPI_Isend(&(sendbuffer[i][0]), sendcount[i], MPI_BYTE, i, 0,
                MPI_COMM_WORLD, &request[i]);
    }
  }
#pragma omp parallel for
  for (int i = 0; i < size; ++i) {
    if (recvcount[i] != 0 && rank != i) {
      MPI_Recv(&(recvbuffer[i][0]), recvcount[i], MPI_BYTE, i, MPI_ANY_TAG,
               MPI_COMM_WORLD, &status[i]);
    }
    if (sendcount[i] != 0) {
      MPI_Wait(&request[i], &status[i]);
    }
  }
} // mpi_alltoallv_p2p

// MIN REDUCTION MPI -------------------

void reduce_min(double &value) {
  MPI_Allreduce(MPI_IN_PLACE, &value, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
}

void reduce_min(float &value) {
  MPI_Allreduce(MPI_IN_PLACE, &value, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
}

void reduce_min(uint64_t &value) {
  MPI_Allreduce(MPI_IN_PLACE, &value, 1, MPI_INT64_T, MPI_MIN, MPI_COMM_WORLD);
}

void reduce_min(int &value) {
  MPI_Allreduce(MPI_IN_PLACE, &value, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
}

// SUM REDUCTION MPI ----------------------

void reduce_sum(double &value) {
  MPI_Allreduce(MPI_IN_PLACE, &value, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

void reduce_sum(float &value) {
  MPI_Allreduce(MPI_IN_PLACE, &value, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
}

void reduce_sum(uint64_t &value) {
  MPI_Allreduce(MPI_IN_PLACE, &value, 1, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);
}

void reduce_sum(int &value) {
  MPI_Allreduce(MPI_IN_PLACE, &value, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
}

void reduce_sum(point_t &value) {
  MPI_Allreduce(MPI_IN_PLACE, &value[0], gdimension, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
}

// MAX REDUCTION MPI -------------------

void reduce_max(double &value) {
  MPI_Allreduce(MPI_IN_PLACE, &value, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
}

void reduce_max(float &value) {
  MPI_Allreduce(MPI_IN_PLACE, &value, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
}

void reduce_max(uint64_t &value) {
  MPI_Allreduce(MPI_IN_PLACE, &value, 1, MPI_INT64_T, MPI_MAX, MPI_COMM_WORLD);
}

void reduce_max(int &value) {
  MPI_Allreduce(MPI_IN_PLACE, &value, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
}

void output_branches_VTK(std::vector<range_t> &recv_branches,
                         std::vector<int> &count, size_t iteration) {
  char filename[255];
  sprintf(filename, "file_%lu.vtk", iteration);
  remove(filename);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  size_t npoints = recv_branches.size() * 4;
  size_t nlines = recv_branches.size() * 4;

  std::ostringstream oss_data;
  oss_data << "# vtk DataFile Version 3.0" << std::endl;
  oss_data << "vtk output" << std::endl;
  oss_data << "ASCII" << std::endl;
  oss_data << "DATASET POLYDATA" << std::endl;
  oss_data << "POINTS " << npoints << " float" << std::endl;
  // Print the points
  for (auto &b : recv_branches) {
    oss_data << b[0][0] << " " << b[0][1] << " 0.0" << std::endl;
    oss_data << b[0][0] << " " << b[1][1] << " 0.0" << std::endl;
    oss_data << b[1][0] << " " << b[1][1] << " 0.0" << std::endl;
    oss_data << b[1][0] << " " << b[0][1] << " 0.0" << std::endl;
  }
  oss_data << "LINES " << nlines << " " << nlines * 3 << std::endl;
  size_t start = 0;
  for (auto &b : recv_branches) {
    oss_data << "2 " << start << " " << start + 1 << std::endl;
    oss_data << "2 " << start + 1 << " " << start + 2 << std::endl;
    oss_data << "2 " << start + 2 << " " << start + 3 << std::endl;
    oss_data << "2 " << start + 3 << " " << start << std::endl;
    start += 4;
  }
  oss_data << "CELL_DATA " << nlines << std::endl
           << "scalars cellvar float" << std::endl;
  oss_data << "LOOKUP_TABLE default" << std::endl;
  // Output process data
  size_t r = 0;
  for (auto &c : count) {
    for (size_t v = 0; v < c; ++v) {
      oss_data << r << " " << r << " " << r << " " << r << std::endl;
    }
    ++r;
  }

  // Open file in append mode
  std::ofstream out(filename, std::ios_base::app);
  out << oss_data.str();
  out.close();
}

/**
 * @brief Communication to one other rank
 * @param [in] rank my rank for this communication
 * @param [in] partner my partner for this communication
 * @param [in] buffer The buffer used to send and store the data
 * @param [in] nsend The number of data to send, avoiding sending the received
 * ones
 * @param [in] last The place of the last data received
 * @return void
 */
template <typename T>
void mpi_one_to_one(const int rank, const int partner, std::vector<T> &buffer,
                    const int nsend, int &last) {
  int size;
  const int sizeofT = sizeof(T);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  MPI_Request request;
  if (partner < size) {
    // I send
    if (rank < partner) {
      // Send size
      MPI_Isend(&(buffer[0]), nsend * sizeofT, MPI_BYTE, partner, 1,
                MPI_COMM_WORLD, &request);
    } else {
      MPI_Status status;
      // Read the size of the message
      MPI_Probe(partner, 1, MPI_COMM_WORLD, &status);
      // Get the size
      int nrecv = 0;
      MPI_Get_count(&status, MPI_BYTE, &nrecv);
      buffer.resize(buffer.size() + nrecv / sizeofT);
      MPI_Recv(&(buffer[last]), nrecv, MPI_BYTE, partner, 1, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
      last = buffer.size();
    }
    // Other rank send
    if (rank > partner) {
      // Send size
      MPI_Isend(&(buffer[0]), nsend * sizeofT, MPI_BYTE, partner, 1,
                MPI_COMM_WORLD, &request);
    } else {
      MPI_Status status;
      // Read the size of the message
      MPI_Probe(partner, 1, MPI_COMM_WORLD, &status);
      // Get the size
      int nrecv = 0;
      MPI_Get_count(&status, MPI_BYTE, &nrecv);
      buffer.resize(buffer.size() + nrecv / sizeofT);
      MPI_Recv(&(buffer[last]), nrecv, MPI_BYTE, partner, 1, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
      last = buffer.size();
    }
    // Wait for request
    MPI_Status status;
    MPI_Wait(&request, &status);
  }
} // mpi_one_to_one

}; // namespace mpi_utils

#endif // _mpisph_utils_
