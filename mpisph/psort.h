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

#pragma once

#include "mpi.h"
#include <numeric>
#include <vector>

/**
 * MPI distributed sort
 * "A Novel Parallel Sorting Algorithm for Contemporary Architectues"
 * David R. Cheng & al.
 * See Combblas
 **/
namespace psort {

class Split
{
public:
  template<typename _Iterator, typename _Compare>
  void split(_Iterator first,
    _Iterator last,
    int * dist,
    _Compare comp,
    std::vector<std::vector<int>> & right_ends,
    MPI_Datatype & MPI_valueType) {
    typedef typename std::iterator_traits<_Iterator>::value_type _ValueType;

    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int n_real = size;
    for(int i = 0; i < size; ++i)
      if(dist[i] == 0) {
        n_real = i;
        break;
      }

    std::copy(dist, dist + size, right_ends[size].begin());

    // union of [0, right_end[i+1]) on each processor produces dist[i] total
    // values
    std::vector<int> targets_v(size - 1);
    int * targets = &targets_v.at(0);
    std::partial_sum(dist, dist + (size - 1), targets);

    // keep a list of ranges, trying to "activate" them at each branch
    std::vector<std::pair<_Iterator, _Iterator>> d_ranges(size - 1);
    std::vector<std::pair<int *, int *>> t_ranges(size - 1);
    d_ranges[0] = std::pair<_Iterator, _Iterator>(first, last);
    t_ranges[0] = std::pair<int *, int *>(targets, targets + (size - 1));

    // invariant: subdist[i][rank] == d_ranges[i].second - d_ranges[i].first
    // amount of data each proc still has in the search
    std::vector<std::vector<int>> subdist(size - 1, std::vector<int>(size));
    std::copy(dist, dist + size, subdist[0].begin());

    // for each processor, d_ranges - first
    std::vector<std::vector<int>> outleft(size - 1, std::vector<int>(size, 0));

    for(int n_act = 1; n_act > 0;) {
      for(int k = 0; k < n_act; ++k) {
        assert(subdist[k][rank] == d_ranges[k].second - d_ranges[k].first);
      }
      //------- generate n_act guesses
      // for the allgather, make a flat array of nproc chunks, each with n_act
      // elts
      _ValueType * mymedians = new _ValueType[n_act];
      _ValueType * medians = new _ValueType[size * n_act];
      for(int k = 0; k < n_act; ++k) {
        _ValueType * ptr = &d_ranges[k].first[0];
        int index = subdist[k][rank] / 2;
        mymedians[k] = ptr[index];
      } // for
      MPI_Allgather(mymedians, n_act, MPI_valueType, medians, n_act,
        MPI_valueType, MPI_COMM_WORLD);
      delete[] mymedians;

      // compute the weighted median of medians
      std::vector<_ValueType> queries(n_act);

      for(int k = 0; k < n_act; ++k) {
        std::vector<int> ms_perm_v(n_real);
        int * ms_perm = &ms_perm_v.at(0);
        for(int i = 0; i < n_real; ++i)
          ms_perm[i] = i * n_act + k;
        std::sort(ms_perm, ms_perm + n_real,
          PermCompare<_ValueType, _Compare>(medians, comp));
        int mid = accumulate(subdist[k].begin(), subdist[k].end(), 0) / 2;
        int query_ind = -1;
        for(int i = 0; i < n_real; ++i) {
          if(subdist[k][ms_perm[i] / n_act] == 0)
            continue;
          mid -= subdist[k][ms_perm[i] / n_act];
          if(mid <= 0) {
            query_ind = ms_perm[i];
            break;
          }
        } // for

        assert(query_ind >= 0);
        queries[k] = medians[query_ind];
      } // for
      delete[] medians;

      //------- find min and max ranks of the guesses
      std::vector<int> ind_local_v(2 * n_act);
      int * ind_local = &ind_local_v.at(0);
      for(int k = 0; k < n_act; ++k) {
        std::pair<_Iterator, _Iterator> ind_local_p = std::equal_range(
          d_ranges[k].first, d_ranges[k].second, queries[k], comp);

        ind_local[2 * k] = ind_local_p.first - first;
        ind_local[2 * k + 1] = ind_local_p.second - first;
      } // for

      std::vector<int> ind_all_v(2 * n_act * size);
      int * ind_all = &ind_all_v.at(0);
      MPI_Allgather(ind_local, 2 * n_act, MPI_INT, ind_all, 2 * n_act, MPI_INT,
        MPI_COMM_WORLD);
      // sum to get the global range of indices
      std::vector<std::pair<int, int>> ind_global(n_act);
      for(int k = 0; k < n_act; ++k) {
        ind_global[k] = std::make_pair(0, 0);
        for(int i = 0; i < size; ++i) {
          ind_global[k].first += ind_all[2 * (i * n_act + k)];
          ind_global[k].second += ind_all[2 * (i * n_act + k) + 1];
        } // for
      } // for

      // state to pass on to next iteration
      std::vector<std::pair<_Iterator, _Iterator>> d_ranges_x(size - 1);
      std::vector<std::pair<int *, int *>> t_ranges_x(size - 1);
      std::vector<std::vector<int>> subdist_x(size - 1, std::vector<int>(size));
      std::vector<std::vector<int>> outleft_x(
        size - 1, std::vector<int>(size, 0));
      int n_act_x = 0;

      for(int k = 0; k < n_act; ++k) {
        int * split_low = std::lower_bound(
          t_ranges[k].first, t_ranges[k].second, ind_global[k].first);
        int * split_high = std::upper_bound(
          t_ranges[k].first, t_ranges[k].second, ind_global[k].second);

        // iterate over targets we hit
        for(int * s = split_low; s != split_high; ++s) {
          assert(*s > 0);
          // a bit sloppy: if more than one target in range, excess won't zero
          // out
          int excess = *s - ind_global[k].first;
          // low procs to high take excess for stability
          for(int i = 0; i < size; ++i) {
            int amount = std::min(ind_all[2 * (i * n_act + k)] + excess,
              ind_all[2 * (i * n_act + k) + 1]);
            right_ends[(s - targets) + 1][i] = amount;
            excess -= amount - ind_all[2 * (i * n_act + k)];
          } // for
        } // for

        if((split_low - t_ranges[k].first) > 0) {
          t_ranges_x[n_act_x] = std::make_pair(t_ranges[k].first, split_low);
          // lop off local_ind_low..end
          d_ranges_x[n_act_x] =
            std::make_pair(d_ranges[k].first, first + ind_local[2 * k]);
          for(int i = 0; i < size; ++i) {
            subdist_x[n_act_x][i] =
              ind_all[2 * (i * n_act + k)] - outleft[k][i];
            outleft_x[n_act_x][i] = outleft[k][i];
          } // for
          ++n_act_x;
        } // if

        if((t_ranges[k].second - split_high) > 0) {
          t_ranges_x[n_act_x] = std::make_pair(split_high, t_ranges[k].second);
          // lop off begin..local_ind_high
          d_ranges_x[n_act_x] =
            std::make_pair(first + ind_local[2 * k + 1], d_ranges[k].second);
          for(int i = 0; i < size; ++i) {
            subdist_x[n_act_x][i] =
              outleft[k][i] + subdist[k][i] - ind_all[2 * (i * n_act + k) + 1];
            outleft_x[n_act_x][i] = ind_all[2 * (i * n_act + k) + 1];
          } // for
          ++n_act_x;
        } // if
      } // for

      t_ranges = t_ranges_x;
      d_ranges = d_ranges_x;
      subdist = subdist_x;
      outleft = outleft_x;
      n_act = n_act_x;
    } // for
  } // split

private:
  template<typename T, typename _Compare>
  class PermCompare
  {
  private:
    T * weights;
    _Compare comp;

  public:
    PermCompare(T * w, _Compare c) : weights(w), comp(c) {}
    bool operator()(int a, int b) {
      return comp(weights[a], weights[b]);
    }
  };
}; // class

template<typename TYPE, typename _Compare>
void
psort(std::vector<TYPE> & vec, _Compare comp, int * dist_in) {

  typename std::vector<TYPE>::iterator first = vec.begin();
  typename std::vector<TYPE>::iterator last = vec.end();

  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  MPI_Datatype MPI_valueType;
  MPI_Type_contiguous(sizeof(TYPE), MPI_CHAR, &MPI_valueType);
  MPI_Type_commit(&MPI_valueType); // ABAB: Any type committed needs to be freed
                                   // to claim storage

  int * dist = new int[size];
  for(int i = 0; i < size; ++i)
    dist[i] = dist_in[i];

  std::sort(first, last, comp);

  // For one rank, no work
  if(size == 1) {
    MPI_Type_free(&MPI_valueType);
    return;
  }

  // Find splitters
  std::vector<std::vector<int>> right_ends(size + 1, std::vector<int>(size, 0));
  Split mysplit;
  mysplit.split(first, last, dist, comp, right_ends, MPI_valueType);

  // Communicate to destination
  int * boundaries = new int[size + 1];

  // Should be _Distance, but MPI wants ints
  char errMsg[] = "32-bit limit for MPI has overflowed";
  int n_loc_ = last - first;
  if(n_loc_ > INT_MAX)
    throw std::overflow_error(errMsg);
  int n_loc = static_cast<int>(n_loc_);
  std::vector<TYPE> trans_data(n_loc);

  // Calculate the counts for redistributing data
  int * send_counts = new int[size];
  int * send_disps = new int[size];
  for(int i = 0; i < size; ++i) {
    int scount = right_ends[i + 1][rank] - right_ends[i][rank];
    if(scount > INT_MAX)
      throw std::overflow_error(errMsg);
    send_counts[i] = static_cast<int>(scount);
  }
  send_disps[0] = 0;
  std::partial_sum(send_counts, send_counts + size - 1, send_disps + 1);

  int * recv_counts = new int[size];
  int * recv_disps = new int[size];
  for(int i = 0; i < size; ++i) {
    int rcount = right_ends[rank + 1][i] - right_ends[rank][i];
    if(rcount > INT_MAX)
      throw std::overflow_error(errMsg);
    recv_counts[i] = static_cast<int>(rcount);
  }

  recv_disps[0] = 0;
  std::partial_sum(recv_counts, recv_counts + size - 1, recv_disps + 1);

  assert(std::accumulate(recv_counts, recv_counts + size, 0) == n_loc);

  // Do the transpose
  MPI_Alltoallv(&first[0], send_counts, send_disps, MPI_valueType,
    &trans_data[0], recv_counts, recv_disps, MPI_valueType, MPI_COMM_WORLD);

  for(int i = 0; i < size; ++i)
    boundaries[i] = (int)recv_disps[i];
  boundaries[size] = (int)n_loc; // for the merging

  delete[] recv_counts;
  delete[] recv_disps;
  delete[] send_counts;
  delete[] send_disps;

  std::sort(trans_data.begin(), trans_data.end(), comp);
  // Merge streams from all processors
  // std::sort(first, last, comp);
  vec = trans_data;

  delete[] boundaries;
  delete[] dist;
  // delete [] trans_data;
  MPI_Type_free(&MPI_valueType);

  // Finish
  return;
}
} // namespace psort