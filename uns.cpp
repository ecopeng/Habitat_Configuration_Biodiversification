/*

Simulations

*/

#include <mpi.h>
#include "Funs.cpp"
#include <stdlib.h>

int main(int argc, char** argv) {
  int rank, n_ranks;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &n_ranks);
  div_t res_0 = div(rank, 37 * 4);
  int r = res_0.quot, rem_ = res_0.rem;
  div_t res_1 = div(rem_, 37);
  int iNet = res_1.quot, iRho = res_1.rem;
  S_uns(iNet, r, iRho);
  MPI_Finalize();
  return 0;
}