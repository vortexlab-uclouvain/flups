#include "mpi.h"
#include <stdlib.h>
#include <stdio.h>

// todo handle nf != 1

inline static size_t localindex(size_t i0, size_t i1, size_t i2, const size_t n[3], int axis) {
  size_t i[3] = {i0,i1,i2};
  int axis1 = (axis+1)%3;
  int axis2 = (axis+2)%3;
  return (i[axis]*n[axis1]+i[axis1])*n[axis2]+i[axis2];
}

void reorder_mpi(size_t n[3], size_t nf, int nproc0[3], int axis0, int nproc1[3], int axis1, double *v) {
  size_t n0[3],n1[3];
  for (int i = 0; i < 3; ++i) {
    n0[i] = n[i]/nproc0[i];
    n1[i] = n[i]/nproc1[i];
  }
  size_t nlocal = n0[0]*n0[1]*n0[2];
  double *bufsend = malloc(sizeof(double)*nf*nlocal);
  double *bufrecv = malloc(sizeof(double)*nf*nlocal);
  int rank, comm_size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&comm_size);
  int rankd0[3] = {
    rank/(nproc0[1]*nproc0[2]),
    (rank%(nproc0[1]*nproc0[2]))/nproc0[2],
    rank%nproc0[2]
  };
  int rankd1[3] = {
    rank/(nproc1[1]*nproc1[2]),
    (rank%(nproc1[1]*nproc1[2]))/nproc1[2],
    rank%nproc1[2]
  };
  int *nsend = malloc((comm_size+1)*sizeof(int));
  int *nrecv = malloc((comm_size+1)*sizeof(int));
  int *ssend = malloc((comm_size+1)*sizeof(int));
  int *srecv = malloc((comm_size+1)*sizeof(int));
  for (int i = 0 ; i < comm_size+1; ++i) {
    nsend[i] = 0;
    nrecv[i] = 0;
    ssend[i] = 0;
    srecv[i] = 0;
  }
  int destrankd[3];
  for (size_t i0 = rankd0[0]*n0[0]; i0 < (rankd0[0]+1)*n0[0]; ++i0) {
    destrankd[0] = i0/n1[0];
    for (size_t i1 = rankd0[1]*n0[1]; i1 < (rankd0[1]+1)*n0[1]; ++i1) {
      destrankd[1] = i1/n1[1];
      for (size_t i2 = rankd0[2]*n0[2]; i2 < (rankd0[2]+1)*n0[2]; ++i2) {
        destrankd[2] = i2/n1[2];
        int destrank = destrankd[0]*(nproc1[1]*nproc1[2])+destrankd[1]*nproc1[2]+destrankd[2];
        nsend[destrank] ++;
      }
    }
  }
  MPI_Alltoall(nsend,1,MPI_INT,nrecv,1,MPI_INT,MPI_COMM_WORLD);
  for (int i = 0; i < comm_size; ++i) {
    ssend[i+1] = ssend[i] + nsend[i];
    srecv[i+1] = srecv[i] + nrecv[i];
    nsend[i] = 0;
  }
  for (size_t i0 = 0; i0 < n0[0]; ++i0) {
    destrankd[0] = (rankd0[0]*n0[0]+i0)/n1[0];
    for (size_t i1 = 0; i1 < n0[1]; ++i1) {
      destrankd[1] = (rankd0[1]*n0[1]+i1)/n1[1];
      for (size_t i2 = 0; i2 < n0[2]; ++i2) {
        destrankd[2] = (rankd0[2]*n0[2]+i2)/n1[2];
        int destrank = destrankd[0]*(nproc1[1]*nproc1[2])+destrankd[1]*nproc1[2]+destrankd[2];
        bufsend[ssend[destrank]+(nsend[destrank]++)] = v[localindex(i0,i1,i2,n0,axis0)];
      }
    }
  }
  MPI_Alltoallv(bufsend,nsend,ssend,MPI_DOUBLE,bufrecv,nrecv,srecv,MPI_DOUBLE,MPI_COMM_WORLD);
  int origrankd[3];
  for (size_t i0 = 0; i0 < n1[0]; ++i0) {
    origrankd[0] = (rankd1[0]*n1[0]+i0)/n0[0];
    for (size_t i1 = 0; i1 < n1[1]; ++i1) {
      origrankd[1] = (rankd1[1]*n1[1]+i1)/n0[1];
      for (size_t i2 = 0; i2 < n1[2]; ++i2) {
        origrankd[2] = (rankd1[2]*n1[2]+i2)/n0[2];
        int origrank = origrankd[0]*(nproc0[1]*nproc0[2])+origrankd[1]*nproc0[2]+origrankd[2];
        v[localindex(i0,i1,i2,n1,axis1)] = bufrecv[srecv[origrank]++];
      }
    }
  }

  free(bufsend);
  free(bufrecv);
  free(nsend);
  free(nrecv);
}

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  size_t n[3] = {8,8,8};
  int nproc[3] = {1,2,2};
  size_t nloc[3];
  for (int i = 0; i < 3; ++i){
    nloc[i] = n[i]/nproc[i];
  }
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  int rankd[3] = {
    rank/(nproc[1]*nproc[2]),
    (rank%(nproc[1]*nproc[2]))/nproc[2],
    rank%nproc[2]
  };
  double *v = malloc(sizeof(double)*nloc[0]*nloc[1]*nloc[2]);
  for (size_t i = 0; i < nloc[0]; i++) {
    for (size_t j = 0; j < nloc[1]; j++) {
      for (size_t k = 0; k < nloc[2]; k++) {
        double x = 1./n[0]*(i+rankd[0]*nloc[0]);
        double y = 1./n[1]*(j+rankd[1]*nloc[1]);
        double z = 1./n[2]*(k+rankd[2]*nloc[2]);
        v[(i*nloc[1]+j)*nloc[2]+k] = z;
      }
    }
  }
  printf("init\n");
  for (size_t i = 0; i < nloc[0]; i++) {
    for (size_t j = 0; j < nloc[1]; j++) {
      for (size_t k = 0; k < nloc[2]; k++) {
        printf("%lu %lu %lu : %g\n", i,j,k,v[(i*nloc[1]+j)*nloc[2]+k]);
      }
    }
  }
  int nproc1[3] = {2,1,2};
  reorder_mpi(n, 1, nproc, 0, nproc1, 1, v);
  printf("y alligned\n");
  for (size_t i = 0; i < nloc[0]; i++) {
    for (size_t j = 0; j < nloc[1]; j++) {
      for (size_t k = 0; k < nloc[2]; k++) {
        printf("%lu %lu %lu : %g\n", i,j,k,v[(i*nloc[1]+j)*nloc[2]+k]);
      }
    }
  }
  int nproc2[3] = {2,2,1};
  reorder_mpi(n, 1, nproc1, 1, nproc2, 2, v);
  printf("z alligned\n");
  for (size_t i = 0; i < nloc[0]; i++) {
    for (size_t j = 0; j < nloc[1]; j++) {
      for (size_t k = 0; k < nloc[2]; k++) {
        printf("%lu %lu %lu : %g\n", i,j,k,v[(i*nloc[1]+j)*nloc[2]+k]);
      }
    }
  }
  free(v);
  MPI_Finalize();
}

