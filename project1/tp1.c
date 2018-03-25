#include <mpi.h>
#include <stdio.h>

int main(int argc, char *argv[])
{
  int r,n,pname;
  char h[MPI_MAX_PROCESSOR_NAME];
  
  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &r);
  MPI_Comm_size (MPI_COMM_WORLD, &n);
  MPI_Get_processor_name (h, &pname);
  printf("hello world from process %d of %d on %s\n", r, n, h);
  MPI_Finalize();
  return 0;
}
