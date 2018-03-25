#include <mpi.h>
#include <math.h>
#include <stdio.h>

#define BLOCK_LOW(id,p,m) ((id)*(m)/(p))
#define BLOCK_HIGH(id,p,m) (BLOCK_LOW((id)+1,p,m)-1)
#define BLOCK_SIZE(id,p,m) (BLOCK_LOW((id)+1,p,m)-BLOCK_LOW(id,p,m))
#define BLOCK_OWNER(id,p,m) (((p)*((id)+1)-1)/(m))
#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))
#define CACHE_LOW(k, c) ((k)*(c)*2+3)
#define CACHE_HIGH(k, c) (CACHE_LOW((k)+1, c)-2)

/* Command line example
Compile: mpicc -o findprime findprime.c
Execute: mpirun -n 5 ./findprime 1000 4
*/

int main(int argc, char *argv[])
{
  int count; /* local prime count */
  double elapsed_time; /* parallel execution time */
  int first; /* index of first multiple */
  int global_count; /* global prime count */
  int high_value; /* highest value on this proc */
  int i, j, k, flag;
  int rem; /* length between low_value and first */
  int id; /* process ID number */
  // int index; /* index of current prime */
  int low_value; /* lowest value on this proc */
  char *marked; /* portion of 2,...,n */
  int *prime_list; /* store all primes */
  int n0, n; /* sieving from 3,...,n0 (n is largest odd number) */
  int m; /* number of odd numbers */
  int p; /* number of processes */
  // int proc0_size; /* size of proc 0's subarray */
  int prime; /* current prime */
  int size; /* elements in 'marked' */
  int cache_size; /* one block's size */
  int prime_number; /* number of prime numbers */
  int cache_low; /* lowest value in cache */
  int cache_high; /* highest value in cache*/

  MPI_Init(&argc, &argv);

  /* start the timer */
  MPI_Barrier(MPI_COMM_WORLD); // block until all processes in the communicator have achieved this routine
  elapsed_time = - MPI_Wtime();

  MPI_Comm_rank(MPI_COMM_WORLD, &id);
  MPI_Comm_size(MPI_COMM_WORLD, &p);

  /* check if assign the upper bound */
  if(argc == 1){
    if (!id) printf("Command line: %s have not assigned the upper bound\n", argv[0]);
    MPI_Finalize();
    exit(1);
  }
  n0 = atoi(argv[1]);
  if (n0%2 == 0)  n = n0-1;
  else  n = n0;
  m = (n-1)/2;

  cache_size = atoi(argv[2]);

  /* figure out this process's share of the array,
    as well as the integers represented by the first
    and last array elements */
  low_value = 3 + BLOCK_LOW(id, p, m) * 2;
  high_value = 3 + BLOCK_HIGH(id, p, m) * 2;
  size = BLOCK_SIZE(id, p, m);
  printf("\n\nProcess %d checks odd numbers from %d to %d\n", id, low_value, high_value);
  printf("Primes found from process %d are", id);

  /* Allocate this process's share of the array */
  marked = (char *) malloc (size);
  if (marked == NULL) {
    printf("Cannot allocate enough memory\n");
    MPI_Finalize();
    exit(1);
  }

  for (i=0; i<size; i++)  marked[i] = 0;
  // if (!id)  index = 0; // root process

  /* calculate prime list */
  prime_list = (int *) malloc (n);
  prime_list[0] = 3;
  k = 1;
  for (i=5; i * i <= n; i+=2) {
    flag = 0;
    for (j=3; j * j < MIN(n, i); j+=2) {
      if (i%j == 0) {
          flag = 1; /* composite number*/
          break;
      }
    }
    if (flag == 0) {
      prime_list[k++] = i;
    }
  }
  prime_number = k;

  /* initialize number block/cache for checking */
  cache_low = low_value;
  cache_high = cache_low + 2 * cache_size;

  /* mark composite numbers as 1 */
  while (cache_low <= high_value) {
  // block loop
    // printf("cache: %d - %d\n", cache_low, cache_high);
    for (i=0; i<prime_number; i++) {
    // prime loop
      prime = prime_list[i];
      // printf("prime: %d; ", prime);
      if (prime * prime > cache_low) {
        first = (prime * prime - low_value) / 2;
      }
      else {
        if (!(cache_low % prime))  first = (cache_low - low_value) / 2;
        else {
          rem = prime - (cache_low % prime);
          first = ((rem % 2) * prime + rem) / 2 + (cache_low - low_value) / 2; // if first is even, rem % 2 = 1
        }
      }

      while (low_value + 2 * first <= MIN(cache_high, high_value) && low_value + 2 * first >= MAX(cache_low, low_value)) {
        // printf("composite: %d; ", low_value + 2 * first);
        marked[first] = 1;
        first += prime;
        // printf("next: %d", low_value + 2 * first);
      }
      // printf("\n");
    }

    // update block
    cache_low = cache_high + 2;
    cache_high = cache_low + 2 * cache_size;

    // printf("\n");
  }

  count = 0;
  if (!id) {
    count = 1; // add 2 in the final prime list
    // printf(" 2 ");
  }
  for (i=0; i<size; i++) {
    if (!marked[i]) {
      count++;
      printf(" %d ", low_value+i*2);
    }
  }
  printf("\n\n");
  MPI_Reduce(&count, &global_count, 1, MPI_INT, MPI_SUM,
           0, MPI_COMM_WORLD);

  /* stop timer */
  elapsed_time += MPI_Wtime();

  /* print the results */
  if (!id) {
    printf("%d primes are less than or equal to %d\n",
          global_count, n0);
    printf("Total elapsed time: %10.6f\n", elapsed_time);
  }

  MPI_Finalize();
  return 0;


}
