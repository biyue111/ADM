#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <mpi.h>
//#include "MyMPI."

#define DEBUG 1
#define SUBMATRIX_HEIGHT(m, p) (m + m%p)/p
#define SUBMATRIX_WIDTH(n, p) (n + n%p)/p

int malloc2ddouble(double ***array, int n, int m)
{
    double *p = (double *)malloc(n*m*sizeof(double));
    if (!p) return -1;

    (*array) = (double **)malloc(n*sizeof(double*));
    if (!(*array)) 
    {
        free(p);
        return -1;
    }

    int i;
    for (i=0; i<n; i++)
        (*array)[i] = &(p[i*m]);
    return 0;
}

int main(int argc, char *argv[])
{
    MPI_Comm cart_comm; //Cartesian topology communicator
    int p; // number of processes
    int id; // process ID number 
    int root; // the root process ID
    int grid_coords[2]; // process coordinates
    int periodic[2];
    int size[2]; // number of rows of blocks (==sqrt(p))
    
    double **A;
    int dimA[2]; //dimension of matrix A
    double **B;
    int dimB[2]; //dimension of matrix B
    double **subA; //submatrix of A
    int dimsubA[2];
    double **subB;
    int dimsubB[2];

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    /* Create cartesien topology communicator */
    size[0] = size[1] = (int)sqrt((double)p);
    periodic[0] = periodic[1] = 0;
    MPI_Dims_create(p, 2, size);
    MPI_Cart_create(MPI_COMM_WORLD, 2, size, periodic, 1, &cart_comm);
    MPI_Cart_coords(cart_comm, id, 2, grid_coords);
    /* find the root ID */
    int root_coord[2] ;
    root_coord[0] = root_coord[1] = 0;
    root = MPI_Cart_rank(cart_comm, root_coord, &root);
    
    /* The first process read the input file*/
    if(id == root)
    {
    
        FILE* file = fopen ("input.txt", "r");

        fscanf(file, "%d", &dimA[0]);
        fscanf(file, "%d", &dimA[1]);
        MPI_Bcast(&dimA[0], 1, MPI_INT, root, cart_comm);
        MPI_Bcast(&dimA[1], 1, MPI_INT, root, cart_comm);

        int i,j;
        malloc2ddouble(&A, dimA[0], dimA[1]);

        for (i=0; i<dimA[0]; i++)
            for (j=0; j<dimA[1]; j++)
            {
                fscanf(file, "%lf", &A[i][j]);
            }
#ifdef DEBUG
        printf("The matrix A is:\n");
        for (i=0; i<dimA[0]; i++)
        {
            for(j=0; j<dimA[1]; j++)
                printf("%f", A[i][j]);
            printf("\n");
        }
#endif
    }
    /* Create subMatrix */
    MPI_Datatype subAtype;
    dimsubA[0] = SUBMATRIX_HEIGHT(dimA[0], size[0]);
    dimsubA[1] = SUBMATRIX_WIDTH(dimA[1], size[1]);
    malloc2ddouble(&subA, dimsubA[0], dimsubB[1]); 
    int starts[2] = {0,0};
    MPI_Type_create_subarray(2, dimA, dimsubA, starts, MPI_ORDER_C, MPI_INT, &subAtype);
    MPI_Type_commit(&subAtype);

    double *Aptr = NULL;
    if (id == 0) Aptr = &(A[0][0]);

    MPI_Scatter(Aptr, 1, subAtype, &(subA[0][0]), dimsubA[0]*dimsubA[1], MPI_INT, root, cart_comm);

#ifdef DEBUG
    /* Each process print received submatrix of A */
    int k = 0;
    for (k=0; k<p; k++)
    {
        if (k = id)
        {
            printf("The submatrix in process %d is:", k);
            for (int i=0; i<dimsubA[0]; i++)
            {
                for (int j=0; j<dimsubA[1]; j++)
                    printf("%f ", subA[i][j]);
                printf("\n");
            }
        }
        MPI_Barrier(cart_comm);

    }
#endif
}
