#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <mpi.h>
//#include "MyMPI."

#define DEBUG 1
#define SUBMATRIX_HEIGHT(m, p) (m +p + m%p)/p
#define SUBMATRIX_WIDTH(n, p) (n + p + n%p)/p

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

int find_counts_displs(int *counts, int *displs, int *size, int *dimsub)
{
    int k = 0;
    for (int i=0; i<size[0]; i++)
    {
       for (int j=0; j<size[1];j++)
       {
           counts[k] = 1;
           displs[k] = i * size[1] * dimsub[0] + j;
           k++;
       }
    }
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
    MPI_Cart_create(MPI_COMM_WORLD, 2, size, periodic, 0, &cart_comm); // Use the old ID
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

        dimB[0] = dimA[1];
        dimB[1] = dimA[0];

#ifdef DEBUG
    printf("dimA: %d %d\n", dimA[0], dimA[1]);
    printf("dimB: %d %d\n", dimB[0], dimB[1]);
#endif

        int i,j;
        malloc2ddouble(&A, dimA[0], dimA[1]);

        for (i=0; i<dimA[0]; i++)
            for (j=0; j<dimA[1]; j++)
            {
                fscanf(file, "%lf", &A[i][j]);
            }

        malloc2ddouble(&B, dimB[0], dimB[1]);
        for (i=0; i<dimB[0];i++)
            for (j=0; j<dimB[1]; j++)
            {
                fscanf(file, "%lf", &B[i][j]);
            }  
        fclose(file);
    }
    MPI_Bcast(&dimA[0], 1, MPI_INT, root, cart_comm);
    MPI_Bcast(&dimA[1], 1, MPI_INT, root, cart_comm);
    MPI_Bcast(&dimB[0], 1, MPI_INT, root, cart_comm);
    MPI_Bcast(&dimB[1], 1, MPI_INT, root, cart_comm);
    /* Create subMatrix */
    MPI_Datatype subAtype, r_subAtype;
    MPI_Datatype subBtype, r_subBtype;
    dimsubA[0] = SUBMATRIX_HEIGHT(dimA[0], size[0]);
    dimsubA[1] = SUBMATRIX_HEIGHT(dimA[1], size[1]);
    dimsubB[0] = SUBMATRIX_HEIGHT(dimB[0], size[0]);
    dimsubB[1] = SUBMATRIX_WIDTH(dimB[1], size[1]);

#ifdef DEBUG
    printf("Dimension of subA in process %d: %d %d\n", id, dimsubA[0], dimsubA[1]);
    printf("Dimension of subB in process %d: %d %d\n", id, dimsubB[0], dimsubB[1]);
#endif

    malloc2ddouble(&subA, dimsubA[0], dimsubA[1]); 
    int starts[2] = {0,0};
    MPI_Type_create_subarray(2, dimA, dimsubA, starts, MPI_ORDER_C, MPI_DOUBLE, &subAtype);
    MPI_Type_create_resized(subAtype, 0, dimsubA[1]*sizeof(double), &r_subAtype);
    MPI_Type_commit(&r_subAtype);

    malloc2ddouble(&subB, dimsubB[0], dimsubB[1]);
    MPI_Type_create_subarray(2, dimB, dimsubB, starts, MPI_ORDER_C, MPI_DOUBLE, &subBtype);
    MPI_Type_create_resized(subBtype, 0, dimsubB[1]*sizeof(double), &r_subBtype);
    MPI_Type_commit(&r_subBtype);

    double *Aptr = NULL;
    double *Bptr = NULL;
    if (id == root) Aptr = &(A[0][0]);
    if (id == root) Bptr = &(B[0][0]);

#ifdef DEBUG
    /* Pring matrix A and matrix B*/
    if (id == root)
    {
        printf("the matrix a is:\n");
        for (int i=0; i<dimA[0]; i++)
        {
            for(int j=0; j<dimA[1]; j++)
                printf("%f ", A[i][j]);
            printf("\n");
        }
    }

    if (id == root)
    {
        printf("the matrix a is:\n");
        for (int i=0; i<dimB[0]; i++)
        {
            for(int j=0; j<dimB[1]; j++)
                printf("%f ", B[i][j]);
            printf("\n");
        }
    }
#endif
    int counts[p];
    int displs[p];
    find_counts_displs(counts, displs, size, dimsubA);
    MPI_Scatterv(Aptr, counts, displs, r_subAtype, &(subA[0][0]), dimsubA[0]*dimsubA[1], MPI_DOUBLE, root, cart_comm);

    find_counts_displs(counts, displs, size, dimsubB);
    MPI_Scatterv(Bptr, counts, displs, r_subBtype, &(subB[0][0]), dimsubB[0]*dimsubB[1], MPI_DOUBLE, root, cart_comm);

#ifdef DEBUG
    /* Each process print received submatrix of A */
    MPI_Barrier(cart_comm);
    for (int k=0; k<p; k++)
    {
        if (k == id)
        {
            printf("The submatrix A in process (%d %d)is:\n", grid_coords[0], grid_coords[1]);
            for (int i=0; i<dimsubA[0]; i++)
            {
                for (int j=0; j<dimsubA[1]; j++)
                    printf("%f ", subA[i][j]);
                printf("\n");
            }
        }
        MPI_Barrier(cart_comm);
    }

    /* Each process print received submatrix of B */
    MPI_Barrier(cart_comm);
    for (int k=0; k<p; k++)
    {
        if (k == id)
        {
            printf("The submatrix B in process (%d %d)is:\n", grid_coords[0], grid_coords[1]);
            for (int i=0; i<dimsubB[0]; i++)
            {
                for (int j=0; j<dimsubB[1]; j++)
                    printf("%f ", subB[i][j]);
                printf("\n");
            }
        }
        MPI_Barrier(cart_comm);
    }

#endif
    MPI_Finalize();
    return 0;
}
