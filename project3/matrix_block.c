#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <mpi.h>
//#include "MyMPI."

#define DEBUG 1
#define INPUT_FILE "input.txt"
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

    for (int i=0; i<n; i++)
        (*array)[i] = &(p[i*m]);
    return 0;
}
int clear_matrix(double **matrix, int *dim)
{
    for (int i=0; i<dim[0];i++)
        for(int j=0; j<dim[1];j++)
            matrix[i][j] = 0;
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

int matrix_multip(double **A, double **B, double **C, int m, int n)
{
    /* matrix A has m rows and n columns
     * matrix B has n rows and m columns
     * matrix C is m*m
     */
    for (int i=0; i<m; i++)
        for(int j=0; j<m; j++)
            for(int k=0; k<n; k++)
            {
                C[i][j] = C[i][j] + A[i][k]*B[k][j];
            }
    
    return 0;
}

int print_matrix(int process_num, double **matrix, int *dim, int id, int *grid_coords, char *msg, MPI_Comm WORLD_COMM)
{
    MPI_Barrier(WORLD_COMM);
    for (int k=0; k<process_num; k++)
    {
        if (k == id)
        {
            printf("%s in process (%d %d) is:\n", msg, grid_coords[0], grid_coords[1]);
            for (int i=0; i<dim[0]; i++)
            {
                for (int j=0; j<dim[1]; j++)
                    printf("%f ", matrix[i][j]);
                printf("\n");
            }
        }
        MPI_Barrier(WORLD_COMM);
    }
    return 0;
}

int main(int argc, char *argv[])
{
    MPI_Comm cart_comm; //Cartesian topology communicator
    double elapsed_time; //parallel execution time 
    int p; // number of processes
    int id; // process ID number 
    int root; // the root process ID
    int grid_coords[2]; // process coordinates
    int periodic[2];
    int size[2]; // number of rows of blocks (==sqrt(p))

    
    double **A;
    int dimA[2]; //dimension of matrix A
    int r_dimA[2];
    double **B;
    int dimB[2]; //dimension of matrix B
    int r_dimB[2];
    double **subA; //submatrix of A
    double **subA_buf;
    int dimsubA[2];
    double **subB;
    double **subB_buf;
    int dimsubB[2];
    double **C;
    int dimC[2];
    int r_dimC[2];
    double **subC;
    int dimsubC[2];

    MPI_Init(&argc, &argv);
    /* start the timer */
    
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
    
        FILE* file = fopen (INPUT_FILE, "r");

        fscanf(file, "%d", &dimA[0]);
        fscanf(file, "%d", &dimA[1]);

        dimB[0] = dimA[1];
        dimB[1] = dimA[0];

#ifdef DEBUG
    printf("dimA: %d %d\n", dimA[0], dimA[1]);
    printf("dimB: %d %d\n", dimB[0], dimB[1]);
#endif

        int i,j;
        // Assign memory for A and zero-filling
        // Make sure that A can be divided equaly to each process
        malloc2ddouble(&A, dimA[0] + dimA[0]%size[0], dimA[1] + dimA[1]%size[1]);

        for (i=0; i<dimA[0]; i++)
            for (j=0; j<dimA[1]; j++)
            {
                fscanf(file, "%lf", &A[i][j]);
            }

        malloc2ddouble(&B, dimB[0] + dimB[0]%size[0], dimB[1] + dimB[1]%size[1]);
        for (i=0; i<dimB[0];i++)
            for (j=0; j<dimB[1]; j++)
            {
                fscanf(file, "%lf", &B[i][j]);
            }  
        fclose(file);
    }
    
#ifdef DEBUG
    /* Pring matrix A and matrix B*/
    if (id == root)
    {
        printf("the matrix A is:\n");
        for (int i=0; i<dimA[0]; i++)
        {
            for(int j=0; j<dimA[1]; j++)
                printf("%f ", A[i][j]);
            printf("\n");
        }
    }

    if (id == root)
    {
        printf("the matrix B is:\n");
        for (int i=0; i<dimB[0]; i++)
        {
            for(int j=0; j<dimB[1]; j++)
                printf("%f ", B[i][j]);
            printf("\n");
        }
    }
#endif
    /* Launch the timer */
    MPI_Barrier(MPI_COMM_WORLD); // block until all processes in the communicator have achieved this routine
    elapsed_time = - MPI_Wtime();
    
    /* Broadcast the dimA and create matrix C */
    MPI_Bcast(&dimA[0], 1, MPI_INT, root, cart_comm);
    MPI_Bcast(&dimA[1], 1, MPI_INT, root, cart_comm);
    MPI_Bcast(&dimB[0], 1, MPI_INT, root, cart_comm);
    MPI_Bcast(&dimB[1], 1, MPI_INT, root, cart_comm);
    r_dimA[0] = dimA[0] + dimA[0]%size[0];
    r_dimA[1] = dimA[1] + dimA[1]%size[1];
    r_dimB[0] = dimB[0] + dimB[0]%size[0];
    r_dimB[1] = dimB[1] + dimB[1]%size[1];

    dimC[0] = dimA[0];
    dimC[1] = dimB[1];
    r_dimC[0] = r_dimA[0];
    r_dimC[1] = r_dimB[1];
    malloc2ddouble(&C, r_dimC[0], r_dimC[1]);
    clear_matrix(C, r_dimC);
    /* Only one process */
    if (p == 1)
    {
        matrix_multip(A, B, C, dimA[0], dimA[1]);
        printf("The result: \n");
        for(int i=0; i<dimC[0]; i++)
        {
            for(int j=0; j<dimC[1]; j++)
                printf("%.1lf ", C[i][j]);
            printf("\n");
        }
        return 0;
    }

#if DEBUG
    printf("dimC: %d %d\n", dimC[0], dimC[1]);
#endif
    
    /* Create subMatrixes */
    MPI_Datatype subAtype, r_subAtype;
    MPI_Datatype subBtype, r_subBtype;
    dimsubA[0] = SUBMATRIX_HEIGHT(dimA[0], size[0]);
    dimsubA[1] = SUBMATRIX_HEIGHT(dimA[1], size[1]);
    dimsubB[0] = SUBMATRIX_HEIGHT(dimB[0], size[0]);
    dimsubB[1] = SUBMATRIX_WIDTH(dimB[1], size[1]);
    dimsubC[0] = dimsubA[0];
    dimsubC[1] = dimsubB[1];

#ifdef DEBUG
    printf("Dimension of subA in process %d: %d %d\n", id, dimsubA[0], dimsubA[1]);
    printf("Dimension of subB in process %d: %d %d\n", id, dimsubB[0], dimsubB[1]);
#endif

    malloc2ddouble(&subA, dimsubA[0], dimsubA[1]); 
    malloc2ddouble(&subA_buf, dimsubA[0], dimsubA[1]);
    int starts[2] = {0,0};
    MPI_Type_create_subarray(2, r_dimA, dimsubA, starts, MPI_ORDER_C, MPI_DOUBLE, &subAtype);
    MPI_Type_create_resized(subAtype, 0, dimsubA[1]*sizeof(double), &r_subAtype);
    MPI_Type_commit(&r_subAtype);

    malloc2ddouble(&subB, dimsubB[0], dimsubB[1]);
    malloc2ddouble(&subB_buf, dimsubB[0], dimsubB[1]);
    MPI_Type_create_subarray(2, r_dimB, dimsubB, starts, MPI_ORDER_C, MPI_DOUBLE, &subBtype);
    MPI_Type_create_resized(subBtype, 0, dimsubB[1]*sizeof(double), &r_subBtype);
    MPI_Type_commit(&r_subBtype);

    malloc2ddouble(&subC, dimsubC[0], dimsubC[1]);

    double *Aptr = NULL;
    double *Bptr = NULL;
    if (id == root) Aptr = &(A[0][0]);
    if (id == root) Bptr = &(B[0][0]);

    int counts[p];
    int displs[p];
    find_counts_displs(counts, displs, size, dimsubA);
    MPI_Scatterv(Aptr, counts, displs, r_subAtype, &(subA[0][0]), dimsubA[0]*dimsubA[1], MPI_DOUBLE, root, cart_comm);

    find_counts_displs(counts, displs, size, dimsubB);
    MPI_Scatterv(Bptr, counts, displs, r_subBtype, &(subB[0][0]), dimsubB[0]*dimsubB[1], MPI_DOUBLE, root, cart_comm);

#ifdef DEBUG
    /* Each process print received submatrix of A */
/*    MPI_Barrier(cart_comm);
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
*/
    print_matrix(p, subA, dimsubA, id, grid_coords, "The submatrix A", cart_comm);

    /* Each process print received submatrix of B */
    //MPI_Barrier(cart_comm);
    //for (int k=0; k<p; k++)
    //{
    //    if (k == id)
    //    {
    //        printf("The submatrix B in process (%d %d)is:\n", grid_coords[0], grid_coords[1]);
    //        for (int i=0; i<dimsubB[0]; i++)
    //        {
    //            for (int j=0; j<dimsubB[1]; j++)
    //                printf("%f ", subB[i][j]);
    //            printf("\n");
    //        }
    //    }
    //    MPI_Barrier(cart_comm);
    //}
    print_matrix(p, subB, dimsubB, id, grid_coords, "The submatrix B", cart_comm);
#endif
    
    /* Rearrange Blocks */
    int cntA_out[p];
    int cntA_in[p];
    int cntB_out[p];
    int cntB_in[p];
    int dispA_out[p];
    int dispA_in[p];
    int dispB_out[p];
    int dispB_in[p];

    int destA_coords[2];
    int destA_rank;
    int receA_coords[2];
    int receA_rank;
    int destB_coords[2];
    int destB_rank;
    int receB_coords[2];
    int receB_rank;
    MPI_Status status;
    for (int i=0; i<p; i++)
    {
        cntA_out[i] = 0;
        cntA_in[i] = 0;
        cntB_out[i] = 0;
        cntB_in[i] = 0;
        dispA_out[i] = 0;
        dispA_in[i] =dimsubA[0] * dimsubA[1];
        dispB_out[i] = 0;
        dispB_in[i] = dimsubB[0] * dimsubB[1];
    }
    destA_coords[0] = grid_coords[0];
    destA_coords[1] = (grid_coords[1] - grid_coords[0] + size[1])%size[1];
    MPI_Cart_rank(cart_comm, destA_coords, &destA_rank);
    cntA_out[destA_rank] = dimsubA[0] * dimsubA[1]; 
    receA_coords[0] = grid_coords[0];
    receA_coords[1] = (grid_coords[1] + grid_coords[0] + size[1])%size[1];
    MPI_Cart_rank(cart_comm, receA_coords, &receA_rank);
    cntA_in[receA_rank] = dimsubA[0] * dimsubA[1]; 
    dispA_in[receA_rank] = 0;
    MPI_Sendrecv_replace(&(subA[0][0]), dimsubA[0]*dimsubA[1], MPI_DOUBLE, destA_rank, 0, receA_rank, 0, cart_comm, &status);

#ifdef DEBUG
    printf("Process %d %d for A, dest_coords: %d %d, rece_coords: %d %d\n", 
            grid_coords[0],grid_coords[1], destA_coords[0], destA_coords[1],
            receA_coords[0], receA_coords[1]);
#endif

    destB_coords[0] = (grid_coords[0] - grid_coords[1] + size[0])%size[0];
    destB_coords[1] = grid_coords[1];
    MPI_Cart_rank(cart_comm, destB_coords, &destB_rank);
    cntB_out[destB_rank] = dimsubB[0] * dimsubB[1];
    receB_coords[0] = (grid_coords[0] + grid_coords[1] + size[0])%size[0];
    receB_coords[1] = grid_coords[1];
    MPI_Cart_rank(cart_comm, receB_coords, &receB_rank);
    cntB_in[receB_rank] = dimsubB[0] * dimsubB[1];
    dispB_in[receB_rank] = 0;
    MPI_Sendrecv_replace(&(subB[0][0]), dimsubB[0]*dimsubB[1], MPI_DOUBLE, destB_rank, 0, receB_rank, 0, cart_comm, &status);

#ifdef DEBUG
    MPI_Barrier(cart_comm);
    printf("Test point 1\n");
#endif

#ifdef DEBUG
    /* Each process print received submatrix of A */
    //MPI_Barrier(cart_comm);
    //for (int k=0; k<p; k++)
    //{
    //    if (k == id)
    //    {
    //        printf("Rearranged: The submatrix A in process (%d %d)is:\n", grid_coords[0], grid_coords[1]);
    //        for (int i=0; i<dimsubA[0]; i++)
    //        {
    //            for (int j=0; j<dimsubA[1]; j++)
    //                printf("%f ", subA[i][j]);
    //            printf("\n");
    //        }
    //    }
    //    MPI_Barrier(cart_comm);
    //}
    print_matrix(p, subA, dimsubA, id, grid_coords, "Rearranged: The submatrix A", cart_comm);

    ///* Each process print received submatrix of B */
    //MPI_Barrier(cart_comm);
    //for (int k=0; k<p; k++)
    //{
    //    if (k == id)
    //    {
    //        printf("Rearranged: The submatrix B in process (%d %d)is:\n", grid_coords[0], grid_coords[1]);
    //        for (int i=0; i<dimsubB[0]; i++)
    //        {
    //            for (int j=0; j<dimsubB[1]; j++)
    //                printf("%f ", subB[i][j]);
    //            printf("\n");
    //        }
    //    }
    //    MPI_Barrier(cart_comm);
    //}
    print_matrix(p, subB, dimsubB, id, grid_coords, "Rearranged: The submatrix B", cart_comm);
#endif

    matrix_multip(subA, subB, subC, dimsubA[0], dimsubA[1]);

    destA_coords[0] = grid_coords[0];
    destA_coords[1] = (grid_coords[1] - 1 + size[1])%size[1];
    MPI_Cart_rank(cart_comm, destA_coords, &destA_rank);
    receA_coords[0] = grid_coords[0];
    receA_coords[1] = (grid_coords[1] + 1 + size[1])%size[1];
    MPI_Cart_rank(cart_comm, receA_coords, &receA_rank);

    destB_coords[0] = (grid_coords[0] - 1 + size[0])%size[0];
    destB_coords[1] = grid_coords[1];
    MPI_Cart_rank(cart_comm, destB_coords, &destB_rank);
    receB_coords[0] = (grid_coords[0] + 1 + size[0])%size[0];
    receB_coords[1] = grid_coords[1];
    MPI_Cart_rank(cart_comm, receB_coords, &receB_rank);
    for (int i=0; i<size[0] - 1; i++)
    {
        MPI_Sendrecv_replace(&(subA[0][0]), dimsubA[0]*dimsubA[1], MPI_DOUBLE, destA_rank, 0, receA_rank, 0, cart_comm, &status);
        MPI_Sendrecv_replace(&(subB[0][0]), dimsubB[0]*dimsubB[1], MPI_DOUBLE, destB_rank, 0, receB_rank, 0, cart_comm, &status);
        matrix_multip(subA, subB, subC, dimsubA[0], dimsubA[1]);

 #ifdef DEBUG
    /* Each process print received submatrix of A */
    char msg[100];    
    sprintf(msg,"Loop %d: The submatrix A", i);
    print_matrix(p, subA, dimsubA, id, grid_coords, msg, cart_comm);
        //MPI_Barrier(cart_comm);
        //for (int k=0; k<p; k++)
        //{
        //    if (k == id)
        //    {
        //        printf("Loop %d:The submatrix A in process (%d %d)is:\n",i, grid_coords[0], grid_coords[1]);
        //        for (int i=0; i<dimsubA[0]; i++)
        //        {
        //            for (int j=0; j<dimsubA[1]; j++)
        //                printf("%f ", subA[i][j]);
        //            printf("\n");
        //        }
        //    }
        //    MPI_Barrier(cart_comm);
        //}

        /* Each process print received submatrix of B */
        sprintf(msg,"Loop %d: The submatrixB", i);
        print_matrix(p, subB, dimsubB, id, grid_coords, msg, cart_comm);
        //MPI_Barrier(cart_comm);
        //for (int k=0; k<p; k++)
        //{
        //    if (k == id)
        //    {
        //        printf("The submatrix B in process (%d %d)is:\n", grid_coords[0], grid_coords[1]);
        //        for (int i=0; i<dimsubB[0]; i++)
        //        {
        //            for (int j=0; j<dimsubB[1]; j++)
        //                printf("%f ", subB[i][j]);
        //            printf("\n");
        //        }
        //    }
        //    MPI_Barrier(cart_comm);
        //}

#endif
#ifdef DEBUG
        /* Print each sub matrix of C */
        sprintf(msg,"Loop %d: The submatrixC", i);
        print_matrix(p, subC, dimsubC, id, grid_coords, msg, cart_comm);
        //MPI_Barrier(cart_comm);
        //for (int k=0; k<p; k++)
        //{
        //    if (k == id)
        //    {
        //        printf("Loop %d:The submatrix C in process (%d %d)is:\n", i, grid_coords[0], grid_coords[1]);
        //        for (int i=0; i<dimsubC[0]; i++)
        //        {
        //            for (int j=0; j<dimsubC[1]; j++)
        //                printf("%f ", subC[i][j]);
        //            printf("\n");
        //        }
        //    }
        //    MPI_Barrier(cart_comm);
        //}
#endif
    }

    /* Gather matrix C */
    MPI_Datatype subCtype, r_subCtype;
#ifdef DEBUG
    MPI_Barrier(cart_comm);
    printf("dimsubC: %d %d\n", dimsubC[0], dimsubC[1]);
#endif
    MPI_Type_create_subarray(2, r_dimC, dimsubC, starts, MPI_ORDER_C, MPI_DOUBLE, &subCtype);
    MPI_Type_create_resized(subCtype, 0, dimsubC[1]*sizeof(double), &r_subCtype);
    MPI_Type_commit(&r_subCtype);
    find_counts_displs(counts, displs, size, dimsubC);

    MPI_Gatherv(&(subC[0][0]), dimsubC[0]*dimsubC[1], MPI_DOUBLE, &(C[0][0]), counts, displs, r_subCtype, root, cart_comm);
    
    /* stop timer */
    elapsed_time += MPI_Wtime();
    
    /* print result */
    if (id == root)
    {
        printf("The result: \n");
        for(int i=0; i<dimC[0]; i++)
        {
            for(int j=0; j<dimC[1]; j++)
                printf("%.1lf ", C[i][j]);
            printf("\n");
        }
        printf("Total elapsed time: %10.6f\n", elapsed_time);
    }

    MPI_Finalize();
    return 0;
}
