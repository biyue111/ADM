#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
//#define DEBUG 
#define MAX(A, B) ((A) > (B) ? (A) : (B))

int find_counts_displs(int *counts, int *displs, int size, int *dimsub)
{
    int k = 0;
    for (int i=0; i<size; i++)
    {
       
        counts[k] = 1;
        displs[k] = i * dimsub[0];
        k++;
       
    }
}  

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

int main(int argc, char *argv[])
{
    int size; // number of process
    int id; // process ID number
    int root; //root process ID
    double elapsed_time; // parallel execution time

    double **A; //matrixA
    double **B; //matrixB
    double **r_C; //matrixC with the addition of 0
    double **C; //matrixC
    double **subA;
    double **subB;
    double **subC;
    int dimA[2]; //dim of A
    int dimB[2]; //dim of B
    int dimC[2];
    int dimsubA[2];
    int dimsubB[2];
    int dimsubC[2];
    int r_dimA[2]; //real dim of A (after addition of 0)
    int r_dimB[2]; //real dim of B (after addition of 0)
    int r_dimC[2]; //for gather subC



    double *row_A; //
    double *col_B;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    root = 0;

    /* Read two matrices from the files */
    if(id == root){
        printf("\n-----Multiply two matrices with row in parallel-----\n\n");
        if(argc != 3){
            printf("Command line: %s has not enough matrices!\n", argv[0]);
            MPI_Finalize();
            exit(1);
        }


        FILE* file1 = fopen (argv[1], "r");
        fscanf(file1, "%d", &dimA[0]);
        fscanf(file1, "%d", &dimA[1]);

        FILE* file2 = fopen (argv[2], "r");
        fscanf(file2, "%d", &dimB[0]);
        fscanf(file2, "%d", &dimB[1]);

        if(dimA[1] != dimB[0]){
            printf("Command line: Two matrices can not be multiplied!\n");
            MPI_Finalize();
            exit(1);
        }

        dimC[0] = dimA[0];
        dimC[1] = dimB[1];

        r_dimA[1] = dimA[1];
        r_dimB[1] = dimB[1];

        if(size >= dimA[0])
            r_dimA[0] = size;
        else{
            int temp = ceil((double)dimA[0] / size);
            r_dimA[0] = size * temp ;
        }

        if(size >= dimB[0])
            r_dimB[0] = size;
        else{
            int temp = ceil((double)dimB[0] / size);
            r_dimB[0] = size * temp;
        }

        malloc2ddouble(&A, r_dimA[0], r_dimA[1]);

        for(int i = 0; i < dimA[0]; i++){
            for(int j = 0; j < r_dimA[1]; j++){
                fscanf(file1, "%lf", &A[i][j]);
            }
        }
        for(int i = dimA[0]; i < r_dimA[0]; i++){ // Add 0 to ensure the matrix can be divided into process size
            for(int j = 0; j < r_dimA[1]; j++){
                A[i][j] = 0;
            }
        }

        malloc2ddouble(&B, r_dimB[0], r_dimB[1]);
        for(int i = 0; i < dimB[0]; i++){
            for(int j = 0; j < r_dimB[1]; j++){
                fscanf(file2, "%lf", &B[i][j]);
            }
        }
        for(int i = dimB[0]; i < r_dimB[0]; i++){
            for(int j = 0; j < r_dimB[1]; j++){
                B[i][j] = 0;
            }
        }

        fclose(file1);
        fclose(file2);

        printf("Matrix A with row: %d, col: %d\n", dimA[0], dimA[1]);
        //for(int i = 0; i < dimA[0]; i++){
        //    for(int j = 0; j < dimA[1]; j++){
        //        printf("%f ", A[i][j]);
        //    }
        //    printf("\n");
        //}
        printf("\n");

        printf("Matrix B with row: %d, col: %d\n", dimB[0], dimB[1]);
        //for(int i = 0; i < dimB[0]; i++){
        //    for(int j = 0; j < dimB[1]; j++){
        //        printf("%f ", B[i][j]);
        //    }
        //    printf("\n");
        //}
        printf("\n");

#ifdef DEBUG
    printf("process num %d\n", size);
    printf("A %d %d\n", r_dimA[0], r_dimA[1]);
    for(int i = 0; i < r_dimA[0]; i++){
        for(int j = 0; j < r_dimA[1]; j++){
            printf("%f ", A[i][j]);
        }
        printf("\n");
    }
    printf("B %d %d\n", r_dimB[0], r_dimB[1]);
    for(int i = 0; i < r_dimB[0]; i++){
        for(int j = 0; j < r_dimB[1]; j++){
            printf("%f ", B[i][j]);
        }
        printf("\n");
    }
#endif

        malloc2ddouble(&C, dimC[0], dimC[1]);
    }

    /* Begin to record the time */
    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time = - MPI_Wtime(); 


    MPI_Bcast(&r_dimA[0], 1, MPI_INT, root, MPI_COMM_WORLD);
    MPI_Bcast(&r_dimA[1], 1, MPI_INT, root, MPI_COMM_WORLD);
    MPI_Bcast(&r_dimB[0], 1, MPI_INT, root, MPI_COMM_WORLD);
    MPI_Bcast(&r_dimB[1], 1, MPI_INT, root, MPI_COMM_WORLD);

    dimsubA[0] = r_dimA[0] / size;
    dimsubA[1] = r_dimA[1];
    dimsubB[0] = r_dimB[0] / size;
    dimsubB[1] = r_dimB[1];
    dimsubC[0] = dimsubA[0];
    dimsubC[1] = dimsubB[1];
    r_dimC[0] = dimsubC[0] * size;
    r_dimC[1] = dimsubC[1];
    malloc2ddouble(&r_C, r_dimC[0], r_dimC[1]);

#ifdef DEBUG
    printf("subdimA: %d %d\n", dimsubA[0], dimsubA[1]);
#endif
    /* Create Sub Matrix*/
    int starts[2] = {0, 0};
    
    malloc2ddouble(&subA, dimsubA[0], dimsubA[1]);

    MPI_Datatype subAtype, r_subAtype;
    MPI_Type_create_subarray(2, r_dimA, dimsubA, starts, MPI_ORDER_C, MPI_DOUBLE, &subAtype);
    MPI_Type_create_resized(subAtype, 0, dimsubA[1]*sizeof(double), &r_subAtype);
    MPI_Type_commit(&r_subAtype);

    malloc2ddouble(&subB, dimsubB[0], dimsubB[1]);
    MPI_Datatype subBtype, r_subBtype;
    MPI_Type_create_subarray(2, r_dimB, dimsubB, starts, MPI_ORDER_C, MPI_DOUBLE, &subBtype);
    MPI_Type_create_resized(subBtype, 0, dimsubB[1]*sizeof(double), &r_subBtype);
    MPI_Type_commit(&r_subBtype);

    malloc2ddouble(&subC, dimsubC[0], dimsubC[1]);
    for(int i = 0; i < dimsubC[0]; i++){
        for(int j = 0; j < dimsubC[1]; j++){
            subC[i][j] = 0;
        }
    }

    /* Scatter the matrix by several rows to each process */
    int sendcountA[size];
    int senddispA[size];

    double *Aptr = NULL;
    if(id == root) Aptr = &(A[0][0]);
    find_counts_displs(sendcountA, senddispA, size, dimsubA);
    MPI_Scatterv(Aptr, sendcountA, senddispA, r_subAtype, &(subA[0][0]), dimsubA[0]*dimsubA[1], MPI_DOUBLE, root, MPI_COMM_WORLD);

    int sendcountB[size];
    int senddispB[size];
    double *Bptr = NULL;
    if(id == root) Bptr = &(B[0][0]);
    find_counts_displs(sendcountB, senddispB, size, dimsubB);
    MPI_Scatterv(Bptr, sendcountB, senddispB, r_subBtype, &(subB[0][0]), dimsubB[0]*dimsubB[1], MPI_DOUBLE, root, MPI_COMM_WORLD);
    
#ifdef DEBUG
    MPI_Barrier(MPI_COMM_WORLD);
    
    for (int k=0; k<size; k++){
        if (k == id){
            printf("process num %d subA %d %d\n", id, dimsubA[0], dimsubA[1]);
            for(int i = 0; i < dimsubA[0]; i++){
                for(int j = 0; j < dimsubA[1]; j++){
                    printf("%f ",subA[i][j]);
                }
                printf("\n");
            }
            
            printf("process num %d subB %d %d\n", id, dimsubA[0], dimsubA[1]);
            for(int i = 0; i < dimsubB[0]; i++){
                for(int j = 0; j < dimsubB[1]; j++){
                    printf("%f ",subB[i][j]);
                }
                printf("\n");
            }
            
        }
        MPI_Barrier(MPI_COMM_WORLD);
        printf("\n");
    }
#endif


    for(int p = 0; p < size; p++){

    	/* Do the multiplication of several rows in two matrices */
        for(int i = 0; i < dimsubA[0]; i++){
            for(int m = 0; m < dimsubB[0]; m++){
                for(int n = 0; n < dimsubC[1]; n++){
                    subC[i][n] += subA[i][m +  ((int)(id + p) % size) * dimsubB[0]] * subB[m][n];
                }
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);

#ifdef DEBUG
    MPI_Barrier(MPI_COMM_WORLD);
    
    for (int k=0; k<size; k++){
        if (k == id){
            printf("p: %d process num %d subA %d %d\n",p , id, dimsubA[0], dimsubA[1]);
            for(int i = 0; i < dimsubA[0]; i++){
                for(int j = 0; j < dimsubA[1]; j++){
                    printf("%f ",subA[i][j]);
                }
                printf("\n");
            }
            printf("p: %d process num %d subB %d %d\n",p , id, dimsubB[0], dimsubB[1]);
            for(int i = 0; i < dimsubB[0]; i++){
                for(int j = 0; j < dimsubB[1]; j++){
                    printf("%f ",subB[i][j]);
                }
                printf("\n");
            }
            printf("p: %d process num %d subC %d %d\n",p , id, dimsubC[0], dimsubC[1]);
            for(int i = 0; i < dimsubC[0]; i++){
                for(int j = 0; j < dimsubC[1]; j++){
                    printf("%f ",subC[i][j]);
                }
                printf("\n");
            }
            printf("\n");
            
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
#endif

    	/* Do the communication between processes and send the task as a ring */
        int dest_id, recv_id;
        MPI_Status status;
        if(id == (size - 1))
            dest_id = root;
        else
            dest_id = id + 1;

        if(id == root)
            recv_id = (size - 1);
        else
            recv_id = id - 1;

        MPI_Sendrecv_replace(&(subB[0][0]), dimsubB[0] * dimsubB[1], MPI_DOUBLE, recv_id, 0, dest_id, 0, MPI_COMM_WORLD, &status);
        MPI_Barrier(MPI_COMM_WORLD);
    }
    
#ifdef DEBUG
    MPI_Barrier(MPI_COMM_WORLD);
    
    for (int k=0; k<size; k++){
        if (k == id){
            printf("process num %d subC %d %d\n", id, dimsubC[0], dimsubC[1]);
            for(int i = 0; i < dimsubC[0]; i++){
                for(int j = 0; j < dimsubC[1]; j++){
                    printf("%f ",subC[i][j]);
                }
                printf("\n");
            }
            
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
#endif

    /* Gather all sub-matrices to a whole matrix */
    MPI_Datatype subCtype, r_subCtype;
    MPI_Type_create_subarray(2, r_dimC, dimsubC, starts, MPI_ORDER_C, MPI_DOUBLE, &subCtype);
    MPI_Type_create_resized(subCtype, 0, dimsubC[1]*sizeof(double), &r_subCtype);
    MPI_Type_commit(&r_subCtype);

    int sendcountC[size];
    int senddispC[size];
    find_counts_displs(sendcountC, senddispC, size, dimsubC);

    MPI_Gatherv(&(subC[0][0]), dimsubC[0]*dimsubC[1], MPI_DOUBLE, &(r_C[0][0]),sendcountC, senddispC, r_subCtype, root, MPI_COMM_WORLD);
    

#ifdef DEBUG
    if(id == root){
        
        printf("process num %d\n", size);
        printf("r_C %d %d\n", r_dimC[0], r_dimC[1]);
        for(int i = 0; i < r_dimC[0]; i++){
            for(int j = 0; j < r_dimC[1]; j++){
              printf("%f ", r_C[i][j]);
            }
            printf("\n");
        }
    }
#endif

    /* Stop time */
    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time += MPI_Wtime();
    

    if(id == root){
        
        for(int i = 0; i < dimC[0]; i++){
            for(int j = 0; j < dimC[1]; j++){
                C[i][j] = r_C[i][j];
            }
        }

        printf("Success to multiply matrix A and matrix B\n");

        /*
        printf("Result for A * B :\n");
        for(int i = 0; i < dimC[0]; i++){
            for(int j = 0; j < dimC[1]; j++){
                printf("%f ", C[i][j]);
            }
            printf("\n");
        }
        */
        
    }


    if(id == root) printf("\nTotal elapsed time: %10.6f\n", elapsed_time);
    MPI_Finalize();
    return 0;
}
