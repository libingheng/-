#include <stdio.h>
#include <mpi.h>
#include <malloc.h>
#include <algorithm>
#include <iostream>
#define M 1000
#define N 1000

bool check(double *a, double *b, int n)
{
    for (int i = 0; i < n; i++)
    {
        if (a[i] != b[i])
            return false;
    }
    return true;
}
int main(int argc, char **argv)
// int main()
{
    // printf("initialize done1\n");

    for (int loop = 0; loop < 5; loop++)
    {
        int my_rank;
        int comm_sz;
        double local_M;
        int i, j, k;
        double start, finish;
        double tem;
        MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
        // printf("initialize done\n");

        //平均每个线程获得的matrixA的行数
        local_M = M / comm_sz;

        double *local_Matrix_A = (double *)malloc(local_M * N * sizeof(double));

        double *Matrix_A = NULL;
        double *Matrix_B = (double *)malloc(M * N * sizeof(double));

        double *local_result = (double *)malloc(local_M * N * sizeof(double));

        double *result_Matrix = NULL;
        if (my_rank == 0)
        {
            // printf("process %d of %d\n",my_rank,comm_sz);

            Matrix_A = (double *)malloc(M * N * sizeof(double));
            for (i = 0; i < M; i++)
            {
                for (j = 0; j < N; j++)
                {
                    Matrix_A[1 * N + j] = rand() / double(RAND_MAX);
                    Matrix_B[1 * N + j] = rand() / double(RAND_MAX);
                }
            }

            start = MPI_Wtime();

            MPI_Scatter(Matrix_A, local_M * N, MPI_DOUBLE, local_Matrix_A, local_M * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

            MPI_Bcast(Matrix_B, M * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

            for (i = 0; i < local_M; i++)
                for (j = 0; j < M; j++)
                {
                    tem = 0;
                    for (k = 0; k < N; k++)
                        tem += local_Matrix_A[i * M + k] * Matrix_B[j * M + k];
                    local_result[i * M + j] = tem;
                }
            free(local_Matrix_A);

            result_Matrix = (double *)malloc(M * N * sizeof(double));

            MPI_Gather(local_result, local_M * N, MPI_DOUBLE, result_Matrix, local_M * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            // the rest
            int rest = M % comm_sz;
            if (rest != 0)
                for (i = M - rest - 1; i < M; i++)
                    for (j = 0; j < M; j++)
                    {
                        tem = 0;
                        for (k = 0; k < N; k++)
                            tem += Matrix_A[i * M + k] * Matrix_B[j * M + k];
                        result_Matrix[i * M + j] = tem;
                    }
            MPI_Barrier(MPI_COMM_WORLD);
            finish = MPI_Wtime();

            printf("\nProcess %d used time = %lf ms\n", my_rank, double(finish - start) * 1000);

            start = MPI_Wtime();
            double *currect_answer = (double *)malloc(M * N * sizeof(double));
            for (i = 0; i < M; i++)
                for (j = 0; j < M; j++)
                {
                    tem = 0;
                    for (k = 0; k < N; k++)
                        tem += Matrix_A[i * M + k] * Matrix_B[j * M + k];
                    currect_answer[i * M + j] = tem;
                }
            finish = MPI_Wtime();

            printf("串行计算用时 = %lf ms\n", double(finish - start) * 1000);

            if (check(currect_answer, result_Matrix, M * N))
            {
                printf("correct\n");
            }

            free(Matrix_A);
            free(Matrix_B);
            free(local_result);
            free(currect_answer);
        }
        else
        {
            // printf("process %d of %d\n", my_rank, comm_sz);

            MPI_Scatter(Matrix_A, local_M * N, MPI_DOUBLE, local_Matrix_A, local_M * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Bcast(Matrix_B, M * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            for (i = 0; i < local_M; i++)
                for (j = 0; j < M; j++)
                {
                    tem = 0;
                    for (k = 0; k < N; k++)
                        tem += local_Matrix_A[i * M + k] * Matrix_B[j * M + k];
                    local_result[i * M + j] = tem;
                }
            free(local_Matrix_A);
            free(Matrix_B);

            MPI_Gather(local_result, local_M * N, MPI_DOUBLE, result_Matrix, local_M * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);

            free(local_result);
        }
        MPI_Finalize();
    }

    return 0;
}
