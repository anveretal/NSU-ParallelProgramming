#include <iostream>
#include <unistd.h>
#include <mpi.h>

using namespace std;

int main() {
    int size;
    int rank;

    MPI_Init(NULL, NULL);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int matrix[] = { // 1, 5, 9, 13,   2, 6, 10, 14,   3, 7, 11, 15,   4, 8, 12, 16
        1, 5, 9,  13,
        2, 6, 10, 14,
        3, 7, 11, 15,
        4, 8, 12, 16
    };

    MPI_Datatype temp;
    // MPI_Type_vector(4, 2, 4, MPI_INT, &temp);
    // MPI_Type_commit(&temp);

    // MPI_Datatype two_cols;
    // MPI_Type_create_resized(temp, 0, 2 * sizeof(int), &two_cols);
    // MPI_Type_commit(&two_cols);
    // MPI_Type_free(&temp);

    // int *sub = (int *) malloc(8 * sizeof(int));

    // MPI_Scatter(matrix, 1, two_cols, sub, 8, MPI_INT, 0, MPI_COMM_WORLD);

    // for (int r = 0; r < size; r++) {
    //     MPI_Barrier(MPI_COMM_WORLD);
    //     if (rank == r) {
    //         for (int i = 0; i < 4; i++) {
    //             for (int j = 0; j < 2; j++) {
    //                 printf("%d ", sub[i * 2 + j]);
    //             }
    //             putchar('\n');
    //         }
    //         putchar('\n');
    //     }
    // }
    // sleep(1);

    // if (rank == 0) {
    //     printf("\n--------------------\n\n");
    // }
    // sleep(1);

    MPI_Type_vector(2, 2, 4, MPI_INT, &temp);
    MPI_Type_commit(&temp);

    MPI_Datatype block;
    MPI_Type_create_resized(temp, 0, 2 * sizeof(int), &block);
    MPI_Type_commit(&block);
    MPI_Type_free(&temp);

    int *sub_block = (int *) malloc(4 * sizeof(int));

    int sendcounts[] = {1, 1, 1, 1};
    int displs[] = {0, 1, 4, 5};

    MPI_Scatterv(matrix, sendcounts, displs, block, sub_block, 4, MPI_INT, 0, MPI_COMM_WORLD);
    //MPI_Scatter(matrix, 1, block, sub_block, 4, MPI_INT, 0, MPI_COMM_WORLD);

    for (int r = 0; r < size; r++) {
        MPI_Barrier(MPI_COMM_WORLD);
        if (rank == r) {
            for (int i = 0; i < 2; i++) {
                for (int j = 0; j < 2; j++) {
                    printf("%d ", sub_block[i * 2 + j]);
                }
                putchar('\n');
            }
            putchar('\n');
        }
    }

    MPI_Finalize();

    return 0;
}