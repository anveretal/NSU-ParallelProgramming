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
        1, 7,  13, 19, 25, 31,
        2, 8,  14, 20, 26, 32,
        3, 9,  15, 21, 27, 33,
        4, 10, 16, 22, 28, 34,
        5, 11, 17, 23, 29, 35,
        6, 12, 18, 24, 30, 36
    };

    MPI_Datatype temp;
    MPI_Type_vector(3, 2, 6, MPI_INT, &temp);
    MPI_Type_commit(&temp);

    MPI_Datatype block;
    MPI_Type_create_resized(temp, 0, 2 * sizeof(int), &block);
    MPI_Type_commit(&block);
    MPI_Type_free(&temp);

    int *sub_block = (int *) malloc(6 * sizeof(int));

    int sendcounts[] = {1, 1, 1, 1, 1, 1};
    int displs[] = {0, 1, 2, 9, 10, 11}; // i * (K / subCols) * (M / subRows) + j

    MPI_Scatterv(matrix, sendcounts, displs, block, sub_block, 6, MPI_INT, 0, MPI_COMM_WORLD);
    //MPI_Scatter(matrix, 1, block, sub_block, 4, MPI_INT, 0, MPI_COMM_WORLD);

    for (int r = 0; r < size; r++) {
        MPI_Barrier(MPI_COMM_WORLD);
        if (rank == r) {
            for (int i = 0; i < 3; i++) {
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