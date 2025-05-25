#include <iostream>
#include <mpi.h>

int main() {
    int size;
    int rank;

    MPI_Init(NULL, NULL);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int local_sum[2] = {1, 1};
    int *sum;

    if (rank == 0) {
        sum = (int *) malloc(2 * sizeof(int));
    }

    MPI_Reduce(local_sum, sum, 2, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    for (int i = 0; i < size; i++) {
        MPI_Barrier(MPI_COMM_WORLD);
        if (rank == i) {
            printf("%d %d\n", sum[0], sum[1]);
        }
    }

    MPI_Finalize();

    return 0;
}