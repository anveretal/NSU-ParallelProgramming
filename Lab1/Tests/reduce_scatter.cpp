#include <iostream>
#include <mpi.h>

int main() {
    int size;
    int rank;

    MPI_Init(NULL, NULL);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int local_sum = 1;
    int sum;

    int recvcounts[] = {2};

    //MPI_Reduce(local_sum, sum, 2, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce_scatter(&local_sum, &sum, recvcounts, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    for (int i = 0; i < size; i++) {
        MPI_Barrier(MPI_COMM_WORLD);
        if (rank == i) {
            printf("%d\n", sum);
        }
    }

    MPI_Finalize();

    return 0;
}