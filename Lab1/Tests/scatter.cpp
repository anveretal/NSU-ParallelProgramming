#include <iostream>
#include <stdlib.h>
#include <mpi.h>

int main() {
    int rank;
    int size_received_values = 2;
    int *received_values = (int *) calloc(size_received_values, sizeof(int));
    int *buf;

    MPI_Init(NULL, NULL);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int size_buf = 10;
    if (rank == 0) {
        buf = (int *) malloc(size_buf * sizeof(int));
        for (int i = 1; i <= size_buf; i++) {
            buf[i - 1] = i;
        }
    }
    MPI_Scatter(buf, 2, MPI_INT, received_values, 2, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        for (int i = 0; i < size_buf; i++) {
            printf("%d ", buf[i]);
        }
    }

    printf("\n\nReceived value in %d: %d, %d\n", rank, received_values[0], received_values[1]);

    MPI_Finalize();

    return 0;
}