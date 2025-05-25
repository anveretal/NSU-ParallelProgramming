#include <iostream>
#include <stdlib.h>
#include <mpi.h>

int main() {
    int size;
    int rank;

    MPI_Init(NULL, NULL);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int sendbuf;
    int *recvbuf;

    if (rank == 0) {
        recvbuf = (int *) malloc(size * sizeof(int));
        sendbuf = 0;
    }
    for (int i = 1; i < size; i++) {
        if (rank == i) {
            sendbuf = i;
        }
    }

    MPI_Gather(&sendbuf, 1, MPI_INT, recvbuf, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        for (int i = 0; i < size; i++) {
            printf("%d ", recvbuf[i]);
        }
    }

    MPI_Finalize();

    return 0;
}