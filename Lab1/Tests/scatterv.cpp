#include <iostream>
#include <stdlib.h>
#include <mpi.h>

int main() {
    int root = 0;
    int size;
    int rank;

    MPI_Init(NULL, NULL);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int *sendbuf;
    int *recvbuf = (int *) calloc(10, sizeof(int));

    // Подразумевается 4 процесса
    int *sendcounts;
    int *displs;
    int recvcount; // он не должен быть общий для всех, он должен быть уникальный для каждого процесса
    if (rank == root) {
        sendbuf = (int *) malloc(10 * sizeof(int));
        for (int i = 0; i < 10; i++) {
            sendbuf[i] = i;
        }
        sendcounts = (int *) malloc(4 * sizeof(int));
        sendcounts[0] = 2;
        sendcounts[1] = 2;
        sendcounts[2] = 1;
        sendcounts[3] = 3;
        displs = (int *) malloc(4 * sizeof(int));
        displs[0] = 0;
        displs[1] = 0;
        displs[2] = 0;
        displs[3] = 0;
    }
    MPI_Scatter(sendcounts, 1, MPI_INT, &recvcount, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatterv(sendbuf, sendcounts, displs, MPI_INT, recvbuf, recvcount, MPI_INT, root, MPI_COMM_WORLD);

    for (int i = 0; i < size; i++) {
        MPI_Barrier(MPI_COMM_WORLD);
        if (rank == i) {
            for (int i = 0; i < 10; i++) {
                printf("%d ", recvbuf[i]);
            }
            putchar('\n');
        }
    }

    MPI_Finalize();

    return 0;
}