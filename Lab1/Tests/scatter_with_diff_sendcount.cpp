#include <iostream>
#include <mpi.h>

int main() {
    int size;
    int rank;

    MPI_Init(NULL, NULL);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int *sendbuf;
    int  recvbuf[4];
    int  counts[4] = {2, 1, 2, 3}; // подразумевается 4 процесса
    int  count = 1;

    if (rank == 0) {
        sendbuf = (int *) malloc(8 * sizeof(int));
        for (int i = 0; i < 8; i++) {
            sendbuf[i] = i;
        }
    }

    for (int i = 0; i < 4; i++) {
        if (rank == i) {
            count = counts[i];
            // break;
        }
    }

    // Вызовется ошибка, т.к.:
    // This message is split into n equal segments,
    // the ith segment is sent to the ith process in the group,
    // and each process receives this message as above.
    MPI_Scatter(sendbuf, count, MPI_INT, recvbuf, count, MPI_INT, 0, MPI_COMM_WORLD);

    for (int i = 0; i < 4; i++) {
        MPI_Barrier(MPI_COMM_WORLD);
        if (rank == i) {
            for (int j = 0; j < 4; j++) {
                printf("%d ", recvbuf[j]);
            }
            putchar('\n');
        }
    }

    MPI_Finalize();

    return 0;
}