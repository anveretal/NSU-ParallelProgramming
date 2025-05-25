#include <iostream>
#include <mpi.h>

int main() {
    int rank;

    MPI_Init(NULL, NULL);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        char buf[6] = "Hello";
        MPI_Send(buf, 6, MPI_CHAR, 1, 1234, MPI_COMM_WORLD);
    }
    else if (rank == 1) {
        char buf[6];
        printf("Before sending: %s\n", buf);
        MPI_Recv(buf, 6, MPI_CHAR, 0, 1234, MPI_COMM_WORLD, NULL);
        printf("After  sending: %s\n", buf);
    }

    MPI_Finalize();

    return 0;
}