#include <iostream>
#include <unistd.h>
#include <mpi.h>

int main() {
    int gsize;
    int myrank;
    
    MPI_Init(NULL, NULL);

    MPI_Comm_size(MPI_COMM_WORLD, &gsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    int buffer;
    for (int i = 0; i < gsize; i++) {
        if (myrank == i) {
            buffer = i;
        }
    }

    for (int i = 0; i < gsize; i++) {
        MPI_Barrier(MPI_COMM_WORLD);
        if (myrank == i) {
            std::cout << myrank << ' ' << buffer << std::endl;
        }
    }
    
    usleep(100000);

    MPI_Sendrecv_replace(&buffer, 1, MPI_INT, (myrank + 1) % gsize, 12345, (myrank - 1) % gsize, 12345, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    for (int i = 0; i < gsize; i++) {
        MPI_Barrier(MPI_COMM_WORLD);
        if (myrank == i) {
            std::cout << myrank << ' ' << buffer << std::endl;
        }
    }

    MPI_Finalize();
    return 0;
}