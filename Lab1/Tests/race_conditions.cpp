#include <iostream>
#include <mpi.h>

int main() {
    int size;
    int rank;

    MPI_Init(NULL, NULL);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    for (long long int i = 0; i < 10000; i++) {
        std::cout << rank;
    }
    std::cout << std::endl;

    MPI_Finalize();

    return 0;
}