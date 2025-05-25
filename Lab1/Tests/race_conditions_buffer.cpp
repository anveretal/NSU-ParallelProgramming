#include <iostream>
#include <mpi.h>

int main() {
    int size;
    int rank;

    MPI_Init(NULL, NULL);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int *buf = (int *) malloc(size * sizeof(int));
    for (int i = 0; i < size; i++) {
        buf[i] = i;
    }

    printf("Rank: %d\n", rank);
    for (int i = 0; i < size; i++) {
        printf("%d ", buf[i]);
    }
    printf("\n\n");

    MPI_Finalize();

    return 0;
}