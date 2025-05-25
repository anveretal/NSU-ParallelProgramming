#include <iostream>
#include <mpi.h>

using namespace std;

int main() {
    int size;
    int rank;

    MPI_Init(NULL, NULL);
    
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // для 4-х процессов
    int recvcounts[4] = {2, 1, 2, 3};
    int displs[4] = {0, 0, 0, 0};
    for (int i = 1; i < 4; i++) {
        displs[i] = displs[i - 1] + recvcounts[i - 1];
    }

    int *sendbuf;
    int sendcount;
    for (int i = 0; i < 4; i++) {
        if (rank == i) {
            sendbuf = (int *) malloc(recvcounts[i] * sizeof(int));
            for (int j = 0; j < recvcounts[i]; j++) {
                sendbuf[j] = i;
            }
            sendcount = recvcounts[i];
        }
    }

    int *recvbuf = (int *) malloc(8 * sizeof(int));

    MPI_Allgatherv(sendbuf, sendcount, MPI_INT, recvbuf, recvcounts, displs, MPI_INT, MPI_COMM_WORLD);

    for (int i = 0; i < 4; i++) {
        MPI_Barrier(MPI_COMM_WORLD);
        if (rank == i) {
            for (int j = 0; j < 8; j++) {
                cout << recvbuf[j] << ' ';
            }
            cout << endl;
        }
    }

    MPI_Finalize();

    return 0;
}