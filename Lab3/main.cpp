#include <iostream>
#include <mpi.h>
#include <vector>

constexpr int M = 4;
constexpr int N = 4;
constexpr int K = 4;

void print_matrix(std::vector<double> &matrix, int rows, int cols) {
    for (int r = 0; r < rows; r++) {
        for (int c = 0; c < cols; c++) {
            printf("%lf ", matrix[r * cols + c]);
        }
        putchar('\n');
    }
}

void create_matrix(std::vector<double>& matrix, int rows, int cols, double val) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            matrix[i * cols + j] = val;
        }
    }
}

void mul(double* A, double* B, double* C, int* matrix_sizes, int* grid_sizes, MPI_Comm comm) {
    MPI_Bcast(matrix_sizes, 3, MPI_INT, 0, comm);
    MPI_Bcast(grid_sizes, 2, MPI_INT, 0, comm);

    MPI_Comm comm2d;
    int periods[2] = {0, 0};
    MPI_Cart_create(comm, 2, grid_sizes, periods, 0, &comm2d);

    int size;
    int rank;
    MPI_Comm_size(comm2d, &size);
    MPI_Comm_rank(comm2d, &rank);
    int coords[2];
    MPI_Cart_coords(comm2d, rank, 2, coords);

    MPI_Comm comm1d[2];
    for (int i = 0; i < 2; ++i) {
        int remainDims[2] = {i == 0, i == 1};
        MPI_Cart_sub(comm2d, remainDims, &comm1d[i]);
    }

    int subRows = matrix_sizes[0] / grid_sizes[0];
    int subCols = matrix_sizes[2] / grid_sizes[1];

    MPI_Datatype typeB, typeC;
    int* sendcountsB = nullptr;
    int* displsScatterB = nullptr;
    int* recvcountsC = nullptr;
    int* displsGatherC = nullptr;

    if (rank == 0) {
        MPI_Datatype tempVec;
        MPI_Type_vector(matrix_sizes[1], subCols, matrix_sizes[2], MPI_DOUBLE, &tempVec);
        MPI_Type_commit(&tempVec);

        MPI_Type_create_resized(tempVec, 0, subCols * sizeof(double), &typeB);
        MPI_Type_commit(&typeB);
        MPI_Type_free(&tempVec);

        sendcountsB = new int[grid_sizes[1]]();
        displsScatterB = new int[grid_sizes[1]]();
        for (int i = 0; i < grid_sizes[1]; ++i) {
            sendcountsB[i] = 1;
            displsScatterB[i] = i;
        }

        MPI_Type_vector(subRows, subCols, matrix_sizes[2], MPI_DOUBLE, &tempVec);
        MPI_Type_commit(&tempVec);

        MPI_Type_create_resized(tempVec, 0, subCols * sizeof(double), &typeC);
        MPI_Type_commit(&typeC);
        MPI_Type_free(&tempVec);

        recvcountsC = new int[grid_sizes[0] * grid_sizes[1]]();
        displsGatherC = new int[grid_sizes[0] * grid_sizes[1]]();
        for (int i = 0; i < grid_sizes[0]; ++i) {
            for (int j = 0; j < grid_sizes[1]; ++j) {
                recvcountsC[i * grid_sizes[1] + j] = 1;
                displsGatherC[i * grid_sizes[1] + j] = i * subRows * (matrix_sizes[2] / subCols) + j;
                printf("%d ", displsGatherC[i * grid_sizes[1] + j]);
            }
        }
        printf("\n");
    }

    std::vector<double> subA(subRows * matrix_sizes[1]);
    std::vector<double> subB(matrix_sizes[1] * subCols);
    std::vector<double> subC(subRows * subCols, 0.0);

    // 1. Раздача матрицы A вдоль оси, направленной вниз.
    if (coords[1] == 0) {
        MPI_Scatter(A, subRows * matrix_sizes[1], MPI_DOUBLE, subA.data(), subRows * matrix_sizes[1], MPI_DOUBLE, 0, comm1d[0]);
    }
    // 2. Раздача матрицы B вдоль оси, направленной вправо.
    if (coords[0] == 0) {
        MPI_Scatterv(B, sendcountsB, displsScatterB, typeB, subB.data(), matrix_sizes[1] * subCols, MPI_DOUBLE, 0, comm1d[1]);
    }
    
    // 3. & 4. Трансляция для каждого процесса в сетке 2D.
    MPI_Bcast(subA.data(), subRows * matrix_sizes[1], MPI_DOUBLE, 0, comm1d[1]);
    MPI_Bcast(subB.data(), matrix_sizes[1] * subCols, MPI_DOUBLE, 0, comm1d[0]);

    // 5. Вычисление матрицы C по кусочкам
    for (int i = 0; i < subRows; ++i) {
        for (int k = 0; k < matrix_sizes[1]; ++k) {
            for (int j = 0; j < subCols; ++j) {
                subC[i * subCols + j] += subA[i * matrix_sizes[1] + k] * subB[k * subCols + j];
            }
        }
    }

    for (int r = 0; r < size; r++) {
        MPI_Barrier(comm2d);
        if (rank == r) {
            print_matrix(subC, subRows, subCols);
            putchar('\n');
        }
    }

    // 6. Сборка результирующей матрицы.
    MPI_Gatherv(subC.data(), subRows * subCols, MPI_DOUBLE, C, recvcountsC, displsGatherC, typeC, 0, comm2d);

    if (rank == 0) {
        delete[] sendcountsB;
        delete[] displsScatterB;
        delete[] recvcountsC;
        delete[] displsGatherC;
        MPI_Type_free(&typeB);
        MPI_Type_free(&typeC);
    }

    MPI_Comm_free(&comm2d);
    MPI_Comm_free(&comm1d[0]);
    MPI_Comm_free(&comm1d[1]);
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (argc != 3) {
        if (rank == 0) std::cerr << "Usage: mpirun -np <p1*p2> ./program <p1> <p2>\n";
        MPI_Finalize();
        return 1;
    }

    int p1 = std::atoi(argv[1]);
    int p2 = std::atoi(argv[2]);

    if (p1 <= 0 || p2 <= 0 || p1 * p2 != size) {
        if (rank == 0) std::cerr << "Invalid grid dimensions\n";
        MPI_Finalize();
        return 1;
    }

    int matrix_sizes[3] = {M, N, K};
    int grid_sizes[2] = {p1, p2};

    std::vector<double> A, B, C;
    if (rank == 0) {
        A.resize(M * N);
        B.resize(N * K);
        C.resize(M * K);
        create_matrix(A, M, N, 1.0);
        create_matrix(B, N, K, 1.0);
    }

    double start = MPI_Wtime();
    mul(A.data(), B.data(), C.data(), matrix_sizes, grid_sizes, MPI_COMM_WORLD);
    double end = MPI_Wtime();

    if (rank == 0) {
        std::cout << "Time taken: " << (end - start) << " seconds\n";
        print_matrix(C, M, K);
    }

    MPI_Finalize();
    return 0;
}