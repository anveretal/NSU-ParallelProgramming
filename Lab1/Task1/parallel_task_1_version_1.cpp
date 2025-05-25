#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <chrono>
#include <unistd.h>
#include <mpi.h>

#define EPS 0.00001

using namespace std;



typedef struct double_matrix_part {
    int     rows_count;
    int     cols_count;
    double *rows;
} Double_matrix_part;



double *create_double_vector(int size, double value) {
    double *v;
    if (value == 0) {
        v = (double *) calloc(size, sizeof(double));
        return v;
    }
    v = (double *) malloc(size * sizeof(double));
    for (int i = 0; i < size; i++) {
        v[i] = value;
    }
    return v;
}

double *create_double_matrix(int size, double value, double main_diagonal_value) {
    double *m = (double *) malloc(size * size * sizeof(double));
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (i == j) {
                m[i * size + j] = main_diagonal_value;
            }
            else {
                m[i * size + j] = value;
            }
        }
    }
    return m;
}

void mul_vector_part_by_scalar(double *vec_part, int vec_part_size, const double scalar) {
    for (int i = 0; i < vec_part_size; i++) {
        vec_part[i] *= scalar;
    }
}

void sub_vector_from_vector_part(double *vec_part, int size_vec_part, const double *vec) {
    for (int i = 0; i < size_vec_part; i++) {
        vec_part[i] -= vec[i];
    }
}

void mul_matrix_part_by_vector(double *res, const Double_matrix_part &matrix_part, const double *vec) {
    for (int i = 0; i < matrix_part.rows_count; i++) {
        for (int j = 0; j < matrix_part.cols_count; j++) {
            res[i] += matrix_part.rows[i * matrix_part.cols_count + j] * vec[j];
        }
    }
}

// x_(n+1) = x_n - TAU * (Ax_n - b)
void calc_next_x(const Double_matrix_part &A_part, const double *b, double *x, const double TAU, const int myrank, const int gsize) {
    double *Ax_part = create_double_vector(A_part.rows_count, 0);
    mul_matrix_part_by_vector(Ax_part, A_part, x);

    double *Ax_sub_b_part = Ax_part;
    sub_vector_from_vector_part(Ax_sub_b_part, A_part.rows_count, b);

    double *tau_mul__Ax_sub_b__part = Ax_sub_b_part;
    mul_vector_part_by_scalar(tau_mul__Ax_sub_b__part, A_part.rows_count, TAU);

    int *recvcounts = (int *) malloc(gsize * sizeof(int));
    MPI_Allgather(
        &A_part.rows_count, 1, MPI_INT, recvcounts, 1, MPI_INT, MPI_COMM_WORLD
    );

    int *displs = (int *) calloc(gsize, sizeof(int));
    for (int i = 1; i < gsize; i++) {
        displs[i] = displs[i - 1] + recvcounts[i - 1];
    }

    double *tau_mul__Ax_sub_b__ = create_double_vector(A_part.cols_count, 0);
    MPI_Allgatherv(
        tau_mul__Ax_sub_b__part, A_part.rows_count, MPI_DOUBLE, tau_mul__Ax_sub_b__, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD
    );

    // for (int i = 0; i < gsize; i++) {
    //     if (myrank == i) {
    //         for (int j = 0; j < A_part.cols_count; j++) {
    //             printf("%lf ", tau_mul__Ax_sub_b__[j]);
    //         }
    //         printf("\n\n");
    //     }
    //     MPI_Barrier(MPI_COMM_WORLD);
    // }
    // usleep(100000);

    for (int i = 0; i < A_part.cols_count; i++) {
        x[i] -= tau_mul__Ax_sub_b__[i];
    }

    free(tau_mul__Ax_sub_b__);
    free(displs);
    free(recvcounts);
    free(tau_mul__Ax_sub_b__part);
}

double calc_norm_vector(const double *vec, const int size) {
    double sum = 0;
    for (int i = 0; i < size; i++) {
        sum += vec[i] * vec[i];
    }
    return sqrt(sum);
}

// sqrt(sum_i(vec[i]^2))
double calc_norm_vector_part(const double *vec_part, const int size_vec_part) {
    double local_sum = 0;
    for (int i = 0; i < size_vec_part; i++) {
        local_sum += vec_part[i] * vec_part[i];
    }
    double sum;
    MPI_Allreduce(&local_sum, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return sqrt(sum);
}

// norm(Ax - b) / norm(b)
double calc_accuracy(const Double_matrix_part &A_part, const double *b, const double *x, const int myrank, const int gsize) {
    double *Ax_part = create_double_vector(A_part.rows_count, 0);
    mul_matrix_part_by_vector(Ax_part, A_part, x);

    // for (int i = 0; i < gsize; i++) {
    //     if (myrank == i) {
    //         for (int j = 0; j < A_part.rows_count; j++) {
    //             printf("%lf ", Ax_part[j]);
    //         }
    //         printf("\n\n");
    //     }
    //     MPI_Barrier(MPI_COMM_WORLD);
    // }
    // usleep(100000);
    // exit(0);

    double *Ax_sub_b = Ax_part;
    sub_vector_from_vector_part(Ax_sub_b, A_part.rows_count, b);

    // for (int i = 0; i < gsize; i++) {
    //     if (myrank == i) {
    //         for (int j = 0; j < A_part.rows_count; j++) {
    //             printf("%lf ", Ax_part[j]);
    //         }
    //         printf("\n\n");
    //     }
    //     MPI_Barrier(MPI_COMM_WORLD);
    // }
    // usleep(100000);

    double norm__Ax_sub_b__ = calc_norm_vector_part(Ax_sub_b, A_part.rows_count);
    double norm__b__ = calc_norm_vector(b, A_part.cols_count);

    // for (int i = 0; i < gsize; i++) {
    //     if (myrank == i) {
    //         printf("%f %f\n\n", norm__Ax_sub_b__, norm__b__);
    //     }
    //     MPI_Barrier(MPI_COMM_WORLD);
    // }
    // usleep(100000);

    free(Ax_sub_b);

    return norm__Ax_sub_b__ / norm__b__;
}

bool is_vectors_equal(const double *vec1, const double *vec2, int size) {
    for (int i = 0; i < size; ++i) {
        if (isnan(vec1[i]) || isnan(vec2[i])) {
            cout << "Nan values!" << endl;
            return false;
        }
        if (abs(vec1[i] - vec2[i]) > EPS) {
            return false;
        }
    }
    return true;
}

void print_vector(const double *vec, const int size) {
    cout << "[";
    for (int i = 0; i < size - 1; ++i) {
        cout << vec[i] << ", ";
    }
    cout << vec[size - 1] << "]" << endl;
}

// Ax = b
int main(int argc, char *argv[]) {
    int root = 0;
    MPI_Init(&argc, &argv);

    int gsize;
    int myrank;
    MPI_Comm_size(MPI_COMM_WORLD, &gsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    if (argc < 2) {
        if (myrank == root) {
            cerr << "You need to enter the arguments!" << endl;
        }
        MPI_Finalize();
        return 1;
    }
    const int    N   = stoi(string(argv[1]));
    const double TAU = stod(string(argv[2]));

    double *x = create_double_vector(N, 0);
    double *b = create_double_vector(N, N + 1);

    double *A;
    int    *sendcounts;
    int    *displs;

    if (myrank == root) {
        A = create_double_matrix(N, 1, 2);

        sendcounts = (int *) malloc(gsize * sizeof(int));
        int base = (int)(N / gsize);
        int addition = N % gsize;
        for (int i = 0; i < gsize; i++) {
            sendcounts[i] = base;
            if (i < addition) {
                sendcounts[i]++;
            }
        }

        displs = (int *) calloc(gsize, sizeof(int));
        for (int i = 1; i < gsize; i++) {
            displs[i] = displs[i - 1] + N * sendcounts[i - 1];
        }
    }

    int rows_count;
    MPI_Scatter(
        sendcounts, 1, MPI_INT, &rows_count, 1, MPI_INT, root, MPI_COMM_WORLD
    );

    if (myrank == root) {
        for (int i = 0; i < gsize; i++) {
            sendcounts[i] *= N;
        }
    }
    int doubles_count;
    MPI_Scatter(
        sendcounts, 1, MPI_INT, &doubles_count, 1, MPI_INT, root, MPI_COMM_WORLD
    );

    //vector<double> rows(rows_count * N);
    double *rows;
    if (rows_count) {
        rows = (double *) malloc(doubles_count * sizeof(double));
    }

    Double_matrix_part A_part;
    A_part.cols_count = N;
    A_part.rows_count = rows_count;
    A_part.rows = rows;

    MPI_Scatterv(
        A, sendcounts, displs, MPI_DOUBLE, rows, doubles_count, MPI_DOUBLE, root, MPI_COMM_WORLD
    );

    // for (int i = 0; i < gsize; i++) {
    //     if (myrank == i) {
    //         printf("%d\n", A_part.rows_count);
    //         if (A_part.rows) {
    //             for (int i = 0; i < A_part.rows_count; i++) {
    //                 for (int j = 0; j < N; j++) {
    //                     printf("%lf ", A_part.rows[i * N + j]);
    //                 }
    //                 putchar('\n');
    //             }
    //             putchar('\n');
    //         }
    //     }
    //     MPI_Barrier(MPI_COMM_WORLD);
    // }

    // До этого момента все работает.

    auto start_time = chrono::high_resolution_clock::now();
    int it = 0;
    while (calc_accuracy(A_part, b, x, myrank, gsize) >= EPS) {
        calc_next_x(A_part, b, x, TAU, myrank, gsize);
        it++;
    }
    if (myrank == root) {
        printf("%d\n", it);
    }
    auto end_time = chrono::high_resolution_clock::now();

    usleep(100000);
    if (myrank == root) {
        double *expected_vector = create_double_vector(N, 1);
        if (is_vectors_equal(expected_vector, x, N)) {
            printf("Vectors are equal.\n");
        }
        else {
            printf("Vectors are NOT equal.\n");
        }
        cout << chrono::duration_cast<chrono::milliseconds>(end_time - start_time).count() << endl;
        free(expected_vector);
    }

    if (myrank == root) {
        delete[] A;
        delete[] sendcounts;
        delete[] displs;
    }
    if (rows) {
        free(rows);
    }
    free(x);
    free(b);

    MPI_Finalize();

    return 0;
}