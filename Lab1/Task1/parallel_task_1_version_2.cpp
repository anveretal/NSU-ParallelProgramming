#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <chrono>
#include <unistd.h>
#include <mpi.h>

using namespace std;



class Context {
public:
    int root;
    int gsize;
    int myrank;

    int    N;
    double TAU;
    double EPS;

    int *part_sizes;
    int *part_displs;

    double *temp;

    int *counts;
    int *displs;

    Context() {
        part_sizes = NULL;
        temp = NULL;
    }

    ~Context() {
        if (part_sizes) {
            free(part_sizes);
        }
        if (part_displs) {
            free(part_displs);
        }
        if (temp) {
            free(temp);
        }
        if (counts) {
            free(counts);
        }
        if (displs) {
            free(displs);
        }
    }

    void clean_temp() {
        for (int i = 0; i < N; i++) {
            temp[i] = 0;
        }
    }

    void clean_counts() {
        for (int i = 0; i < gsize; i++) {
            counts[i] = 0;
        }
    }

    void clean_displs() {
        for (int i = 0; i < gsize; i++) {
            displs[i] = 0;
        }
    }
};



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

void mul_vector_part_by_scalar(const Context &context, double *vec_part, const double scalar) {
    for (int i = 0; i < context.part_sizes[context.myrank]; i++) {
        vec_part[i] *= scalar;
    }
}

void sub_vector_part_from_vector_part(const Context &context, double *vec1_part, const double *vec2_part) {
    for (int i = 0; i < context.part_sizes[context.myrank]; i++) {
        vec1_part[i] -= vec2_part[i];
    }
}

void mul_matrix_part_by_vector(const Context &context, double *res, const double *matrix_part, double *vec_part) {
    int cur_rank = context.myrank;
    for (int k = 0; k < context.gsize; k++) {
        for (int i = 0; i < context.part_sizes[context.myrank]; i++) {
            // usleep(100000);
            // printf("[%d %d %d %d]\n", context.myrank, cur_rank, context.part_displs[cur_rank], context.part_sizes[cur_rank]);
            // usleep(100000);
            for (int j = 0; j < context.part_sizes[cur_rank]; j++) {
                res[i] += matrix_part[i * context.N + context.part_displs[cur_rank] + j] * vec_part[j];
            }
        }
        MPI_Sendrecv_replace(
            vec_part, context.part_sizes[0], MPI_DOUBLE,
            (context.myrank + 1 + context.gsize) % context.gsize, 12345,
            (context.myrank - 1 + context.gsize) % context.gsize, 12345,
            MPI_COMM_WORLD, MPI_STATUS_IGNORE
        );
        // for (int i = 0; i < context.gsize; i++) {
        //     MPI_Barrier(MPI_COMM_WORLD);
        //     if (context.myrank == i) {
        //         for (int i1 = 0; i1 < context.part_sizes[context.myrank]; i1++) {
        //             std::cout << res[i1] << ' ';
        //         }
        //         std::cout << std::endl;
        //     }
        // }
        // std::cout << std::endl;

        // usleep(100000);
        // printf("[%d, ", cur_rank);
        cur_rank = (cur_rank - 1 + context.gsize) % context.gsize;
        // printf("%d]\n", cur_rank);
        // usleep(100000); 
    }
}

// x_(n+1) = x_n - TAU * (Ax_n - b)
void calc_next_x(Context &context, const double *A_part, double *x_part, const double *b_part) {
    context.clean_temp();
    mul_matrix_part_by_vector(context, context.temp, A_part, x_part);
    sub_vector_part_from_vector_part(context, context.temp, b_part);
    mul_vector_part_by_scalar(context, context.temp, context.TAU);
    sub_vector_part_from_vector_part(context, x_part, context.temp);
}

// sqrt(sum_i(vec[i]^2))
double calc_norm_vector_part(const Context &context, const double *vec_part) {
    double local_sum = 0;
    for (int i = 0; i < context.part_sizes[context.myrank]; i++) {
        local_sum += vec_part[i] * vec_part[i];
    }
    double sum;
    MPI_Allreduce(&local_sum, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return sqrt(sum);
}

// norm(Ax - b) / norm(b)
double calc_accuracy(Context &context, const double *A_part, double *x_part, const double *b_part) {
    context.clean_temp();
    mul_matrix_part_by_vector(context, context.temp, A_part, x_part);
    
    // usleep(100000);
    // for (int i = 0; i < context.gsize; i++) {
    //     MPI_Barrier(MPI_COMM_WORLD);
    //     if (context.myrank == i) {
    //         for (int i1 = 0; i1 < context.part_sizes[context.myrank]; i1++) {
    //             std::cout << context.temp[i1] << ' ';
    //         }
    //         std::cout << std::endl << std::endl;
    //     }
    // }
    // usleep(10000);
    // MPI_Finalize();
    // exit(0);

    sub_vector_part_from_vector_part(context, context.temp, b_part);

    double norm__Ax_sub_b_part__ = calc_norm_vector_part(context, context.temp);
    double norm__b_part__ = calc_norm_vector_part(context, b_part);

    return norm__Ax_sub_b_part__ / norm__b_part__;
}

bool is_vectors_equal(const Context &context, const double *vec1, const double *vec2) {
    for (int i = 0; i < context.N; ++i) {
        if (isnan(vec1[i]) || isnan(vec2[i])) {
            cout << "Nan values!" << endl;
            return false;
        }
        if (abs(vec1[i] - vec2[i]) > context.EPS) {
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
    Context context;
    context.root = 0;

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &context.gsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &context.myrank);

    if (argc < 2) {
        if (context.myrank == context.root) {
            cerr << "You need to enter the arguments!" << endl;
        }
        MPI_Finalize();
        return 1;
    }

    context.N   = stoi(string(argv[1]));
    context.TAU = stod(string(argv[2]));
    context.EPS = 0.00001;
    
    context.temp   = (double *) malloc(context.N     * sizeof(double));
    context.counts = (int    *) malloc(context.gsize * sizeof(int   ));
    context.displs = (int    *) malloc(context.gsize * sizeof(int   ));

    double *x;
    double *b;

    double *A;
    int    *sendcounts;
    int    *displs;

    if (context.myrank == context.root) {
        A = create_double_matrix(context.N, 1, 2);
        x = create_double_vector(context.N, 0);
        b = create_double_vector(context.N, context.N + 1);
    }

    context.part_sizes = (int *) malloc(context.gsize * sizeof(int));
    int base = (int)(context.N / context.gsize);
    int addition = context.N % context.gsize;
    for (int i = 0; i < context.gsize; i++) {
        context.part_sizes[i] = base;
        if (i < addition) {
            context.part_sizes[i]++;
        }
    }

    sendcounts = (int *) malloc(context.gsize * sizeof(int));
    for (int i = 0; i < context.gsize; i++) {
        sendcounts[i] = context.part_sizes[i] * context.N;
    }

    displs = (int *) calloc(context.gsize, sizeof(int));
    context.part_displs = (int *) calloc(context.gsize, sizeof(int));
    for (int i = 1; i < context.gsize; i++) {
        displs[i] = displs[i - 1] + context.part_sizes[i - 1];
        context.part_displs[i] = context.part_displs[i - 1] + context.part_sizes[i - 1];
    }

    double *x_part = (double *) malloc(context.part_sizes[0] * sizeof(double));
    MPI_Scatterv(
        x, context.part_sizes, displs, MPI_DOUBLE, x_part, context.part_sizes[context.myrank], MPI_DOUBLE, context.root, MPI_COMM_WORLD
    );

    double *b_part = (double *) malloc(context.part_sizes[0] * sizeof(double));
    MPI_Scatterv(
        b, context.part_sizes, displs, MPI_DOUBLE, b_part, context.part_sizes[context.myrank], MPI_DOUBLE, context.root, MPI_COMM_WORLD
    );

    int sendcount; // Amount of doubles that will be send to each process
    MPI_Scatter(
        sendcounts, 1, MPI_INT, &sendcount, 1, MPI_INT, context.root, MPI_COMM_WORLD
    );

    double *A_part; // (N / gsize) rows for each process
    if (context.part_sizes[context.myrank]) {
        A_part = (double *) malloc(sendcount * sizeof(double));
    }

    if (context.myrank == context.root) {
        for (int i = 1; i < context.gsize; i++) {
            displs[i] = displs[i - 1] + context.N * context.part_sizes[i - 1];
        }
    }

    MPI_Scatterv(
        A, sendcounts, displs, MPI_DOUBLE, A_part, sendcount, MPI_DOUBLE, context.root, MPI_COMM_WORLD
    );

    // for (int i = 0; i < context.gsize; i++) {
    //     MPI_Barrier(MPI_COMM_WORLD);
    //     if (context.myrank == i) {
    //         for (int i1 = 0; i1 < context.part_sizes[context.myrank]; i1++) {
    //             for (int i2 = 0; i2 < context.N; i2++) {
    //                 std::cout << A_part[i1 * context.N + i2] << ' ';
    //             }
    //             std::cout << std::endl;
    //         }
    //         std::cout << std::endl;
    //     }
    // }

    auto start_time = chrono::high_resolution_clock::now();
    int it = 0;
    while (calc_accuracy(context, A_part, x_part, b_part) >= context.EPS) {
        calc_next_x(context, A_part, x_part, b_part);
        it++;
    }
    if (context.myrank == context.root) {
        printf("%d\n", it);
    }
    auto end_time = chrono::high_resolution_clock::now();

    MPI_Gatherv(
        x_part, context.part_sizes[context.myrank], MPI_DOUBLE, x, context.part_sizes, context.part_displs,
        MPI_DOUBLE, context.root, MPI_COMM_WORLD
    );

    usleep(100000);
    if (context.myrank == context.root) {
        double *expected_vector = create_double_vector(context.N, 1);
        printf("%f\n", x[0]);
        printf("%f\n", x[3]);
        printf("%f\n", x[4]);
        if (is_vectors_equal(context, expected_vector, x)) {
            printf("Vectors are equal.\n");
        }
        else {
            printf("Vectors are NOT equal.\n");
        }
        cout << chrono::duration_cast<chrono::milliseconds>(end_time - start_time).count() << endl;
        free(expected_vector);
    }

    if (context.myrank == context.root) {
        delete[] A;
        delete[] x;
        delete[] b;
        delete[] sendcounts;
        delete[] displs;
    }
    delete[] x_part;
    delete[] b_part;
    if (A_part) {
        delete[] A_part;
    }

    MPI_Finalize();

    return 0;
}