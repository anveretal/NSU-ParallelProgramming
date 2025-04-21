#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <chrono>
#include <omp.h>

#define EPS 0.00001
#define PI acos(-1.0)

using namespace std;



class Matrix {
public:
    Matrix(int size) : _size(size), _table(vector<double>(size * size)) {}

    Matrix(int size, double value) : _size(size), _table(vector<double>(size * size, 1)) {}

    const int get_size() const {
        return _size;
    }

    const double &operator()(int i, int j) const {
        return _table[i * _size + j];
    }

    double &operator()(int i, int j) {
        return _table[i * _size + j];
    }

private:
    int _size;
    vector<double> _table;
};

vector<double> mul_matrix_by_vector(const Matrix &matrix, const vector<double> &vec) {
    vector<double> res(matrix.get_size(), 0);
    #pragma omp parallel for
    for (int i = 0; i < matrix.get_size(); ++i) {
        for (int j = 0; j < vec.size(); ++j) {
            res[i] += matrix(i, j) * vec[j];
        }
    }
    return res;
}

vector<double> sub_vectors(const vector<double> &vec1, const vector<double> &vec2) {
    vector<double> res(vec1.size());
    #pragma omp parallel for
    for (int i = 0; i < vec1.size(); ++i) {
        res[i] = vec1[i] - vec2[i];
    }
    return res;
}

vector<double> mul_vector_by_scalar(const vector<double> &vec, const double scalar) {
    vector<double> res(vec.size());
    #pragma omp parallel for
    for (int i = 0; i < vec.size(); ++i) {
        res[i] = vec[i] * scalar;
    }
    return res;
}

// x_(n+1) = x_n - TAU * (Ax_n - b)
vector<double> calc_next_x(const Matrix &A, const vector<double> &b, const vector<double> &x, const double TAU) {
    vector<double> Ax = mul_matrix_by_vector(A, x);
    vector<double> Ax_sub_b = sub_vectors(Ax, b);
    vector<double> tau_mul__Ax_sub_b = mul_vector_by_scalar(Ax_sub_b, TAU);

    return sub_vectors(x, tau_mul__Ax_sub_b);
}

// sqrt(sum_i(vec[i]^2))
double calc_norm_vector(const vector<double> &vec) {
    double sum = 0;
    #pragma omp parallel for reduction(+:sum)
    for (int i = 0; i < vec.size(); ++i) {
        sum += vec[i] * vec[i];
    }
    return sqrt(sum);
}

// norm(Ax - b) / norm(b)
double calc_accuracy(const Matrix &A, const vector<double> &b, const vector<double> &x) {
    vector<double> Ax = mul_matrix_by_vector(A, x);
    vector<double> Ax_sub_b = sub_vectors(Ax, b);
    double norm_Ax_sub_b = calc_norm_vector(Ax_sub_b);
    double norm_b = calc_norm_vector(b);

    return norm_Ax_sub_b / norm_b;
}

bool is_vectors_equal(const vector<double> &vec1, const vector<double> &vec2) {
    for (int i = 0; i < vec1.size(); ++i) {
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

void fill_matrix_2_on_main_diagonal(Matrix &matrix, int size) {
    #pragma omp parllel for
    for (int i = 0; i < size; ++i) {
        matrix(i, i) = 2;
    }
}

vector<double> create_vector_u(int size) {
    vector<double> u(size);
    #pragma omp parallel for
    for (int i = 0; i < size; ++i) {
        u[i] = sinf((2 * PI * i) / size);
    }
    return u;
}

vector<double> create_vector_b(Matrix &A, vector<double> &u) {
    return mul_matrix_by_vector(A, u);
}

void print_vector(const vector<double> &vec) {
    cout << "[";
    for (int i = 0; i < vec.size() - 1; ++i) {
        cout << vec[i] << ", ";
    }
    cout << vec[vec.size() - 1] << "]" << endl;
}



// Ax = b
int main(int argc, char *argv[]) {
    if (argc < 2) {
        cerr << "You need to enter the arguments!" << endl;
        return 1;
    }
    const int N = stoi(string(argv[1]));
    const double TAU = stod(string(argv[2]));

    // Ax = b, need to find x.
    Matrix A(N, 1);
    fill_matrix_2_on_main_diagonal(A, N);
    vector<double> u = create_vector_u(N);
    vector<double> b = create_vector_b(A, u);
    vector<double> x(N, 0); // Initial vector

    auto start_time = chrono::high_resolution_clock::now();

    vector<double> Ax(A.get_size());

    double sum__Ax_sub_b = 0;
    double sum_b = 0;
    #pragma omp parallel
    {
        while (true) {
            // Ax
            #pragma omp for
            for (int i = 0; i < A.get_size(); ++i) {
                Ax[i] = 0;
                for (int j = 0; j < x.size(); ++j) {
                    Ax[i] += A(i, j) * x[j];
                }
            }
            // Calculating norms
            #pragma omp single
            {
                sum__Ax_sub_b = 0;
                sum_b = 0;
            }
            #pragma omp for reduction(+:sum__Ax_sub_b, sum_b)
            for (int i = 0; i < A.get_size(); ++i) {
                double Ax_sub_b = (Ax[i] - b[i]);
                sum__Ax_sub_b += Ax_sub_b * Ax_sub_b;
                sum_b += b[i] * b[i];
            }
            //sumb = sqrt(sumb); ? Потому что пока один поток перезаписывает, другой тоже ее перезаписывает
            double norm__Ax_sub_b = sqrt(sum__Ax_sub_b); // тут все хорошо потому, что нормы — локальные для каждого потока
            double norm_b = sqrt(sum_b);                 // а суммы — общая переменная только для чтения

            // Check accuracy:  norm(Ax - b) / norm(b) < eps
            if (norm__Ax_sub_b / norm_b < EPS) {
                break;
            }

            // Calculating next vector x
            // x_(n+1) = x_n - TAU * (Ax_n - b)
            #pragma omp for
            for (int i = 0; i < A.get_size(); ++i) {
                x[i] = x[i] - TAU * (Ax[i] - b[i]);
            }
        }
    }
    auto end_time = chrono::high_resolution_clock::now();

    if (is_vectors_equal(u, x)) {
        cout << "Vectors are equal." << endl;
        cout << chrono::duration_cast<chrono::milliseconds>(end_time - start_time).count() << endl;
    }
    else {
        cout << "Vectors are NOT equal!" << endl;
    }

    return 0;
}