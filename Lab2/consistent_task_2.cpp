#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <chrono>

#define TAU 0.001
#define EPS 0.00001
#define PI acos(-1.0)

using namespace std;



class Matrix {
public:
    Matrix(int size) : _table(vector<vector<double>>(size, vector<double>(size))) {
        _size = size;
    }

    Matrix(int size, double value) : _table(vector<vector<double>>(size, vector<double>(size, value))) {
        _size = size;
    }

    const int get_size() const {
        return _size;
    }

    vector<double> &operator[](int i) {
        return _table[i];
    }

private:
    int _size;
    vector<vector<double>> _table;
};

vector<double> mul_matrix_by_vector(Matrix &matrix, const vector<double> &vec) {
    vector<double> res(matrix.get_size(), 0);
    for (int i = 0; i < matrix.get_size(); ++i) {
        for (int j = 0; j < vec.size(); ++j) {
            res[i] += matrix[i][j] * vec[j];
        }
    }
    return res;
}

vector<double> sub_vectors(const vector<double> &vec1, const vector<double> &vec2) {
    vector<double> res(vec1.size());
    for (int i = 0; i < vec1.size(); ++i) {
        res[i] = vec1[i] - vec2[i];
    }
    return res;
}

vector<double> mul_vector_by_scalar(const vector<double> &vec, const double scalar) {
    vector<double> res(vec.size());
    for (int i = 0; i < res.size(); ++i) {
        res[i] = vec[i] * scalar;
    }
    return res;
}

// [x]_(n+1) = [x]_n - tau * (A[x]_n - b)
vector<double> calc_next_x(Matrix &matrix, vector<double> &b, vector<double> &x_n) {
    return sub_vectors(
        x_n,
        mul_vector_by_scalar(
            sub_vectors(
                mul_matrix_by_vector(matrix, x_n), b
            ),
            TAU
        )
    );
}

double calc_norm_vector(vector<double> vec) {
    double sum = 0;
    for (int i = 0; i < vec.size(); ++i) {
        sum += vec[i] * vec[i];
    }
    return sqrt(sum);
}

double calc_accuracy(Matrix &matrix, vector<double> &b, vector<double> &x) {
    return calc_norm_vector(sub_vectors(mul_matrix_by_vector(matrix, x), b)) / calc_norm_vector(b);
}

bool is_vectors_equal(vector<double> &vec1, vector<double> &vec2) {
    for (int i = 0; i < vec1.size(); ++i) {
        if (isnan(vec1[i]) || isnan(vec2[i])) {
            cout << "Nan values! " << endl;
            return false;
        }
        if (abs(vec1[i] - vec2[i]) > EPS) {
            return false;
        }
    }
    return true;
}

Matrix create_matrix_filled_1_and_2_on_main_diagonal(int size) {
    Matrix matrix(size, 1);
    for (int i = 0; i < size; ++i) {
        matrix[i][i] = 2;
    }
    return matrix;
}

vector<double> create_vector_u(int size) {
    vector<double> u(size);
    for (int i = 0; i < size; ++i) {
        u[i] = sinf((2 * PI * i) / size);
    }
    return u;
}

vector<double> create_vector_b(Matrix &A, vector<double> &u) {
    return mul_matrix_by_vector(A, u);
}

void print_vector(vector<double> &vec) {
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
    int N = stoi(string(argv[1]));

    // Ax = b, need to find x.
    Matrix A = create_matrix_filled_1_and_2_on_main_diagonal(N);
    vector<double> u = create_vector_u(N);
    vector<double> b = create_vector_b(A, u);

    auto start_time = chrono::high_resolution_clock::now();
    vector<double> x(N, 0); // Initial vector
    while (calc_accuracy(A, b, x) > EPS) {
        x = calc_next_x(A, b, x);
    }
    auto end_time = chrono::high_resolution_clock::now();

    if (is_vectors_equal(u, x)) {
        cout << "Vectors are equal." << endl;
        cout << chrono::duration_cast<chrono::milliseconds>(end_time - start_time).count() << endl;
    }

    return 0;
}