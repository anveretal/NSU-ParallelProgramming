#include <iostream>
#include <string>
#include <vector>
#include <cmath>

#define TAU 0.001
#define EPS 0.00001
#define PI 3.141592653589793

using namespace std;



class Matrix {
public:
    Matrix(int size) : _table(vector<vector<double>>(size, vector<double>(size, 1))) {
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
    vector<double> res(vec1.capacity());

    for (int i = 0; i < vec1.capacity(); ++i) {
        res[i] = vec1[i] - vec2[i];
    }
    return res;
}

vector<double> mul_vector_by_scalar(const vector<double> vec, const double scalar) {
    vector<double> res(vec.capacity());

    for (int i = 0; i < res.capacity(); ++i) {
        res[i] = vec[i] * scalar;
    }
    return res;
}

// [x]_(n+1) = [x]_n - tau * (A[x]_n - b)
vector<double> calc_next(Matrix &matrix, vector<double> &b, vector<double> &cur_vec) {
    return sub_vectors(
        cur_vec,
        mul_vector_by_scalar(
            sub_vectors(
                mul_matrix_by_vector(matrix, cur_vec), b
            ),
            TAU
        )
    );
}

double calc_norm_vector(vector<double> vec) {
    double sum = 0;
    
    for (int i = 0; i < vec.capacity(); ++i) {
        sum += vec[i] * vec[i];
    }
    return sqrt(sum);
}

double calc_accuracy(Matrix matrix, vector<double> b, vector<double> cur_vec) {
    return calc_norm_vector(sub_vectors(mul_matrix_by_vector(matrix, cur_vec), b)) / calc_norm_vector(b);
}

vector<double> get_u_vector(int size) {
    vector<double> u(size);

    for (int i = 0; i < size; ++i) {
        u[i] = sin((2 * PI * i) / size);
    }

    return u;
}

bool is_vectors_equal(vector<double> &vec1, vector<double> &vec2) {
    for (int i = 0; i < vec1.capacity(); ++i) {
        if (abs(vec1[i] - vec2[i]) > EPS) {
            return false;
        }
    }
    return true;
}

// Ax = b
int main(int argc, char *argv[]) {

    if (argc < 2) {
        cerr << "You need to enter the arguments!" << endl;
        return 1;
    }
    int matrix_size = stoi(string(argv[1]));

    Matrix matrix(matrix_size);
    for (int i = 0; i < matrix_size; ++i) {
        matrix[i][i] = 2;
    }
    //vector<double> b(matrix_size, matrix_size + 1);
    vector<double> u = get_u_vector(matrix_size);
    vector<double> b = mul_matrix_by_vector(matrix, u);
    vector<double> vec(matrix_size);

    while (calc_accuracy(matrix, b, vec) > EPS) {
        vec = calc_next(matrix, b, vec);
    }

    cout << '[' << endl;
    for (int i = 0; i < vec.capacity(); ++i) {
        cout << vec[i] << endl;
    }
    cout << ']' << endl;

    cout << (is_vectors_equal(vec, u) ? "yes" : "no") << endl;

    vector<double> a1(matrix_size, 0);
    vector<double> a2(matrix_size, 1);

    cout << (is_vectors_equal(a1, a2) ? "yes" : "no") << endl;

    return 0;
}