#include <iostream>
#include <vector>

using namespace std;



class Matrix {
public:
    Matrix(int size) : _size(size), _table(vector<double>(size * size, 1)) {}

    Matrix(int size, double value) : _size(size), _table(vector<double>(size * size, value)) {}

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

int main() {

    Matrix *p_matrix;
    p_matrix = new Matrix(3);

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            printf("%f ", (*p_matrix)(i, j));
        }
        putchar('\n');
    }

    return 0;
}