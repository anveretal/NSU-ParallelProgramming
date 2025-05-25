#include<cmath>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <mpi.h>

int* dataStarts;
int* dataLengths;
int proc_rank;
int num_of_proc;
int N;

using namespace std;

void mult(double* a, const int size, const double tau) {
    for (int i = 0; i < size; i++) {
        a[i] *= tau;
    }
}

void clear_mult_res(double* a)
{
    for (int i = 0;i < dataLengths[0];i++)
    {
        a[i] = 0;
    }
}

void shiftData(double* x)
{
    MPI_Sendrecv_replace(x, dataLengths[0], MPI_DOUBLE, (proc_rank + num_of_proc + 1) % num_of_proc, 3,
    (proc_rank + num_of_proc - 1) % num_of_proc, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

}

void mult(double* matrix, double* vect, double* res, const int size) {
    int currRank = proc_rank;
    int offset = 0;
    clear_mult_res(res);
    for (int iter = 0;iter < num_of_proc;iter++)
    {
        offset = dataStarts[currRank];
        for (int i = 0;i < dataLengths[proc_rank];i++)
        {
            double sum = 0;
            for (int j = dataStarts[currRank];j < dataLengths[currRank] + dataStarts[currRank];j++)
            {
                sum += matrix[i * size + j] * vect[j - offset];
            }

            res[i] += sum;
        }
        shiftData(vect);
        currRank = (currRank + num_of_proc - 1) % num_of_proc;
    }

}

void sub(const double* a, const double* b, double* c, const int size) {
    for (int i = 0; i < size; i++) {
        c[i] = a[i] - b[i];
    }
}

struct Context {
    double* A;
    double* x;
    double* b;
    double tau;
    double* multRes;
    double* subRes;
    const int size;
    const double epsilon;

    Context(const int n, const double epsilon) : size(n), epsilon(epsilon) {
        A = new double[size * dataLengths[proc_rank]];
        x = new double[dataLengths[0]];
        b = new double[dataLengths[0]];
        tau = 0.000001;
        multRes = new double[dataLengths[0]];
        subRes = new double[dataLengths[0]];
        for (int i = 0; i < dataLengths[proc_rank]; i++) {
            x[i] = 0;
            b[i] = size + 1;
        }
        for (int i = 0; i < dataLengths[proc_rank]; i++) {
            for (int j = 0; j < size; j++) {
                if (dataStarts[proc_rank] + i == j) {
                    A[i * size + j] = 2.0;
                }
                else {
                    A[i * size + j] = 1.0;
                }
            }

        }
    }
    ~Context() {
        delete[] A;
        delete[] x;
        delete[] b;
        delete[] multRes;
        delete[] subRes;
    }
};



double Norm(double* A)
{
    double sum = 0;
    for (int i = 0;i < dataLengths[proc_rank];i++)
    {
        sum += A[i] * A[i];
    }
    double result = 0;
    MPI_Allreduce(&sum, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);


    result = sqrt(result);

    return result;
}

bool isCloseEnough(Context& cont) {
    mult(cont.A, cont.x, cont.multRes, cont.size);
    sub(cont.multRes, cont.b, cont.subRes, dataLengths[proc_rank]);
    double res = Norm(cont.subRes) / Norm(cont.b);


    return res < cont.epsilon;

}
void next(Context& cont) {
    mult(cont.A, cont.x, cont.multRes, cont.size);
    sub(cont.multRes, cont.b, cont.subRes, dataLengths[proc_rank]);
    mult(cont.subRes, dataLengths[proc_rank], cont.tau);
    sub(cont.x, cont.subRes, cont.x, dataLengths[proc_rank]);
}
void InitRowsPerProcess(int n) {
    dataStarts = new int[num_of_proc];
    dataLengths = new int[num_of_proc];
    int perProcess = n / num_of_proc;
    int ost = n % num_of_proc;
    int start = 0;
    for (int i = 0; i < ost; i++) { // 0 - 0, 1 - perProcess + 1, 2 - 2*PerProcess + 2
        dataStarts[i] = start;
        dataLengths[i] = perProcess + 1;
        start += perProcess + 1;
    }
    for (int i = ost; i < num_of_proc; i++) {
        dataStarts[i] = start;
        dataLengths[i] = perProcess;
        start += perProcess;
    }

}

int main(int argc, char* argv[]) {
    double begin, end;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_of_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);
    const double epsilon = pow(10, -5);
    N = atoi(argv[1]);
    InitRowsPerProcess(N);
    Context cont = Context(N, epsilon);

    begin = MPI_Wtime();
    while (!isCloseEnough(cont)) {
        next(cont);
    }
    end = MPI_Wtime();
    printf("%f\n", cont.x[0]);
    printf("%f\n", cont.x[3]);
    printf("%f\n", cont.x[4]);

    cout << "Time diff = " << (end - begin) << "[s]" << endl;
    delete[] dataStarts;
    delete[] dataLengths;
    MPI_Finalize();
}