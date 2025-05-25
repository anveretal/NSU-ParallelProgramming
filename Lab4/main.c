#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define X0 -1
#define Y0 -1
#define Z0 -1
#define Dx 2.0
#define Dy 2.0
#define Dz 2.0
#define eps 1e-8
#define a 1e5
#define Nx 400
#define Ny 400
#define Nz 400

#define IDX(x, y, z) ((z) * Nx * Ny + (y) * Nx + (x))

typedef struct {
    int rank, size;
    int layerHeight;
    double hx, hy, hz;
    double* phi, * prevPhi;
    double* upLayer, * downLayer;
    double multiplier;
} Grid;

double phi_analytic(double x, double y, double z) {
    return x * x + y * y + z * z;
}

double ro(double x, double y, double z) {
    return 6.0 - a * phi_analytic(x, y, z);
}

void initializeGrid(Grid* g) {

    int baseLayerHeight = Nz / g->size;
    int extraLayers = Nz % g->size;

    g->layerHeight = baseLayerHeight + (g->rank < extraLayers ? 1 : 0);

    g->hx = Dx / (Nx - 1);
    g->hy = Dy / (Ny - 1);
    g->hz = Dz / (Nz - 1);

    size_t layerSize = Nx * Ny * g->layerHeight;
    g->phi = calloc(layerSize, sizeof(double));
    g->prevPhi = calloc(layerSize, sizeof(double));
    g->downLayer = malloc(sizeof(double) * Nx * Ny);
    g->upLayer = malloc(sizeof(double) * Nx * Ny);

    g->multiplier = 1.0 / (2.0 / (g->hx * g->hx) + 2.0 / (g->hy * g->hy) + 2.0 / (g->hz * g->hz) + a);


    for (int z = 0; z < g->layerHeight; ++z) {
        for (int y = 0; y < Ny; ++y) {
            for (int x = 0; x < Nx; ++x) {
                double val = 0;
                if (x == 0 || x == Nx - 1 || y == 0 || y == Ny - 1 ||
                    (z == 0 && g->rank == 0) || (z == g->layerHeight - 1 && g->rank == g->size - 1)) {
                    double gx = X0 + x * g->hx;
                    double gy = Y0 + y * g->hy;
                    double gz = Z0 + (z + g->rank * baseLayerHeight) * g->hz;
                    val = phi_analytic(gx, gy, gz);
                }
                g->phi[IDX(x, y, z)] = val;
                g->prevPhi[IDX(x, y, z)] = val;
            }
        }
    }
}

// Обмен соседей границами
void bordersExchange(Grid* g, MPI_Request reqs[4]) {
    // Обмен границами с нижним соседом
    if (g->rank != 0) {
        MPI_Isend(&g->prevPhi[0], Nx * Ny, MPI_DOUBLE, g->rank - 1, 0, MPI_COMM_WORLD, &reqs[0]);
        MPI_Irecv(g->downLayer, Nx * Ny, MPI_DOUBLE, g->rank - 1, 1, MPI_COMM_WORLD, &reqs[1]);
    }
    else {
        reqs[0] = MPI_REQUEST_NULL;
        reqs[1] = MPI_REQUEST_NULL;
    }
    // Обмен границами с верхним соседом
    if (g->rank != g->size - 1) {
        MPI_Isend(&g->prevPhi[IDX(0, 0, g->layerHeight - 1)], Nx * Ny, MPI_DOUBLE, g->rank + 1, 1, MPI_COMM_WORLD, &reqs[2]);
        MPI_Irecv(g->upLayer, Nx * Ny, MPI_DOUBLE, g->rank + 1, 0, MPI_COMM_WORLD, &reqs[3]);
    }
    else {
        reqs[2] = MPI_REQUEST_NULL;
        reqs[3] = MPI_REQUEST_NULL;
    }
}

// Вычисление внутренних точек
void calculateInternal(Grid* g, char* flag) {
    for (int z = 1; z < g->layerHeight - 1; ++z) {
        for (int y = 1; y < Ny - 1; ++y) {
            for (int x = 1; x < Nx - 1; ++x) {
                double gx = X0 + x * g->hx;
                double gy = Y0 + y * g->hy;
                double gz = Z0 + (z + g->rank * g->layerHeight) * g->hz;

                double lap = (
                    g->prevPhi[IDX(x - 1, y, z)] + g->prevPhi[IDX(x + 1, y, z)]) / (g->hx * g->hx) +
                    (g->prevPhi[IDX(x, y - 1, z)] + g->prevPhi[IDX(x, y + 1, z)]) / (g->hy * g->hy) +
                    (g->prevPhi[IDX(x, y, z - 1)] + g->prevPhi[IDX(x, y, z + 1)]) / (g->hz * g->hz);

                double newVal = g->multiplier * (lap - ro(gx, gy, gz));
                if (fabs(newVal - g->prevPhi[IDX(x, y, z)]) > eps) *flag = 0;
                g->phi[IDX(x, y, z)] = newVal;
            }
        }
    }
}

// Довычисление внутринних точек на нижней границе слоя (каждого процесса) после обмена границами
void calculateLowerBoundary(Grid* g, char* flag) {
    if (g->rank == 0) return;
    int z = 0;
    for (int y = 1; y < Ny - 1; ++y) {
        for (int x = 1; x < Nx - 1; ++x) {
            double gx = X0 + x * g->hx;
            double gy = Y0 + y * g->hy;
            double gz = Z0 + (z + g->rank * g->layerHeight) * g->hz;

            double lap = (
                g->prevPhi[IDX(x - 1, y, z)] + g->prevPhi[IDX(x + 1, y, z)]) / (g->hx * g->hx) +
                (g->prevPhi[IDX(x, y - 1, z)] + g->prevPhi[IDX(x, y + 1, z)]) / (g->hy * g->hy) +
                (g->downLayer[y * Nx + x] + g->prevPhi[IDX(x, y, z + 1)]) / (g->hz * g->hz);

            double newVal = g->multiplier * (lap - ro(gx, gy, gz));
            if (fabs(newVal - g->prevPhi[IDX(x, y, z)]) > eps) *flag = 0;
            g->phi[IDX(x, y, z)] = newVal;
        }
    }
}

// Довычисление внутринних точек на верхней границе слоя (каждого процесса) после обмена границами
void calculateUpperBoundary(Grid* g, char* flag) {
    if (g->rank == g->size - 1) return;
    int z = g->layerHeight - 1;
    for (int y = 1; y < Ny - 1; ++y) {
        for (int x = 1; x < Nx - 1; ++x) {
            double gx = X0 + x * g->hx;
            double gy = Y0 + y * g->hy;
            double gz = Z0 + (z + g->rank * g->layerHeight) * g->hz;

            double lap = (
                g->prevPhi[IDX(x - 1, y, z)] + g->prevPhi[IDX(x + 1, y, z)]) / (g->hx * g->hx) +
                (g->prevPhi[IDX(x, y - 1, z)] + g->prevPhi[IDX(x, y + 1, z)]) / (g->hy * g->hy) +
                (g->prevPhi[IDX(x, y, z - 1)] + g->upLayer[y * Nx + x]) / (g->hz * g->hz);

            double newVal = g->multiplier * (lap - ro(gx, gy, gz));
            if (fabs(newVal - g->prevPhi[IDX(x, y, z)]) > eps) *flag = 0;
            g->phi[IDX(x, y, z)] = newVal;
        }
    }
}


double calculateMaxDiff(Grid* g) {
    double max = 0.0;
    for (int z = 0; z < g->layerHeight; ++z) {
        for (int y = 0; y < Ny; ++y) {
            for (int x = 0; x < Nx; ++x) {
                double gx = X0 + x * g->hx;
                double gy = Y0 + y * g->hy;
                double gz = Z0 + (z + g->rank * g->layerHeight) * g->hz;
                double diff = fabs(g->phi[IDX(x, y, z)] - phi_analytic(gx, gy, gz));
                if (diff > max) max = diff;
            }
        }
    }
    double globalMax = 0.0;
    MPI_Allreduce(&max, &globalMax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    return globalMax;
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    Grid g;
    MPI_Comm_rank(MPI_COMM_WORLD, &g.rank);
    MPI_Comm_size(MPI_COMM_WORLD, &g.size);

    double startTime = 0;
    if (g.rank == 0) startTime = MPI_Wtime();

    initializeGrid(&g);

    int iter = 0;
    char flag = 0;
    do {
        ++iter;
        flag = 1;
        double* tmp = g.prevPhi;
        g.prevPhi = g.phi;
        g.phi = tmp;

        MPI_Request reqs[4];
        bordersExchange(&g, reqs);
        calculateInternal(&g, &flag);
        MPI_Waitall(4, reqs, MPI_STATUSES_IGNORE); // Если не ждать Send, то как только поток получит свои Recv, то пойдет считать границы, которые не факт, что еще все получены у других процессов.
        calculateLowerBoundary(&g, &flag);
        calculateUpperBoundary(&g, &flag);

        char globalFlag;
        MPI_Allreduce(&flag, &globalFlag, 1, MPI_CHAR, MPI_LAND, MPI_COMM_WORLD);
        flag = globalFlag;

    } while (!flag);

    double maxDiff = calculateMaxDiff(&g);
    if (g.rank == 0) {
        double endTime = MPI_Wtime();
        printf("Iterations: %d\n", iter);
        printf("Time: %lf sec\n", endTime - startTime);
        printf("Max diff: %e\n", maxDiff);
    }

    MPI_Finalize();
    return 0;
}