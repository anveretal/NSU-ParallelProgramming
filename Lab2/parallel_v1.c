#include <stdio.h>
#include <omp.h>

int main() {
    int sum = 0;
    #pragma omp parallel
    {
        int local_sum = 0;
        #pragma omp for
        for (int i = 1; i <= 16; i++) {
            local_sum += i;
            #pragma omp critical
            printf("Thread %d: i = %d\n", omp_get_thread_num(), i);
        }

        #pragma omp atomic
        sum += local_sum;
    }
    printf("Sum = %d\n", sum);

    return 0;
}
