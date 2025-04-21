#include <iostream>
#include <vector>
#include <omp.h>

using namespace std;

// void break_in_omp_for_test() {
//     int flag = 0;
//     #pragma omp parallel for
//     for (int i = 0; i < 8; i++) {
//         if (i == 5) {
//             flag = 1;
//             break;
//         }
//     }
//     cout << flag << endl;
// }

// void local_variables_test() {
//     int a = 0;
//     #pragma omp parallel
//     {
//         #pragma omp for reduction(+:a)
//         for (int i = 0; i < 8; i++) {
//             a += 1;
//         }

//         #pragma omp critical
//         {
//             printf("Thread %d: %d\n", omp_get_thread_num(), a);
//         }
//     }
//     cout << a << endl;
// }

int main() {
    //local_variables_test();
    
    #pragma omp parallel
    {
        printf("Hello, World!!!\n");
    }
    // int sum = 0;
    // #pragma omp parallel
    // {
    //     int local_sum = 0;
    //     #pragma omp for
    //     for (int i = 1; i <= 16; i++) {
    //         local_sum += i;
    //         #pragma omp critical
    //         printf("Thread %d: i = %d\n", omp_get_thread_num(), i);
    //     }
    //     #pragma omp atomic
    //     sum += local_sum;
    // }
    // printf("Sum = %d\n", sum);

    // vector<int> vec(8, 0);
    // int cnt = 0;

    // #pragma omp parallel
    // {
    //     #pragma omp single
    //     {
    //         for (int z = 0; z < 2; z++) {
    //             #pragma omp for
    //             for (int i = 0; i < vec.size(); ++i) {
    //                 //vec[i]++;

    //                 #pragma omp atomic
    //                 vec[0]++;
    //             }
    //         }
    //     }
    // }
    // cout << cnt;

    // #pragma omp parallel
    // {
    //     #pragma omp critical
    //     {
    //         cout << "Hello, World!" << endl;
    //     }
    // }



    return 0;
}
