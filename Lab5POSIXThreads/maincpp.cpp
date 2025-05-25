#include <iostream>
#include <vector>
#include <pthread.h>
#include <mpi.h>
#include <cmath>
#include <mutex>
#include <queue>
#include <cstdlib>
#include <ctime>

#define TAG_REQUEST 0
#define TAG_TASK    1
#define TAG_STOP    2

struct Task {
    int repeatNum;
};

int rank, size;
bool stop_flag = false;

std::vector<Task> taskList;
std::mutex task_mutex;

int tasks_processed = 0;

// === РАБОЧИЙ ПОТОК ===
void* worker_thread(void*) {
    while (!stop_flag) {
        task_mutex.lock();
        if (taskList.empty()) {
            task_mutex.unlock();
            std::this_thread::sleep_for(std::chrono::milliseconds(50));
            continue;
        }

        Task t = taskList.back();
        taskList.pop_back();
        task_mutex.unlock();

        // Выполнение задания
        double res = 0;
        for (int i = 0; i < t.repeatNum; ++i)
            res += sqrt(i);
        tasks_processed++;
    }
    return nullptr;
}

// === ПОТОК-ЗАПРОСЧИК ===
void* requester_thread(void*) {
    while (!stop_flag) {
        task_mutex.lock();
        bool need_more_tasks = taskList.size() < 5; // Если мало задач — запрашиваем
        task_mutex.unlock();

        if (need_more_tasks) {
            int other_rank = rand() % size;
            if (other_rank == rank) continue;

            int request = 1;
            MPI_Send(&request, 1, MPI_INT, other_rank, TAG_REQUEST, MPI_COMM_WORLD);
        }

        std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }
    return nullptr;
}

// === ПОТОК-СЛУШАТЕЛЬ ===
void* listener_thread(void*) {
    while (!stop_flag) {
        MPI_Status status;
        int flag;
        MPI_Iprobe(MPI_ANY_SOURCE, TAG_REQUEST, MPI_COMM_WORLD, &flag, &status);

        if (flag) {
            int request;
            MPI_Recv(&request, 1, MPI_INT, status.MPI_SOURCE, TAG_REQUEST, MPI_COMM_WORLD, &status);

            task_mutex.lock();
            std::vector<Task> send_tasks;
            int count = std::min((int)taskList.size(), 5); // Отправляем до 5 задач
            for (int i = 0; i < count; ++i) {
                send_tasks.push_back(taskList.back());
                taskList.pop_back();
            }
            task_mutex.unlock();

            int num_tasks = send_tasks.size();
            MPI_Send(&num_tasks, 1, MPI_INT, status.MPI_SOURCE, TAG_TASK, MPI_COMM_WORLD);

            if (num_tasks > 0) {
                int *data = new int[num_tasks];
                for (int i = 0; i < num_tasks; ++i)
                    data[i] = send_tasks[i].repeatNum;

                MPI_Send(data, num_tasks, MPI_INT, status.MPI_SOURCE, TAG_TASK + 1, MPI_COMM_WORLD);
                delete[] data;
            }
        }

        // Приём задач
        MPI_Iprobe(MPI_ANY_SOURCE, TAG_TASK, MPI_COMM_WORLD, &flag, &status);
        if (flag) {
            int num_tasks;
            MPI_Recv(&num_tasks, 1, MPI_INT, status.MPI_SOURCE, TAG_TASK, MPI_COMM_WORLD, &status);

            if (num_tasks > 0) {
                int *data = new int[num_tasks];
                MPI_Recv(data, num_tasks, MPI_INT, status.MPI_SOURCE, TAG_TASK + 1, MPI_COMM_WORLD, &status);

                task_mutex.lock();
                for (int i = 0; i < num_tasks; ++i)
                    taskList.emplace_back(Task{data[i]});
                task_mutex.unlock();

                delete[] data;
            }
        }

        std::this_thread::sleep_for(std::chrono::milliseconds(50));
    }
    return nullptr;
}

// === ГЕНЕРАЦИЯ ЗАДАНИЙ ===
void generate_tasks(int iterCounter) {
    taskList.clear();
    srand(time(nullptr) + rank + iterCounter * 1000);

    for (int i = 0; i < 20; ++i) {
        int wave = abs(rank - (iterCounter % size));
        int base = 1000000;
        int repeatNum = base + wave * 1000000;
        taskList.push_back(Task{repeatNum});
    }
}

// === Точка входа ===
int main(int argc, char** argv) {
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int provided;
    if (provided != MPI_THREAD_MULTIPLE) {
        std::cerr << "This MPI implementation does not support MPI_THREAD_MULTIPLE mode." << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    const int NUM_ITERATIONS = 3;

    pthread_t worker_tid, requester_tid, listener_tid;

    for (int iter = 0; iter < NUM_ITERATIONS; ++iter) {
        if (rank == 0) std::cout << "=== Итерация " << iter << " ===" << std::endl;

        generate_tasks(iter);

        // Запуск потоков
        pthread_create(&worker_tid, NULL, worker_thread, NULL);
        pthread_create(&requester_tid, NULL, requester_thread, NULL);
        pthread_create(&listener_tid, NULL, listener_thread, NULL);

        // Ждём окончания обработки
        while (true) {
            task_mutex.lock();
            bool done = taskList.empty();
            task_mutex.unlock();

            if (done) break;

            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }

        // Остановка потоков
        stop_flag = true;
        pthread_join(worker_tid, NULL);
        pthread_join(requester_tid, NULL);
        pthread_join(listener_tid, NULL);

        stop_flag = false;

        // Синхронизация перед следующей итерацией
        MPI_Barrier(MPI_COMM_WORLD);
    }

    // Вывод результата
    std::cout << "Процесс " << rank << " выполнил " << tasks_processed << " заданий." << std::endl;

    MPI_Finalize();
    return 0;
}