#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <mpi.h>
#include <unistd.h>
#include <stdatomic.h>

#define TASKS_COUNT 10
#define TASK_LIST_CAPACITY 100


typedef struct {
    int repeat_count;  // Количество повторений (вес задания)
    int task_id;
} Task;

typedef struct {
    Task* tasks;
    int size;
    int capacity;
    pthread_mutex_t lock;
} TaskList;


double global_result = 0;
int iteration_number = 0;

int rank, size;

TaskList local_tasks;
atomic_int tasksProcessed = 0;
atomic_int is_all_local_tasks_done = 0;


void init_task_list(TaskList* tl, int capacity) {
    tl->tasks = (Task*)malloc(capacity * sizeof(Task));
    tl->size = 0;
    tl->capacity = capacity;
    pthread_mutex_init(&tl->lock, NULL);
}

void freeTaskList(TaskList* tl) {
    free(tl->tasks);
    pthread_mutex_destroy(&tl->lock);
}

void addTask(TaskList* tl, Task task) {
    pthread_mutex_lock(&tl->lock);
    if (tl->size < tl->capacity) {
        tl->tasks[tl->size++] = task;
    }
    pthread_mutex_unlock(&tl->lock);
}

int getTask(TaskList* tl, Task* task) {
    pthread_mutex_lock(&tl->lock);
    if (tl->size == 0) {
        pthread_mutex_unlock(&tl->lock);
        return 0;  // Нет заданий
    }
    *task = tl->tasks[--tl->size];
    pthread_mutex_unlock(&tl->lock);
    return 1;  // Задание получено
}

void executeTask(Task task) {
    double localRes = 0;
    for (int i = 0; i < task.repeat_count; i++) {
        localRes += sqrt(i);
    }
    pthread_mutex_lock(&local_tasks.lock);
    global_result += localRes;
    tasksProcessed++;
    pthread_mutex_unlock(&local_tasks.lock);
}

void* worker_thread(void* arg) {
    while (!is_all_local_tasks_done) {
        Task task;
        if (getTask(&local_tasks, &task)) {
            executeTask(task);
        } else {
            // Попробовать получить задания от других процессов
            int req = 1;  // Запрос на задания
            for (int i = 0; i < size; i++) {
                if (i != rank) {
                    MPI_Send(&req, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                }
            }
            usleep(1000);
        }
    }
    return NULL;
}

// Поток для обработки входящих запросов
void* requestHandlerThread(void* arg) {
    while (!is_all_local_tasks_done) {
        MPI_Status status;
        int req;
        // Проверяем есть ли входящие запросы (неблокирующий прием)
        int flag;
        MPI_Iprobe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &flag, &status);
        if (flag) {
            MPI_Recv(&req, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
            
            pthread_mutex_lock(&local_tasks.lock);
            int tasksToSend = local_tasks.size / 2;  // Отправляем половину заданий
            if (tasksToSend > 0) {
                MPI_Send(&tasksToSend, 1, MPI_INT, status.MPI_SOURCE, 1, MPI_COMM_WORLD);
                MPI_Send(&local_tasks.tasks[local_tasks.size - tasksToSend], 
                         tasksToSend * sizeof(Task), MPI_BYTE, 
                         status.MPI_SOURCE, 1, MPI_COMM_WORLD);
                local_tasks.size -= tasksToSend;
            } else {
                tasksToSend = 0;
                MPI_Send(&tasksToSend, 1, MPI_INT, status.MPI_SOURCE, 1, MPI_COMM_WORLD);
            }
            pthread_mutex_unlock(&local_tasks.lock);
        } else {
            usleep(1000);
        }
    }
    return NULL;
}

// Поток для приема заданий от других процессов
void* taskReceiverThread(void* arg) {
    while (!is_all_local_tasks_done) {
        MPI_Status status;
        int tasksAvailable;
        // Проверяем есть ли входящие задания (неблокирующий прием)
        int flag;
        MPI_Iprobe(MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &flag, &status);
        if (flag) {
            MPI_Recv(&tasksAvailable, 1, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);
            if (tasksAvailable > 0) {
                Task* receivedTasks = (Task*)malloc(tasksAvailable * sizeof(Task));
                MPI_Recv(receivedTasks, tasksAvailable * sizeof(Task), MPI_BYTE, 
                        status.MPI_SOURCE, 1, MPI_COMM_WORLD, &status);
                
                pthread_mutex_lock(&local_tasks.lock);
                for (int i = 0; i < tasksAvailable; i++) {
                    if (local_tasks.size < local_tasks.capacity) {
                        local_tasks.tasks[local_tasks.size++] = receivedTasks[i];
                    }
                }
                pthread_mutex_unlock(&local_tasks.lock);
                
                free(receivedTasks);
            }
        } else {
            usleep(1000);
        }
    }
    return NULL;
}

void generate_tasks(int i) {
    pthread_mutex_lock(&local_tasks.lock);
    local_tasks.size = 0;
    
    for (int i = 0; i < TASKS_COUNT; i++) {
        Task task;
        task.repeat_count = 10000 * (abs(rank - (i % size)) + 1000);
        task.task_id = rank * TASKS_COUNT + i;

        if (local_tasks.size < local_tasks.capacity) {
            local_tasks.tasks[local_tasks.size++] = task;
        }
    }
    tasksProcessed = 0;
    pthread_mutex_unlock(&local_tasks.lock);
}

int main(int argc, char** argv) {
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    if (provided < MPI_THREAD_MULTIPLE) {
        printf("Warning: MPI_THREAD_MULTIPLE not supported!\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    init_task_list(&local_tasks, TASK_LIST_CAPACITY);

    pthread_t worker, handler, receiver;
    pthread_create(&worker, NULL, worker_thread, NULL);
    pthread_create(&handler, NULL, requestHandlerThread, NULL);
    pthread_create(&receiver, NULL, taskReceiverThread, NULL);

    for (iteration_number = 1; iteration_number <= 10; iteration_number++) {
        if (rank == 0) {
            printf("Starting iteration %d\n", iteration_number);
        }

        generate_tasks(iteration_number);
        is_all_local_tasks_done = 0;

        int is_global_tasks_done;
        do {
            usleep(1000);

            pthread_mutex_lock(&local_tasks.lock);
            int local_size = local_tasks.size;
            pthread_mutex_unlock(&local_tasks.lock);
            
            int local_done = (local_size == 0);
            MPI_Allreduce(&local_done, &is_global_tasks_done, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
        } while (!is_global_tasks_done);

        MPI_Barrier(MPI_COMM_WORLD);
    }

    is_all_local_tasks_done = 1;
    pthread_join(worker, NULL);
    pthread_join(handler, NULL);
    pthread_join(receiver, NULL);

    freeTaskList(&local_tasks);
    MPI_Finalize();
    return 0;
}