#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <mpi.h>
#include <unistd.h>
#include <stdatomic.h>

#define TASKS_COUNT 10
#define TASK_LIST_CAPACITY 100
#define REQUEST_DELAY_US 1000

typedef struct {
    int repeat_count;  // Количество повторений (вес задания)
    int task_id;
} Task;

typedef struct {
    Task* tasks;
    int size;
    int capacity;
    pthread_mutex_t lock;
    pthread_cond_t cond;
} TaskList;

typedef struct {
    double global_result;
    atomic_int tasks_processed;
    atomic_int is_done;
    int iteration;
} SharedData;

SharedData shared;
int rank, size;
TaskList local_tasks;

void init_task_list(TaskList* tl, int capacity) {
    tl->tasks = (Task*)malloc(capacity * sizeof(Task));
    if (!tl->tasks) {
        perror("Failed to allocate task list");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    tl->size = 0;
    tl->capacity = capacity;
    pthread_mutex_init(&tl->lock, NULL);
    pthread_cond_init(&tl->cond, NULL);
}

void free_task_list(TaskList* tl) {
    free(tl->tasks);
    pthread_mutex_destroy(&tl->lock);
    pthread_cond_destroy(&tl->cond);
}

void add_task(TaskList* tl, Task task) {
    pthread_mutex_lock(&tl->lock);
    if (tl->size < tl->capacity) {
        tl->tasks[tl->size++] = task;
        pthread_cond_signal(&tl->cond);
    }
    pthread_mutex_unlock(&tl->lock);
}

int get_task(TaskList* tl, Task* task) {
    pthread_mutex_lock(&tl->lock);
    while (tl->size == 0 && !shared.is_done) {
        pthread_cond_wait(&tl->cond, &tl->lock);
    }
    
    if (tl->size == 0) {
        pthread_mutex_unlock(&tl->lock);
        return 0;
    }
    
    *task = tl->tasks[--tl->size];
    pthread_mutex_unlock(&tl->lock);
    return 1;
}


void execute_task(Task task) {
    double local_res = 0.0;
    for (int i = 0; i < task.repeat_count; i++) {
        local_res += sqrt(i);
    }
    pthread_mutex_lock(&local_tasks.lock);
    shared.global_result += local_res;
    shared.tasks_processed++;
    pthread_mutex_unlock(&local_tasks.lock);
}

void* worker_thread(void* arg) {
    while (!shared.is_done) {
        Task task;
        if (get_task(&local_tasks, &task)) {
            execute_task(task);
        } else {
            int req = 1;
            for (int i = 0; i < size && !shared.is_done; i++) {
                if (i != rank) {
                    MPI_Send(&req, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                }
            }
            usleep(REQUEST_DELAY_US);
        }
    }
    return NULL;
}

void* request_handler_thread(void* arg) {
    while (!shared.is_done) {
        MPI_Status status;
        int flag;
        
        MPI_Iprobe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &flag, &status);
        if (flag) {
            int req;
            MPI_Recv(&req, 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD, &status);
            
            pthread_mutex_lock(&local_tasks.lock);
            int tasks_to_send = local_tasks.size / 2;
            
            if (tasks_to_send > 0) {
                MPI_Send(&tasks_to_send, 1, MPI_INT, status.MPI_SOURCE, 1, MPI_COMM_WORLD);
                MPI_Send(&local_tasks.tasks[local_tasks.size - tasks_to_send], 
                        tasks_to_send * sizeof(Task), MPI_BYTE, 
                        status.MPI_SOURCE, 1, MPI_COMM_WORLD);
                local_tasks.size -= tasks_to_send;
            } else {
                tasks_to_send = 0;
                MPI_Send(&tasks_to_send, 1, MPI_INT, status.MPI_SOURCE, 1, MPI_COMM_WORLD);
            }
            pthread_mutex_unlock(&local_tasks.lock);
        } else {
            usleep(REQUEST_DELAY_US);
        }
    }
    return NULL;
}

void* task_receiver_thread(void* arg) {
    while (!shared.is_done) {
        MPI_Status status;
        int flag;
        
        MPI_Iprobe(MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &flag, &status);
        if (flag) {
            int tasks_available;
            MPI_Recv(&tasks_available, 1, MPI_INT, status.MPI_SOURCE, 1, MPI_COMM_WORLD, &status);
            
            if (tasks_available > 0) {
                Task* received_tasks = (Task*)malloc(tasks_available * sizeof(Task));
                if (!received_tasks) {
                    perror("Failed to allocate received tasks");
                    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
                }
                
                MPI_Recv(received_tasks, tasks_available * sizeof(Task), MPI_BYTE, 
                        status.MPI_SOURCE, 1, MPI_COMM_WORLD, &status);
                
                pthread_mutex_lock(&local_tasks.lock);
                for (int i = 0; i < tasks_available; i++) {
                    if (local_tasks.size < local_tasks.capacity) {
                        local_tasks.tasks[local_tasks.size++] = received_tasks[i];
                    }
                }
                pthread_cond_signal(&local_tasks.cond);
                pthread_mutex_unlock(&local_tasks.lock);
                
                free(received_tasks);
            }
        } else {
            usleep(REQUEST_DELAY_US);
        }
    }
    return NULL;
}

void generate_tasks(int iteration) {
    pthread_mutex_lock(&local_tasks.lock);
    local_tasks.size = 0;
    
    for (int i = 0; i < TASKS_COUNT; i++) {
        Task task;
        task.repeat_count = 10000 * (abs(rank - (iteration % size)) + 1000);
        task.task_id = rank * TASKS_COUNT + i + iteration * TASKS_COUNT * size;

        if (local_tasks.size < local_tasks.capacity) {
            local_tasks.tasks[local_tasks.size++] = task;
        }
    }
    
    shared.tasks_processed = 0;
    pthread_cond_broadcast(&local_tasks.cond);
    pthread_mutex_unlock(&local_tasks.lock);
}

int main(int argc, char** argv) {
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    if (provided < MPI_THREAD_MULTIPLE) {
        fprintf(stderr, "Error: MPI_THREAD_MULTIPLE not supported!\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    init_task_list(&local_tasks, TASK_LIST_CAPACITY);
    shared.global_result = 0.0;
    shared.tasks_processed = 0;
    shared.is_done = 0;

    pthread_t worker, handler, receiver;
    if (pthread_create(&worker, NULL, worker_thread, NULL) ||
        pthread_create(&handler, NULL, request_handler_thread, NULL) ||
        pthread_create(&receiver, NULL, task_receiver_thread, NULL)) {
            
        fprintf(stderr, "Error creating threads\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    for (int iter = 1; iter <= 10; iter++) {
        shared.iteration = iter;
        if (rank == 0) {
            printf("Starting iteration %d\n", iter);
        }

        generate_tasks(iter);
        shared.is_done = 0;

        int local_done, global_done;
        do {
            usleep(REQUEST_DELAY_US);
            
            pthread_mutex_lock(&local_tasks.lock);
            local_done = (local_tasks.size == 0);
            pthread_mutex_unlock(&local_tasks.lock);
            
            MPI_Allreduce(&local_done, &global_done, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
        } while (!global_done);

        MPI_Barrier(MPI_COMM_WORLD);
        shared.is_done = 1;
        pthread_cond_broadcast(&local_tasks.cond);
        
        if (rank == 0) {
            printf("Iteration %d completed. Global result: %f\n", iter, shared.global_result);
        }
    }

    pthread_join(worker, NULL);
    pthread_join(handler, NULL);
    pthread_join(receiver, NULL);

    free_task_list(&local_tasks);
    MPI_Finalize();

    return 0;
}