#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <mpi.h>
#include <unistd.h>

#define DEFAULT_TASK_LIST_CAPACITY 10000
#define DELAY 100000

#define REQ_TO_RECV      0
#define SEND_RECV_TASKS  1
#define ACK              2


typedef struct {
    int id;
    int repeat_count;
    int is_completed;
} Task;

typedef struct {
    int size;
    int capacity;

    Task* tasks;

    pthread_mutex_t lock;
    pthread_cond_t  cond;
} TaskList;



int rank;
int size;

int result;

int processed_tasks_count;
int is_all_processes_done;

long long exp_weight; // expected weight
long long act_weight; // actual   weight

TaskList* task_list;



TaskList* create_task_list_by_default(void) {
    TaskList* tl = (TaskList*) malloc(sizeof(TaskList));
    tl->size = 0;
    tl->capacity = DEFAULT_TASK_LIST_CAPACITY;
    tl->tasks = (Task*) malloc(tl->capacity * sizeof(Task));
    if (!tl->tasks) {
        perror("Failed to allocate task list");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    pthread_mutex_init(&tl->lock, NULL);
    pthread_cond_init(&tl->cond, NULL);
    return tl;
}

void free_task_list(TaskList* tl) {
    free(tl->tasks);
}

void destroy_task_list(TaskList* tl) {
    free_task_list(tl);
    pthread_mutex_destroy(&tl->lock);
    pthread_cond_destroy(&tl->cond);
    free(tl);
}

// int put_task_to_task_list(Task* task) {

// }

Task* get_task_from_task_list(TaskList* tl) {
    Task* task;
    pthread_mutex_lock(&tl->lock);
    while (tl->size == 0 && !is_all_processes_done) {
        pthread_cond_wait(&tl->cond, &tl->lock);
    }
    
    if (tl->size == 0) {
        pthread_mutex_unlock(&tl->lock);
        return NULL;
    }
    
    task = &tl->tasks[--tl->size];
    pthread_mutex_unlock(&tl->lock);
    
    return task;
}


void generate_tasks(int iteration) {
    pthread_mutex_lock(&task_list->lock);
    task_list->size = 0;
    
    for (int i = 0; i < 100 * rank + 100; i++) {
        Task task;
        task.id = rank;
        task.repeat_count = 1000000 * (rank) + 1000000;
        task.is_completed = 0;

        exp_weight += task.repeat_count;

        if (task_list->size < task_list->capacity) {
            task_list->tasks[task_list->size++] = task;
        }
        else {
            perror("Task list capacity <= size");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
    }  
    pthread_cond_broadcast(&task_list->cond);
    pthread_mutex_unlock(&task_list->lock);
}

void execute_task(Task* task) {
    double local_result = 0.0;
    for (int i = 0; i < task->repeat_count; i++) {
        local_result += sqrt(i);
    }
    result += local_result;
    act_weight += task->repeat_count;

    processed_tasks_count++;
}

void* worker_thread(void* ignore) {
    while (!is_all_processes_done) {
        Task* task = get_task_from_task_list(task_list);
        if (task != NULL) {
            execute_task(task);
        }
    }
    return NULL;
}

void* handler_thread(void* ignore) {
    while (!is_all_processes_done) {
        if (task_list->size <= 100 / 2) {
            int req = 1;
            for (int i = 0; i < size && !is_all_processes_done; i++) {
                if (i != rank) {
                    MPI_Send(&req, 1, MPI_INT, i, REQ_TO_RECV, MPI_COMM_WORLD);
                }
            }
        }

        MPI_Status status;
        int is_have_req;
        
        MPI_Iprobe(MPI_ANY_SOURCE, REQ_TO_RECV, MPI_COMM_WORLD, &is_have_req, &status);
        if (is_have_req) {
            int req;
            MPI_Recv(&req, 1, MPI_INT, status.MPI_SOURCE, REQ_TO_RECV, MPI_COMM_WORLD, &status);
            
            pthread_mutex_lock(&task_list->lock);
            int send_tasks_count = task_list->size / 2;
            
            if (send_tasks_count > 0) {
                MPI_Send(&send_tasks_count, 1, MPI_INT, status.MPI_SOURCE, SEND_RECV_TASKS, MPI_COMM_WORLD);
                MPI_Send(&task_list->tasks[task_list->size - send_tasks_count], send_tasks_count * sizeof(Task), MPI_BYTE,
                        status.MPI_SOURCE, SEND_RECV_TASKS, MPI_COMM_WORLD);
                
                int ack;
                MPI_Recv(&ack, 1, MPI_INT, status.MPI_SOURCE, ACK, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                if (ack == 1) {
                    task_list->size--;
                }
            } else {
                send_tasks_count = 0;
                MPI_Send(&send_tasks_count, 1, MPI_INT, status.MPI_SOURCE, SEND_RECV_TASKS, MPI_COMM_WORLD);
            }
            pthread_mutex_unlock(&task_list->lock);
        }

        usleep(DELAY);
    }
    return NULL;
}



void* task_receiver_thread(void* arg) {
    while (!is_all_processes_done) {
        MPI_Status status;
        int is_have_req;
        
        MPI_Iprobe(MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &is_have_req, &status);
        if (is_have_req) {
            int recv_tasks_count;
            MPI_Recv(&recv_tasks_count, 1, MPI_INT, status.MPI_SOURCE, SEND_RECV_TASKS, MPI_COMM_WORLD, &status);
            
            if (recv_tasks_count > 0) {
                Task* received_tasks = (Task*)malloc(recv_tasks_count * sizeof(Task));
                MPI_Recv(received_tasks, recv_tasks_count * sizeof(Task), MPI_BYTE, status.MPI_SOURCE,
                        SEND_RECV_TASKS, MPI_COMM_WORLD, &status);
                
                int ack = 1;
                MPI_Send(&ack, 1, MPI_INT, status.MPI_SOURCE, ACK, MPI_COMM_WORLD);
                
                pthread_mutex_lock(&task_list->lock);
                for (int i = 0; i < recv_tasks_count; i++) {
                    if (task_list->size < task_list->capacity) {
                        //received_tasks[i].is_completed = 0;
                        task_list->tasks[task_list->size++] = received_tasks[i];
                    }
                    else {
                        perror("receiver: size >= capacity");
                        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
                    }
                }
                pthread_cond_signal(&task_list->cond);
                pthread_mutex_unlock(&task_list->lock);
                free(received_tasks);
            }
            else {
                int ack = 0;
                MPI_Send(&ack, 1, MPI_INT, status.MPI_SOURCE, ACK, MPI_COMM_WORLD);
            }
        }
    }
    return NULL;
}



int main(int argc, char* argv[]) {
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    if (provided < MPI_THREAD_MULTIPLE) {
        fprintf(stderr, "Error: MPI_THREAD_MULTIPLE not supported!\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    task_list = create_task_list_by_default();

    pthread_t worker, handler, receiver;
    if (pthread_create(&worker, NULL, worker_thread, NULL) ||
        pthread_create(&handler, NULL, handler_thread, NULL) ||
        pthread_create(&receiver, NULL, task_receiver_thread, NULL)) {
            
        fprintf(stderr, "Error creating threads\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    for (int iteration = 1; iteration <= 1; iteration++) {
        is_all_processes_done = 0;

        processed_tasks_count = 0;

        long long global_exp_weight;
        long long global_act_weight;
        exp_weight = 0;
        act_weight = 0;

        if (rank == 0) {
            printf("Starting iteration %d #===================#\n\n", iteration);
        }

        generate_tasks(iteration);

        MPI_Allreduce(&exp_weight, &global_exp_weight, 1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);

        int is_current_process_done;
        do {
            usleep(DELAY);

            pthread_mutex_lock(&task_list->lock);
            is_current_process_done = (task_list->size == 0);
            MPI_Allreduce(&is_current_process_done, &is_all_processes_done, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
            pthread_mutex_unlock(&task_list->lock);

            for (int i = 0; i < size; i++) {
                MPI_Barrier(MPI_COMM_WORLD);
                if (rank == i) {
                    printf("Rank = %d :: %d\n", rank, task_list->size);
                }
            }
            if (rank == 0) putchar('\n');
        } while (!is_all_processes_done);
        MPI_Barrier(MPI_COMM_WORLD);

        pthread_cond_broadcast(&task_list->cond);

        MPI_Allreduce(&act_weight, &global_act_weight, 1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);

        if (rank == 0) {
            printf("EXPECTED = %lld :: ACTUAL = %lld\n", global_exp_weight, global_act_weight);
        }

        for (int i = 0; i < size; i++) {
            MPI_Barrier(MPI_COMM_WORLD);
            if (rank == i) {
                printf("Rank = %d\n", rank);
                printf("Expected weight = %lld :: Actual weight = %lld\n", exp_weight, act_weight);
                printf("Tasks processed = %d\n", processed_tasks_count);
                printf("\n");
            }
        }
    }

    pthread_join(worker, NULL);
    pthread_join(handler, NULL);
    pthread_join(receiver, NULL);

    destroy_task_list(task_list);
    MPI_Finalize();

    return 0;
}