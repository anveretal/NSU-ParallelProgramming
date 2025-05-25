#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <mpi.h>
#include <unistd.h>
#include <stdbool.h>

#define MAX_TASKS 100
#define THRESHOLD 0.5

typedef struct {
    int repeatNum;
    int completed; // Используем int вместо bool для MPI
} Task;

typedef struct {
    Task tasks[MAX_TASKS];
    int count;
    pthread_mutex_t mutex;
} TaskList;

typedef struct {
    int rank;
    int size;
    TaskList *taskList;
    double *globalRes;
    int *iterCounter;
    int *tasksCompleted;
    int *allTasksDone; // Используем int вместо bool для MPI
} ThreadData;

pthread_mutex_t listMutex = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t condVar = PTHREAD_COND_INITIALIZER;
int all_processes_done = 0;

void executeTask(Task *task, double *globalRes) {
    for(int i = 0; i < task->repeatNum; i++) {
        *globalRes += sqrt(i);
    }
    task->completed = 1;
}

void* workerThread(void* arg) {
    ThreadData* data = (ThreadData*)arg;
    
    while (1) {
        pthread_mutex_lock(&data->taskList->mutex);
        
        int hasTasks = 0;
        for (int i = 0; i < data->taskList->count; i++) {
            if (!data->taskList->tasks[i].completed) {
                hasTasks = 1;
                break;
            }
        }
        
        if (!hasTasks && *data->allTasksDone) {
            pthread_mutex_unlock(&data->taskList->mutex);
            break;
        }
        
        if (!hasTasks) {
            pthread_mutex_unlock(&data->taskList->mutex);
            usleep(1000);
            continue;
        }
        
        Task* taskToExecute = NULL;
        for (int i = 0; i < data->taskList->count; i++) {
            if (!data->taskList->tasks[i].completed) {
                taskToExecute = &data->taskList->tasks[i];
                break;
            }
        }
        
        pthread_mutex_unlock(&data->taskList->mutex);
        
        if (taskToExecute) {
            executeTask(taskToExecute, data->globalRes);
            (*data->tasksCompleted)++;
        }
    }
    
    return NULL;
}

void* requesterThread(void* arg) {
    ThreadData* data = (ThreadData*)arg;
    
    while (!*data->allTasksDone) {
        pthread_mutex_lock(&data->taskList->mutex);
        
        int remainingTasks = 0;
        for (int i = 0; i < data->taskList->count; i++) {
            if (!data->taskList->tasks[i].completed) {
                remainingTasks++;
            }
        }
        
        double remainingRatio = (double)remainingTasks / data->taskList->count;
        
        if (remainingRatio < THRESHOLD) {
            for (int i = 0; i < data->size; i++) {
                if (i != data->rank) {
                    MPI_Send(&data->rank, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                }
            }
        }
        
        pthread_mutex_unlock(&data->taskList->mutex);
        
        int flag;
        MPI_Status status;
        MPI_Iprobe(MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &flag, &status);
        
        if (flag) {
            Task receivedTask;
            MPI_Recv(&receivedTask, sizeof(Task), MPI_BYTE, status.MPI_SOURCE, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
            pthread_mutex_lock(&data->taskList->mutex);
            if (data->taskList->count < MAX_TASKS) {
                data->taskList->tasks[data->taskList->count] = receivedTask;
                data->taskList->count++;
            }
            pthread_mutex_unlock(&data->taskList->mutex);
        }
        
        usleep(100000);
    }
    
    return NULL;
}

void* responderThread(void* arg) {
    ThreadData* data = (ThreadData*)arg;
    
    while (!*data->allTasksDone) {
        int flag;
        MPI_Status status;
        MPI_Iprobe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &flag, &status);
        
        if (flag) {
            int requestingRank;
            MPI_Recv(&requestingRank, 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
            pthread_mutex_lock(&data->taskList->mutex);
            
            Task taskToSend;
            int taskFound = 0;
            
            for (int i = 0; i < data->taskList->count; i++) {
                if (!data->taskList->tasks[i].completed) {
                    taskToSend = data->taskList->tasks[i];
                    data->taskList->tasks[i].completed = 1;
                    taskFound = 1;
                    break;
                }
            }
            
            pthread_mutex_unlock(&data->taskList->mutex);
            
            if (taskFound) {
                MPI_Send(&taskToSend, sizeof(Task), MPI_BYTE, requestingRank, 1, MPI_COMM_WORLD);
            }
        }
        
        usleep(100000);
    }
    
    return NULL;
}

void generateTasks(TaskList *taskList, int rank, int iterCounter, int size) {
    pthread_mutex_lock(&taskList->mutex);
    
    taskList->count = MAX_TASKS;
    for (int i = 0; i < MAX_TASKS; i++) {
        int weight = abs(rank - (iterCounter % size)) + 1;
        taskList->tasks[i].repeatNum = 1000 * weight * (1 + (rand() % 5));
        taskList->tasks[i].completed = 0;
    }
    
    pthread_mutex_unlock(&taskList->mutex);
}

int allTasksCompleted(TaskList *taskList) {
    pthread_mutex_lock(&taskList->mutex);
    
    int allDone = 1;
    for (int i = 0; i < taskList->count; i++) {
        if (!taskList->tasks[i].completed) {
            allDone = 0;
            break;
        }
    }
    
    pthread_mutex_unlock(&taskList->mutex);
    return allDone;
}

int main(int argc, char** argv) {
    int rank, size, provided;
    
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    if (provided != MPI_THREAD_MULTIPLE) {
        printf("MPI_THREAD_MULTIPLE not supported\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    TaskList taskList;
    pthread_mutex_init(&taskList.mutex, NULL);
    taskList.count = 0;
    
    double globalRes = 0;
    int iterCounter = 0;
    int tasksCompleted = 0;
    int allTasksDone = 0;
    
    ThreadData threadData = {
        .rank = rank,
        .size = size,
        .taskList = &taskList,
        .globalRes = &globalRes,
        .iterCounter = &iterCounter,
        .tasksCompleted = &tasksCompleted,
        .allTasksDone = &allTasksDone
    };
    
    pthread_t worker, requester, responder;
    pthread_create(&worker, NULL, workerThread, &threadData);
    pthread_create(&requester, NULL, requesterThread, &threadData);
    pthread_create(&responder, NULL, responderThread, &threadData);
    
    while (iterCounter < 10) {
        generateTasks(&taskList, rank, iterCounter, size);
        allTasksDone = 0;
        tasksCompleted = 0;
        
        while (!allTasksCompleted(&taskList)) {
            usleep(100000);
        }
        
        allTasksDone = 1;
        
        MPI_Barrier(MPI_COMM_WORLD);
        
        if (rank == 0) {
            printf("Iteration %d completed\n", iterCounter);
        }
        
        iterCounter++;
    }
    
    pthread_join(worker, NULL);
    pthread_join(requester, NULL);
    pthread_join(responder, NULL);
    
    int totalTasksCompleted;
    MPI_Reduce(&tasksCompleted, &totalTasksCompleted, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    
    if (rank == 0) {
        printf("\nTotal tasks completed by all processes: %d\n", totalTasksCompleted);
    }
    
    printf("Process %d completed %d tasks\n", rank, tasksCompleted);
    
    pthread_mutex_destroy(&taskList.mutex);
    MPI_Finalize();
    
    return 0;
}