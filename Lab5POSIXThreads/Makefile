.PHONY: all run clean

all: lab5

lab5: main.c
	mpicc -o lab5 main3.c -pthread -lm

run:
	mpiexec -n 4 ./lab5

clean:
	rm lab5