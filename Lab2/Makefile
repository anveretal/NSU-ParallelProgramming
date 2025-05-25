# INFO #---------------------------------#
# -march=native
# Включает оптимизации под процессор, на котором выполняется программа.
#----------------------------------------#

OMP_DIR = -I/opt/homebrew/opt/libomp/include
COM_FLG = -O3 -march=native

# TASK_1 #-------------------------------#
# Consistent
ct1:
	g++-14 $(COM_FLG) ./task_1/consistent_task_1.cpp -o ct1

# Parallel version 1
pt1v1:
	g++-14 -fopenmp $(OMP_DIR) ./task_1/parallel_task_1_version_1.cpp -o pt1v1

# Parallel version 2
pt1v2:
	g++-14 -fopenmp $(COM_FLG) $(OMP_DIR) ./task_1/parallel_task_1_version_2.cpp -o pt1v2
#----------------------------------------#



# TASK_2 #-------------------------------#
# Consistent// x_(n+1) = x_n - TAU * (Ax_n - b)
ct2:
	g++-14 $(COM_FLG) ./task_2/consistent_task_2.cpp -o ct2

# Parallel version 1
pt2v1:
	g++-14 -fopenmp $(OMP_DIR) ./task_2/parallel_task_2_version_1.cpp -o pt2v1

# Parallel version 2
pt2v2:
	g++-14 -fopenmp $(OMP_DIR) ./task_2/parallel_task_2_version_2.cpp -o pt2v2
#----------------------------------------#



# TEST #---------------------------------#
test:
	g++-14 -fopenmp $(OMP_DIR) test.cpp -o test
#----------------------------------------#