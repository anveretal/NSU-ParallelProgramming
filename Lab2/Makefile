# -march=native
# Включает оптимизации под процессор, на котором выполняется программа.

# Consistent task 1
ct1:
	g++-14 -O3 -march=native consistent_task_1.cpp -o ct1

# Consistent task 2
ct2:
	g++-14 -O3 -march=native consistent_task_2.cpp -o ct2

# Parallel taks 1 version 1
pt1v1:
	g++-14 -fopenmp -O3 -march=native -I/opt/homebrew/opt/libomp/include parallel_task_1_version_1.cpp -o pt1v1

# Parallel taks 1 version 2
pt1v2:
	g++-14 -fopenmp -I/opt/homebrew/opt/libomp/include parallel_task_1_version_2.cpp -o pt1v2

# Parallel task 2 version 1
pt2v1:
	g++-14 -fopenmp -I/opt/homebrew/opt/libomp/include parallel_task_2_version_1.cpp -o pt2v1

# Parallel task 2 version 2
pt2v2:
	g++-14 -fopenmp -I/opt/homebrew/opt/libomp/include parallel_task_2_version_2.cpp -o pt2v2