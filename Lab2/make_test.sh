OMP_DIR="-I/opt/homebrew/opt/libomp/include"

g++-14 -fopenmp $OMP_DIR test.cpp -o test