#include "matrices.h"
#include <iostream>
#include <cstdlib>
#include <random>
#include <cstring>
int main(int argc, char *argv[]) {
  LLA* matrix = MMA_read(argv[1]);
  LLA* matrix2 = MMA_read(argv[2]);
  if(matrix->n != matrix2->n) {
    std::cerr << "matrices must be same size" << std::endl;
    matrix->freeMemory();
    matrix2->freeMemory();
    delete matrix;
    delete matrix2;
    return 1;
  }
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(0,matrix->n - 1);
  //int row=dis(gen);
  int row=1;
  std::cout << "matrix 1" << std::endl;
  print_matrix(matrix);
  std::cout << "matrix 2" << std::endl;
  print_matrix(matrix2);

  std::cout << "swapping both " << row << std::endl;
  swap_both(matrix,matrix2,row);
  std::cout << "matrix 1" <<std::endl;
  print_matrix(matrix);
  std::cout << "matrix 2" <<std::endl;
  print_matrix(matrix2);

  
  matrix->freeMemory();
  matrix2->freeMemory();
  delete matrix;
  delete matrix2;
  return 0;
}
