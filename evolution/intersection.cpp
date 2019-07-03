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
  
  std::cout << "matrix 1" << std::endl;
  print_matrix(matrix);
  std::cout << "matrix 2" << std::endl;
  print_matrix(matrix2);

  std::cout << "finding intersection" << std::endl;
  LLA* intersect = intersection(matrix,matrix2);
  print_matrix(intersect);
  MMA_write(intersect,argv[3]);
  matrix->freeMemory();
  matrix2->freeMemory();
  intersect->freeMemory();
  delete matrix;
  delete matrix2;
  delete intersect;
  return 0;
}
