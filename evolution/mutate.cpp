#include "matrices.h"
#include <iostream>
#include <cstdlib>
#include <random>
#include <cstring>
int main(int argc, char *argv[]) {
  LLA* matrix = MMA_read(argv[1]);
  
  print_matrix(matrix);
  std::cout << "create mutation" << std::endl;
  LLA_mod mutate(matrix,5);
  int add=10;
  mutate.add_edges(add,1,1000);
  print_matrix(matrix);
  std::cout << "restore" << std::endl;
  //mutate.replace_edges();
  mutate.reverse_edges();
  print_matrix(matrix);
  mutate.free_memory();
  matrix->freeMemory();
  //second_mat->freeMemory();
  delete matrix;
  
  return 0;
}
