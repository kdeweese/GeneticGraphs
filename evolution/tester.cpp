#include "matrices.h"
#include <iostream>
int main() {
  int n = 10;
  LLA matrix(n);
  mat_add_edge(&matrix, 0, 1, 1.5);
  mat_add_edge(&matrix, 0, 4, 4);
  mat_add_edge(&matrix, 0, 8, 2);
  mat_add_edge(&matrix, 0, 6, 3);
  mat_add_edge(&matrix, 0, 8, 4);
  mat_add_edge(&matrix, 0, 4, 4);
  mat_add_edge(&matrix, 9, 5, 1);
  mat_add_edge(&matrix, 9, 1, 3);
  mat_add_edge(&matrix, 4, 4, 4);
  mat_add_edge(&matrix, 5, 10, 4);
  
  print_matrix(&matrix);
  std::cout << "done" << std::endl;
  matrix.freeMemory();
  
  return 0;
}
