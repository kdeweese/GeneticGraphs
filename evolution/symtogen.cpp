#include "matrices.h"
#include <iostream>
#include <cstdlib>

int main(int argc, char* argv[]) {
  LLA* matrix = MMA_read(argv[1]);
  MMA_write(matrix, argv[2]);
  delete matrix;
  return 0;
}
