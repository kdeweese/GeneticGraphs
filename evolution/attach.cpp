#include "matrices.h"
#include <iostream>
#include <cstdlib>
#include <random>
#include <cstring>
#include <sys/stat.h>

int main(int argc, char *argv[]) {
  LLA* matrix = MMA_read(argv[1]);
  LLA* matrix2 = MMA_read(argv[2]);

  attach(matrix,matrix2);
  MMA_write(matrix,argv[3]);
  
  matrix->freeMemory();
  //matrix2->freeMemory();
  
  delete matrix;
  delete matrix2;
  
  return 0;
}
