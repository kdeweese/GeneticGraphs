#include "matrices.h"
#include <iostream>
#include <cstdlib>
#include <random>
int main(int argc, char *argv[]) {
  int n = atoi(argv[1]);
  LLA matrix(n);
  int desiredm=atoi(argv[2]);
  int m=0;
  int startconnected=atoi(argv[4]);
  double minweight=atof(argv[5]);
  double maxweight=atof(argv[6]);
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(0,n-1);
  std::uniform_real_distribution<> dis2(minweight, maxweight);
 
  if(startconnected==1) {
    std::cout << "start connected" << std::endl;
    for(int i=0; i < n-1; ++i) {
      m+=mat_add_edge(&matrix,i,i+1,dis2(gen));
    }
    std::cout << "already added " << m << std::endl;
  }
  while(m < desiredm) {
    m+=mat_add_edge(&matrix,dis(gen),dis(gen),dis2(gen));
  }
  find_in_neighbors(&matrix);
  if(connected(&matrix)) {
    MMA_write(&matrix,argv[3]);
  }
  matrix.freeMemory();
  
  return 0;
}
