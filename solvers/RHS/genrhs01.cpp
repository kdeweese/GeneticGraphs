#include <iostream>
#include <cstdlib>
int main(int argc, char* argv[]) {
  int n = atoi(argv[1]);
  int multi=1;
  if(argc > 2) {
    multi=atoi(argv[2]);
  }
  double step = 2./(n-1);
  double* rhs = new double[n];
  double total=-1;
  rhs[0]=-1;
  for(int i=1; i < n; ++i) {
    rhs[i]=-1+step*i;
    total+=rhs[i];
  }
  double avg=total/n;
  
  std::cout << "%%MatrixMarket matrix array real general" << std::endl;
  std::cout << "%" << std::endl;
  std::cout << n+(n-1)*(multi-1) << " " << 1 << std::endl;
  std::cout << -1-avg << std::endl;
  while(multi > 0) {
    for(int i=1; i < n; ++i) {
      if(multi>1 && i==n-1) {
	std::cout << 0 << std::endl;
      }
      else {
	std::cout << -1+step*i-avg<< std::endl;
      }
    }
    multi-=1;
  }
  delete[] rhs;
}
