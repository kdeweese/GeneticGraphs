#include <iostream>
#include <cstdlib>
#include <random>
#include <iomanip>
int main(int argc, char* argv[]) {
  int n = atoi(argv[1]);
  /*int multi=1;
  if(argc > 2) {
    multi=atoi(argv[2]);
    }*/
  //double step = 2./(n-1);
  double* rhs = new double[n];
  double total=0;
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(0, 1);
  double temp;
  int precdigits=16;
  for(int i=0; i < n; ++i) {
    temp=dis(gen);
    rhs[i]=round(pow(10.,precdigits)*temp)/pow(10.,precdigits);
    total+=rhs[i];
  }
  double avg=total/n;
  for(int i=0; i < n; ++i) {
    rhs[i]=rhs[i]-avg;
  }
  
  std::cout << "%%MatrixMarket matrix array real general" << std::endl;
  std::cout << "%" << std::endl;
  std::cout << n << " " << 1 << std::endl;
  for(int i=0; i < n; ++i) {
    std::cout << std::setprecision(precdigits+1) << rhs[i] << std::endl;
  }
  delete[] rhs;
}
