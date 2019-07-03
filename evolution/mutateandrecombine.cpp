#include "matrices.h"
#include <iostream>
#include <cstdlib>
#include <random>
#include <cstring>
#include <sys/stat.h>
std::string random_string( size_t length )
{
  auto randchar = []() -> char
    {
        const char charset[] =
        "0123456789"
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        "abcdefghijklmnopqrstuvwxyz";
        const size_t max_index = (sizeof(charset) - 1);
        return charset[ rand() % max_index ];
    };
  std::string str(length,0);
  std::generate_n( str.begin(), length, randchar );
  return str;
}
inline bool exists(const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}


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
  int fixedm=matrix->m;
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(0,matrix->n - 1);
  int row=dis(gen);
  swap_rows(matrix,matrix2,row);
  int off1=matrix->m-fixedm;
  int off2=matrix2->m-fixedm;
  int mutatect=atoi(argv[3]);
  int exchange=10;
  LLA_mod mutate(matrix,exchange+abs(off1));
  for(int i=0; i < mutatect; ++i) {
    mutate.remove_edges(exchange);
    mutate.add_edges(exchange-off1,1,1000);
    std::string outfile="testfiles/"+random_string(32)+".mtx";
    while(exists(outfile)) {
      outfile="testfiles/"+random_string(32)+".mtx";
    }
    find_in_neighbors(matrix);
    if(connected(matrix))
      MMA_write(matrix,outfile);
    mutate.reverse_edges();
    mutate.replace_edges();
  }
  mutate.free_memory();
  LLA_mod mutate2(matrix2,exchange+abs(off2));
  for(int i=0; i < mutatect; ++i) {
    mutate2.remove_edges(exchange);
    mutate2.add_edges(exchange-off2,1,1000);
    std::string outfile="testfiles/"+random_string(32)+".mtx";
    while(exists(outfile)) {
      outfile="testfiles/"+random_string(32)+".mtx";
    }
    find_in_neighbors(matrix2);
    if(connected(matrix2))
      MMA_write(matrix2,outfile);

    mutate2.reverse_edges();
    mutate2.replace_edges();
  }
  mutate2.free_memory();
  
  matrix->freeMemory();
  matrix2->freeMemory();
  
  delete matrix;
  delete matrix2;
  
  return 0;
}
