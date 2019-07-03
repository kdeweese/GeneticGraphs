#include <experimental/filesystem>
#include <iostream>
#include <fstream>
#include "/home/kdeweese/geneticbase/processing/matrices.h"
#include <sys/stat.h>
#include <cstring>
int main(int argc, char *argv[]) {
  namespace stdfs = std::experimental::filesystem;
  stdfs::path path=argv[1];
  const stdfs::directory_iterator end{};
  std::vector<std::string> filenames;
  for(stdfs::directory_iterator iter{path}; iter != end; ++iter) {
    filenames.push_back(iter->path().string());
  }
  for(int i=0; i < filenames.size(); ++i) {
    int n,m,test;
    bool symmetric=false;

    std::ifstream mat_file(filenames[i].c_str());

    std::string header;
    getline(mat_file,header);

    if(header.find("symmetric") != std::string::npos) {
      symmetric=true;
    }

    if(!symmetric) {
      std::cout << "only works with symmetric matrices" << std::endl;
      return -1;
    }

    while(mat_file.peek() == '%') {
      mat_file.ignore(2048, '\n');
    }

    mat_file >> n >> test >> m;
    if(m != atoi(argv[3]) || n !=atoi(argv[2])) {
      std::cout << filenames[i] << std::endl;
    }
  }
}
