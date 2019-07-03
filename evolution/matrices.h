#include <cstddef>
#include <iostream>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cstring>
#include <random>
#include <algorithm>
#include <iomanip>
#include <queue>

struct LL_node {
  LL_node* next_ptr;
  LL_node* prev_ptr;
  double weight;
  int idx;
};

struct LLA {
  int n;
  int m;
  
  int *num_entries;
  int *in_num_entries;
  LL_node **neighbor_list;
  LL_node **in_neighbor_list;

  LLA(int size) {
    n=size;
    m=0;
    num_entries = new int[n];
    in_num_entries = new int[n];
    neighbor_list = new LL_node*[n];
    in_neighbor_list = new LL_node*[n];
    for(int i=0; i < n; ++i) {
      neighbor_list[i]=NULL;
      in_neighbor_list[i]=NULL;
      num_entries[i]=0;
      in_num_entries[i]=0;
    }
  }


  void freeMemory() const
  {
    
    for(int i=0; i < n; ++i) {
      LL_node *temp_node = neighbor_list[i];
      LL_node *node;
      for(int j=0; j < num_entries[i]; ++j) {
        node = temp_node;
        temp_node=temp_node->next_ptr;
        delete node;
      }
    }
    
    for(int i=0; i < n; ++i) {
      LL_node *temp_node = in_neighbor_list[i];
      if(temp_node != NULL) {
	LL_node *node;
	for(int j=0; j < in_num_entries[i]; ++j) {
	  node = temp_node;
	  temp_node=temp_node->next_ptr;
	  delete node;
	}
      }
    }
    
    delete[] neighbor_list;
    delete[] num_entries;
    delete[] in_num_entries;
    delete[] in_neighbor_list;
  }
};



int mat_add_edge(LLA* mat, int source, int dest, double weight) {
  if(source >= mat->n || dest >= mat->n || source == dest) {
    return 0;
  }
  if(dest < source) {
    std::swap(dest,source);
  }
  LL_node *temp_node = mat->neighbor_list[source];
  //list was empty
  if(temp_node==NULL) {
    LL_node* node = new LL_node;
    node->weight=weight;
    node->idx=dest;
    node->next_ptr=NULL;
    node->prev_ptr=NULL;
    mat->neighbor_list[source]=node;
    mat->m++;
    mat->num_entries[source]++;
    return 1;
  }
  else if(temp_node->idx == dest) {
    return 0;
  }
  //add new first element
  else if(temp_node->idx > dest) {
    LL_node* node = new LL_node;
    node->weight=weight;
    node->idx=dest;
    node->next_ptr=temp_node;
    node->prev_ptr=NULL;
    mat->neighbor_list[source]=node;
    temp_node->prev_ptr=node;
    mat->m++;
    mat->num_entries[source]++;
    return 1;
  }
  //add some other element
  else {
    while(temp_node->next_ptr != NULL && temp_node->next_ptr->idx < dest) {
      temp_node=temp_node->next_ptr;
    }
    
    
    if(temp_node->next_ptr == NULL) {
      LL_node* node = new LL_node;
      node->idx=dest;
      node->prev_ptr=temp_node;
      node->next_ptr=NULL;
      node->weight=weight;
      temp_node->next_ptr=node;
      mat->m++;
      mat->num_entries[source]++;
      return 1;
    }
    else if(temp_node->next_ptr->idx == dest) {
      return 0;
    }
    else {
      LL_node* node = new LL_node;
      node->idx=dest;
      node->next_ptr=temp_node->next_ptr;
      node->prev_ptr=temp_node;
      node->weight=weight;
      temp_node->next_ptr=node;
      node->next_ptr->prev_ptr=node;
      mat->m++;
      mat->num_entries[source]++;
      return 1;
    }
  }
    
}

int mat_add_edge_in(LLA* mat, int source, int dest, double weight) {
  if(source >= mat->n || dest >= mat->n || source == dest) {
    return 0;
  }
  if(dest > source) {
    std::swap(dest,source);
  }
  LL_node *temp_node = mat->in_neighbor_list[source];
  //list was empty
  if(temp_node==NULL) {
    LL_node* node = new LL_node;
    node->weight=weight;
    node->idx=dest;
    node->next_ptr=NULL;
    node->prev_ptr=NULL;
    mat->in_neighbor_list[source]=node;
    mat->in_num_entries[source]++;
    return 1;
  }
  else if(temp_node->idx == dest) {
    return 0;
  }
  //add new first element
  else if(temp_node->idx > dest) {
    LL_node* node = new LL_node;
    node->weight=weight;
    node->idx=dest;
    node->next_ptr=temp_node;
    node->prev_ptr=NULL;
    mat->in_neighbor_list[source]=node;
    temp_node->prev_ptr=node;
    mat->in_num_entries[source]++;
    return 1;
  }
  //add some other element
  else {
    while(temp_node->next_ptr != NULL && temp_node->next_ptr->idx < dest) {
      temp_node=temp_node->next_ptr;
    }
    
    
    if(temp_node->next_ptr == NULL) {
      LL_node* node = new LL_node;
      node->idx=dest;
      node->prev_ptr=temp_node;
      node->next_ptr=NULL;
      node->weight=weight;
      temp_node->next_ptr=node;
      mat->in_num_entries[source]++;
      return 1;
    }
    else if(temp_node->next_ptr->idx == dest) {
      return 0;
    }
    else {
      LL_node* node = new LL_node;
      node->idx=dest;
      node->next_ptr=temp_node->next_ptr;
      node->prev_ptr=temp_node;
      node->weight=weight;
      temp_node->next_ptr=node;
      node->next_ptr->prev_ptr=node;
      mat->in_num_entries[source]++;
      return 1;
    }
  }
    
}


int mat_remove_edge(LLA* mat, int source, int dest) {
  if(source >= mat->n || dest >= mat->n || source == dest) {
    return 0;
  }
  if(dest < source) {
    std::swap(dest,source);
  }
  LL_node *temp_node = mat->neighbor_list[source];
  while(temp_node != NULL) {
    if(temp_node->idx == dest) {
      if(temp_node->prev_ptr !=NULL) {
        temp_node->prev_ptr->next_ptr=temp_node->next_ptr;
      }
      else {
        mat->neighbor_list[source]=temp_node->next_ptr;
      }
      if(temp_node->next_ptr!=NULL) {
        temp_node->next_ptr->prev_ptr=temp_node->prev_ptr;
      }
      (mat->num_entries[source])--;
      (mat->m)--;
      delete temp_node;
      return 1;
    }
    else {
      temp_node = temp_node->next_ptr;
    }
  }
  return 0;
}

struct LLA_mod {
  LLA* mat;
  int *removed_i;
  int *removed_j;
  double *removed_w;
  int removed_ct;

  int *added_i;
  int *added_j;
  int added_ct;

  LLA_mod(LLA *A, int target) {
    mat=A;
    removed_ct=0;
    added_ct=0;
    removed_i = new int[target];
    removed_j = new int[target];
    removed_w = new double[target];
    added_i = new int[target];
    added_j = new int[target];

  }

  void remove_edges(int m) {
    removed_ct=0;
    std::random_device rd;
    std::mt19937 gen(rd());
    
    for(int k=0; k < m; ++k) {

      std::uniform_int_distribution<> dis(0,mat->m - 1);
      int e = dis(gen);
      int i=0;
      while(e >= mat->num_entries[i]) {
        e=e-(mat->num_entries[i]);
        ++i;
      }
      LL_node *temp_node = mat->neighbor_list[i];
      for(int j=0; j < e; ++j) {
        temp_node=temp_node->next_ptr;
      }
      removed_i[removed_ct]=i;
      removed_j[removed_ct]=temp_node->idx;
      removed_w[removed_ct]=temp_node->weight;
      
      if(e==0) {
        mat->neighbor_list[i]=temp_node->next_ptr;
      }
      else {
        temp_node->prev_ptr->next_ptr=temp_node->next_ptr;
      }
      if(temp_node->next_ptr !=NULL)
        temp_node->next_ptr->prev_ptr=temp_node->prev_ptr;
      
      delete temp_node;
      (mat->m)--;
      (mat->num_entries[i])--;
      removed_ct++;
    }
  }

  void add_edges(int m, double minweight, double maxweight) {
    added_ct=0;
    int success;
    int i,j;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(minweight, maxweight);
    std::uniform_int_distribution<> dis2(0,mat->n - 1);
    while(added_ct < m) {
      success=0;
      i=dis2(gen);
      j=dis2(gen);
      if(j < i)
        std::swap(i,j);

      success=mat_add_edge(mat,i,j,dis(gen));
      if(success==1) {
	added_i[added_ct]=i;
        added_j[added_ct]=j;
        added_ct++;
      }
    }

  }

  void replace_edges() {
    for(int i=0; i < removed_ct; ++i) {
      if(mat_add_edge(mat,removed_i[i],removed_j[i],removed_w[i])==0) {
	std::cerr << "complain " << removed_i[i] << " " << removed_j[i] << std::endl;
      }
    }
    removed_ct=0;
  }

  void reverse_edges() {
    for(int i=0; i < added_ct; ++i) {
      mat_remove_edge(mat,added_i[i],added_j[i]);
    }
    added_ct=0;
  }

  void free_memory() {
    delete[] removed_i;
    delete[] removed_j;
    delete[] removed_w;
    delete[] added_i;
    delete[] added_j;
  }
};


void swap_rows(LLA *mat, LLA *mat2, int row) {
  LL_node *temp_node = mat->neighbor_list[row];
  int temp_ct = mat->num_entries[row];

  
  mat->neighbor_list[row]=mat2->neighbor_list[row];
  mat->num_entries[row]=mat2->num_entries[row];
  mat->m+=(mat2->num_entries[row]-temp_ct);

  mat2->neighbor_list[row]=temp_node;
  mat2->m+=(temp_ct-mat2->num_entries[row]);
  mat2->num_entries[row]=temp_ct;
  
}

void swap_both(LLA *mat, LLA *mat2, int vtx) {
  LL_node *temp_node = mat->neighbor_list[vtx];
  int temp_ct = mat->num_entries[vtx];

  
  mat->neighbor_list[vtx]=mat2->neighbor_list[vtx];
  mat->num_entries[vtx]=mat2->num_entries[vtx];
  mat->m+=(mat2->num_entries[vtx]-temp_ct);

  mat2->neighbor_list[vtx]=temp_node;
  mat2->m+=(temp_ct-mat2->num_entries[vtx]);
  mat2->num_entries[vtx]=temp_ct;

  for(int i=0; i < vtx; ++i) {
    int found1=0;
    int found2=0;
    LL_node *temp_node = mat->neighbor_list[i];
    LL_node *temp_node2 = mat2->neighbor_list[i];
    while(temp_node != NULL) {
      if(temp_node->idx == vtx) {
	found1=1;
	break;
      }
      else if(temp_node->idx > vtx) {
	found1=2;
	break;
      }
      if(temp_node->next_ptr == NULL) {
	found1=3;
	break;
      }
      temp_node=temp_node->next_ptr;
    }
    while(temp_node2 != NULL) {
      if(temp_node2->idx == vtx) {
	found2=1;
	break;
      }
      else if(temp_node2->idx > vtx) {
	found2=2;
	break;
      }
      if(temp_node2->next_ptr == NULL) {
	found2=3;
	break;
      }
      temp_node2=temp_node2->next_ptr;
    }

    if(found1==1 && found2==1) {
      std::swap(temp_node->weight,temp_node2->weight);
    }
    else if(found1==1) {
      if(temp_node->prev_ptr != NULL) {
	temp_node->prev_ptr->next_ptr=temp_node->next_ptr;
      }
      else {
	mat->neighbor_list[i]=temp_node->next_ptr;
      }
      if(temp_node->next_ptr != NULL) {
	temp_node->next_ptr->prev_ptr=temp_node->prev_ptr;
      }

      if(found2==2) {

	if(temp_node2->prev_ptr != NULL) {
	  temp_node2->prev_ptr->next_ptr=temp_node;
	  temp_node->prev_ptr=temp_node2->prev_ptr;
	  temp_node2->prev_ptr=temp_node;
	  temp_node->next_ptr=temp_node2;
	}
	else {
	  mat2->neighbor_list[i]=temp_node;
	  temp_node2->prev_ptr=temp_node;
	  temp_node->next_ptr=temp_node2;
	  temp_node->prev_ptr=NULL;
	}
      }
      else if(found2==3) {
      
	if(temp_node2->next_ptr != NULL) {
	  temp_node2->next_ptr->prev_ptr=temp_node;
	  temp_node->next_ptr=temp_node2->next_ptr;
	  temp_node2->next_ptr=temp_node;
	  temp_node->prev_ptr=temp_node2;
	}
	else {
	  temp_node2->next_ptr=temp_node;
	  temp_node->prev_ptr=temp_node2;
	  temp_node->next_ptr=NULL;
	}

      }
      else if(found2==0) {
	mat2->neighbor_list[i]=temp_node;
	temp_node->prev_ptr=NULL;
	temp_node->next_ptr=NULL;
      }
      mat2->m++;
      mat->m--;
      mat2->num_entries[i]++;
      mat->num_entries[i]--;
    }
   else if(found2==1) {
      if(temp_node2->prev_ptr != NULL) {
	temp_node2->prev_ptr->next_ptr=temp_node2->next_ptr;
      }
      else {
	mat2->neighbor_list[i]=temp_node2->next_ptr;
      }
      if(temp_node2->next_ptr != NULL) {
	temp_node2->next_ptr->prev_ptr=temp_node2->prev_ptr;
      }

      if(found1==2) {

	if(temp_node->prev_ptr != NULL) {
	  temp_node->prev_ptr->next_ptr=temp_node2;
	  temp_node2->prev_ptr=temp_node->prev_ptr;
	  temp_node->prev_ptr=temp_node2;
	  temp_node2->next_ptr=temp_node;
	}
	else {
	  mat->neighbor_list[i]=temp_node2;
	  temp_node->prev_ptr=temp_node2;
	  temp_node2->next_ptr=temp_node;
	  temp_node2->prev_ptr=NULL;
	}
      }
      else if(found1==3) {
      
	if(temp_node->next_ptr != NULL) {
	  temp_node->next_ptr->prev_ptr=temp_node2;
	  temp_node2->next_ptr=temp_node->next_ptr;
	  temp_node->next_ptr=temp_node2;
	  temp_node2->prev_ptr=temp_node;
	}
	else {
	  temp_node->next_ptr=temp_node2;
	  temp_node2->prev_ptr=temp_node;
	  temp_node2->next_ptr=NULL;
	}

      }
      else if(found1==0) {
	mat->neighbor_list[i]=temp_node2;
	temp_node2->prev_ptr=NULL;
	temp_node2->next_ptr=NULL;
      }
      mat->m++;
      mat2->m--;
      mat->num_entries[i]++;
      mat2->num_entries[i]--;
    } 
  }
}

LLA* intersection(LLA *mat, LLA*mat2) {
  LLA* mat_intersect = new LLA(mat->n);
  for(int i=0; i < mat->n; ++i) {
    LL_node *temp_node = mat->neighbor_list[i];
    LL_node *temp_node2 = mat2->neighbor_list[i];
    while(temp_node != NULL && temp_node2 !=NULL) {
      if(temp_node->idx == temp_node2->idx) {
        mat_add_edge(mat_intersect,i,temp_node->idx, std::min(temp_node->weight,temp_node2->weight));
        temp_node=temp_node->next_ptr;
        temp_node2=temp_node2->next_ptr;
      }
      else if(temp_node->idx < temp_node2->idx) {
        temp_node=temp_node->next_ptr;
      }
      else if(temp_node2->idx < temp_node->idx) {
        temp_node2=temp_node2->next_ptr;
      }
    }
  }

  return mat_intersect;
}

void find_in_neighbors(LLA *mat) {
  for(int i=0; i < mat->n; ++i) {
    LL_node *temp_node = mat->in_neighbor_list[i];
    LL_node *node;

    for(int j=0; j < mat->in_num_entries[i]; ++j) {
      node = temp_node;
      temp_node=temp_node->next_ptr;
      delete node;
    }
    mat->in_num_entries[i]=0;
    mat->in_neighbor_list[i]=NULL;
  }

  
  for(int i=0; i < mat->n; ++i) {
    LL_node *temp_node = mat->neighbor_list[i];
    while(temp_node != NULL) {
      if(mat_add_edge_in(mat, temp_node->idx, i, temp_node->weight)==0)
	std::cout << "complain" << std::endl;
      temp_node=temp_node->next_ptr;
    }
  }
}

bool connected(LLA *mat) {
  int *visited = new int[mat->n];
  for(int i=0; i < mat->n; ++i) {
    visited[i]=0;
  }
  std::queue<int> myqueue;
  myqueue.push(0);
  int count=1;
  visited[0]=1;
  while(!myqueue.empty()) {
    int current = myqueue.front();
    myqueue.pop();
    LL_node *temp_node=mat->neighbor_list[current];
    while(temp_node != NULL) {
      if(visited[temp_node->idx]==0) {
	visited[temp_node->idx]=1;
	myqueue.push(temp_node->idx);
	count++;
      }
      temp_node=temp_node->next_ptr;
    }
    temp_node=mat->in_neighbor_list[current];
    while(temp_node != NULL) {
      if(visited[temp_node->idx]==0) {
	visited[temp_node->idx]=1;
	myqueue.push(temp_node->idx);
	count++;
      }
      temp_node=temp_node->next_ptr;
    }
  }
  
  delete[] visited;

  
  return(count==mat->n);

}

void attach(LLA *mat, LLA *mat2) {
  int n1=mat->n;
  int n2=mat2->n;
  int newn = n1 +  n2 - 1;
  mat->n=newn;
  mat->m+=mat2->m;
  
  int *new_num_entries = new int[newn];
  for(int i=0; i < n1-1; ++i) {
    new_num_entries[i]=mat->num_entries[i];
  }
  new_num_entries[n1-1]=mat->num_entries[n1-1] + mat2->num_entries[0];
  for(int i=0; i < n2-1; ++i) {
    new_num_entries[i+n1] = mat2->num_entries[i+1];
  }

  LL_node **new_neighbor_list = new LL_node*[newn];
  for(int i=0; i < n1-1; ++i) {
    new_neighbor_list[i]=mat->neighbor_list[i];
  }
  for(int i=-1; i < n2-1; ++i) {
    new_neighbor_list[i+n1] = mat2->neighbor_list[i+1];
    LL_node *temp_node = new_neighbor_list[i+n1];
    while(temp_node != NULL) {
      temp_node->idx+=(n1-1);
      temp_node=temp_node->next_ptr;
    }
  }

  int *new_in_num_entries = new int[newn];
  for(int i=0; i < n1-1; ++i) {
    new_in_num_entries[i]=mat->in_num_entries[i];
  }
  new_in_num_entries[n1-1]=mat->in_num_entries[n1-1] + mat2->in_num_entries[0];
  for(int i=0; i < n2-1; ++i) {
    new_in_num_entries[i+n1] = mat2->in_num_entries[i+1];
  }

  
  LL_node **new_in_neighbor_list = new LL_node*[newn];
  for(int i=0; i < n1-1; ++i) {
    new_in_neighbor_list[i]=mat->in_neighbor_list[i];
  }
  for(int i=-1; i < n2-1; ++i) {
    new_in_neighbor_list[i+n1] = mat2->in_neighbor_list[i+1];
    LL_node *temp_node = new_in_neighbor_list[i+n1];
    while(temp_node != NULL) {
      temp_node->idx+=(n1-1);
      temp_node=temp_node->next_ptr;
    }
  }
  
  delete[] (mat->num_entries);
  delete[] (mat->neighbor_list);
  delete[] (mat->in_neighbor_list);
  delete[] (mat->in_num_entries);
  delete[] (mat2->num_entries);
  delete[] (mat2->neighbor_list);
  delete[] (mat2->in_neighbor_list);
  delete[] (mat2->in_num_entries);
  mat->num_entries = new_num_entries;
  mat->neighbor_list = new_neighbor_list;
  mat->in_neighbor_list = new_in_neighbor_list;
  mat->in_num_entries = new_in_num_entries;
}

void MMA_write(LLA *mat, std::string filename) {
  std::ofstream mat_file(filename.c_str());
  //int precdigits=6;
  mat_file << "%%MatrixMarket matrix coordinate real symmetric" << std::endl;
  mat_file << "%%" << std::endl;
  mat_file << mat->n << ' ' << mat->n  << ' ' << (mat->m) << std::endl;
  for(int i=0; i < mat->n; ++i) {
    LL_node *temp_node = mat->neighbor_list[i];
    for(int j=0; j < mat->num_entries[i]; ++j) {
      mat_file << i+1 << " " << temp_node->idx+1 << " "  << temp_node->weight << std::endl;
      //mat_file << temp_node->idx+1 << " " << i+1 << " " << std::setprecision(precdigits+1) << temp_node->weight << std::endl;
      temp_node=temp_node->next_ptr;
    }
  }
}
void MML_write(LLA *mat, std::string filename) {
  std::ofstream mat_file(filename.c_str());
  
  mat_file << "%%MatrixMarket matrix coordinate real general" << std::endl;
  mat_file << "%%" << std::endl;
  mat_file << mat->n << ' ' << mat->n  << ' ' << 2*mat->m + mat->n<< std::endl;
  double *diagonal = new double[mat->n];
  for(int i=0; i < mat->n; ++i) {
    LL_node *temp_node = mat->neighbor_list[i];
    for(int j=0; j < mat->num_entries[i]; ++j) {
      mat_file << i+1 << " " << temp_node->idx+1 << " " << -temp_node->weight << std::endl;
      mat_file << temp_node->idx+1 << " " << i+1 << " " << -temp_node->weight << std::endl;
      diagonal[i]+=temp_node->weight;
      diagonal[temp_node->idx]+=temp_node->weight;
      temp_node=temp_node->next_ptr;
    }
  }
  for(int i=0; i < mat->n;++i) {
    mat_file << i+1 << " " << i+1 << " " << diagonal[i] << std::endl;
  }
  delete[] diagonal;
}


LLA* MMA_read(std::string filename) {
  int n,m,test;
  bool symmetric=false;
  
  std::ifstream mat_file(filename.c_str());

  std::string header;
  getline(mat_file,header);

  if(header.find("symmetric") != std::string::npos) {
    symmetric=true;
  }

  if(!symmetric) {
    std::cout << "only works with symmetric matrices" << std::endl;
    return NULL;
  }

  while(mat_file.peek() == '%') {
    mat_file.ignore(2048, '\n');
  }

  mat_file >> n >> test >> m;

  if(n != test) {
    std::cout << "only works with square matrices" << std::endl;
    return NULL;
  }

  LLA *mat = new LLA(n);

  std::string linetext;
  getline(mat_file,linetext);

  for(int line=0; line<m; ++line) {
    int i,j;
    double w;
    getline(mat_file,linetext);
    std::string delimeter=" ";
    std::string token=linetext.substr(0,linetext.find(delimeter));
    i=atoi(token.c_str());
    linetext.erase(0,linetext.find(delimeter)+delimeter.length());
    token=linetext.substr(0,linetext.find(delimeter));
    j=atoi(token.c_str());
    linetext.erase(0,linetext.find(delimeter)+delimeter.length());
    token=linetext.substr(0,linetext.find(delimeter));
    w=atof(token.c_str());

    mat_add_edge(mat, i-1, j-1, w);
    
  }

  return mat;
}

void print_matrix(LLA *mat) {
  std::cout << "Vertices: " << mat->n << ", Edges: " << mat->m << std::endl;
  for(int i=0; i < mat->n; ++i) {
    LL_node *temp_node = mat->neighbor_list[i];
    for(int j=0; j < mat->num_entries[i]; ++j) {
      std::cout << "Edge from " << i << " to " << temp_node->idx << " with weight " << temp_node->weight << std::endl;
      temp_node=temp_node->next_ptr;
    }
  }
}
