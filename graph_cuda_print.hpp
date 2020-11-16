
template <typename T>
T * cuda_copy_vector(const T * cuda_v , int size)
{
  T*host_v = (T*)malloc(sizeof(T) * (size));
  cudaMemcpy(host_v, cuda_v, sizeof(T)  * (size), cudaMemcpyDeviceToHost);
 return host_v ;
}


double ** cuda_copy_vector_ptr(double ** cuda_v , int size)
{
  double**host_v = (double**)malloc(sizeof(double*) * (size));
  cudaMemcpy(host_v, cuda_v, sizeof(double*)  * (size), cudaMemcpyDeviceToHost);
 return host_v ;
}


template <typename T>
void cuda_print_vector_none0(const T * cuda_v , int size)
{
  T *host_v = (T*)malloc(sizeof(T) * size);
  cudaMemcpy(host_v, cuda_v, sizeof(T)  * size, cudaMemcpyDeviceToHost);
  std::cout << "cuda_print_vector:\n";
  int count = 0;
  for(int i=0; i < size; ++i ){
   if(host_v[i] > 0){
    std::cout << "i:" << i << " " << count << " " << host_v[i]  << "\n";
    count++;
   }
  }
  free(host_v);
}

template <typename T>
void cuda_print_vector(int i1,int i2,const T * cuda_v)
{
  T *host_v = (T*)malloc(sizeof(T) * (i2-i1));
  cudaMemcpy(host_v, cuda_v+i1, sizeof(T)  * (i2-i1), cudaMemcpyDeviceToHost);
  std::cout << "cuda_print_vector:\n";
  for(int i=0; i < (i2-i1); ++i ){std::cout << "i:" << i1 + i << " " << host_v[i]  << "\n";}
  free(host_v);
}

template <typename T>
void cuda_print_vector_AB(int size,const T * cuda_A,const T * cuda_B)
{
  T *host_B = (T*)malloc(sizeof(T) * (size));
  T *host_A = (T*)malloc(sizeof(T) * (size));
  cudaMemcpy(host_A, cuda_A+size, sizeof(T)  * (size), cudaMemcpyDeviceToHost);
  cudaMemcpy(host_B, cuda_B+size, sizeof(T)  * (size), cudaMemcpyDeviceToHost);
  std::cout << "cuda_print_vector:\n";
  for(int i=0; i < size; ++i ){
  if(host_A[i] < 1 ){ std::cout << "A i:" << i << " " << host_A[i] << " " << host_B[i]  << "\n";}
  if(host_B[i] > 2 ){ std::cout << "B i:" << i << " " << host_A[i] << " " << host_B[i]  << "\n";}
  }
  
  free(host_A);
  free(host_B);
}


__global__ void cuda_get_node(
int k,int topK , 
double ** cuda_Edges_d_by_Nodes, int * cuda_Nodes, int * cuda_NodesStat, 
double * cuda_Edges_d, double * cuda_Edges_D,
int *cuda_Edges_v1, int *cuda_Edges_v2 ,

double * cuda_node_d, double * cuda_node_D, int *cuda_node_v1, int *cuda_node_v2, int *cuda_edge_k

 )
{
     
     int tid = blockIdx.x * blockDim.x + threadIdx.x;
     if (tid < 1){
     
       double ** ptr = cuda_Edges_d_by_Nodes + cuda_Nodes[k]; 
     
       for(int i=0; i < topK; ++i ){
        double * e = ptr[i];
        int indE  = e - cuda_Edges_d;
        
        cuda_node_d[i]  = cuda_Edges_d[indE];
        cuda_node_D[i]  = cuda_Edges_D[indE];
        cuda_node_v1[i] = cuda_Edges_v1[indE];
        cuda_node_v2[i] = cuda_Edges_v2[indE];
        cuda_edge_k[i]  = indE;
        
       }
     }
     
}

void cuda_print_node
(
int k,int topK , 
double ** cuda_Edges_d_by_Nodes, int * cuda_Nodes, int * cuda_NodesStat, 
double * cuda_Edges_d, double * cuda_Edges_D,
int *cuda_Edges_v1, int *cuda_Edges_v2
)
{
  double * cuda_node_d;
  cudaMalloc((void**)&cuda_node_d, sizeof(double) * topK);
  double * cuda_node_D;
  cudaMalloc((void**)&cuda_node_D, sizeof(double) * topK);
  int *cuda_node_v1;
  cudaMalloc((void**)&cuda_node_v1, sizeof(int) * topK);
  int *cuda_node_v2;
  cudaMalloc((void**)&cuda_node_v2, sizeof(int) * topK);
  int *cuda_edge_k;
  cudaMalloc((void**)&cuda_edge_k, sizeof(int) * topK);
   
  //###############
  cuda_get_node<<<1,1>>>(k,topK ,cuda_Edges_d_by_Nodes, cuda_Nodes, cuda_NodesStat ,  cuda_Edges_d, cuda_Edges_D, cuda_Edges_v1, cuda_Edges_v2, cuda_node_d,  cuda_node_D, cuda_node_v1 , cuda_node_v2 , cuda_edge_k  );
                                                    
  //############### 
                                                      
  double * node_d = (double*)malloc(sizeof(double) * topK);
  double * node_D = (double*)malloc(sizeof(double) * topK);
  int *node_v1 = (int*)malloc(sizeof(int) * topK);
  int *node_v2 = (int*)malloc(sizeof(int) * topK);
  int *edge_k  = (int*)malloc(sizeof(int) * topK);
  
  cudaMemcpy(node_d, cuda_node_d, sizeof(double)  * topK, cudaMemcpyDeviceToHost);
  cudaMemcpy(node_D, cuda_node_D, sizeof(double)  * topK, cudaMemcpyDeviceToHost);
  cudaMemcpy(node_v1, cuda_node_v1, sizeof(int)  * topK, cudaMemcpyDeviceToHost);
  cudaMemcpy(node_v2, cuda_node_v2, sizeof(int)  * topK, cudaMemcpyDeviceToHost);
  cudaMemcpy(edge_k, cuda_edge_k, sizeof(int)  * topK, cudaMemcpyDeviceToHost);
  
  std::cout << "print_nodeZ Node:" << k << "\n";
  
  for(int i=0; i < topK; ++i ){
   //if(node_d[i] < EF_porog){
    std::cout << "Edges:" << i << " " << node_d[i]  << " " << node_D[i] << " " << node_v1[i] << " " << node_v2[i] << " " << edge_k[i] << "\n";
   //}
  }
  
  cudaFree(cuda_node_d);
  cudaFree(cuda_node_D);
  cudaFree(cuda_node_v1);
  cudaFree(cuda_node_v2);
  cudaFree(cuda_edge_k);
   free(node_d);
   free(node_D);
   free(node_v1);
   free(node_v2);
   free(edge_k);
}



void cuda_print_node_porog
(
double  porog,
int k,int topK , 
double ** cuda_Edges_d_by_Nodes, int * cuda_Nodes, int * cuda_NodesStat, 
double * cuda_Edges_d, double * cuda_Edges_D,
int *cuda_Edges_v1, int *cuda_Edges_v2,
double *  cuda_Nodes_grad
)
{
  double * cuda_node_d;
  cudaMalloc((void**)&cuda_node_d, sizeof(double) * topK);
  double * cuda_node_D;
  cudaMalloc((void**)&cuda_node_D, sizeof(double) * topK);
  int *cuda_node_v1;
  cudaMalloc((void**)&cuda_node_v1, sizeof(int) * topK);
  int *cuda_node_v2;
  cudaMalloc((void**)&cuda_node_v2, sizeof(int) * topK);
  int *cuda_edge_k;
  cudaMalloc((void**)&cuda_edge_k, sizeof(int) * topK);
   
  //###############
  cuda_get_node<<<1,1>>>(k,topK ,cuda_Edges_d_by_Nodes, cuda_Nodes, cuda_NodesStat ,  cuda_Edges_d, cuda_Edges_D, cuda_Edges_v1, cuda_Edges_v2, cuda_node_d,  cuda_node_D, cuda_node_v1 , cuda_node_v2 , cuda_edge_k  );
                                                    
  //############### 
                                                      
  double * node_d = (double*)malloc(sizeof(double) * topK);
  double * node_D = (double*)malloc(sizeof(double) * topK);
  int *node_v1 = (int*)malloc(sizeof(int) * topK);
  int *node_v2 = (int*)malloc(sizeof(int) * topK);
  int *edge_k  = (int*)malloc(sizeof(int) * topK);
  
  cudaMemcpy(node_d, cuda_node_d, sizeof(double)  * topK, cudaMemcpyDeviceToHost);
  cudaMemcpy(node_D, cuda_node_D, sizeof(double)  * topK, cudaMemcpyDeviceToHost);
  cudaMemcpy(node_v1, cuda_node_v1, sizeof(int)  * topK, cudaMemcpyDeviceToHost);
  cudaMemcpy(node_v2, cuda_node_v2, sizeof(int)  * topK, cudaMemcpyDeviceToHost);
  cudaMemcpy(edge_k, cuda_edge_k, sizeof(int)  * topK, cudaMemcpyDeviceToHost);
  
  
  double max_grad = 0;
  cudaMemcpy(&max_grad, cuda_Nodes_grad+k, sizeof(double) , cudaMemcpyDeviceToHost);
  std::cout << "print_nodeZ Node:" << k << " grad " << max_grad << "\n";
  
  for(int i=0; i < topK; ++i ){
   if(node_d[i] < porog){
   
    std::cout << "Edges:" << i << " " << node_d[i]  << " " << node_D[i] << " " << node_v1[i] << " " << node_v2[i] << " " << edge_k[i] << "\n";
   }
  }
  
  cudaFree(cuda_node_d);
  cudaFree(cuda_node_D);
  cudaFree(cuda_node_v1);
  cudaFree(cuda_node_v2);
  cudaFree(cuda_edge_k);
   free(node_d);
   free(node_D);
   free(node_v1);
   free(node_v2);
   free(edge_k);
}






