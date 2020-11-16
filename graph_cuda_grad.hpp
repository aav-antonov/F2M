__global__  void cuda_recompute_Edges(double * cuda_Edges_d, double * cuda_Edges_D, int * cuda_Edges_v1,  int * cuda_Edges_v2, double * cuda_Nodes_F, int n ){
  
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  
  if (tid < n){
  
    int v1     =  cuda_Edges_v1[tid];
    int v2     =  cuda_Edges_v2[tid];
    cuda_Edges_d[tid] = cuda_Edges_D[tid] + 0.5*cuda_Nodes_F[v1] + 0.5*cuda_Nodes_F[v2] ;
    
  }
}





__global__  void cuda_get_Nodes_StatAB( double comb_porog, double ** cuda_Edges_d_by_Nodes, int    *  cuda_Nodes, int    *  cuda_NodesStat, int * cuda_Nodes_statA, int * cuda_Nodes_statB, int n )
{
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  
  if (tid < n){
  double ** ptr = cuda_Edges_d_by_Nodes + cuda_Nodes[tid];
  int A = 0, B = 0;
  
   for(int i=0; i < cuda_NodesStat[tid]; ++i ){
    double v = *ptr[i];
        
     if(v < - comb_porog ){A++;}
else if(v <   comb_porog ){B++;} 
    
    
   }
   
   cuda_Nodes_statA[tid] = A;
   cuda_Nodes_statB[tid] = B;
   
  }
}






__global__  void cuda_update_Nodes_F_BLOCK(double * cuda_Nodes_grad, double * cuda_Nodes_F, int* cuda_block_id,int n )
{
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  
  if (tid < n){  if(cuda_block_id[tid] == 0){ cuda_Nodes_F[tid] += cuda_Nodes_grad[tid]; }}
}


__global__  void cuda_update_Nodes_F(double * cuda_Nodes_grad, double * cuda_Nodes_F, int n )
{
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  
  if (tid < n){
   cuda_Nodes_F[tid]  +=  cuda_Nodes_grad[tid];
  }
}



__global__  void cuda_get_Nodes_grad_R( double **cuda_Edges_d_by_Nodes, int *cuda_Nodes, int *cuda_NodesStat, double *cuda_Nodes_grad, int n )
{
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  
  if (tid < n){
  double ** ptr = cuda_Edges_d_by_Nodes + cuda_Nodes[tid];
  double v0 = *ptr[0] , v1 = 1e6, v2 = 1e6;
    
   for(int i=1; i < cuda_NodesStat[tid]; ++i ){
    double v = *ptr[i];
    
    if(v < v0){v2 = v1; v1 = v0; v0 = v;continue;}
    if(v < v1){v2 = v1;          v1 = v;continue;}
    if(v < v2){                  v2 = v;continue;}
   }

   cuda_Nodes_grad[tid]  =  -(v1 + v2) / 2;
 
  }
}


__global__  void get_cuda_EF_nodes_NORM(int  *  cuda_EF_NodesStat, int  *  cuda_EF_Nodes , int n ){
   
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  
  if (tid < n){
   cuda_EF_Nodes[tid] -= cuda_EF_NodesStat[tid];
  }
}

__global__  void get_cuda_EF(
double ** cuda_Edges_d_by_Nodes, int    *    cuda_Nodes,   
int * cuda_EF_NodesStat,
double ** cuda_EF,               int    * cuda_EF_Nodes ,    
int n 
)
{
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  
  if (tid < n){
   for(int i=0; i < cuda_EF_NodesStat[tid]; ++i ){
     cuda_EF[cuda_EF_Nodes[tid] + i] = cuda_Edges_d_by_Nodes[cuda_Nodes[tid]+i]; 
    }
  }
  
}












