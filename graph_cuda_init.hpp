

void graph_c::init(){
  
  unsigned long long dtime; 
  dtime = dtime_usec(0);
  
  A->Ac = this;
  
  Vn = A->Vn;
  En = A->En;
   
  Edges_D   = (double*)malloc(sizeof(double) * En);
  Edges_d   = (double*)malloc(sizeof(double) * En);
    
  Edges_v1  = (int*)malloc(sizeof(int) * En);
  Edges_v2  = (int*)malloc(sizeof(int) * En);
  
  
  Nodes       = (int*)malloc(sizeof(int) * Vn);
  NodesStat   = (int*)malloc(sizeof(int) * Vn);
  Nodes_F     = (double*)malloc(sizeof(double) * Vn);
  
  
  for(int i=0; i < Vn; ++i ){NodesStat[i]=0;Nodes_F[i]= 0;}
  for(int i=0; i < En; ++i ){
   
   Edges_D[i] = A->E[i]->D;
   Edges_d[i] = A->E[i]->d;
   
   int v1,v2;
   v1 = A->E[i]->v1->id;
   if(v1 < A->E[i]->v2->id){v2 = A->E[i]->v2->id;}else{v1 = A->E[i]->v2->id;v2 = A->E[i]->v1->id;}
  
   Edges_v1[i] = v1;
   Edges_v2[i] = v2;
   NodesStat[v1]++;
   NodesStat[v2]++;  
  }
  
  
  std::unordered_map<int,int> node_count;
  int n1 = 0;
  for(int i=0; i < Vn; ++i ){
   
   node_count[i]=0;
   Nodes[i] =  n1;
   n1 +=  NodesStat[i];
  }
  

  i_v1      = (int*)malloc(sizeof(int) * En);
  i_v2      = (int*)malloc(sizeof(int) * En);
  
  for(int i=0; i < En; ++i ){
   int id1 = Edges_v1[i];
   int id2 = Edges_v2[i];
  
   i_v1[i] =  node_count[id1];
   i_v2[i] =  node_count[id2];
   node_count[id1]++;
   node_count[id2]++;
  }
  
   
  dtime = dtime_usec(dtime);
  std::cout << "graph_c init() time: " << dtime/(double)USECPSEC  << std::endl;
  
  
}

//#################################################
//#################################################
//#################################################

template <typename T>
__global__  void setvalue(T    *  v, T val, int n )
{
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid < n){v[tid]    =   val;}
}


__global__ void cuda_Edges_d_by_Nodes_fill(double *cuda_Edges_d, int *cuda_Edges_v1, int *cuda_Edges_v2 , int *cuda_i_v1, int *cuda_i_v2, double ** cuda_Edges_d_by_Nodes, int * cuda_Nodes, int * cuda_NodesStat ,int n)
{
     int tid = blockIdx.x * blockDim.x + threadIdx.x;
     if (tid < n){
       int v1 = cuda_Edges_v1[tid];
       int v2 = cuda_Edges_v2[tid];
       int i1 = cuda_Nodes[v1] + cuda_i_v1[tid];
       int i2 = cuda_Nodes[v2] + cuda_i_v2[tid];
       cuda_Edges_d_by_Nodes[i1] = &cuda_Edges_d[tid];
       cuda_Edges_d_by_Nodes[i2] = &cuda_Edges_d[tid];  
       
     }
}





void graph_c::init_cuda(){
   
  unsigned long long dtime; 
  dtime = dtime_usec(0);
 
    
  grid_size_Vn  = ((  Vn + block_size)   / block_size);
  grid_size_En  = ((  En + block_size)   / block_size);
  
  //std::cout << "init_cuda() " << Vn << " " << grid_size_Vn << "\n"; 
  // Allocate device memory 
  cudaMalloc((void**)&cuda_Edges_D, sizeof(double) * En);
  cudaMalloc((void**)&cuda_Edges_d, sizeof(double) * En);
  cudaMalloc((void**)&cuda_Edges_v1, sizeof(int) * En);
  cudaMalloc((void**)&cuda_Edges_v2, sizeof(int) * En);
  
  cudaMalloc((void**)&cuda_Edges_d_by_Nodes, sizeof(double*)  * 2 * En);
  cudaMalloc((void**)&cuda_Nodes,  sizeof(int) * Vn);
  cudaMalloc((void**)&cuda_NodesStat,  sizeof(int) * Vn);
   
  cudaMalloc((void**)&cuda_Nodes_grad, sizeof(double) * Vn);
  cudaMalloc((void**)&cuda_Nodes_F,    sizeof(double) * Vn);

  cudaMalloc((void**)&cuda_EF, sizeof(double*)  * 1);
  cudaMalloc((void**)&cuda_EF_Nodes,      sizeof(int) * Vn);
  cudaMalloc((void**)&cuda_EF_NodesStat,  sizeof(int) * Vn);
  
  cudaMalloc((void**)&cuda_block_id, sizeof(int) * Vn);
  
  cudaMalloc((void**)&cuda_COMB_NodesStat,  sizeof(int) * Vn);
  cudaMalloc((void**)&cuda_Nodes_A,  sizeof(int) * Vn);
  cudaMalloc((void**)&cuda_Nodes_B,  sizeof(int) * Vn);
  
  cudaMalloc((void**)&cuda_comb_sizeV,  sizeof(int) * Vn);
  cudaMalloc((void**)&cuda_N,   sizeof(int) );
  cudaMalloc((void**)&cuda_NN,  sizeof(int) );
  
  cudaMalloc((void**)&cuda_comb_sizeE,  sizeof(int) * Vn);
  cudaMalloc((void**)&cuda_E,   sizeof(int) );
  cudaMalloc((void**)&cuda_EE,  sizeof(int) );
  
  cudaMalloc((void**)&cuda_stagnation_control,  sizeof(int) * Vn);
  
  cudaCheckErrors("cudamalloc fail");
  //checkGpuMem();
  
  
  // Transfer data from host to device memory
  cudaMemcpy(cuda_Edges_D, Edges_D, sizeof(double) * En, cudaMemcpyHostToDevice);
  cudaMemcpy(cuda_Edges_d, Edges_d, sizeof(double) * En, cudaMemcpyHostToDevice);
    
  cudaMemcpy(cuda_Edges_v1, Edges_v1, sizeof(int) * En, cudaMemcpyHostToDevice);
  cudaMemcpy(cuda_Edges_v2, Edges_v2, sizeof(int) * En, cudaMemcpyHostToDevice);
  
  cudaMemcpy(cuda_Nodes, Nodes, sizeof(int) * Vn, cudaMemcpyHostToDevice);
  cudaMemcpy(cuda_NodesStat, NodesStat, sizeof(int) * Vn, cudaMemcpyHostToDevice);
  
  cudaCheckErrors("cudaMemcpy fail last");
  
  
  /* *************************************************************************************************************************** */
  int *cuda_i_v1, *cuda_i_v2;
  cudaMalloc((void**)&cuda_i_v1, sizeof(int) * En);
  cudaMalloc((void**)&cuda_i_v2, sizeof(int) * En);
  
  cudaMemcpy(cuda_i_v1, i_v1, sizeof(int) * En, cudaMemcpyHostToDevice);
  cudaMemcpy(cuda_i_v2, i_v2, sizeof(int) * En, cudaMemcpyHostToDevice);
  
  cuda_Edges_d_by_Nodes_fill<<<grid_size_En,block_size>>>(cuda_Edges_d, cuda_Edges_v1, cuda_Edges_v2 , cuda_i_v1,cuda_i_v2, cuda_Edges_d_by_Nodes, cuda_Nodes, cuda_NodesStat ,En);
  cudaDeviceSynchronize();
  
  cudaFree(cuda_i_v1);
  cudaFree(cuda_i_v2);
  
  free(i_v1);
  free(i_v2);
  
  cudaCheckErrors("cudaMemcpy fail cuda_i_v2");
  /* *************************************************************************************************************************** */
    
  setvalue<<<grid_size_Vn,block_size>>>(cuda_Nodes_grad, 0.00, Vn );
  setvalue<<<grid_size_Vn,block_size>>>(cuda_Nodes_F, 0.00 , Vn );
  setvalue<<<grid_size_Vn,block_size>>>(cuda_block_id, 0 , Vn );
  setvalue<<<grid_size_Vn,block_size>>>(cuda_stagnation_control, 0 , Vn );
  
  
  cudaCheckErrors("cudaMemcpy fail setvalue");
  
  dtime = dtime_usec(dtime);
  std::cout << "cuda init() time: " << dtime/(double)USECPSEC  << std::endl;
  
}  


graph_c::~graph_c(){

free(Edges_d);
free(Edges_D);

cudaFree(Edges_d);
cudaFree(Edges_D);
cudaFree(Edges_v1);
cudaFree(Edges_v2);
cudaFree(cuda_Edges_d_by_Nodes);
cudaFree(cuda_Nodes);
cudaFree(cuda_NodesStat);
    
cudaFree(cuda_EF_Nodes); 
cudaFree(cuda_EF_NodesStat);
cudaFree(cuda_EF);
    
cudaFree(cuda_stagnation_control);
    
cudaFree(cuda_Nodes_A);
cudaFree(cuda_Nodes_B);
cudaFree(cuda_block_id);
    
cudaFree(cuda_COMB_NodesStat);


cudaFree(cuda_Nodes_grad);
cudaFree(cuda_Nodes_F);
 
 
 
cudaFree(cuda_N);
cudaFree(cuda_NN);

cudaFree(comb_id);

cudaFree(cuda_comb_sizeV) ;
cudaFree(comb_start);
cudaFree(comb_stat);
cudaFree(comb_V);
    

cudaFree(cuda_comb_sizeE) ;
cudaFree(comb_start_E);
cudaFree(comb_stat_E);
cudaFree(comb_E);
 
cudaFree(comb_id_IRREGULAR);
   
    
}










