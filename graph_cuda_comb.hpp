
int graph::host_combs_check(){
  
  unsigned long long dtimeBG; 
  dtimeBG = dtime_usec(0);
  
  EF_thread(MAX_THREADS , 1);
  V_thread(MAX_THREADS , 1 , error_B  );
  V_thread(MAX_THREADS , 2 , error_B  );
  V_thread(MAX_THREADS , 3 , error_B * error_B_comb  );
  
  COMB_thread(MAX_THREADS );
  COMB_solve_thread(MAX_THREADS );
  get_comb_IRREGULAR_A();
  
  if(comb_IRREGULAR_A.size() > 0){
   
   COMB_solve_AxB_thread(MAX_THREADS );
   get_comb_IRREGULAR_B();
   
 }
   
   dtimeBG = dtime_usec(dtimeBG);
   std::cout << "host_combs_check() time: " << dtimeBG/(double)USECPSEC  << std::endl;
  
   
return comb_IRREGULAR_B.size();   
}



void graph_c::get_combs_cuda(int id){
  
  
  std::cout << "get_combs_cuda() error_A:" << error_A << " \n";
  
  get_combs(error_B * error_B_comb);
  
  cudaCheckErrors("cudaMemcpy fail get_combs_cuda() 1");
  
  if(id > 0){
  
   free(A->comb_V_start);
   free(A->comb_V_stat);
   free(A->comb_V);
   
   free(A->comb_E_start);
   free(A->comb_E_stat);
   free(A->comb_E);
   
  }
  
  copyBackToHost();
  
  cudaCheckErrors("cudaMemcpy fail get_combs_cuda() 2");
  
}  


void graph_c::copyBackToHost(){ 
   
   cudaMemcpy(Nodes_F,   cuda_Nodes_F, sizeof(double)*Vn , cudaMemcpyDeviceToHost);
   A->F = Nodes_F;
   
   A->combs_size = N;
      
   A->comb_V_start  = (int*)malloc(sizeof(int) * N);
   A->comb_V_stat   = (int*)malloc(sizeof(int) * N);
   A->comb_V        = (int*)malloc(sizeof(int) * NN);
   std::cout << "graph_c::copyBackToHost()\n";
   cudaMemcpy(A->comb_V_start,   comb_start, sizeof(int)*N  , cudaMemcpyDeviceToHost);
   cudaMemcpy(A->comb_V_stat,     comb_stat, sizeof(int)*N  , cudaMemcpyDeviceToHost);
   cudaMemcpy(A->comb_V,            comb_V,  sizeof(int)*NN , cudaMemcpyDeviceToHost);
   
   
   A->comb_E_start  = (int*)malloc(sizeof(int) * E);
   A->comb_E_stat   = (int*)malloc(sizeof(int) * E);
   A->comb_E        = (int*)malloc(sizeof(int) * EE);
      
   cudaMemcpy(A->comb_E_start,   comb_start_E, sizeof(int)*E  , cudaMemcpyDeviceToHost);
   cudaMemcpy(A->comb_E_stat,     comb_stat_E, sizeof(int)*E  , cudaMemcpyDeviceToHost);
   cudaMemcpy(A->comb_E,               comb_E, sizeof(int)*EE , cudaMemcpyDeviceToHost);
   
   std::cout << "graph_c::copyBackToHost()\n";
   
}   







///#################################################################################
///#################################################################################
///#################################################################################
///#################################################################################
///#################################################################################

__global__  void get_cuda_comb_v2(int  *  cuda_comb_size, int  * comb_id , int * comb_start, int * comb_stat, int n )
{
 int tid = blockIdx.x * blockDim.x + threadIdx.x;
  
  if (tid <1){
   int k = 0;
   for(int i=0; i < n; ++i ){
     if(cuda_comb_size[i]> 0){
     comb_id[k] = i;
     comb_stat[k] = cuda_comb_size[i];
     if(k>0){comb_start[k] = comb_start[k-1] + comb_stat[k-1];}else{comb_start[k] = 0;}  
     k++; 
     }
   }
  
  }
}

__global__  void get_cuda_comb_v1(int  *  cuda_comb_size, int  * comb_N , int  * comb_NN , int n )
{
 int tid = blockIdx.x * blockDim.x + threadIdx.x;
  
  if (tid <1){
   *comb_N  = 0;
   *comb_NN = 0;
   for(int i=0; i < n; ++i ){
     if(cuda_comb_size[i]> 0){*comb_N += 1;*comb_NN += cuda_comb_size[i]; }
   }
   
   
  }
}







__global__  void get_cuda_comb_fill(
double comb_porog, int *comb_id, int *comb_start, int *comb_stat, int *comb_V, int *comb_start_E, int *comb_stat_E, int *comb_E,
double ** cuda_Edges_d_by_Nodes, int *cuda_Nodes, int *cuda_NodesStat, double * cuda_Edges_d, int * cuda_Edges_v1, int * cuda_Edges_v2,
int n 
)
{
  
  int bid = blockIdx.x * blockDim.x + threadIdx.x;
  int tid = comb_id[bid];
  
  if(bid < n){
       
   int filter_edge[comb_size_MAX];int size_e =0;
   int filter_node[comb_size_MAX];int size_f =0; 
   int list_A[comb_size_list];int size_A = 0;
  
  filter_node[0] = tid;size_f++;
  list_A[0]      = tid;size_A++;

 while(size_A > 0){
   int list_B[comb_size_list];int size_B = 0;
  
  for(int iA=0; iA < size_A; ++iA ){
   int lA = list_A[iA];
   double ** ptr = cuda_Edges_d_by_Nodes + cuda_Nodes[lA];
  
   for(int i=0; i < cuda_NodesStat[lA]; ++i ){
   
     if( (*ptr[i] >= -comb_porog)  ){
      if(  (*ptr[i] <= comb_porog) ){
       
       
       double * e = ptr[i];
       int indE  = e - cuda_Edges_d;
        
       int v1   =  cuda_Edges_v1[indE];
       int v2   =  cuda_Edges_v2[indE];
        
       int de = array_filtr(filter_edge,  size_e , indE );
       if(de == -1){  filter_edge[size_e] = indE;size_e++;} 
        
       int v = v1;
       if(v1 == lA){v =v2;}
       
       int da = array_filtr(filter_node,  size_f , v  );
       if(da == -1){ list_B[size_B] = v; size_B++;filter_node[size_f] = v;size_f++;}
       
     }}
    }
   }
   
    for(int iB=0; iB < size_B; ++iB ){list_A[iB] = list_B[iB];}
    size_A = size_B; size_B = 0;
   
  }
   
   for(int i=0; i < size_f; ++i ){comb_V[comb_start[bid] + i] = filter_node[i];}
   for(int i=0; i < size_e; ++i ){comb_E[comb_start_E[bid] + i] = filter_edge[i];}
   
 }
}

__global__  void get_cuda_comb_size(  
double comb_porog, 
double ** cuda_Edges_d_by_Nodes, int *cuda_Nodes, int *cuda_NodesStat, 
int *cuda_comb_size_v,int *cuda_comb_size_e,
double * cuda_Edges_d, int * cuda_Edges_v1, int * cuda_Edges_v2,
int n )
{
  
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  
  if(tid < n){
  int ind_break = 0;
  
   int filter_edge[comb_size_MAX];int size_e =0;
   int filter_node[comb_size_MAX];int size_f =0;
   int list_A[comb_size_list];int size_A = 0;
  
  filter_node[0] = tid;size_f++;
  list_A[0]      = tid;size_A++;

 while(size_A > 0){
   int list_B[comb_size_list];int size_B = 0;
 
  for(int iA=0; iA < size_A; ++iA ){
   int lA = list_A[iA];
   double ** ptr = cuda_Edges_d_by_Nodes + cuda_Nodes[lA];
  
   for(int i=0; i < cuda_NodesStat[lA]; ++i ){
   
    if( (*ptr[i] >= -comb_porog)  ){
     if(  (*ptr[i] <= comb_porog) ){
       
       ind_break = -1;;
       
       double * e = ptr[i];
       int indE  = e - cuda_Edges_d;
        
       int v1   =  cuda_Edges_v1[indE];
       int v2   =  cuda_Edges_v2[indE];
       
       if(v1 < tid){ind_break = 1;break;}
       if(v2 < tid){ind_break = 1;break;}
       
       int de = array_filtr(filter_edge,  size_e , indE );
       if(de == -1){  filter_edge[size_e] = indE;size_e++;}
       
              
       int v = v1;
       if(v1 == lA){v =v2;}
       
       int da = array_filtr(filter_node,  size_f , v  );
       if(da == -1){ list_B[size_B] = v; size_B++;filter_node[size_f] = v;size_f++;}
       
    }}
    
    }
    
    if(ind_break == 1){break;}
    
   }
   
   if(ind_break == 1){break;}
   
    for(int iB=0; iB < size_B; ++iB ){list_A[iB] = list_B[iB];}
    
    size_A = size_B; size_B = 0;
   
   }
   
  if(ind_break >= 0){cuda_comb_size_v[tid] = 0;cuda_comb_size_e[tid] = 0;}
  if(ind_break == -1){cuda_comb_size_v[tid] = size_f;cuda_comb_size_e[tid] = size_e; }
 }
}




/* ********************************************************************* */
/* ********************************************************************* */
/* ********************************************************************* */
/* ********************************************************************* */

void graph_c::free_combs_cuda(){

  cudaFree(comb_id);
  cudaFree(comb_start);
  cudaFree(comb_stat);
  cudaFree(comb_V);
  cudaFree(comb_start_E);
  cudaFree(comb_stat_E);
  cudaFree(comb_E);
  
  cudaCheckErrors("cudaMemcpy fail get_combs_cuda() B");
}  

void graph_c::get_combs(double comb_porog){  
  
  unsigned long long dtimeBG; 
  dtimeBG = dtime_usec(0); 
 
  cuda_partition_v3<<<grid_size_Vn,block_size>>>(cuda_Edges_d_by_Nodes, cuda_Nodes, cuda_EF_NodesStat, cuda_COMB_NodesStat, 0.25, Vn);
  cudaDeviceSynchronize();
 
  std::cout << "comb_porog: " << comb_porog  << std::endl;
  
  get_cuda_comb_size<<<grid_size_Vn,block_size>>>(comb_porog,cuda_Edges_d_by_Nodes, cuda_Nodes, cuda_COMB_NodesStat, cuda_comb_sizeV,cuda_comb_sizeE,  cuda_Edges_d, cuda_Edges_v1, cuda_Edges_v2,Vn );
  cudaDeviceSynchronize();
  
  
  get_cuda_comb_v1<<<1,1>>>(cuda_comb_sizeV, cuda_N , cuda_NN , Vn );
  cudaDeviceSynchronize();
  get_cuda_comb_v1<<<1,1>>>(cuda_comb_sizeE, cuda_E , cuda_EE , Vn );
  cudaDeviceSynchronize();
  
  cudaMemcpy(&N,   cuda_N, sizeof(int) , cudaMemcpyDeviceToHost);
  cudaMemcpy(&NN, cuda_NN, sizeof(int) , cudaMemcpyDeviceToHost);
  cudaMemcpy(&E,   cuda_E, sizeof(int) , cudaMemcpyDeviceToHost);
  cudaMemcpy(&EE, cuda_EE, sizeof(int) , cudaMemcpyDeviceToHost);
  
  std::cout << "comb_N:" << N << " " << NN  << "\n";
  std::cout << "comb_E:" << E << " " << EE  << "\n";
  
  //free_combs_cuda();
  
  cudaMalloc((void**)&comb_id,  sizeof(int)*N );
  cudaMalloc((void**)&comb_start,  sizeof(int)*N );
  cudaMalloc((void**)&comb_stat,  sizeof(int)*N );
  cudaMalloc((void**)&comb_V,  sizeof(int)*NN );
  cudaMalloc((void**)&comb_start_E,  sizeof(int)*E );
  cudaMalloc((void**)&comb_stat_E,  sizeof(int)*E );
  cudaMalloc((void**)&comb_E,  sizeof(int)*EE );
  
  
  get_cuda_comb_v2<<<1,1>>>(cuda_comb_sizeV, comb_id , comb_start, comb_stat, Vn );
  cudaDeviceSynchronize();

  get_cuda_comb_v2<<<1,1>>>(cuda_comb_sizeE, comb_id , comb_start_E, comb_stat_E, Vn );
  cudaDeviceSynchronize();
  
  int grid_size_Cn  = ((  N + block_size)   / block_size);
    
  get_cuda_comb_fill<<<grid_size_Cn,block_size>>>(comb_porog, comb_id, comb_start, comb_stat, comb_V, comb_start_E, comb_stat_E, comb_E, cuda_Edges_d_by_Nodes, cuda_Nodes, cuda_COMB_NodesStat,  cuda_Edges_d, cuda_Edges_v1, cuda_Edges_v2,  N );
  cudaDeviceSynchronize();
    
  dtimeBG = dtime_usec(dtimeBG);
  
  cuda_get_Nodes_StatAB<<<grid_size_Vn,block_size>>>(comb_porog,cuda_Edges_d_by_Nodes, cuda_Nodes, cuda_COMB_NodesStat, cuda_Nodes_A,cuda_Nodes_B,Vn);
  cudaDeviceSynchronize();
  
  cudaCheckErrors("cudaMemcpy fail get_combs_cuda() A");
  
  std::cout << "cuda_comb RUN dtime: " << dtimeBG/(double)USECPSEC  << std::endl;


  
}

//############################################################################################
//############################################################################################
//############################################################################################
//############################################################################################



__global__  void shift_comb_IRREGULAR_B
(
double shiftV, int type,
int *comb_id_IRREGULAR ,  int * cuda_block_id, double *cuda_F,
int *comb_id, int *comb_start, int *comb_stat, int *comb_V, 
int * cuda_Nodes_A, int * cuda_Nodes_B,
int n 
)
{
  
  int Tid = blockIdx.x * blockDim.x + threadIdx.x;
  
  if(Tid < n){
   int tid  = comb_id_IRREGULAR[Tid] ;
   int cmb_s    = comb_start[tid];
   int cmb_size = comb_stat[tid];
   
    int shift = 0 ;
    for(int i=0; i < cmb_size; ++i ){
    
     int nid =  comb_V[cmb_s +i];
     int A   = cuda_Nodes_A[nid];
     int B   = cuda_Nodes_B[nid];
     
     int v = 4 - 2*A - B;
      
     if(v != 0){
      if(v > 0){shift = 1;}else{shift = -1;}
     
      double a_shift = shiftV * shift ;
      //double b_shift = 1e-2 * shift * (-1) ;   
      
      cuda_F[nid] += a_shift;
      cuda_block_id[nid] = 1;
      if(type == 1){break;}
     }
   }
 }
}

void graph_c::copy_IRREGULAR_B_ToDevice(){ 
   
   setvalue<<<grid_size_Vn,block_size>>>(cuda_block_id, 0 , Vn );
   
   comb_id_IRREGULAR_size = A->comb_IRREGULAR_B.size();
   int * tmp_comb_IRREGULAR_B  = (int*)malloc(sizeof(int) *comb_id_IRREGULAR_size );
   for(int i=0;i<comb_id_IRREGULAR_size;i++){tmp_comb_IRREGULAR_B[i] = A->comb_IRREGULAR_B[i];}
   cudaMalloc((void**)&comb_id_IRREGULAR,  sizeof(int) * comb_id_IRREGULAR_size);
   cudaMemcpy(comb_id_IRREGULAR,   tmp_comb_IRREGULAR_B, sizeof(int)*comb_id_IRREGULAR_size , cudaMemcpyHostToDevice);
   
   
   
   for(int i=0;i<comb_id_IRREGULAR_size;i++){std::cout << "tmp_comb_IRREGULAR_B[i]: " << i << " " << tmp_comb_IRREGULAR_B[i]  << std::endl;}
   free(tmp_comb_IRREGULAR_B);
   cudaCheckErrors("cudaCheckErrors copy_IRREGULAR_B_ToDevice()");
}   
   
   
void graph_c::cuda_IRREGULAR_B_recompute(double shiftV, int type){ 
   
   
   int IRRn = comb_id_IRREGULAR_size;
   int grid_size_IRRn  = ((  IRRn + block_size)   / block_size);
   shift_comb_IRREGULAR_B<<<grid_size_IRRn,block_size>>>(shiftV,type,comb_id_IRREGULAR,cuda_block_id, cuda_Nodes_F, comb_id, comb_start, comb_stat, comb_V, cuda_Nodes_A,cuda_Nodes_B,IRRn);
   cudaDeviceSynchronize();
   
   cudaFree(comb_id_IRREGULAR);
   
   cudaCheckErrors("cuda_IRREGULAR_B_recompute A");
   
   cuda_recompute_Edges<<<grid_size_En,block_size>>>(cuda_Edges_d, cuda_Edges_D, cuda_Edges_v1,  cuda_Edges_v2,  cuda_Nodes_F, En );
   cudaDeviceSynchronize();
   for(int go = 0 ; go < 20; go++){ block_EF_iter_BLOCK();}
   
  
   setvalue<<<grid_size_Vn,block_size>>>(cuda_block_id, 0 , Vn );
   cudaDeviceSynchronize();
   
   
   run_fix_cuda();
    
   std::cout << "cuda_IRREGULAR_B_recompute error_A " << error_A << " \n"; 
    
   cudaCheckErrors("cuda_IRREGULAR_B_recompute B");
   //cuda_print_vector(12466, 12467, cuda_Nodes_F);
   //cuda_print_node(12466,6 ,cuda_Edges_d_by_Nodes,  cuda_Nodes, cuda_NodesStat, cuda_Edges_d, cuda_Edges_D, cuda_Edges_v1,cuda_Edges_v2);
   
   
   
   
}   


///#################################################################################
///#################################################################################
///#################################################################################
///#################################################################################
///#################################################################################


__global__  void get_cuda_comb_one( 
int tid, 
double comb_porog, 
int *comb_V, 
int *comb_V_size, 
int *comb_E,  int *comb_E_size, 
double ** cuda_Edges_d_by_Nodes, int *cuda_Nodes, int *cuda_NodesStat, double * cuda_Edges_d, int * cuda_Edges_v1, int * cuda_Edges_v2)
{
  
  
   int filter_edge[comb_size_MAX];int size_e =0;
   int filter_node[comb_size_MAX];int size_f =0; 
   int list_A[comb_size_list];int size_A = 0;
  
  filter_node[0] = tid;size_f++;
  list_A[0]      = tid;size_A++;

 while(size_A > 0){
   int list_B[comb_size_list];int size_B = 0;
  
  for(int iA=0; iA < size_A; ++iA ){
   int lA = list_A[iA];
   double ** ptr = cuda_Edges_d_by_Nodes + cuda_Nodes[lA];
  
   for(int i=0; i < cuda_NodesStat[lA]; ++i ){
   
     if( (*ptr[i] >= -comb_porog)  ){
      if(  (*ptr[i] <= comb_porog) ){
       
       
       double * e = ptr[i];
       int indE  = e - cuda_Edges_d;
        
       int v1   =  cuda_Edges_v1[indE];
       int v2   =  cuda_Edges_v2[indE];
        
       int de = array_filtr(filter_edge,  size_e , indE );
       if(de == -1){  filter_edge[size_e] = indE;size_e++;} 
        
       int v = v1;
       if(v1 == lA){v =v2;}
       
       int da = array_filtr(filter_node,  size_f , v  );
       if(da == -1){ list_B[size_B] = v; size_B++;filter_node[size_f] = v;size_f++;}
       
     }}
    }
   }
   
    for(int iB=0; iB < size_B; ++iB ){list_A[iB] = list_B[iB];}
    size_A = size_B; size_B = 0;
   
  }
   
   for(int i=0; i < size_f; ++i ){comb_V[i] = filter_node[i];}
   for(int i=0; i < size_e; ++i ){comb_E[i] = filter_edge[i];}
   
   *comb_V_size = size_f;
   *comb_E_size = size_e;
}




void graph_c::get_comb_vid(int vid, double comb_porog){  
  
  cuda_partition_v3<<<grid_size_Vn,block_size>>>(cuda_Edges_d_by_Nodes, cuda_Nodes, cuda_EF_NodesStat, cuda_COMB_NodesStat, 0.1, Vn);
  cudaDeviceSynchronize();
  
  
  int *cuda_comb_VID, *cuda_comb_VID_size,  *cuda_comb_EID, *cuda_comb_EID_size;
  cudaMalloc((void**)&cuda_comb_VID,  sizeof(int) * 300);
  cudaMalloc((void**)&cuda_comb_VID_size,  sizeof(int) * 1);
  cudaMalloc((void**)&cuda_comb_EID,  sizeof(int) * 600);
  cudaMalloc((void**)&cuda_comb_EID_size,  sizeof(int) * 1);
  
  std::cout << "comb_porog: " << comb_porog  << std::endl;
  
  get_cuda_comb_one<<<1,1>>>(vid, comb_porog,  cuda_comb_VID, cuda_comb_VID_size,  cuda_comb_EID, cuda_comb_EID_size,
  cuda_Edges_d_by_Nodes, cuda_Nodes, cuda_COMB_NodesStat,  cuda_Edges_d, cuda_Edges_v1, cuda_Edges_v2);
  cudaDeviceSynchronize();
  
  int *comb_VID, comb_VID_size,  *comb_EID, comb_EID_size;
  
  cudaMemcpy(&comb_VID_size, cuda_comb_VID_size, sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(&comb_EID_size, cuda_comb_EID_size, sizeof(int), cudaMemcpyDeviceToHost);
  
  comb_VID  = (int*)malloc(sizeof(int) * comb_VID_size);
  comb_EID  = (int*)malloc(sizeof(int) * comb_EID_size);
  
  cudaMemcpy(comb_VID, cuda_comb_VID, sizeof(int)*comb_VID_size, cudaMemcpyDeviceToHost);
  cudaMemcpy(comb_EID, cuda_comb_EID, sizeof(int)*comb_EID_size, cudaMemcpyDeviceToHost);
  
  for(int i=0; i < comb_VID_size; ++i ){
   std::cout << "comb_VID: " << comb_VID[i]  << std::endl;
   cuda_print_node_porog(10*comb_porog, comb_VID[i],4 ,cuda_Edges_d_by_Nodes,  cuda_Nodes, cuda_NodesStat, cuda_Edges_d, cuda_Edges_D,    cuda_Edges_v1,cuda_Edges_v2, cuda_Nodes_grad);
  }
  
  
  
}




