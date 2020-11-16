

template <typename T>
__device__ int qpart_ptr_B(T ** list,  int size , double avg  )
{  

    T * temp;
    int pivot, i, j;
    int low  = 0;
    int high = size-1; 
    
    
        i = low;
        j = high;
        
        while (i < j) 
        {   while (*list[i] <= avg && i <= high){i++;}
            while (*list[j] > avg && j  >= low ){j--;}
            if (i < j){temp = list[i]; list[i] = list[j]; list[j] = temp;}
        }
        
        pivot = j+1;
        
return pivot;        
}

template <typename T>
__device__ int qpart_ptr_EF(T ** list,  int size , double avg  )
{  

    T * temp;
    int pivot, i, j;
    int low  = 0;
    int high = size-1; 
    
    int go = 1;
    while(go > 0){
    go = 0;
       
        i = low;
        j = high;
        
        while (i < j) 
        {   while (*list[i] <= avg && i <= high){i++;}
            while (*list[j] > avg && j  >= low ){j--;}
            if (i < j){temp = list[i]; list[i] = list[j]; list[j] = temp;}
        }
    
    pivot = j+1;
    if(pivot < 3){avg += 1;go =1;}
  
    }    
return pivot;        
}


__global__  void cuda_partition_v2(double **cuda_Edges_d_by_Nodes, int *cuda_Nodes, int *cuda_NodesStat, int *cuda_pivots, double avg, int n)
{
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  
  if (tid < n){
    if(cuda_NodesStat[tid] < 10){cuda_pivots[tid] = cuda_NodesStat[tid];}
    else{
    cuda_pivots[tid] = qpart_ptr_EF(cuda_Edges_d_by_Nodes + cuda_Nodes[tid], cuda_NodesStat[tid],avg);
    } 
    //cuda_pivots[tid] = pivot;
    //if(pivot > 5){ int pivotB = qpart_ptr_B(cuda_Edges_d_by_Nodes + cuda_Nodes[tid], pivot,1e-3);}
   }
}


//#################################################
//#################################################

__global__  void cuda_partition_v3(double **cuda_Edges_d_by_Nodes, int *cuda_Nodes, int *cuda_NodesStat, int *cuda_pivots, double avg, int n)
{
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  
  if (tid < n){
    cuda_pivots[tid] = qpart_ptr_B(cuda_Edges_d_by_Nodes + cuda_Nodes[tid], cuda_NodesStat[tid],avg);
   }
}



//#################################################
//#################################################

__global__  void cuda_partition_v4(double **cuda_Edges_d_by_Nodes, int *cuda_Nodes, int *cuda_NodesStat, int *cuda_pivots, int avg_k, int n)
{
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  
  if (tid < n){
  
   double ** ptr = cuda_Edges_d_by_Nodes + cuda_Nodes[tid];
    
  double avg = *ptr[0]; 
  if(cuda_NodesStat[tid] < avg_k){cuda_pivots[tid] = cuda_NodesStat[tid];}
   
  else{
   
   for(int i=1; i < avg_k; ++i ){
    double v = *ptr[i];
    if(v > avg){avg = v;}
   }
   
   cuda_pivots[tid] = qpart_ptr_EF(cuda_Edges_d_by_Nodes + cuda_Nodes[tid], cuda_NodesStat[tid],avg);
   
  }

   
   
   
  }
}




