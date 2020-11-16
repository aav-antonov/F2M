#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <stdlib.h>
#include <map>
#include <algorithm>    
#include <math.h>
#include <unordered_map>
#include <random>
#include <cmath>
#include <iomanip>      // std::setprecision
#include <cstring> 
#include <boost/algorithm/string/predicate.hpp>
//#################################################
#include <time.h>
#include <sys/time.h>
//#################################################
#include <cuda.h>
#include <cuda_runtime.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/sequence.h> 
#include <thrust/extrema.h>
//##########################################
//##########################################
#include "F2M.h"
#include "graph_cuda_define.hpp"
#include "graph_cuda_maxmin.hpp"
#include "graph_cuda_print.hpp"

#include "FXY.hpp"
#include "graph_tour.hpp"
#include "graph_cuda_init.hpp"

#include "graph_cuda_partition.hpp"
#include "graph_cuda_grad.hpp"
#include "graph_cuda_comb.hpp"



//#################################################
//#################################################
//#################################################

void graph::RUN(){


//#######################
dtimeRUN = dtime_usec(0);
//#######################
std::cout << "graph::RUN() " <<  std::endl;


error_B      = 1e-7; 
error_B_comb = 1e2;

graph_c * Pc = Ac;

Pc->error_B       = error_B;
Pc->error_B_comb  = error_B_comb ;

Pc->run_start_cuda();
Pc->get_combs_cuda(0 );
int IRR_B = host_combs_check();


int goshift = 0;
if(IRR_B > 0){goshift = 1;}

 double shiftV = 0.5;
 
 while(goshift == 1){ 
     
     Pc->copy_IRREGULAR_B_ToDevice();
     Pc->cuda_IRREGULAR_B_recompute(shiftV,0);
     
     Pc->free_combs_cuda();
     Pc->get_combs_cuda(1);  
     
     IRR_B = host_combs_check();
     
     if(IRR_B == 0){break;}
     
     Pc->copy_IRREGULAR_B_ToDevice();
     Pc->cuda_IRREGULAR_B_recompute(shiftV,1);
     
     Pc->free_combs_cuda();
     Pc->get_combs_cuda(1);  
     
     IRR_B = host_combs_check();
     
          
     if(IRR_B == 0){break;}
     
     shiftV = - shiftV;
     
     Pc->copy_IRREGULAR_B_ToDevice();
     Pc->cuda_IRREGULAR_B_recompute(shiftV,0);
     
     Pc->free_combs_cuda();
     Pc->get_combs_cuda(1);  
     
     IRR_B = host_combs_check();
     
     if(IRR_B == 0){break;}
     
     shiftV = - shiftV;
     
     Pc->copy_IRREGULAR_B_ToDevice();
     Pc->cuda_IRREGULAR_B_recompute(shiftV,1);
     
     Pc->free_combs_cuda();
     Pc->get_combs_cuda(1);  
     
     IRR_B = host_combs_check();
     
     if(IRR_B == 0){break;}
     
     shiftV = - shiftV;
     
     shiftV *= 2;
     
     if(fabs(shiftV) > 32){break;}  
 }



//#######################
dtimeRUN = dtime_usec(dtimeRUN);
std::cout << "dtimeRUN: " << dtimeRUN/(double)USECPSEC  << std::endl;
//#######################

 if(IRR_B == 0){
  for(int i=0; i < Vn; ++i ){V[i]->get_eA_X( this );}
  
  for(int i=0; i < Combs.size(); ++i ){
   for(auto jf : Combs[i]->XB){   
    X[jf.first] = Combs[i]->XB[jf.first];
    }
  }
     
  check_solution();
  print_file_solution(IRR_B);
  print_file_report(IRR_B);
  
  
 }else{
 
 print_file_report(IRR_B);
 
 }
 
}



//#################################################
//#################################################
//#################################################

void graph_c::run_start_cuda(){ 
 
  
  unsigned long long dtimeBG; 
  dtimeBG = dtime_usec(0);
 
  block_GLOB();

  
  std::cout << "block_GLOB ENDS " <<  std::endl;
   
  
  double M_grad =check_OPT();
  update_EF(1);
  
  int factor = 100;
  for(int go =1; go<= 10000; go++){ 
    
   block_EF(factor);
   double M_grad =check_OPT();
   update_EF(1);

   
   if(factor < 5000){factor += 300;}
   if(fabs(M_grad) < error_B){error_A = fabs(M_grad);break;}
  }
  
  
  
  dtimeBG = dtime_usec(dtimeBG);
  std::cout << "RUN dtime: " << dtimeBG/(double)USECPSEC  << std::endl;
  
  
  cudaCheckErrors("cudaMemcpy fail run_start_cuda()");
  
  
}


  

void graph_c::run_fix_cuda(){ 
  
  unsigned long long dtimeBG; 
  
  dtimeBG = dtime_usec(0);
  
  double M_grad =check_OPT();
  update_EF(1);
  
  
  for(int go =1; go<= 100; go++){ 
    
   block_EF(500);
   double M_grad =check_OPT();
   update_EF(1);
   if(fabs(M_grad) < error_B){error_A = fabs(M_grad);break;}
  }
  
  
  
  dtimeBG = dtime_usec(dtimeBG);
  
} 


void graph_c::update_EF(int ind_to_partition ){ 
  
 
  if(ind_to_partition == -1){
   cuda_partition_v2<<<grid_size_Vn,block_size>>>(cuda_Edges_d_by_Nodes,  cuda_Nodes, cuda_NodesStat, cuda_EF_NodesStat, 1.00,Vn);
   cudaDeviceSynchronize();
  }
  
    
  thrust::inclusive_scan(thrust::device,cuda_EF_NodesStat, cuda_EF_NodesStat + Vn, cuda_EF_Nodes);
  cudaDeviceSynchronize();
  
  cudaMemcpy(&size_EF, cuda_EF_Nodes + Vn-1, sizeof(int) , cudaMemcpyDeviceToHost);
  //std::cout << "size_EF: " << size_EF  << std::endl;   
  
  get_cuda_EF_nodes_NORM<<<grid_size_Vn,block_size>>>( cuda_EF_NodesStat, cuda_EF_Nodes, Vn ); 
  cudaDeviceSynchronize();
   
  cudaFree(cuda_EF);
  cudaMalloc((void**)&cuda_EF, sizeof(double**) * size_EF);
  
  get_cuda_EF<<<grid_size_Vn,block_size>>>(cuda_Edges_d_by_Nodes,  cuda_Nodes, cuda_EF_NodesStat, cuda_EF, cuda_EF_Nodes , Vn);
   
}

double graph_c::check_OPT(){ 
   
  cuda_recompute_Edges<<<grid_size_En,block_size>>>( cuda_Edges_d,  cuda_Edges_D,  cuda_Edges_v1,   cuda_Edges_v2, cuda_Nodes_F, En);
  cudaDeviceSynchronize();
    
  cuda_partition_v2<<<grid_size_Vn,block_size>>>(cuda_Edges_d_by_Nodes,  cuda_Nodes, cuda_NodesStat, cuda_EF_NodesStat, 1.00,Vn);
  cudaDeviceSynchronize();
  
  cuda_get_Nodes_grad_R<<<grid_size_Vn,block_size>>>(cuda_Edges_d_by_Nodes, cuda_Nodes, cuda_EF_NodesStat, cuda_Nodes_grad ,Vn) ;
  cudaDeviceSynchronize();
  
  int iM_grad   = max_node_grad_S(cuda_Nodes_grad,  Vn);
  double M_grad = max_node_grad(cuda_Nodes_grad,  Vn);
  std::cout << "check_OPT() M_grad: " << M_grad << " " << iM_grad << std::endl;
  
  
return M_grad ;
} 
 
  

   


void graph_c::block_EF_iter_BLOCK(){ 
   
    cuda_get_Nodes_grad_R<<<grid_size_Vn,block_size>>>(cuda_EF, cuda_EF_Nodes, cuda_EF_NodesStat, cuda_Nodes_grad,Vn) ;
    cudaDeviceSynchronize();
       
    cuda_update_Nodes_F_BLOCK<<<grid_size_Vn,block_size>>>(cuda_Nodes_grad,  cuda_Nodes_F ,  cuda_block_id,Vn);
    cudaDeviceSynchronize();
    
    cuda_recompute_Edges<<<grid_size_En,block_size>>>(cuda_Edges_d,cuda_Edges_D,cuda_Edges_v1,cuda_Edges_v2,cuda_Nodes_F,En);
    cudaDeviceSynchronize();
    
}      



void graph_c::block_EF_iter(){ 
    
    cuda_get_Nodes_grad_R<<<grid_size_Vn,block_size>>>(cuda_EF, cuda_EF_Nodes, cuda_EF_NodesStat, cuda_Nodes_grad,Vn);
    cudaDeviceSynchronize();
    
    cuda_update_Nodes_F<<<grid_size_Vn,block_size>>>(cuda_Nodes_grad,  cuda_Nodes_F , Vn);
    cudaDeviceSynchronize();
  
    cuda_recompute_Edges<<<grid_size_En,block_size>>>(cuda_Edges_d,cuda_Edges_D,cuda_Edges_v1,cuda_Edges_v2,cuda_Nodes_F,En);
    cudaDeviceSynchronize();
    
}      



void graph_c::block_EF(int N){  
    
  cuda_get_Nodes_grad_R<<<grid_size_Vn,block_size>>>(cuda_EF, cuda_EF_Nodes, cuda_EF_NodesStat, cuda_Nodes_grad,Vn) ;
  cudaDeviceSynchronize();
  
  double M_gradA = max_node_grad(cuda_Nodes_grad,  Vn);
  std::cout << "block_EF: M_gradA " << M_gradA  << "\n"; 
   
  int count_M = 0 ;//, count_stagnation = 0;
  
  for(int go = 0 ; go < N; go++){
    
    block_EF_iter();
   
    count_M++;
    if(count_M >= 50){
    
     
     int max_id = max_node_grad_S(cuda_Nodes_grad,  Vn);
     double M_gradB = max_node_grad(cuda_Nodes_grad,  Vn);
               
     if(fabs(M_gradA) < error_B){ std::cout << "block_EF: " << N << " go: " << go << " break\n";break;}
     
     count_M = 0;
    }
    
   }
   
  
  
}  

//##############################################################################
//##############################################################################
//##############################################################################



  

 

//##############################################################################
//##############################################################################
//##############################################################################

template <typename T>
__global__  void setvalue2V(T * v1, T * v2, T val, int n )
{
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid < n){
   if(v2[tid] > val){v1[tid]    =   val;}else{v1[tid]    =   v2[tid];}
  }
}

void graph_c::block_GLOB(){ 


  unsigned long long  dtimeBG; 
  
  dtimeBG = dtime_usec(0);
  
  int glob_K = 5;
  
  setvalue2V<<<grid_size_Vn,block_size>>>(cuda_EF_NodesStat, cuda_NodesStat, glob_K, Vn );
    
  
  for(int go = 0 ; go < 200; go++){
   
   block_GLOB_N(50);
   
   cuda_partition_v4<<<grid_size_Vn,block_size>>>(cuda_Edges_d_by_Nodes,  cuda_Nodes, cuda_NodesStat, cuda_EF_NodesStat, glob_K,Vn);
   cudaDeviceSynchronize();
   
   //int max_pivot =  max_stat_v(cuda_EF_NodesStat, Vn);
   //int min_pivot =  min_stat_v(cuda_EF_NodesStat, Vn);
   //std::cout << "pivot: " << min_pivot << " " << max_pivot << std::endl;
   
   cuda_get_Nodes_grad_R<<<grid_size_Vn,block_size>>>(cuda_Edges_d_by_Nodes, cuda_Nodes, cuda_EF_NodesStat, cuda_Nodes_grad ,Vn) ;
   cudaDeviceSynchronize();
   
   
   double M_grad = max_node_grad(cuda_Nodes_grad,  Vn);
   //std::cout << "block_GLOB_N: M_grad " << go << " " << M_grad  << std::endl;
   
   
    if(fabs(M_grad) < 0.01){
     std::cout << "block_GLOB_N: M_grad " << go << " " << M_grad  << std::endl;
     
     
    cuda_get_Nodes_grad_R<<<grid_size_Vn,block_size>>>(cuda_Edges_d_by_Nodes, cuda_Nodes, cuda_EF_NodesStat, cuda_Nodes_grad ,Vn) ;
    cudaDeviceSynchronize();
  
    double M_grad = max_node_grad(cuda_Nodes_grad,  Vn);
    std::cout << "block_GLOB_N: M_grad " << go << " " << M_grad  << std::endl;
     
    break;
    }

   
   }
 
  dtimeBG = dtime_usec(dtimeBG);
  std::cout << "dtimeBG: " << dtimeBG/(double)USECPSEC  << std::endl;
     
}

void graph_c::block_GLOB_N(int N){  
  
    
  for(int go = 0 ; go < N; go++){
   
  cuda_get_Nodes_grad_R<<<grid_size_Vn,block_size>>>(cuda_Edges_d_by_Nodes, cuda_Nodes, cuda_EF_NodesStat, cuda_Nodes_grad ,Vn) ;
  cudaDeviceSynchronize();
  
  cuda_update_Nodes_F<<<grid_size_Vn,block_size>>>(cuda_Nodes_grad,  cuda_Nodes_F ,Vn);
 
  cuda_recompute_Edges<<<grid_size_En,block_size>>>(cuda_Edges_d, cuda_Edges_D, cuda_Edges_v1, cuda_Edges_v2,cuda_Nodes_F, En );
  cudaDeviceSynchronize();
  }
  
 
}









