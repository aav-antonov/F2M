
#define block_size 256
//#define error_B 1e-5
//#define error_B_comb 1e2
//#define MAX_edges_per_node 10

#define comb_size_MAX 550
#define comb_size_list 150



#define cudaCheckErrors(msg) \
    do { \
        cudaError_t __err = cudaGetLastError(); \
        if (__err != cudaSuccess) { \
            fprintf(stderr, "Fatal error: %s (%s at %s:%d)\n", \
                msg, cudaGetErrorString(__err), \
                __FILE__, __LINE__); \
            fprintf(stderr, "*** FAILED - ABORTING\n"); \
            exit(1); \
        } \
    } while (0)


#include <stdlib.h>
#include <time.h>
#include <sys/time.h>

unsigned long long dtime_usec(unsigned long long prev){
#define USECPSEC 1000000ULL
  timeval tv1;
  gettimeofday(&tv1,0);
  return ((tv1.tv_sec * USECPSEC)+tv1.tv_usec) - prev;
}


__device__ int array_filtr(int* filtr,  int size_filtr , int key  )
{  
   for(int i=0; i < size_filtr; ++i ){
    if(filtr[i] == key){return i;}
   }
    
return -1;        
}









