#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <stdlib.h>
#include <cstdlib>
#include <map>
#include <algorithm>    // std::sort  
#include <math.h>
#include <iomanip>
#include "F2M.h"

//##################################################
//##################################################
//## MULTITHREAD
//##################################################

typedef struct {
int thread_id;
int i1;
int i2;
std::vector<edge*> * E;
double * F;
int type ;    
} thread_argument;

void * runthread_EF(void *arg) {
  thread_argument*  tArg =(thread_argument*)arg;
  
  int thread_id                        =   tArg->thread_id;
  int i1                               =   tArg->i1;
  int i2                               =   tArg->i2;
  std::vector<edge*>  * E              =   tArg->E;
  double * F                           =   tArg->F;
  int type                             =   tArg->type;
   
  
  if(type == 1){for(int i=i1; i < i2; ++i ){(*E)[i]->d = (*E)[i]->D + F[(*E)[i]->v1->id]/2 + F[(*E)[i]->v2->id]/2;}}
  

  pthread_exit(NULL);
} 



void graph::EF_thread(int k_THREADS , int type)
{   
   pthread_t threads[k_THREADS];
   int rc;
   int i;
	 
   pthread_attr_t attr;
   void *status;

   // Initialize and set thread joinable
   pthread_attr_init(&attr);
   pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

   int thread_size = 0;
   if(type == 1){ thread_size = (E.size() / k_THREADS) + 1;}
   //if(type == 2){ thread_size = (E.size()  / k_THREADS) + 1;}
   //if(type == 3){ thread_size = (E.size()  / k_THREADS) + 1;}
   std::cout << "thread_size, " << thread_size << " " << E.size()  << "\n";

   for(int i=0; i < k_THREADS; ++i ){
      //std::cout << "main() : creating thread, " << i << "\n";
            
      thread_argument* arg = new thread_argument();
      
      arg->thread_id               = i;
      arg->i1                      = (i  )*thread_size;
      arg->i2                      = (i+1)*thread_size;
      arg->F                       =   F;
      arg->E                       = & E;
      arg->type                    = type;
       
      if(type == 1){if(arg->i2 > E.size() ){arg->i2 = E.size();}}
      //if(type == 2){if(arg->i2 > E.size() ){arg->i2 = E.size();}} 
      //if(type == 3){if(arg->i2 > E.size() ){arg->i2 = E.size();}}

      pthread_create(&threads[i], NULL, runthread_EF, arg);

   }

  // free attribute and wait for the other threads
   
   for(int i=0; i < k_THREADS; ++i ){ 
      rc = pthread_join(threads[i], &status);
      if(rc){std::cout << "Error:unable to join," << rc << "\n";exit(-1);}
   }
}

//##################################################
//##################################################     
//##################################################






