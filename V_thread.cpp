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
//#################################################
//##################################################
//## MULTITHREAD
//##################################################
//##################################################
//##################################################

typedef struct {
int thread_id;
int i1;
int i2;
std::vector<node*> * V;
double V_break;
int type ; 
double grad_max;
int i_max; 
graph * G; 
} thread_argument;

void * runthread_V(void *arg) {
  thread_argument*  tArg =(thread_argument*)arg;
  
  int thread_id                        =   tArg->thread_id;
  int i1                               =   tArg->i1;
  int i2                               =   tArg->i2;
  std::vector<node*>  * V              =   tArg->V;
  int type                             =   tArg->type; 
  double V_break                       =   tArg->V_break; 
  graph * G                            =   tArg->G; 
  
  if(type == 1){for(int i=i1; i < i2; ++i ){ (*V)[i]->resort();}}

  if(type == 2){
   tArg->grad_max = fabs( (*V)[i1]->grad);
   tArg->i_max = i1;
   for(int i=i1; i < i2; ++i ){ 
    
    (*V)[i]->get_grad(V_break);
    if( fabs( (*V)[i]->grad) > tArg->grad_max){tArg->grad_max = fabs( (*V)[i]->grad); tArg->i_max = i;}
   }
  }
  
 if(type == 3){for(int i=i1; i < i2; ++i ){(*V)[i]->get_statusAB(V_break);}}
 
// if(type == 4){for(int i=i1; i < i2; ++i ){(*V)[i]->get_eA_X( G );}}
  
  pthread_exit(NULL);
  
} 


void graph::V_thread(int k_THREADS, int type , double  V_break  )
{   
   
   
   pthread_t threads[k_THREADS];
   int rc;
   int i;
   
   pthread_attr_t attr;
   void *status;

   // Initialize and set thread joinable
   pthread_attr_init(&attr);
   pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
   
   int thread_size = V.size() / k_THREADS + 1;
   
   std::vector<thread_argument*> TA(k_THREADS);
   for(int i=0; i < k_THREADS; ++i ){
      //std::cout << "main() : creating thread, " << i << "\n";
          
      thread_argument* arg = new thread_argument();
      TA[i] = arg;
      
      arg->thread_id               = i;
      arg->i1                      = (i  )*thread_size;
      arg->i2                      = (i+1)*thread_size;
      arg->V                       = &V;
      arg->type                    = type;
      arg->V_break                 = V_break; 
      arg->G                       = this; 
      
      if(arg->i2 > V.size() ){arg->i2 = V.size();} 

      pthread_create(&threads[i], NULL, runthread_V, arg);

   }

  // free attribute and wait for the other threads
  //pthread_attr_destroy(&attr);  
  for(int i=0; i < k_THREADS; ++i ){ 
      rc = pthread_join(threads[i], &status);
      if (rc){
         std::cout << "Error:unable to join," << rc << "\n";
         exit(-1);
      }
  }
   

  
   
   
    if(type == 2){
    
     int     i_max  = -1;
     double grad_max = 0;
     for(int i=0; i < k_THREADS; ++i ){ 
      if(fabs(TA[i]->grad_max) > grad_max){ grad_max = fabs(TA[i]->grad_max);i_max = TA[i]->i_max; }
     }
    
    std::cout << "V_thread i_max grad_max: " << i_max << " " << grad_max<< "\n";
    
    }
}
 
 
//##################################################
//##################################################
//##################################################
/*
typedef struct {
int thread_id;
int i1;
int i2;
int  E_max;
graph * Ga;
graph * Gb; 
} thread_argument_b;

void * runthread_VE(void *arg) {
  thread_argument_b*  tArg =(thread_argument_b*)arg;
  
  int thread_id                        =   tArg->thread_id;
  int i1                               =   tArg->i1;
  int i2                               =   tArg->i2;
  
  int  E_max = =   tArg->E_max;
    
  graph * Ga                            =   tArg->Ga;
  graph * Gb                            =   tArg->Gb; 
  
  for(int i=i1; i < i2; ++i ){ 
   
   for (int ie=0;ie < Ga->V[i]->edges.size();ie++){
    
    int id1 = Ga->V[i]->edges[ie]->v1->id;
    int id2 = Ga->V[i]->edges[ie]->v2->id;
    
    if(id1 < i){continue;}
    if(id2 < i){continue;}
    
    int rankDo1 = Ga->V[i]->edges[ie]->rankDo1;
    int rankDo2 = Ga->V[i]->edges[ie]->rankDo2;
    
    if(rankDo1 <= E_max ){Ga->Ef[Ga->V[i]->edges[ie]->id]= 1;}
    if(rankDo2 <= E_max ){Ga->Ef[Ga->V[i]->edges[ie]->id]= 1;}
   
   }
  }
  
  pthread_exit(NULL);
  
} 


void graph::VE_thread(int k_THREADS, int E_max, graph * Gb  )
{   
   pthread_t threads[k_THREADS];
   int rc;
   int i;
   
   pthread_attr_t attr;
   void *status;

   // Initialize and set thread joinable
   pthread_attr_init(&attr);
   pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
   
   int thread_size = V.size() / k_THREADS + 1;
   
   
   for(int i=0; i < k_THREADS; ++i ){
                
      thread_argument_b* arg = new thread_argument_b();
      
      arg->thread_id               = i;
      arg->i1                      = (i  )*thread_size;
      arg->i2                      = (i+1)*thread_size;
      arg->E_max                   = E_max;
      arg->Ga                      = this; 
      arg->Gb                      = Gb; 
      
      if(arg->i2 > V.size() ){arg->i2 = V.size();} 

      pthread_create(&threads[i], NULL, runthread_VE, arg);

   }

  // free attribute and wait for the other threads
  //pthread_attr_destroy(&attr);  
  for(int i=0; i < k_THREADS; ++i ){ 
      rc = pthread_join(threads[i], &status);
      if (rc){
         std::cout << "Error:unable to join," << rc << "\n";
         exit(-1);
      }
  }
    
    
} 

*/




