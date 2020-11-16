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
#include "Ax_b.h"
//#################################################
//#################################################
std::vector<std::string> split(const std::string &text, char sep);
std::string IntToString ( int number );
std::string DoubleToString ( double number );
std::string CharToString ( char number );
//#################################################
//#################################################
//#################################################
//######################################################
typedef struct {
int thread_id;
int i1;
int i2;
graph * G;
std::vector<comb*> * VCB;;
} thread_argument_C;

 
comb::comb( graph * GA, int iC ){
	G = GA;
	
	//std::cout << "comb::comb iC " <<  iC << "\n";
	
	int v_s    = G->comb_V_start[iC];
    int v_size = G->comb_V_stat[iC];
    
    //std::cout << "v_s   " <<  v_s   << "\n";
    //std::cout << "v_size   " <<  v_size   << "\n";
    
    for(int i=0; i < v_size; ++i ){ nodes[G->comb_V[v_s+i]]++; }
    int e_s    = G->comb_E_start[iC];
    int e_size = G->comb_E_stat[iC];
    
    for(int i=0; i < e_size; ++i ){ edgesB[G->comb_E[e_s+i]]++; }
    
    //std::cout << "comb::comb  end \n";
   
}


void * runthread_COMB(void *arg) {
  thread_argument_C*  tArg =(thread_argument_C*)arg;
  
  int thread_id                         =   tArg->thread_id;
  int i1                                =   tArg->i1;
  int i2                                =   tArg->i2;
  graph * G                             =   tArg->G;
  std::vector<comb*> * VCB              =   tArg->VCB;
  
  for(int i=i1; i < i2; ++i ){
     comb * C = new comb( G, i );
     (*VCB).push_back(C);
  }

  pthread_exit(NULL);
} 


void graph::COMB_thread(int k_THREADS )
{ 

   std::cout << "void graph::COMB_thread Combs.size(), " << Combs.size() << "\n";
   for(int i=0; i < Combs.size(); ++i ){delete Combs[i];}
   Combs.clear();
   std::cout << "after clear Combs.size(), " << Combs.size() << "\n";
   
   pthread_t threads[k_THREADS];
   int rc;
   
   pthread_attr_t attr;
   void *status;

   // Initialize and set thread joinable
   pthread_attr_init(&attr);
   pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
   
    std::cout << "combs_size cuda, " << combs_size << "\n";
   
   int thread_size = combs_size / k_THREADS + 1;
   
   std::vector<std::vector<comb*>>  VCB;
   for(int i=0; i < k_THREADS; ++i ){
     std::vector<comb*> VCBi;
     VCB.push_back(VCBi);
   }
   
   for(int i=0; i < k_THREADS; ++i ){
      thread_argument_C* arg = new thread_argument_C();
      
      arg->thread_id               = i;
      arg->i1                      = (i  )*thread_size;
      arg->i2                      = (i+1)*thread_size;
      arg->G                       = this ;
      arg->VCB                     = &VCB[i] ;
      if(arg->i2 > combs_size ){arg->i2 = combs_size;}
      
      pthread_create(&threads[i], NULL, runthread_COMB, arg);

   }
  
  // free attribute and wait for the other threads
   for(int i=0; i < k_THREADS; ++i ){ 
      rc = pthread_join(threads[i], &status);
      if (rc){
         std::cout << "Error:unable to join," << rc << "\n";
         exit(-1);
      }
     }

  // assemble output
  for(int i=0; i < k_THREADS; ++i ){
   for(int j=0; j < VCB[i].size(); ++j ){
     Combs.push_back(VCB[i][j]);
   }
  }


std::cout << "Combs.size() " <<  Combs.size() << "\n";
}
//########################################################################
//########################################################################


void  comb::print(int q){
 
 if(q > 0){ for ( auto jf : edgesB){G->E[jf.first]->print(XB[jf.first]);}}
 for ( auto jf : nodes){G->V[jf.first]->print(5);}
}


void  comb::print(){
 
 for ( auto jf : nodes){
  int sv =  4 - 2* G->V[jf.first]->statusA -   G->V[jf.first]->statusB;
  std::cout << "comb; " << jf.first << " " << G->V[jf.first]->statusA << " " << G->V[jf.first]->statusB << "  sv "<< sv << "\n";
 }
}


void  comb::printE(){
 
 for ( auto jf : nodes){
  std::cout << "comb v; " << jf.first << " " << G->V[jf.first]->statusA << " " << G->V[jf.first]->statusB << " " << G->F[jf.first] << "\n";
 } 
 
 std::cout << " edges \n";
 
 for ( auto jf : edgesB){
  std::cout << "comb e; " << G->E[jf.first]->v1->id << "  sv "<< G->E[jf.first]->v2->id << "\n";
 }
}

//########################################################################
//########################################################################
void graph::get_comb_IRREGULAR_A(  ){

   comb_IRREGULAR_A.clear();
   
   
   for(int i=0; i < Combs.size(); ++i ){
    if(  Combs[i]->STATUS == -1){comb_IRREGULAR_A.push_back(i);}
   }
   std::cout << "get_comb_IRREGULAR_A; " << comb_IRREGULAR_A.size() << "\n";
}

void comb::solve_simple(  ){

   for (auto i : nodes){
    int sv =  4 - 2* G->V[i.first]->statusA -   G->V[i.first]->statusB;
    if(sv != 0){STATUS = -1;break;}
   }
   
   if(STATUS == 1){
    for (auto i : nodes){
     for (int j=0;j< G->V[i.first]->eB.size(); j++){
      XB[G->V[i.first]->eB[j]->id] = 0.5;
     }
    }
   }
  
}



void * runthread_COMB_solve(void *arg) {
  thread_argument_C*  tArg =(thread_argument_C*)arg;
  
  int thread_id                         =   tArg->thread_id;
  int i1                                =   tArg->i1;
  int i2                                =   tArg->i2;
  graph * G                             =   tArg->G;
  
  for(int i=i1; i < i2; ++i ){ G->Combs[i]->solve_simple();}

  pthread_exit(NULL);
} 


void graph::COMB_solve_thread(int k_THREADS )
{ 
      
   pthread_t threads[k_THREADS];
   int rc;
   
   pthread_attr_t attr;
   void *status;

   // Initialize and set thread joinable
   pthread_attr_init(&attr);
   pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
   
   int thread_size = Combs.size() / k_THREADS + 1;
        
   for(int i=0; i < k_THREADS; ++i ){
      thread_argument_C* arg = new thread_argument_C();
      
      arg->thread_id               = i;
      arg->i1                      = (i  )*thread_size;
      arg->i2                      = (i+1)*thread_size;
      arg->G                       = this ;
      
      if(arg->i2 > Combs.size() ){arg->i2 = Combs.size();}
      
      pthread_create(&threads[i], NULL, runthread_COMB_solve, arg);

   }
  
  // free attribute and wait for the other threads
   for(int i=0; i < k_THREADS; ++i ){ 
      rc = pthread_join(threads[i], &status);
      if (rc){
         std::cout << "Error:unable to join," << rc << "\n";
         exit(-1);
      }
     }

}

//########################################################################
//###############COMB_solve_AxB_thread####################################
//########################################################################

void graph::get_comb_IRREGULAR_B(  ){
    std::cout << "get_comb_IRREGULAR_B; START " << Combs.size() << "\n";
   comb_IRREGULAR_B.clear(); 
   for(int i=0; i < Combs.size(); ++i ){
    if(  Combs[i]->AxB_status == 0){comb_IRREGULAR_B.push_back(i);}
   }
   
   for(int i=0; i < comb_IRREGULAR_B.size() ; ++i ){
     std::cout << "get_comb_IRREGULAR_B;" << i << " " << comb_IRREGULAR_B[i] << "\n";
     
     if(comb_IRREGULAR_B[i] == 37027){Combs[comb_IRREGULAR_B[i]]->printE();}
   }
   std::cout << "get_comb_IRREGULAR_B; " << comb_IRREGULAR_B.size() << "\n";
}

void comb::solve_AxB( ){
    
    std::vector<int>  AxB;
    std::unordered_map<int,int>  mapX, mapXR;
    int countX = 0;
	for ( auto jf : nodes){
	 int j = jf.first;
	 for(int i1=0; i1 < G->V[jf.first]->eB.size() ; ++i1 ){	
       int eBid = G->V[jf.first]->eB[i1]->id;
	   int i,v;
	   if(mapX.count(eBid)){i = mapX[eBid];v = 1;}
                    else{i = countX; v = 1;mapXR[i] = eBid;mapX[eBid] = i;countX++;}

				AxB.push_back(i);
				AxB.push_back(v);
	}
   
    AxB.push_back(-1);
    int b = 2 - G->V[jf.first]->eA.size();
	AxB.push_back(b);
     
   } 


    Ax_b * AB =  new Ax_b(AxB);
	//AB->print();
    AB->solve();
    AxB_status = AB->go_max;
    AxB_error  = AB->d_max_error;
    for(int j = 0; j < AB->x.size(); j++){XB[mapXR[j]] = AB->x[j];}
	delete  AB;

}

void * runthread_COMB_solve_AxB(void *arg) {
  thread_argument_C*  tArg =(thread_argument_C*)arg;
  
  int thread_id                         =   tArg->thread_id;
  int i1                                =   tArg->i1;
  int i2                                =   tArg->i2;
  graph * G                             =   tArg->G;
  
  for(int i=i1; i < i2; ++i ){ G->Combs[G->comb_IRREGULAR_A[i]]->solve_AxB();}

  pthread_exit(NULL);
} 


void graph::COMB_solve_AxB_thread(int k_THREADS )
{ 
      
   pthread_t threads[k_THREADS];
   int rc;
   
   pthread_attr_t attr;
   void *status;

   // Initialize and set thread joinable
   pthread_attr_init(&attr);
   pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

   int thread_size = comb_IRREGULAR_A.size() / k_THREADS + 1;
        
   for(int i=0; i < k_THREADS; ++i ){
      thread_argument_C* arg = new thread_argument_C();
      
      arg->thread_id               = i;
      arg->i1                      = (i  )*thread_size;
      arg->i2                      = (i+1)*thread_size;
      arg->G                       = this ;
      
      if(arg->i2 > comb_IRREGULAR_A.size() ){arg->i2 = comb_IRREGULAR_A.size();}
      
      pthread_create(&threads[i], NULL, runthread_COMB_solve_AxB, arg);

   }
  
  // free attribute and wait for the other threads
   for(int i=0; i < k_THREADS; ++i ){ 
      rc = pthread_join(threads[i], &status);
      
           
      if (rc){
         std::cout << "Error:unable to join," << rc << "\n";
         exit(-1);
      }
     }
  
  std::cout << "COMB_solve_AxB_thread END\n";
  
}
//########################################################################


















