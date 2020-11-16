#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <stdlib.h>
#include <map>
#include <algorithm>    // std::sort  
#include <math.h>
#include "F2M.h"
#include <unordered_map>
#include <random>

bool compare_edge_Do(edge * a, edge * b) {return (a->Do < b->Do);}
bool compare_edge(edge * a, edge * b) {return (a->d < b->d);}

void node::resortDo(){std::sort(edges.begin(), edges.end(), compare_edge_Do);}
void node::rankEDo(){for(int i = 0; i < edges.size(); i++){if(ID == edges[i]->v1->ID){edges[i]->rankDo1 = i ;}else{edges[i]->rankDo2 = i;}}}
void node::resort(){std::sort(edges.begin(), edges.end(), compare_edge);}
void node::resort(int top_k){std::sort(edges.begin(), edges.begin() + top_k, compare_edge);}


void node::get_grad(double V_break){ 
  grad = -(edges[1]->d + edges[2]->d)/2.0;
  if( (edges[1]->d < -1e3*V_break) && (edges[2]->d > 1e3*V_break) ){grad = 0;}
}



void node::print(int k , std::ofstream & myfile){
   myfile << "NODE " << ID << "\n";
   for (int i=0;i<edges.size();i++){
     myfile << "edges: " << edges[i]->id << " " << edges[i]->d << " " << edges[i]->v1->ID << " " << edges[i]->v2->ID << "\n"; 
     if(i == k){break;}
   }
   myfile << "ENDNODE\n";	 
}

void node::print(int k ){
   std::cout << "NODE " << ID << " " << id << " " << statusA << " " << statusB << "\n";
   for (int i=0;i<edges.size();i++){
     std::cout << "edges: " << edges[i]->id << " " << edges[i]->d << " " << edges[i]->D << " " << edges[i]->v1->id << " " << edges[i]->v2->id <<  " " << edges[i]->rankDo1 << " " << edges[i]->rankDo2 << "\n"; 
     if(i == k){break;}
   }
   std::cout << "ENDNODE\n";	 
}


void node::get_statusAB(double V_break ){
    statusA =0 ; statusB = 0 ;
    eA.clear();eB.clear();
	for (int i=0;i<edges.size();i++){
	   if((edges[i]->d ) <  V_break){
        if((edges[i]->d ) > -V_break){statusB++;eB.push_back(edges[i]);}else{statusA++;eA.push_back(edges[i]);} //
	   }else{break;}
    }
}



void node::get_eA_X(graph * G ){ 

 for (int i=0;i<eA.size();i++){
  int v0 = ID;
  int v1 = G->E[eA[i]->id]->v1->ID;
  int v2 = G->E[eA[i]->id]->v2->ID;
  if(v0 < v1){  G->X[eA[i]->id] = 1;}
  if(v0 < v2){  G->X[eA[i]->id] = 1;}
 }
 
}






























