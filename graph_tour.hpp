
graph::graph(std::string file, int  E_max)
{


//#######################
unsigned long long dtime; 
dtime = dtime_usec(0);
//#######################

fileOUT_report   = file + ".top" + std::to_string(E_max) + ".F2M_report"; 
fileOUT_solution = file + ".top" + std::to_string(E_max) + ".F2M_solution"; 
double d_min = 1e16 , d_max = -1e16; 


XYDSRT  = new FXY  (file, E_max);
std::string fileDM  = file + ".top" + std::to_string(E_max) + ".DM"; 
XYDSRT->save_to_file(fileDM, E_max);

int Ek = XYDSRT->cuda_MyEdgeSort_size;


std::unordered_map<int,int> stat_node;
std::unordered_map<int,node*> map_node;
int countEMAX = 0 ;
Vn = 0;
En = 0; 
for (int e = 0; e<Ek; e++)
{  
      if(countEMAX == XYDSRT->Vn){
         std::cout << "countEMAX break: " << e  << std::endl;
         break;
      }
      
      //int k1     = XYDSRT->h_MyEdgeSort[e].k1;
      //int k2     = XYDSRT->h_MyEdgeSort[e].k2;
      
      int i     = XYDSRT->h_MyEdgeSort[e].k1;
      int j     = XYDSRT->h_MyEdgeSort[e].k2;
      double  d = XYDSRT->h_MyEdgeSort[e].d;
      
      //if(e < 10){ std::cout << "d: " << i << " " << j << " " << d  << std::endl;}
      
      if(d == 0){
        
        
        std::cout << "ERROR d = 0: " << i << " " << j << " " << d  << std::endl;
        exit(0);
      }
      
      
        
      if(Vn < XYDSRT->Vn){  
      node * v1;
      node * v2;
      if(map_node.count(i)){ v1 = map_node[i];}else{ v1 = new node(i);map_node[i] = v1;Vn++;} //
      if(map_node.count(j)){ v2 = map_node[j];}else{ v2 = new node(j);map_node[j] = v2;Vn++;} //
      }
      
      
      
      if(stat_node[i] > E_max){if(stat_node[j] > E_max){continue;}}
      
      stat_node[i]++;
      stat_node[j]++;
      
      if(stat_node[i] == E_max){countEMAX++;stat_node[i]++;}
      if(stat_node[j] == E_max){countEMAX++;stat_node[j]++;}
      
      if(d < d_min){d_min = d;}
      if(d > d_max){d_max = d;}
      
      
      
      edge * eA = new edge(En,map_node[i],map_node[j], d);
      
      long long key = key_vv(i,j); 
      if(vv2edge.count(key)){
       std::cout << "ERROR key duplicate : " << key << " " << i << " " << j  << " " << d  << std::endl;
        exit(0);
      } 
      
      vv2edge[key]=En;
      
      En++;
      
      map_node[i]->edges.push_back(eA);
      map_node[j]->edges.push_back(eA);
      E.push_back(eA);
}

free(XYDSRT->h_MyEdgeSort);

for(int i = 0; i < E.size(); i++){E[i]->add_distance(d_min);}
V.resize(Vn);
int iV = 0;
for (auto i : map_node){V[iV] = i.second;V[iV]->id = iV; iV++;}
for(int i = 0; i < V.size(); i++){V[i]->rankEDo();}

dtime = dtime_usec(dtime);
std::cout << "creation graph object dtime: " << dtime/(double)USECPSEC  << std::endl;
//####################### 
MAX_THREADS = std::thread::hardware_concurrency();
//####################### 

std::cout << "Edges:" << En << "\n";
std::cout << "Nodes:" << Vn << "\n";


}

//#####################################################//

//#####################################################//
void  graph::getOVAL_XY(){

  OVAL = 0;
   
  for(auto ie : X ){
    
    int v1 = E[ie.first]->v1->ID;
    int v2 = E[ie.first]->v2->ID;
    
    double d = X[ie.first];
        
    double x1 = XYDSRT->X[v1];
    double y1 = XYDSRT->Y[v1];
    
    double x2 = XYDSRT->X[v2];
    double y2 = XYDSRT->Y[v2];
    
    double dX = (x2-x1);
    double dY = (y2-y1);
    
    OVAL += d * sqrt(dX*dX + dY*dY);
    
  }
  
  std::cout << "OVAL XY: " << OVAL   << "\n";
}   

void graph::check_solution(){

  std::unordered_map<int,float> stst;
  for(auto ie : X ){
    
    stst[E[ie.first]->v1->id] += X[ie.first];
    stst[E[ie.first]->v2->id] += X[ie.first]; 
    
  }
    
  int count2 = 0 , ERROR_STAT = 0;
  for(int i=0; i < V.size() ; ++i ){
   
    if(stst.count(i) ){
     if( (stst[i] > 2.0 - 1e-2 )&&(stst[i] < 2.0 + 1e-2 )  ){count2++;}
     else{std::cout << "graph::check_solution() ERROR_STAT i = " << i << "  " << stst[i] << "\n";ERROR_STAT++;V[i]->print(5);}
    
    }
  }
  std::cout << "count2 " << count2  << " Vn " << Vn << "\n";

} 
	
void graph::getOVAL(){
  
  OVAL = 0;
  
  MAXrankDo = 0;
  
  for(auto ie : X ){
    int eid = ie.first;
    double Deid = E[ie.first]->Do;
    int r1 = E[ie.first]->rankDo1;
    int r2 = E[ie.first]->rankDo2;
    
    int r = r1;
    if(r > r2){r = r2;}
    if(r > MAXrankDo){MAXrankDo = r;}
    
    OVAL += X[eid] * Deid;
  }
  
}

void graph::print_file_solution(int IRR_B){

  getOVAL();
  //graph::getOVAL_XY();
  std::ofstream myfile;
  myfile.open (fileOUT_solution);
  
  
  for(auto ie : X ){
    int eid = ie.first;
    
    int k1    = E[eid]->v1->ID;
    int k2    = E[eid]->v2->ID;
    
    //int k1    = XYDSRT->K[E[eid]->v1->ID];
    //int k2    = XYDSRT->K[E[eid]->v2->ID];
    
    //if(E[eid]->v1->ID == 0){std::cout  << "E[eid]->v1->ID == 0: " << E[eid]->v1->id << " " << k1 << "\n";}
    //if(E[eid]->v2->ID == 0){std::cout  << "E[eid]->v2->ID == 0: " << E[eid]->v2->id << " " << k2 << "\n";}
    //if(k1 == 0){std::cout  << "E[eid]->v1->ID == 0: " << E[eid]->v1->id << " " << k1 << "\n";}
    //if(k2 == 0){std::cout  << "E[eid]->v2->ID == 0: " << E[eid]->v2->id << " " << k2 << "\n";}
    
    
    myfile   << k1 << " " <<  k2 << " " << std::setprecision(4) << X[eid]  <<  "\n";
   }
  
  
  myfile.close();
}



void graph::print_file_report(int IRR_B){
  
  std::ofstream myfile;
  myfile.open (fileOUT_report);
  if(IRR_B == 0){
  myfile << "Optimal Found: " << "Yes" << "\n";
  myfile << "Optimal Value: " << std::setprecision(15)  << OVAL << "\n";
  myfile << "Time: " << dtimeRUN/(double)USECPSEC << "\n";
  myfile << "MAX edge rank: " << MAXrankDo << "\n";
  
  std::cout  << "Optimal Found: " << "Yes" << "\n";
  std::cout  << "Optimal Value: " << std::setprecision(15)  << OVAL << "\n";
  std::cout  << "RunTime: " << dtimeRUN/(double)USECPSEC << "\n";
  std::cout  << "MAX edge rank: " << MAXrankDo << "\n";
  std::cout  << "Solution File: " << fileOUT_solution << "\n"; 
  
  
  
  }else{
   myfile << "ERROR: Optimal Solution not Found\n";
  std::cout << "ERROR: Optimal Solution not Found\n";
  }
  myfile.close();
}














