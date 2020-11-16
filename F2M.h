#ifndef F2M_h
#define F2M_h
#include <unordered_map>
#include <thread>


struct MyEdge
{
    int k1;
    int k2;
    
    int v1;
    int v2;
    double d;
    
};


class node;
class edge;
class comb;
class graph;
class graph_c;

class FXY {
   
   friend graph;
   
  private:
   std::string fileXY;
   int Vn;
   int top_k;
   
   int  *K ; double *X, *Y;
   
    
   std::unordered_map<int,std::unordered_map<int,int>> duplicates;
   
   int cuda_MyEdgeSort_size = 0;
   MyEdge *cuda_MyEdgeSort, *h_MyEdgeSort;
   
   double get_KXY(std::string file, int*K , double*X, double*Y);
   
  public:
  void save_to_file(std::string file_save, int max_E);
  void save_to_file_duplicates(std::string file_save);
  FXY(std::string fileXY, int top_k ); 
  ~FXY();
};


class node {


   friend edge;
   friend graph;
   friend graph_c;
   friend comb;
   
   private:
    
    std::vector<edge*> edges;

   public:
    int  id, ID  ;
    double F     = 0 ; //scaling factor
    double grad  = 0 ; //gradient 
    int  statusA = 0 , statusB = 0 ;
    std::vector<edge*> eA, eB;
    
    void resortDo();    
    void rankEDo();
    void resort();
    void resort(int top_k);
    void get_grad(double V_break);
    void print(int k);
    void print(int k , std::ofstream & myfile); 
    void get_statusAB(double V_break );
    void get_eA_X(graph * G );
   
    node(int k){ID = k;};
    ~node(){}
};


class edge
{
	friend node;
	friend graph;
	friend graph_c;
    friend comb;
       
    private:
     
     int rank1, rank2; 
     int rankDo1, rankDo2;
     
     int D_state = 0; 
     void print(){std::cout << "edge " << id << ";" << D << ";"  << ";" << v1->ID << ";" << v2->ID << ";" << rankDo1 << ";" << rankDo2 << "\n"; } 
     void print(double x){std::cout << "edge " << id << " " << D << " " << Do << " " << v1->ID << " " << v2->ID << " " << x << "\n"; }
     void add_distance( double D_minA){Do_norm = Do/D_minA; D = Do_norm;d = D; }
     
     
    public:
    int id;
    node * v1, *v2; 
    double d , D,  Do_norm , Do ;
    
    edge(int idA, node * v1A, node * v2A , double DoA ){ 
            Do = DoA; id = idA;
            if( (*v1A).ID < (*v2A).ID ){v1 = v1A; v2 = v2A;}
                                   else{v2 = v1A; v1 = v2A;}
     }
     ~edge(){}
};




class graph {
   
   friend node;
   friend edge;
   friend graph_c;
   friend comb;
   
   private:
    
    double error_B,error_B_comb , error_A;
       
    std::string fileOUT_solution, fileOUT_report ;
    unsigned long long dtimeRUN;
    graph_c * Ac;
    
    FXY   * XYDSRT;
    
    int MAX_THREADS;  
     
    double  OVAL = 0; 
    int MAXrankDo;
     
    int En, Vn; 
    std::vector<node*> V;
    std::vector<edge*> E;
    
    double * F;
    
    int combs_size;
    int *comb_V_start, *comb_V_stat , *comb_V;
    int *comb_E_start, *comb_E_stat , *comb_E;
    std::unordered_map<long long,int> vv2edge;
       
    void check_solution();
    void getOVAL();
    void getOVAL_XY();
    void print_file_solution(int IRR_B);
    void print_file_report(int IRR_B);
    
   public:
    
    std::unordered_map<int,double>  X;
    
    void EF_thread(int k_THREADS , int type);
    void V_thread(int k_THREADS, int type , double  V_break  );
    
    std::vector<int> comb_IRREGULAR_A , comb_IRREGULAR_B;
    std::vector<comb*> Combs;
    void COMB_thread(int k_THREADS );
    void COMB_solve_thread(int k_THREADS );
    void COMB_solve_AxB_thread(int k_THREADS );
    
    void get_comb_IRREGULAR_A();
    void get_comb_IRREGULAR_B();
    
    void RUN();
    int host_combs_check();
    
    long long key_vv(int i, int j){long long v; ;if(i>j){v = 1e8 * j + i;return v;}else{v= 1e8 * i + j;return v;}} 
    
   
    graph(std::string file, int  E_max);
           
    ~graph(){
     for(int i=0; i < E.size() ; ++i ){delete E[i];}
     for(int i=0; i < V.size() ; ++i ){delete V[i];}
     }

};



class comb {

   friend graph;
   

   private:
    double AxB_error = 0; int AxB_status =1;
    int COMBID , STATUS = 1; 
    graph * G;
    std::unordered_map<int,int> nodes , edgesB ;
    std::unordered_map<int,float>  XB;
    
   
   public:
    void printE();
    void print();
    void print(int q);
    void solve_AxB( );
    void solve_simple();
    
     comb(  graph * GA, int i);
    ~comb(){}

};


class graph_c {

   friend graph;

   private:
    
    
    double error_B,error_B_comb , error_A;
  
    //####################
    //####################
    graph * A;
    int En, Vn; 
    double *Edges_d, *Edges_D;
    int *Edges_v1, *Edges_v2 ;
    int *i_v1,  *i_v2;
   
        
    double  *Nodes_F; 
       
    int *Nodes , *NodesStat;
    
    
    //####################
    //####################
    
    int grid_size_Vn, grid_size_En;
    
    double *cuda_Edges_d, *cuda_Edges_D ;
    int *cuda_Edges_v1, *cuda_Edges_v2 ;
    
    double *cuda_Nodes_grad,  *cuda_Nodes_F;
    int * cuda_Nodes_grad_stagnation_ID , cuda_Nodes_grad_stagnation_ID_size , * cuda_Nodes_grad_stagnation_ID_check;
        
    int *cuda_Nodes, *cuda_NodesStat;
    double **cuda_Edges_d_by_Nodes;
   
    int * cuda_EF_Nodes ,* cuda_EF_NodesStat;
    double ** cuda_EF;
    int       size_EF;
    
    int * cuda_stagnation_control;
    
    int *cuda_Nodes_A, *cuda_Nodes_B;
    int* cuda_block_id;
    
    int *cuda_COMB_NodesStat;
    
    int N  = 0, NN = 0;
         
    int *cuda_N, *cuda_NN;
    int *cuda_comb_sizeV ;
    int *comb_id, *comb_start, *comb_stat, *comb_V;
    
    int E  = 0, EE = 0;      
    int *cuda_E, *cuda_EE;
    int *cuda_comb_sizeE ;
    int *comb_start_E, *comb_stat_E, *comb_E;
    
    int *comb_id_IRREGULAR , comb_id_IRREGULAR_size;

     
    //#####################//
    //#####################//
    
    void copyCombsToHost();
    void copyBackToHost();
    
    void run_fix_cuda();
    void run_start_cuda();
    
    void get_combs_cuda(int id );
    void free_combs_cuda();
    
    void get_comb_vid(int vid, double comb_porog);
    void get_stagnation();
    void get_stagnation_nodes();
    
    double check_OPT();
    
    void block_GLOB();    
    void block_GLOB_N(int N);
    
    void update_EF(int ind_to_partition);
    void block_EF(int N);
    void block_EF_iter();
    int  block_EF_iter_stagnation(int VSID_SIZE, int * cuda_Nodes_grad_VSID_IND,int * cuda_Nodes_grad_VSID, double M_gradB);
    //int  block_EF_iter_stagnation(double M_gradS );
    void block_EF_iter_BLOCK();
    int block_EF_stagnation_check(double * M_gradA);
    
    void get_combs(double comb_porog);
    
    void copy_IRREGULAR_B_ToDevice();
    void cuda_IRREGULAR_B_recompute(double shiftV, int type);
    
    
    void init(); 
    void init_cuda();
    
   public:
   
    graph_c(graph * A1){A = A1;init();init_cuda();};
    ~graph_c();

};






#endif

