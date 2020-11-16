

__global__  void cuda_D2_cutN(MyEdge * cuda_MyEdge, int *cutN, int n);
__global__  void cuda_D2_fill(MyEdge * cuda_MyEdge,int * cuda_start,int * cuda_stat,double * cuda_dD_cutoff,double * cuda_X, double * cuda_Y, int * cuda_K, int n);


__global__  void cuda_D2_stat(int max_E_per_V,double size_per_node,int * cuda_stat, double * cuda_dD_cutoff,double * cuda_X, double * cuda_Y, int * cuda_K, int n);

int get_size(std::string file); 
 
double get_KXY(std::string file, int*K , double*X, double*Y);
int get_sizeB(std::string fileA) ;

__global__  void cuda_MyEdge_copy_on_device(MyEdge * dest, MyEdge * src, int n );


__host__ __device__ bool operator<(const MyEdge &lhs, const MyEdge &rhs) { return (lhs.d < rhs.d); };



FXY::~FXY(){free(h_MyEdgeSort);cudaFree(cuda_MyEdgeSort);}

FXY::FXY(std::string file, int top_K ) // type - 1 no random ,2  - add random
{

unsigned long long dtimeBG; 
dtimeBG = dtime_usec(0);

fileXY = file;
top_k  = top_K;

int N = get_size(fileXY);

std::string fileXY_tmp = fileXY + "_tmp";

K  = (int*)malloc(sizeof(int) * N);
X  = (double*)malloc(sizeof(double) * N);
Y  = (double*)malloc(sizeof(double) * N);
 
double size_per_node =  get_KXY(fileXY_tmp, K , X, Y);
std::remove( fileXY_tmp.c_str() ) ;

std::cout << "N " << N << "\n";
std::cout << "Vn " << Vn << "\n";


if(duplicates.size() > 0){
 std::string fileXY_W = fileXY + "_duplicates";
 save_to_file_duplicates(fileXY_W);
}

N = Vn;
std::cout << "N " << N << "\n";

int *cuda_K; double  *cuda_X, *cuda_Y;    
cudaMalloc((void**)&cuda_X, sizeof(double) * N);
cudaMalloc((void**)&cuda_Y, sizeof(double) * N);
cudaMalloc((void**)&cuda_K, sizeof(int) * N);

cudaMemcpy(cuda_K, K, sizeof(int)    * N, cudaMemcpyHostToDevice);
cudaMemcpy(cuda_X, X, sizeof(double) * N, cudaMemcpyHostToDevice);
cudaMemcpy(cuda_Y, Y, sizeof(double) * N, cudaMemcpyHostToDevice);

cudaCheckErrors("cudaMemcpy fail last");    

//exit(0);

int * cuda_stat;
cudaMalloc((void**)&cuda_stat, sizeof(int) * N);
int * cuda_start;
cudaMalloc((void**)&cuda_start, sizeof(int) * N);
double * cuda_dD_cutoff;
cudaMalloc((void**)&cuda_dD_cutoff, sizeof(double) * N);


int grid_size_N = (( N + block_size)   / block_size);
if(top_k > N){top_k = N-1;}
cuda_D2_stat<<<grid_size_N,block_size>>>(top_k, 2*size_per_node,  cuda_stat ,cuda_dD_cutoff, cuda_X, cuda_Y, cuda_K,N); 
cudaDeviceSynchronize();

thrust::inclusive_scan(thrust::device,cuda_stat, cuda_stat + N, cuda_start);


int min_cuda_stat =  min_stat_v(cuda_stat,  N);


int MyEdge_size;
cudaMemcpy(&MyEdge_size, cuda_start + N - 1, sizeof(int), cudaMemcpyDeviceToHost);
std::cout << "MyEdge_size: " << MyEdge_size << "\n";


MyEdge * cuda_MyEdge; 
cudaMalloc((void**)&cuda_MyEdge, sizeof(MyEdge) * MyEdge_size);

cuda_D2_fill<<<grid_size_N,block_size>>>(cuda_MyEdge ,cuda_start, cuda_stat, cuda_dD_cutoff, cuda_X, cuda_Y, cuda_K, N);
cudaDeviceSynchronize();
cudaCheckErrors("cudaMalloc fail MyEdge\n");


thrust::device_ptr<MyEdge> g_ptr(cuda_MyEdge);
thrust::sort(g_ptr, g_ptr + MyEdge_size);

int *cME;
cudaMalloc((void**)&cME, sizeof(int) );
int grid_size_sN = ((  MyEdge_size + block_size)   / block_size);
cuda_D2_cutN<<<grid_size_sN,block_size>>>(cuda_MyEdge, cME ,MyEdge_size);
cudaMemcpy(&cuda_MyEdgeSort_size, cME, sizeof(int), cudaMemcpyDeviceToHost);
std::cout << "cuda_MyEdgeSort_size: " << cuda_MyEdgeSort_size << "\n";


cudaMalloc((void**)&cuda_MyEdgeSort, sizeof(MyEdge) * cuda_MyEdgeSort_size);
cudaCheckErrors("cudaMalloc fail MyEdge 1\n");

int grid_size_NN = ((  cuda_MyEdgeSort_size + block_size)   / block_size);
cuda_MyEdge_copy_on_device<<<grid_size_NN,block_size>>>(cuda_MyEdgeSort, cuda_MyEdge, cuda_MyEdgeSort_size );
cudaCheckErrors("cudaMalloc fail MyEdge 2\n");


h_MyEdgeSort = (MyEdge*)malloc(sizeof(MyEdge) * cuda_MyEdgeSort_size);
cudaMemcpy(h_MyEdgeSort, cuda_MyEdgeSort, sizeof(MyEdge)* cuda_MyEdgeSort_size, cudaMemcpyDeviceToHost);
cudaCheckErrors("FXY::save_to_file 2\n");

/*
 for(int i = 0;i< 100;i++ ){
 
  std::cout << "h_MyEdgeSort: " << i << " " << h_MyEdgeSort[i].d << " " << h_MyEdgeSort[i].v1 << " " << " " << h_MyEdgeSort[i].v2 << " " << h_MyEdgeSort[i].k1 << " " << h_MyEdgeSort[i].k2   << std::endl; 
 
 }
*/


cudaFree(cuda_X);
cudaFree(cuda_Y);
cudaFree(cuda_K);
cudaFree(cuda_stat);
cudaFree(cuda_start);
cudaFree(cuda_dD_cutoff);
cudaFree(cuda_MyEdge);
cudaFree(cuda_MyEdgeSort);

dtimeBG = dtime_usec(dtimeBG);
std::cout << "FXY total time: " << dtimeBG/(double)USECPSEC  << std::endl;

}

void FXY::save_to_file(std::string file_save, int E_max){
  
  std::unordered_map<int,int> stat_node;
  
  std::ofstream myfile;
  myfile.open (file_save);
  for(int i = 0;i< cuda_MyEdgeSort_size;i++ ){
   int v1 = h_MyEdgeSort[i].v1;
   int v2 = h_MyEdgeSort[i].v2;
   
   if(stat_node[v1] > E_max){if(stat_node[v2] > E_max){continue;}}
      
   stat_node[v1]++;
   stat_node[v2]++;
   
   myfile   << v1 << " " <<  v2 << " " << h_MyEdgeSort[i].d  <<  "\n";
  }
  myfile.close();
 
}


void FXY::save_to_file_duplicates(std::string file_save){
  
  
  
  std::ofstream myfile;
  myfile.open (file_save);
  for(auto iv1 : duplicates){
   myfile   << "node used: " << iv1.first << " duplicates: ";
    for(auto iv2 : duplicates[iv1.first]){
     myfile   << iv2.first << " " ;
    }
   myfile   << "\n" ; 
  }
  myfile.close();
 
}


int get_size(std::string fileA) {

unsigned long long dtimeBG; 
dtimeBG = dtime_usec(0);

int start_data = 0, count_data = 0;
std::string fileB = fileA + "_tmp";
std::ofstream myfile (fileB);
  
 
std::ifstream F(fileA);
if (F.is_open()) {
    
    std::string line;
    while (std::getline(F, line)) {
        
        if (boost::starts_with(line, "EOF")){break;}        
        
        if(start_data == 1){
         count_data++;
         myfile << line << " ";
         continue;
        }
        
        if (boost::starts_with(line, "NODE_COORD_SECTION")){start_data = 1;}
    }
    F.close();
}

myfile.close();

 
dtimeBG = dtime_usec(dtimeBG);
std::cout << "read file time 1: " << dtimeBG/(double)USECPSEC  << std::endl;

if(count_data == 0){
 std::cout << "ERROR MESSAGE:"  << std::endl;
 std::cout << "Your file : " << fileA << " seems to be in the wrong format"  << std::endl;
 std::cout << "You must provide file in TSPLIB format, please see examples at http://www.math.uwaterloo.ca/tsp/vlsi/ " << std::endl;
 exit(0);
}
         
return count_data; 
}

double FXY::get_KXY(std::string file, int*K , double*X, double*Y){

std::unordered_map<double,std::unordered_map<double,int>> filtr_x_y;


unsigned long long dtimeBG; 
dtimeBG = dtime_usec(0);      

double min_x, min_y, max_x, max_y;  
   
	 FILE * in;
	 in = fopen(file.c_str(), "r");
	 int i = 0;
     int read = 1;
     while ( read != EOF) {
      int k= -1;
      double x,y;std::unordered_map<int,int> stat_node;
      read = fscanf (in,"%d",&k);
      read = fscanf (in,"%lf",&x);
      read = fscanf (in,"%lf",&y);
      
      x = std::round(x * 1000);
      y = std::round(y * 1000);
      
      x /= 1000;
      y /= 1000;
      
      if(k > -1){
      
      if(filtr_x_y.count(x)){
       if(filtr_x_y[x].count(y)){
        int k0 = filtr_x_y[x][y];
        duplicates[k0][k]++;
        continue; 
      }}
      
      filtr_x_y[x][y] = k;
      
      
      //###########################
       if(i == 0){
         min_x = x; min_y = y; max_x = x; max_y = y;
       } 
       else{
            
         if(min_x > x){min_x = x;}
         if(min_y > y){min_y = y;}
         if(max_x < x){max_x = x;}
         if(max_y < y){max_y = y;}
            
       }
       //###############################
       
       K[i] = k;
       X[i] = x;
       Y[i] = y;
       i++;
       
      }
     }
     fclose(in);
    
dtimeBG = dtime_usec(dtimeBG);
std::cout << "read file time 2: " << dtimeBG/(double)USECPSEC  << std::endl; 
   
//std::cout << "min_x " << min_x << "\n";
//std::cout << "min_y " << min_y << "\n"; 
//std::cout << "max_x " << max_x << "\n"; 
//std::cout << "max_y " << max_y << "\n";  



Vn = i;

int area =  (max_x - min_x) * (max_x - min_x);
double size_per_node = sqrt((double) area / (double)i); 
std::cout << "size_per_node " << size_per_node << "\n";  
return size_per_node;
}


//##############################################################################
__global__  void cuda_MyEdge_copy_on_device(MyEdge * dest, MyEdge * src, int n )
{
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid < n){
     
     dest[tid].k1 = src[tid].k1; 
     dest[tid].k2 = src[tid].k2; 
     dest[tid].v1 = src[tid].v1; 
     dest[tid].v2 = src[tid].v2; 
     dest[tid].d = src[tid].d; 
  }
}

__global__  void cuda_D2_cutN
( 
MyEdge * cuda_MyEdge,
int *cutN,
int n 
)
{
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  
 if (tid < n){
  if (tid > 0){  
   if(cuda_MyEdge[tid-1].d < 1e10){
    if(cuda_MyEdge[tid].d == 1e10){
     *cutN = tid;
   }}
  }
 }
 
}


__global__  void cuda_D2_fill
( 
MyEdge * cuda_MyEdge,
int * cuda_start, 
int * cuda_stat, 
double * cuda_dD_cutoff, 
double * cuda_X, double * cuda_Y, int * cuda_K, 
int n 
)
{
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  
  if (tid < n){
   
   double x1 = cuda_X[tid];
   double y1 = cuda_Y[tid];
   int k1    = cuda_K[tid];
   int i1    = tid; 
    
    
   double max_D =   cuda_dD_cutoff[tid];
   MyEdge * Q = cuda_MyEdge + cuda_start[tid] - cuda_stat[tid] ;
   
   int j = 0; 
   for(int i=0; i < n; ++i ){
    if(i == tid){continue;}
    
    double x2 = cuda_X[i];
    double y2 = cuda_Y[i];
    
    
    double dX = fabs(x1 - x2);
    if(dX >max_D){continue;}
    double dY = fabs(y1 - y2);
    if(dY > max_D){continue;}
    double dD = sqrt(dY*dY + dX*dX); // (int)
    if( dD  > max_D){continue;}
    
    int k2          =   cuda_K[i];
    double max_D_k2 =   cuda_dD_cutoff[i];
    int i2          = i; 
    if(dD < max_D_k2 ){
     
     if(i1 < i2){
     
     Q[j].k1 = k1; 
     Q[j].k2 = k2;
     
     Q[j].v1 = i1; 
     Q[j].v2 = i2;
     Q[j].d = dD;
     
     }else{
     
     
     Q[j].k1 = k2; 
     Q[j].k2 = k1;
     
     Q[j].v1 = i2; 
     Q[j].v2 = i1;
     Q[j].d = 1e10;
     
     }
     
     j++;
     
    }
    else
    {
    
    Q[j].k1 = k1; 
    Q[j].k2 = k2;
    
    Q[j].v1 = i1; 
    Q[j].v2 = i2; 
    if(i1 > i2){
     Q[j].v1 = i2; 
     Q[j].v2 = i1;
     
     Q[j].k1 = k2; 
     Q[j].k2 = k1;
     
    }
    Q[j].d = dD;
    j++;
    
   }
   
  }
 }
}



__device__  void cuda_D2_stat_SUB
( 
int tid,
int l,
int L,
int max_E_per_V,
double size_per_node,   
int * cuda_stat, double * cuda_dD_cutoff, 
double * cuda_X, double * cuda_Y, int * cuda_K, 
int n 
)
{
  
   double x1 = cuda_X[tid];
   double y1 = cuda_Y[tid];
      
     
   double max_D =   size_per_node * l + size_per_node *L;
   
   int * stat  = (int*)malloc(sizeof(int) * l);
   
   for(int j=l-1; j >=0 ; --j ){stat[j] = 0;}
   
   
   for(int i=0; i < n; ++i ){
    if(i == tid){continue;}
    
    double x2 = cuda_X[i];
    double y2 = cuda_Y[i];
    
    double dX = fabs(x1 - x2);
    if(dX >max_D){continue;}
    double dY = fabs(y1 - y2);
    if(dY > max_D){continue;}
    double dD = sqrt(dY*dY + dX*dX); //(int)
    if( dD  > max_D){continue;}
    
    
    for(int j=l-1; j >=1 ; --j ){
     double dD_limit = size_per_node*j + size_per_node *L;
     if(dD_limit > dD){stat[j]++;}else{break;}
    }
    
   }
   
  cuda_stat[tid]      = stat[0];
  cuda_dD_cutoff[tid] =  size_per_node*0 +  size_per_node *L;
  for(int j=1; j < l ; ++j ){ 
  
   if( stat[j] > max_E_per_V){
   
   cuda_stat[tid]      = stat[j];
   cuda_dD_cutoff[tid] =  size_per_node * j +  size_per_node *L;
   break;
   }else{
   
   cuda_stat[tid]      = stat[j];
   cuda_dD_cutoff[tid] =  size_per_node * j +  size_per_node *L;
   
   }
   
  }
  free(stat);
 
}


__global__  void cuda_D2_stat
( 
int max_E_per_V,
double size_per_node,   
int * cuda_stat, double * cuda_dD_cutoff, 
double * cuda_X, double * cuda_Y, int * cuda_K, 
int n 
)
{
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  
  if (tid < n){
      
   int l = 5;
   
   
   cuda_D2_stat_SUB(tid,l, 0, max_E_per_V, size_per_node,  cuda_stat ,cuda_dD_cutoff, cuda_X, cuda_Y, cuda_K,n);  
   
   int L = l;
   
   while(cuda_stat[tid] < max_E_per_V){
     cuda_D2_stat_SUB(tid,l, L, max_E_per_V, size_per_node,  cuda_stat ,cuda_dD_cutoff, cuda_X, cuda_Y, cuda_K,n); 
     L *= 2 ;
   } 
   
 }
}




 
