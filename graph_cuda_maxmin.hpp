
struct compare_abs
{
  __host__ __device__
  bool operator()(const double  lhs, const double  rhs) { return fabs(lhs) < fabs(rhs);}
};

double  min_node_grad(double *cuda_v,  int V){
  thrust::device_ptr<double> d_ptr = thrust::device_pointer_cast(cuda_v);
  thrust::device_vector<double>::iterator d_it = thrust::min_element(d_ptr, d_ptr + V, compare_abs());
  int min_index = d_it - (thrust::device_vector<double>::iterator)d_ptr;
  double min_grad = 0;
  cudaMemcpy(&min_grad, cuda_v+min_index, sizeof(double) , cudaMemcpyDeviceToHost);
    
  return min_grad;
}

int min_node_grad_S(double *cuda_v,  int V){
  thrust::device_ptr<double> d_ptr = thrust::device_pointer_cast(cuda_v);
  thrust::device_vector<double>::iterator d_it = thrust::min_element(d_ptr, d_ptr + V, compare_abs());
  int min_index = d_it - (thrust::device_vector<double>::iterator)d_ptr;
  
  return min_index;
}




double  max_node_grad(double *cuda_v,  int V){
  thrust::device_ptr<double> d_ptr = thrust::device_pointer_cast(cuda_v);
  thrust::device_vector<double>::iterator d_it = thrust::max_element(d_ptr, d_ptr + V, compare_abs());
  int max_index = d_it - (thrust::device_vector<double>::iterator)d_ptr;
  double max_grad = 0;
  cudaMemcpy(&max_grad, cuda_v+max_index, sizeof(double) , cudaMemcpyDeviceToHost);
    
  return max_grad;
}


int max_node_grad_S(double *cuda_v,  int V){
  thrust::device_ptr<double> d_ptr = thrust::device_pointer_cast(cuda_v);
  thrust::device_vector<double>::iterator d_it = thrust::max_element(d_ptr, d_ptr + V, compare_abs());
  int max_index = d_it - (thrust::device_vector<double>::iterator)d_ptr;
  
  return max_index;
}


int max_stat_v(int *cuda_v,  int Vn){
  thrust::device_ptr<int> d_ptr = thrust::device_pointer_cast(cuda_v);
  thrust::device_vector<int>::iterator d_it = thrust::max_element(d_ptr, d_ptr + Vn);
  int max_index = d_it - (thrust::device_vector<int>::iterator)d_ptr;
  int max_stat =0;
  cudaMemcpy(&max_stat, cuda_v+max_index, sizeof(int) , cudaMemcpyDeviceToHost);
  //std::cout << "max_stat_v " <<  max_index << " " << max_stat  << "\n";
return   max_stat;
}

int min_stat_v(int *cuda_v,  int Vn){
  thrust::device_ptr<int> d_ptr = thrust::device_pointer_cast(cuda_v);
  thrust::device_vector<int>::iterator d_it = thrust::min_element(d_ptr, d_ptr + Vn);
  int max_index = d_it - (thrust::device_vector<int>::iterator)d_ptr;
  int max_stat =0;
  cudaMemcpy(&max_stat, cuda_v+max_index, sizeof(int) , cudaMemcpyDeviceToHost);
  //std::cout << "min_stat_v " <<  max_index << " " << max_stat  << "\n";
return   max_stat;
}


