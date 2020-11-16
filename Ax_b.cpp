#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <stdlib.h>
#include <map>
#include <algorithm>    // std::sort  
#include <math.h>
#include <unordered_map>
#include <random>
#include "Ax_b.h"

//#################################################
//#################################################

//#################################################
//#################################################
void Ax_b::parse_Ax_b(std::vector<int> & A){

	int j       = 0;
	int Aj_norm = 0;
	std::vector<int> Aij, Avj;

    for (int k =0;k < A.size();k=k+2){

    	int i= A[k];
    	int v= A[k+1];
    	
      if(i >= -1){}else{break;}
      if(i > i_count){i_count = i;}

      if(i == -1){
    	 b.push_back(v);
    	 Aj.push_back(sqrt(Aj_norm));
    	 Aj_norm = 0;
    	 j++;

    	 Ai.push_back(Aij);
    	 Av.push_back(Avj);
    	 Aij.clear();
    	 Avj.clear();
      }
      else{

       Aij.push_back(i);
       Avj.push_back(v);
       Aj_norm += v*v;

       x_stat[i]++;

      }

     }



j_count = j;
x.resize(i_count+1);
x_grad.resize(i_count+1);
d.resize(j_count);
b_grad.resize(j_count);
for (int i =0;i<i_count;i++){if(x_stat[i] == 0){std::cout << "i:" << i << "; xi-> error;\n";}}
}

Ax_b::Ax_b(std::string file){
     FILE * in;
     in = fopen(file.c_str(), "r");
     std::vector<int>  A;

     int read = 1;
     int i;

     while ( read != EOF) {
      read = fscanf (in,"%d",&i);
      if( read != EOF){A.push_back(i);}
     }

fclose(in);

parse_Ax_b(A);
}

//#################################################
//#################################################

double Ax_b::ax(int j , std::vector<double> &X){
	double v = 0;
	for(int i = 0; i < Ai[j].size(); i++){v += Av[j][i] * X[Ai[j][i]];}
	return v;
}

void Ax_b::d_j(int j){
	double v = ax(j , x);
	d[j] = (v - b[j]) / Aj[j];
}

void Ax_b::d_update(){for(int j = 0; j < Ai.size(); j++){d_j(j);}}

void Ax_b::x_grad_update(){
	vector2null(x_grad);
	d_max_error = 0;
	for(int j = 0; j < Ai.size(); j++){
		d_j(j);
		if(fabs(d[j]) > d_max_error){d_max_error = fabs(d[j]);}
		for(int i = 0; i < Ai[j].size(); i++){x_grad[Ai[j][i]] += Av[j][i] * d[j] / Aj[j];}
	}
}


void Ax_b::b_grad_update(){
	vector2null(b_grad);
	for(int j = 0; j < Ai.size(); j++){b_grad[j] = ax(j , x_grad);	}
}


double Ax_b::product(std::vector<double> &X1, std::vector<double> &X2){
	double v = 0;
	for(int j = 0; j < X1.size(); j++){	v += X1[j]*X2[j];	}
	return v;
}

double Ax_b::product(std::vector<int> &X1, std::vector<double> &X2){
	double v = 0;
	for(int j = 0; j < X1.size(); j++){	v += X1[j]*X2[j];	}
	return v;
}

void Ax_b::x_update(){

	double gg  = product(b_grad, b_grad);
	double gd  = product(d, b_grad);
	double alfa =  - gd/gg;
	for(int j = 0; j < x.size(); j++){ x[j] += x_grad[j] * alfa;}

}


void Ax_b::solve(){
	//Timer timer;
	//timer.start();
    go_max = 0;
	for(int go = 1;go < 1000;go++){
	 x_grad_update();
	 if( d_max_error < d_max_error_STOP){go_max = go;break;}
	 b_grad_update();
	 x_update();
	}

    if(go_max == 0){std::cout << "Ax_b::solve() d_max_error " << d_max_error << "\n";}
	//print(x);
	//timer.stop();
	//std::cout << "d_max_error " << d_max_error << "\n";
	//std::cout << "go_max " << go_max << "\n";
	//std::cout << "Seconds: " << timer.elapsedSeconds() << std::endl;
	//std::cout << "Milliseconds: " << timer.elapsedMilliseconds() << std::endl;

}

void Ax_b::vector2null(std::vector<double> &X1){for(int j = 0; j < X1.size(); j++){X1[j]=0;}}
void Ax_b::print(std::vector<double> &X1){for(int j = 0; j < X1.size(); j++){std::cout << j << ";print;" << X1[j] <<  "\n";}}
void Ax_b::print(){
	for(int j = 0; j < Ai.size(); j++){
		for(int i = 0; i< Ai[j].size(); i++){
			std::cout << Av[j][i] << "*x_" << Ai[j][i] <<  " ";
		}
		std::cout  << " = " << b[j] <<  "\n";
	}
}

//#################################################
//#################################################






























