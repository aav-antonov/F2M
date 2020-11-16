#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <stdlib.h>
#include "F2M.h"




int main (int argc, char * argv[]){

if( argc != 3){
 std::cout << "usage: ./F2M.RUN.exe fileXY topK\nfileXY == file in TSPLIB format\ntopK == the number of top K closest neighbors for each node"  << std::endl;
 exit(0);
}

std::string fileXY                = std::string(argv[1]);
std::string topK                  = std::string(argv[2]);

int E_max = std::stoi(topK);

graph   * P  = new graph (fileXY,E_max);
graph_c * Pc = new graph_c (P);
P->RUN();
delete P;
delete Pc; 
 
return(0);
}















