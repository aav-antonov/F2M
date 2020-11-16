# remove!!!! .o files before run
# *****************************************************
# Variables to control Makefile operation
NVCC = nvcc
CXX = g++ -std=c++11 


# ****************************************************
# Targets needed to bring the executable up to date
all: F2M.solver.exe

clean:
	rm RUN.exe  


F2M.solver.exe:    graph_cuda.o   node.o  E_thread.o V_thread.o comb_thread.o Ax_b.o   
	$(NVCC)  -o F2M.solver.exe   graph_cuda.o  node.o E_thread.o V_thread.o comb_thread.o Ax_b.o RUN.cpp
    

graph_cuda.o: graph_cuda.cu 
	$(NVCC)  -c  graph_cuda.cu

node.o:
	$(CXX) -c -g node.cpp
	
E_thread.o:
	$(CXX) -c -g E_thread.cpp	


V_thread.o:
	$(CXX) -c -g V_thread.cpp
	
comb_thread.o:
	$(CXX) -c -g comb_thread.cpp
	
Ax_b.o:
	$(CXX) -c -g Ax_b.cpp
				




