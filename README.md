## Fractional 2-matching GPU solver
Fractional 2-matching (F2M) is a linear programming relaxation of the 2-matching problem which, in turn, is a relaxation of traveling salesman problem. Here we present a fractional 2-matching problem GPU solver based on gradient descent algorithm which exploits the special structure of F2M. The algorithm  can be run in parallel  mode and was implemented as CUDA C/C++ program.


## Requirements:

You must (CUDA):

* have a cuda device with [Compute Capability at least 3.x.](https://developer.nvidia.com/cuda-gpus)
 
* have installed the [CUDA Toolkit](https://developer.nvidia.com/cuda-toolkit)
 
* have installed [Thrust](https://developer.nvidia.com/thrust) (the CUDA C++ template library) 

You must (C++):
 
 * have installed [Boost](https://www.boost.org/users/download/)

## Compilation

To compile the code just run "make" command. 

This would create a file:  F2M.solver.exe. 




## Input/Output

To run program : 

./F2M.solver.exe file N

You must provide 2 arguments: 
* file  is file in [TSPLIB](http://www.math.uwaterloo.ca/tsp/vlsi/) format 
* N  the number N (the top N closest neighbors for each node). Please note that for large scale problems (> 20,000 nodes )  the N should not exceed 200-500 or you risk to run out of memory (cuda device memory is commonly a limiting factor).

To run test example: 

 * ./F2M.solver.exe examples/vsli/xqf131.tsp 20"

As output, 2 files are provided:

 * file with optimal solution:  examples/vsli/xqf131.tsp.top20.F2M_solution
 * file with optimal value   :  examples/vsli/xqf131.tsp.top20.F2M_report
 
The optimal solution for F2M is a set of edges with value 1 and a set of edges with value ½. 

You can find a number of exapmles in the folders "examples/vsli/" and "examples/art/" 
