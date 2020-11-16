import sys 
from os import listdir
from os.path import isfile, join

import os
import re



def parse_file(filepath, data_time , data_oval, data_error):
 with open(filepath) as fp:
       
       for line in fp:
        line = line.replace('\n','')
        line = line.replace('\t','')
        if "Time:" in line:
         data_time[f] = line
        if "ERROR:" in line:
         data_error[f] = line 
        if "Optimal Value:" in line:
         data_oval[f] = line 
        
         
         

mypath   = sys.argv[1];#//train with all files
top_k     = sys.argv[2];#//train with all files

onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]

top_k = int(top_k)

command_cod = "./F2M.solver.exe "

data_time , data_oval , data_error = {} , {} , {}


for f in onlyfiles:
 if(f.endswith('.tsp')):
  
  x = re.search("\d+", f)
  size = 1e12
  if (x):
   print (x.group(0))
   size = int(x.group(0))
   print (size, f)
   
  if(size > 100000):
   continue 
  
  
  command = command_cod + mypath + "/" + f + " " + str(top_k)
  file_out_report = mypath + "/" + f + ".top" + str(top_k) + ".F2M_report";
  
  print(command )
  os.system(command)
  parse_file(file_out_report, data_time , data_oval, data_error)
  #break
  
 
file_out = open("report_dir.txt", "w")
  
for f in data_time:
 f_tsp = f
 if f.endswith(".tsp"): 
    f_tsp = f.replace(".tsp", '') 
 print(f_tsp," ", data_time[f]," ", data_oval[f])
 fstr =  f_tsp + " " + data_time[f] + " " + data_oval[f] + "\n"
 file_out.write(fstr) 
  #sys.exit()
  
file_out.close()  

file_error = open("report_dir_error.txt", "w")
  
for f in data_error:
 f_tsp = f
 if f.endswith(".tsp"): 
    f_tsp = f.replace(".tsp", '') 
 print(f_tsp," ", data_error[f])
 fstr =  f_tsp + " " + data_error[f] + "\n"
 file_error.write(fstr) 
  #sys.exit()
  
file_error.close()  
  
  
