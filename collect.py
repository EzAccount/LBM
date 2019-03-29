import os
import matplotlib.pyplot as plt
import numpy as np
Kn = "./Kn_pump";
for i in np.arange(1, 30,1):
    for j in np.array( [0.01, 0.1, 0.2, 0.4]):
      os.system(Kn +" "+str(10*i) +" "+str(j/0.84) );
      wd = str(j) +"_"+  str(10*i);
      os.system("mkdir data/"+ wd);
      os.system("mv tmp/* data/"+ wd +"/");
