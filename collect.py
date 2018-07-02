import os
import matplotlib.pyplot as plt
import numpy as np
Kn = "./Kn_pump.o";
for i in np.arange(2, 3,1):
    for j in [0.1, 0.3, 0.8, 1, 4, 10]:
      os.system(Kn +" "+str(5**i) +" "+str(j) );
      wd = str(j) +"_"+  str(5**i);
      os.system("mkdir data/"+ wd);
      os.system("mv tmp/* data/"+ wd +"/");
