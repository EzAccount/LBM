import os
import matplotlib.pyplot as plt
import numpy as np
Kn = "cmake-build-debug/Kn_pump"
for i in np.arange(50, 550, 50):
    for j in [0.0395, 0.079, 0.158]:
      os.system(Kn +" "+str(i) +" "+str(j) )
      wd = "kn" + str(j) +"_"+  str(i)
      os.system("mv tmp/Tecplot.dat data/"+ wd)
