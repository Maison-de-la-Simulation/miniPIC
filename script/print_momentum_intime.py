
import pandas as pd
import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt

import sys
import os


filename = sys.argv[1]
ndt=0
pX=[]
pY=[]
pZ=[]
EX=[]
EY=[]
EZ=[]

if os.path.isdir(filename):
    f_list=os.listdir(filename)
    f_list.sort()
    for f in f_list:
        ndt+=1
        data = pd.read_csv(filename+f)
        pX[ndt:]=np.array(data['momentx'].to_numpy())
        pY[ndt:]=np.array(data['momenty'].to_numpy())
        pZ[ndt:]=np.array(data['momentz'].to_numpy())
        EX[ndt:]=np.array(data['Ex'].to_numpy())
        EY[ndt:]=np.array(data['Ey'].to_numpy())
        EZ[ndt:]=np.array(data['Ez'].to_numpy())


iterations=np.arange(0,ndt)

dt=0.1

####################
#     THEORIC
####################

q=1.
m=1836.125

p0=np.array([pX[0], pY[0], pZ[0]])
print("p0: {}".format(p0))
E0=np.array([EX[0], EY[0], EZ[0]])
print("E0: {}".format(E0))

Fl0=(q/m)*E0
print("Lorentz force: {}".format(Fl0))



pth=np.zeros((3,ndt))

for t,it in enumerate(iterations):
    pth[:,it] = Fl0*dt*t+p0
    
print("pth0: {}".format(pth[:,0]))




fig, axs = plt.subplots(1,3)
axs[0].plot(iterations, pX, marker='+', label='Compute')
axs[0].plot(iterations, pth[0], label='Theoric')
axs[0].set_xlabel("t")
axs[0].set_ylabel("px")
axs[1].plot(iterations, pY, marker='+', label='Compute')
axs[1].plot(iterations, pth[1], label='Theoric')
axs[1].set_xlabel("t")
axs[1].set_ylabel("py")
axs[2].plot(iterations, pZ, marker='+', label='Compute')
axs[2].plot(iterations, pth[2], label='Theoric')
axs[2].set_xlabel("t")
axs[2].set_ylabel("pz")

axs[0].legend()
axs[1].legend()
axs[2].legend()




plt.show()
