
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import sys
import os


filename = sys.argv[1]
ndt=0
X=[]
Y=[]
Z=[]
pX=[]
pY=[]
pZ=[]

if os.path.isdir(filename):
    f_list=os.listdir(filename)
    f_list.sort()
    for f in f_list:
        ndt+=1
        data = pd.read_csv(filename+f)
        X[ndt:]=np.array(data['posx'].to_numpy())
        Y[ndt:]=np.array(data['posy'].to_numpy())
        Z[ndt:]=np.array(data['posz'].to_numpy())
        pX[ndt:]=np.array(data['momentx'].to_numpy())
        pY[ndt:]=np.array(data['momenty'].to_numpy())
        pZ[ndt:]=np.array(data['momentz'].to_numpy())
        
iterations=np.arange(0,ndt)

dt=0.1

####################
#     THEORIC
####################

x0=np.array([X[0], Y[0], Z[0]])
print("x0: {}".format(x0))
p0=np.array([pX[0], pY[0], pZ[0]])
print("p0: {}".format(p0))
gamma=np.sqrt(1+np.sum(p0**2))
print("gamma: {}".format(gamma))

v0=p0/gamma
print("v0: {}".format(v0))



xth=np.zeros((3,ndt))

for t,it in enumerate(iterations):
    xth[:,it] = v0*t*dt+x0

print("xth0: {}".format(xth[0]))



fig, axs = plt.subplots(1,3)
axs[0].plot(iterations, X, marker='+', label='Compute')
axs[0].plot(iterations, xth[0], label='Theoric')
axs[1].plot(iterations, Y, marker='+', label='Compute')
axs[1].plot(iterations, xth[1], label='Theoric')
axs[2].plot(iterations, Z, marker='+', label='Compute')
axs[2].plot(iterations, xth[2], label='Theoric')

axs[0].legend()
axs[1].legend()
axs[2].legend()


plt.show()




