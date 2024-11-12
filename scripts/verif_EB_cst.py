
import sys
import os
import pandas as pd

import numpy as np
from numpy import linalg as LA
import math

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D



filename = sys.argv[1]
ndt=0
X=[]
Y=[]
Z=[]
pX=[]
pY=[]
pZ=[]
EX=[]
EY=[]
EZ=[]
BX=[]
BY=[]
BZ=[]

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
        EX[ndt:]=np.array(data['Ex'].to_numpy())
        EY[ndt:]=np.array(data['Ey'].to_numpy())
        EZ[ndt:]=np.array(data['Ez'].to_numpy())
        BX[ndt:]=np.array(data['Bx'].to_numpy())
        BY[ndt:]=np.array(data['By'].to_numpy())
        BZ[ndt:]=np.array(data['Bz'].to_numpy())

iterations=np.arange(0,ndt)
dt=0.1

####################
#     THEORIC
####################

q=-1.
#m=1836.125
m=1

X0=np.array([X[0], Y[0], Z[0]])
normX0=LA.norm(X0)
print("X0: {}".format(X0))
print("Norme de X0: {}".format(normX0))

p0=np.array([pX[0], pY[0], pZ[0]])
normp0=LA.norm(p0)
normp0xy = np.sqrt(pY[0]**2 + pX[0]**2)
print("p0: {}".format(p0))
print("Norme de p0: {}".format(normp0))

E0=np.array([EX[0], EY[0], EZ[0]])
normE0=LA.norm(E0)
print("E0: {}".format(E0))
print("Norme de E0: {}".format(normE0))

B0=np.array([BX[0], BY[0], BZ[0]])
normB0=LA.norm(B0)
print("B0: {}".format(B0))
print("Norme de B0: {}".format(normB0))

Rth = m*normp0/(abs(q)*normB0)
print("R theoric0: {}".format(Rth))

xth0 = X[0] - Rth
yth0 = Y[0]
zth0 = Z[0]

gamma0=np.sqrt(1+np.sum((p0)**2))
print("Gamma0, Lorentz: {}".format(gamma0))
w=(abs(q)/m)*(normB0/gamma0)
print("Omega, angular frequency: {}".format(w))
period=(2*math.pi)/w
print("Periode: {}".format(period))


Fl0=q*E0[2]


normp=np.zeros(ndt)
R=np.zeros(ndt)
pth=np.zeros([3,ndt])
xth=np.zeros(ndt)
yth=np.zeros(ndt)

for t,it in enumerate(iterations):

    pth[:,it] = [p0[0]-normp0xy*np.sin(w*t*dt),
               p0[1]+normp0xy*np.cos(w*t*dt),
               p0[2] + Fl0*dt*t]

    gamma = np.sqrt(1 + np.sum(pth[:,it]**2))

    w = (abs(q)/m)*(normB0/gamma) 
    
    normp[it] = LA.norm(np.array([pX[it],pY[it],pZ[it]]))
    
    xth[it] = xth0 + Rth * math.cos(w*t*dt)
    yth[it] = yth0 + Rth * math.sin(w*t*dt)

error_abs = np.sqrt((X - xth)**2 + (Y - yth)**2)

###########

xlim = [xth0 - 1.5*Rth, xth0 + 1.5*Rth]
ylim = [yth0 - 1.5*Rth, yth0 + 1.5*Rth]
zlim = [zth0 - 1.5*Rth, zth0 + 1.5*Rth]



### Moment
fig1, axs = plt.subplots(1,3)
axs[0].plot(iterations, pX, marker='+', label='Compute')
axs[0].plot(iterations, pth[0], label='Theoric')
axs[0].set_xlabel("t")
axs[0].set_ylabel("px")
axs[0].legend()
axs[1].plot(iterations, pY, marker='+', label='Compute')
axs[1].plot(iterations, pth[1], label='Theoric')
axs[1].set_xlabel("t")
axs[1].set_ylabel("py")
axs[1].legend()
axs[2].plot(iterations, pZ, marker='+', label='Compute')
axs[2].plot(iterations, pth[2], label='Theoric')
axs[2].set_xlabel("t")
axs[2].set_ylabel("pz")
axs[2].legend()


### Position
fig2, axs = plt.subplots(1,3)
axs[0].plot(iterations, X, marker='+', label='Compute')
#axs[0].plot(iterations, xth, label='Theoric')
axs[0].set_xlabel("t")
axs[0].set_ylabel("x")
axs[0].legend()
axs[1].plot(iterations, Y, marker='+', label='Compute')
#axs[1].plot(iterations, yth, label='Theoric')
axs[1].set_xlabel("t")
axs[1].set_ylabel("y")
axs[1].legend()
axs[2].plot(iterations, Z, marker='+', label='Compute')
#axs[2].plot(iterations, zth, label='Theoric')
axs[2].set_xlabel("t")
axs[2].set_ylabel("z")
axs[2].legend()


### Position 3D
fig3 = plt.figure()
ax = fig3.add_subplot(111, projection='3d')
ax.plot(X, Y, Z, marker='+', label='Compute path of particle')
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
ax.legend()


plt.show()
