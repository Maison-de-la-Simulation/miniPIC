
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation

import sys
import os
import glob

if len( sys.argv ) == 1:
    print( "Usage: python3 print_particles.py filename.data" )
    exit()

filename = sys.argv[1]
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ndt=0
#X=[]
Y=[]
Z=[]
li=[]
## If directorie make animation with paticles positions
if os.path.isdir(filename):
    f_list=glob.glob(filename+'/particles.0.*.data')
    #f_list=os.listdir(filename)
    f_list.sort()
    #print(f_list)
    for f in f_list:
        data = pd.read_csv(f)
        df=pd.read_csv(f, index_col=None, header=0)
        li.append(df)
        #X[ndt:][:]=np.array(data['posx'].to_numpy())
        Y[ndt:][:]=np.array(data['posy'].to_numpy())
        Z[ndt:][:]=np.array(data['posz'].to_numpy())
        ndt+=1

    df = pd.concat((pd.read_csv(f) for f in f_list), ignore_index=True)
    X=np.array(df['posx'][0].to_numpy())
    X=np.array(df['posx'][1].to_numpy())
    print(X)
    ax.scatter(X, Y, Z)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.set_xlim([0,1])
    ax.set_ylim([0,1])
    ax.set_zlim([0,1])

        
## If unique file print particles positions
if os.path.isfile(filename):
    data = pd.read_csv(filename)
    
    ax.scatter(data.posx, data.posy, data.posz, c=data.species)
    ax.set_xlim([0,1])
    ax.set_ylim([0,1])
    ax.set_zlim([0,1])


plt.show()
