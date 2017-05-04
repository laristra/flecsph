import sys
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

if len(sys.argv) is not 3:
    print "Error in arguments : <data file name beginning> <number of images>"
    sys.exit(2)

nameFiles = sys.argv[1]
nImages = int(sys.argv[2])

print "File: %s\n nImages: %d" % (nameFiles,nImages)

rangePlot = 4

# Same x vector for each value
x = np.zeros(nImages)

for i in range(0,nImages):
    name = nameFiles+str(i).zfill(5)+".txt"
    posX = np.loadtxt(name,usecols=(0,))
    posY = np.loadtxt(name,usecols=(1,))
    posZ = np.loadtxt(name,usecols=(2,))
    fig = plt.figure(figsize=[40,40])
    ax = fig.add_subplot(1,2,1,projection='3d')
    ax.set_xlim([-rangePlot,rangePlot])
    ax.set_ylim([-rangePlot,rangePlot])
    ax.set_zlim([-rangePlot,rangePlot])
    ax.view_init(elev=90,azim=0)
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')
    ax.plot(posX,posY,posZ,".")
    
    ax = fig.add_subplot(1,2,2,projection='3d')
    ax.set_xlim([-rangePlot,rangePlot])
    ax.set_ylim([-rangePlot,rangePlot])
    ax.set_zlim([-rangePlot,rangePlot])
    ax.view_init(elev=90,azim=0)
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')
    ax.view_init(elev=0,azim=90)
    ax.plot(posX,posY,posZ,".")
    nameOutput = nameFiles+"3D_"+str(i).zfill(5)+".png"
    plt.savefig(nameOutput)
    print "\r"+str(i)+"/"+str(nImages)+" Images Generated",
    sys.stdout.flush()
    plt.close(fig)
