#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.animation as animation
import sys

max_num = int(sys.argv[1])
fig = plt.figure()

img= mpimg.imread("frame_0000.png")
im = plt.imshow(img)

global i
i = 1

global max_num
# max_num = 21 # Your number of png files

def updatefig(*args):
    global i
    if i == max_num:
    i = 0
    name = "frame_%04d.png" % i
    i += 1
    img= mpimg.imread(name)
    im = plt.imshow(img)
    return im,

ani = animation.FuncAnimation(fig, updatefig, interval=50, blit=True)
plt.show()
