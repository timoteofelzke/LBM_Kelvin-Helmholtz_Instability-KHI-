# -*- coding: utf-8 -*-
# Autor: Diogo Nardelli Siebert
# Data: 01/09/2016
#
# Faz figura das redes utilizadas em lattice boltzmann
# para integrar a documentação do código gradlbm

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

from mpl_toolkits.mplot3d import proj3d

class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

cx =  np.array( [  0,  1, -1,  0,  0,  0,  0,  1, -1, -1,  1,  1, -1,  1, -1  ] )
cy =  np.array( [  0,  0,  0,  1, -1,  0,  0,  1, -1,  1, -1, -1,  1,  1, -1  ] )
cz =  np.array( [  0,  0,  0,  0,  0,  1, -1,  1, -1,  1, -1,  1, -1, -1,  1  ] ) 
3
ax.set_xlim3d(-1,1)
ax.set_ylim3d(-1,1)
ax.set_zlim3d(-1,1)
ax.set_xlabel("$x$")
ax.set_ylabel("$y$")
ax.set_zlabel("$z$")
ax.set_xticks([-1,0,1])
ax.set_yticks([-1,0,1])
ax.set_zticks([-1,0,1])
ax.set_aspect(1)

for i in range(1,len(cx)):
    x,y,z = cx[i], cy[i], cz[i]
    a = Arrow3D([0,x],[0,y],[0,z], mutation_scale=20, lw=1, arrowstyle="-|>", color="k")
    ax.add_artist(a)
    ax.text(1.10*x, 1.10*y, 1.10*z, str(i), color='red')

plt.savefig("d3q15.png", dpi = 100)


##for x,y in zip(cx,cy):
#hw = 0.1
#hl = 0.1
#for x,y in zip(cx[1:],cy[1:]):
#    ax.arrow(0, 0, (1-hl)*x, (1-hl)*y  , head_width=0.1, head_length=0.1, fc='k', ec='k')
#plt.show()