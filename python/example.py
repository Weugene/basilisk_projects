from __future__ import print_function
import matplotlib.pyplot as plt
import numpy as np
import stream as bas

N = 256
x = np.linspace(0, 1, N)
y = np.linspace(0, 1, N)
X,Y = np.meshgrid(x,y)

def init(i,t):
    bas.omega.f = bas.noise

def graph(i,t):
    print ("t=",t)
    Z = bas.omega.f(X,Y)
    plt.cla()
    plt.imshow(Z)
    plt.pause(0.0001)

bas.init_grid(N)

bas.event(init, t = 0.)
plt.ion()
plt.figure()
bas.event(graph, t = range(0,1000,10))

bas.run()
