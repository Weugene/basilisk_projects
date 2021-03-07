import matplotlib.pyplot as plt
import numpy as np
import stream as bas
from math import *

N = 256
bas.init_grid(N)

a = bas.scalar()
b = bas.scalar()

a.f = lambda x,y: 0.
b.f = lambda x,y: sin(2.*pi*x)*cos(2.*pi*y)


bas.poisson(a,b)

x = np.linspace(0, 1, N)
y = np.linspace(0, 1, N)
X,Y = np.meshgrid(x,y)
plt.imshow(a.f(X,Y))
plt.show()
