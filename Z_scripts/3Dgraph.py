#!/usr/local/bin/python3
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv("results.csv")
print(df)

#with open('results.csv', newline='') as csvfile:
#     spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')

plt.rcParams["figure.figsize"] = 12.8, 9.6
ax = plt.axes(projection='3d')


N = 100
half_N = N // 2
X2, Y2 = np.meshgrid(range(N), range(N))
Z2 = X2*Y2

#X1 = np.reshape(X2, -1)
#Y1 = np.reshape(Y2, -1)
#Z1 = np.reshape(Z2, -1)

ax = plt.axes(projection='3d')
ax.plot_wireframe(X2, Y2, Z2, color='r')
plt.show()
