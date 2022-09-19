import matplotlib.pyplot as plt
import numpy as np
from sklearn.neighbors import NearestNeighbors
import networkx as nx
List = [(-50.6261, 74.3683), (-63.2489, 75.0038), (-76.0384, 75.6219), (-79.8451, 75.7855), (-30.9626, 168.085), (-27.381, 170.967), (-22.9191, 172.928), (-16.5869, 173.087), (-4.813, 172.505), (-109.056, 92.0063), (-96.0705, 91.4232), (-83.255, 90.8563), (-80.7807, 90.7498), (-54.1694, 89.5087), (-41.6419, 88.9191), (-32.527, 88.7737), (-27.6403, 91.0134), (-22.3035, 95.141), (-18.0168, 100.473), (-15.3918, 105.542), (-13.6401, 112.373), (-13.3475, 118.988), (-14.4509, 125.238), (-17.1246, 131.895), (-21.6766, 139.821), (-28.5735, 149.98), (-33.395, 156.344), (-114.702, 83.9644), (-114.964, 87.4599), (-114.328, 89.8325), (-112.314, 91.6144), (-109.546, 92.0209), (-67.9644, 90.179), (-55.2013, 89.5624), (-34.4271, 158.876), (-34.6987, 161.896), (-33.6055, 164.993), (-87.0365, 75.9683), (-99.8007, 76.0889), (-105.291, 76.5448), (-109.558, 77.3525), (-112.516, 79.2509), (-113.972, 81.3335), (2.30014, 171.635), (4.40918, 169.691), (5.07165, 166.974), (5.34843, 163.817), (5.30879, 161.798), (-29.6746, 73.5082), (-42.5876, 74.0206)]

List = [[x[0], x[1]] for x in List]
def distance(P1, P2):
    """
    This function computes the distance between 2 points defined by
     P1 = (x1,y1) and P2 = (x2,y2)
    """

    return ((P1[0] - P2[0])**2 + (P1[1] - P2[1])**2) ** 0.5


def optimized_path(coords, start=None):
    """
    This function finds the nearest point to a point
    coords should be a list in this format coords = [ [x1, y1], [x2, y2] , ...]

    """
    if start is None:
        start = coords[0]
    pass_by = coords
    path = [start]
    print(pass_by)
    print(start)
    pass_by.remove(start)
    while pass_by:
        nearest = min(pass_by, key=lambda x: distance(path[-1], x))
        path.append(nearest)
        pass_by.remove(nearest)
    return np.asarray(path)

# define a start point
x0, y0 = (-29.6746, 73.5082)
start = [x0, y0]

path = optimized_path(List,start)

plt.plot(path[:,0], path[:,1])
plt.show()