from sklearn.neighbors import KDTree
import numpy as np
import networkx as nx

def find_longest_edge(l):
    e1 = G[l[0]][l[1]]['weight']
    e2 = G[l[0]][l[2]]['weight']
    e3 = G[l[1]][l[2]]['weight']
    if e2 < e1 > e3:
        return (l[0], l[1])
    elif e1 < e2 > e3:
        return (l[0], l[2])
    elif e1 < e3 > e2:
        return (l[1], l[2])

G = nx.Graph()  # A graph to hold the nearest neighbours

#X = [(0, 1), (1, 1), (3, 2), (5, 4)]  # Some list of points in 2D
X = [(0, 1), (0, 0), (2, 1),  (3, 2),  (9, 4), (5, 4)]
tree = KDTree(X, leaf_size=2, metric='euclidean')  # Create a distance tree

# Now loop over your points and find the two nearest neighbours
# If the first and last points are also the start and end points of the line you can use X[1:-1]
for p in X:
    #print(p)
    dist, ind = tree.query(p, k=3)
    print (ind)

    # ind Indexes represent nodes on a graph
    # Two nearest points are at indexes 1 and 2. 
    # Use these to form edges on graph
    # p is the current point in the list
    G.add_node(p)
    n1, l1 = X[ind[0][1]], dist[0][1]  # The next nearest point
    n2, l2 = X[ind[0][2]], dist[0][2]  # The following nearest point  
    G.add_edge(p, n1)
    G.add_edge(p, n2)


print (G.edges())  # A list of all the connections between points
print (nx.shortest_path(G, source=(0,1), target=(5,4)))

end_cliques = [i for i in list(nx.find_cliques(G)) if len(i) == 3]
edge_lengths = [find_longest_edge(i) for i in end_cliques]
G.remove_edges_from(edge_lengths)
edges = G.edges()

start_end = [n for n,nbrs in G.adjacency_iter() if len(nbrs.keys()) == 1]
print (nx.shortest_path(G, source=start_end[0], target=start_end[1]))
