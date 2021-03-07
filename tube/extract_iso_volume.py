# state file generated using paraview version 5.8.0

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# trace generated using paraview version 5.8.0
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *
#from paraview.vtk.util import numpy_support # provides unique()
#from vtkmodules.numpy_interface.algorithms import * # provides volume()

import numpy as np
import glob, os, sys
import logging
import timeit
import argparse
logging.basicConfig(format='%(message)s')
log = logging.getLogger(__name__)
import matplotlib.pyplot as plt
from scipy.signal import butter,filtfilt # for fft filtering
from scipy.ndimage import gaussian_filter1d # gaussian filtering
from scipy.signal import savgol_filter
import json
from scipy.spatial import Delaunay

import functools
#import __builtin__
# from inspect import getmodule
#
# print(getmodule(Show))
vtk_from_pvpython=True # pvpython reads from file, otherwise from paraview GUI
# vtk_from_pvpython=False # pvpython reads from file, otherwise from paraview GUI

from matplotlib.pyplot import *
from scipy import interpolate

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate, signal, ndimage
from sklearn.decomposition import PCA
from sklearn.utils import shuffle
from sklearn.cluster import KMeans #, SpectralClustering, AgglomerativeClustering
from sklearn.neighbors import NearestNeighbors
import networkx as nx
import statsmodels.api as sm
import random

import itertools
import operator
from scipy.signal import find_peaks
from tsmoothie.smoother import *
from fractions import Fraction #convert a decimal number into fraction

def my_custom_timer(func):
    @functools.wraps(func)
    def wrapper_timer(*args, **kwargs):
        start_time = timeit.default_timer()
        s = ''
        for arg_name in kwargs:
                s += arg_name
        print("Started {!r}({!r}) {!r} ".format(func.__name__, s, ' line=' + str(sys._getframe().f_back.f_lineno))),
        value = func(*args, **kwargs)
        end_time = timeit.default_timer()
        run_time = end_time - start_time
        print("Finished in {:.4f} seconds.".format(run_time))
        return value
    return wrapper_timer

# print(sys.modules.keys())



@my_custom_timer
def XMLPartitionedUnstructuredGridReader(*args, **kwargs):
    return paraview.simple.XMLPartitionedUnstructuredGridReader(*args, **kwargs)

@my_custom_timer
def PVDReader(*args, **kwargs):
    return paraview.simple.PVDReader(*args, **kwargs)

@my_custom_timer
def GetActiveSource(*args, **kwargs):
    return paraview.simple.GetActiveSource(*args, **kwargs)

@my_custom_timer
def Slice(*args, **kwargs):
    return paraview.simple.Slice(*args, **kwargs)

@my_custom_timer
def Calculator(*args, **kwargs):
    return paraview.simple.Calculator(*args, **kwargs)

@my_custom_timer
def Clip(*args, **kwargs):
    return paraview.simple.Clip(*args, **kwargs)

@my_custom_timer
def PassArrays(*args, **kwargs):
    return paraview.simple.PassArrays(*args, **kwargs)

@my_custom_timer
def Fetch(*args, **kwargs):
    return paraview.servermanager.Fetch(*args, **kwargs)

@my_custom_timer
def CellDatatoPointData(*args, **kwargs):
    return paraview.simple.CellDatatoPointData(*args, **kwargs)

@my_custom_timer
def ResampleToImage(*args, **kwargs):
    return paraview.simple.ResampleToImage(*args, **kwargs)

@my_custom_timer
def Contour(*args, **kwargs):
    return paraview.simple.Contour(*args, **kwargs)

@my_custom_timer
def IsoVolume(*args, **kwargs):
    return paraview.simple.IsoVolume(*args, **kwargs)

@my_custom_timer
def StreamTracer(*args, **kwargs):
    return paraview.simple.StreamTracer(*args, **kwargs)

@my_custom_timer
def ExtractSelection(*args, **kwargs):
    return paraview.simple.ExtractSelection(*args, **kwargs)

@my_custom_timer
def Show(*args, **kwargs):
    return paraview.simple.Show(*args, **kwargs)

@my_custom_timer
def SaveData(*args, **kwargs):
    return paraview.simple.SaveData(*args, **kwargs)
import plotly.graph_objects as go
from pathlib import Path

def plot_graph(list_x, list_y, names, xtitle, ytitle, image_name, list_x_fill=[], list_y_fill=[], mode=[],\
			   dash=['solid', 'dot', 'dash', 'longdash'],\
			   colors=['blue', 'red', 'hsv(120,100,100)', 'green', 'black' ],\
			   marker_size=15, xrange =[], yrange = [], \
			   marker_style = ['circle', 'triangle-up', 'triangle-down','square', 'diamond', 'cross',  'x-thin', 'cross-thin' ],\
			   width=1000, height=500, path='./', yanchor='center', y0_anchor=0.01, xanchor='left', x0_anchor=0.3):
	if mode == []:
		for i in range(len(list_x)):
			mode.append('lines+markers')

	
	while len(marker_style) < len(list_x):
		marker_style[:] = marker_style[:] + marker_style[:]
	figborderlinesize = 0.7
	legborderlinesize = 0.7
	yaxis = dict(
		tickfont = dict(
			family = 'Times New Roman',
			size = 20,
			color = 'black'
		),
		titlefont = dict(
			family = 'Times New Roman',
			size = 25,
			color = 'black'
		),
	)
	xaxis = dict(
		tickfont = dict(
			family = 'Times New Roman',
			size = 20,
			color = 'black'
		),
		titlefont = dict(
			family = 'Times New Roman',
			size = 25,
			color = 'black'
		)
	)

	axis_style = dict(showline=True, gridwidth=1, gridcolor='lightgrey', linewidth=figborderlinesize, linecolor='black', mirror=True, ticks='outside', tickfont = dict(family = 'Times New Roman', size = 20, color = 'black'))
	bg_style = {'plot_bgcolor': 'rgba(255, 255, 255, 1)', 'paper_bgcolor': 'rgba(255, 255, 255, 1)',}


	fig = go.Figure()
	k = len(list_x)
	n_fill = len(list_x_fill)
	if len(list_x_fill) == 2 and len(list_y_fill) == 2:
		fig.add_trace(go.Scatter(x=list_x_fill[1], y=list_y_fill[1], name=names[k+1], mode='lines', fillcolor='blueviolet', line_color='blueviolet', fill='tozeroy')) # fill to trace0 y
		fig.add_trace(go.Scatter(x=list_x_fill[0], y=list_y_fill[0], name=names[k], mode='lines', fillcolor='lightsteelblue',     line_color='indigo', fill='tozeroy')) # fill down to xaxis
	for i,x in enumerate(list_x):
		print('Plot curve number:', i)
		y = np.asarray(list_y[i])
		fig.add_trace(go.Scatter(x=x, y=y, name=names[i],
								 mode=mode[i],
								 marker=dict(
									 size=marker_size,
									 line=dict(width=1)
								 ),
								 marker_symbol=marker_style[i],
								 line=dict(width=2, dash=dash[i]),
								 textfont=dict(
									 family="Times New Roman",
									 size=18,
									 color="LightSeaGreen")
								 ))
		if colors != []:
			fig['data'][i + n_fill]['marker']['line']['color'] = colors[i]
			fig['data'][i + n_fill]['line']['color'] = colors[i]
	fig.update_layout(
		width = width,
		height = height,
		xaxis_title=xtitle,
		yaxis_title=ytitle,
		yaxis = yaxis,
		xaxis = xaxis,
		showlegend=True
	)
	fig.update_layout(bg_style)
	fig.update_xaxes(axis_style)
	fig.update_yaxes(axis_style)
	fig.update_layout(legend=dict(
		bgcolor="White",
		bordercolor="Black",
		borderwidth=figborderlinesize
	))
	fig.update_layout(font=dict(
        family="Times New Roman",
        size=20,
        color="Black"
    ))
	fig.update_layout(
		autosize=False,
		margin=dict(
			l=0,
			r=0,
			b=0,
			t=0,
			pad=0.1
		),
		#     paper_bgcolor="LightSteelBlue",
	)
	fig.update_layout(legend=dict(
		yanchor=yanchor,
		y=y0_anchor,
		xanchor=xanchor,
		x=x0_anchor
	))
	if len(xrange) == 2:
		fig.update_xaxes(range=xrange)
	if len(yrange) == 2:
		fig.update_yaxes(range=yrange)
	#fig.show()
	fn = path + image_name
	print('Write image to file:', fn)
	fig.write_image(str(Path(fn)), engine="kaleido")
	print("Successfully generated:", fn)




def eprint(var):
    log.warning(var)


#Sample Period - 5 sec (t) number of 
#Sampling Freq - 30 samples / s , i.e 30 Hz (fs)
#Total Samples - (fs x t) = 150
#Signal Freq = 6 signal / 5 sec = 1.2 Hz
#Nyquist Frequency = 0.5 * fs
#order = Polynomial order of the signal
#cutoff - which part of fs will be removed 
def butter_lowpass_filter(new_u, data, cutoff=0.1, order = 2):
	a, b = new_u[0], new_u[-1]
	timestep = (b - a)/(new_u.shape[0] - 1) # dt of discretization
	fs = 1.0/timestep # frequency of discretization
	nyq = 0.5*fs
	normal_cutoff = cutoff / (nyq * timestep)
	# Get the filter coefficients 
	b, a = butter(order, normal_cutoff, btype='low', analog=False)
	y = filtfilt(b, a, data)#, method="pad", padlen=5
	return y




def toX(x, y):
    return np.concatenate((x.reshape(-1,1), y.reshape(-1,1)), axis=1)

def fromX(X):
    return X[:,0], X[:,1]

def PCAsort(x, y):
    # PCA
    xy = toX(x, y)
    pca = PCA(2).fit(xy)
    xypca = pca.transform(xy)
    newx, newy = fromX(xypca)

    #sort
    indexSort = np.argsort(newx)
    return newx[indexSort], newy[indexSort], pca

def interp_data(x, y, N):
    f = interpolate.interp1d(x, y, kind='linear')
    newx = np.linspace(x[0], x[-1], N)
    return newx, f(newx)

def align_data(x, y, direction):
    ddir = np.hstack((x[-1] - x[0], y[-1] - y[0]));
    return (x, y) if np.dot(direction, ddir) > 0 else (x[::-1], y[::-1])


# Solution from https://stackoverflow.com/a/63368162/2531400
def XYclean(x, y, window=None, N=None):
    newx, newy, pca = PCAsort(x, y)
    newx, newy = interp_data(newx, newy, N)
    newy = signal.savgol_filter(newy, window, 2)
    xyclean = pca.inverse_transform(toX(newx, newy))
    return fromX(xyclean)



def alpha_shape(points, alpha, only_outer=True):
    """
    Compute the alpha shape (concave hull) of a set of points.
    :param points: np.array of shape (n,2) points.
    :param alpha: alpha value.
    :param only_outer: boolean value to specify if we keep only the outer border
    or also inner edges.
    :return: set of (i,j) pairs representing edges of the alpha-shape. (i,j) are
    the indices in the points array.
    """
    assert points.shape[0] > 3, "Need at least four points"

    def add_edge(edges, i, j):
        """
        Add an edge between the i-th and j-th points,
        if not in the list already
        """
        if (i, j) in edges or (j, i) in edges:
            # already added
            assert (j, i) in edges, "Can't go twice over same directed edge right?"
            if only_outer:
                # if both neighboring triangles are in shape, it's not a boundary edge
                edges.remove((j, i))
            return
        edges.add((i, j))

    tri = Delaunay(points)
    edges = set()
    # Loop over triangles:
    # ia, ib, ic = indices of corner points of the triangle
    for ia, ib, ic in tri.vertices:
        pa = points[ia]
        pb = points[ib]
        pc = points[ic]
        # Computing radius of triangle circumcircle
        # www.mathalino.com/reviewer/derivation-of-formulas/derivation-of-formula-for-radius-of-circumcircle
        a = np.sqrt((pa[0] - pb[0]) ** 2 + (pa[1] - pb[1]) ** 2)
        b = np.sqrt((pb[0] - pc[0]) ** 2 + (pb[1] - pc[1]) ** 2)
        c = np.sqrt((pc[0] - pa[0]) ** 2 + (pc[1] - pa[1]) ** 2)
        s = (a + b + c) / 2.0
        area = np.sqrt(s * (s - a) * (s - b) * (s - c))
        circum_r = a * b * c / (4.0 * (area + 1e-16)) #corrected by Weugene
        if circum_r < alpha:
            add_edge(edges, ia, ib)
            add_edge(edges, ib, ic)
            add_edge(edges, ic, ia)
    return edges

def find_min_max_curve(points, alpha, p0, pN):
    # Computing the alpha shape
    edges = alpha_shape(points, alpha=alpha, only_outer=True)
    #order edges
    edges = stitch_boundaries(edges)
    
    edge_points = np.zeros((len(edges),2))
    k=0
    for i, j in edges:
        edge_points[k,:] = points[[i, j], 0][0] , points[[i, j], 1][0]
        k += 1
    inodes, jnodes = zip(*edges)
    min_x_ind = np.argmin(np.linalg.norm(edge_points - p0,axis=1))
    max_x_ind = np.argmin(np.linalg.norm(edge_points - pN,axis=1))
    print("min_x_ind={} max_x_ind={}".format(min_x_ind, max_x_ind))
#     min_x_ind = np.argmin(edge_points[:, 0])
#     max_x_ind = np.argmax(edge_points[:, 0])
    if min_x_ind < max_x_ind:
        lower_hull = edge_points[min_x_ind:max_x_ind+1, :]
        upper_hull = np.concatenate([edge_points[max_x_ind:, :], edge_points[:min_x_ind+1, :]])
    else:
        upper_hull = edge_points[max_x_ind:min_x_ind+1, :]
        lower_hull = np.concatenate([edge_points[min_x_ind:, :], edge_points[:max_x_ind+1, :]])
    return lower_hull, upper_hull


def find_edges_with(i, edge_set):
    i_first = [j for (x,j) in edge_set if x==i]
    i_second = [j for (j,x) in edge_set if x==i]
    return i_first,i_second

def stitch_boundaries(edges):
    edge_set = edges.copy()
    boundary_lst = []
    while len(edge_set) > 0:
        boundary = []
        edge0 = edge_set.pop()
        boundary.append(edge0)
        last_edge = edge0
        while len(edge_set) > 0:
            i,j = last_edge
            j_first, j_second = find_edges_with(j, edge_set)
            if j_first:
                edge_set.remove((j, j_first[0]))
                edge_with_j = (j, j_first[0])
                boundary.append(edge_with_j)
                last_edge = edge_with_j
            elif j_second:
                edge_set.remove((j_second[0], j))
                edge_with_j = (j, j_second[0])  # flip edge rep
                boundary.append(edge_with_j)
                last_edge = edge_with_j

            if edge0[0] == last_edge[1]:
                break

        boundary_lst.append(boundary)
    return boundary_lst[0]


def calc_thickness(x, y, x_peak, x_mean, prefix):
	# calculate actual min max thickness and averaged min max values of thickness
	args = (x >= x_peak) & (x <= x_mean)
	y_ripple_slice = y[args]
	delta_min = 0.5 - y_ripple_slice.max()
	delta_max = 0.5 - y_ripple_slice.min()
	delta_avg = 0.5 - y_ripple_slice.mean()
	delta_avg_std = np.std(y_ripple_slice, ddof=1)
	print("delta_{}".format(prefix), 'delta_min=', delta_min, 'delta_max=', delta_max, 'delta_avg=', delta_avg, 'delta_avg_std=', delta_avg_std, 'NOTE: here avg is calsculated differently')
	return delta_min, delta_max, delta_avg, delta_avg_std

from matplotlib.pyplot import *
from scipy import interpolate
def find_first_peak(x_fil, y_fil, x0, xmin, xmax, x_mean, time):
	if x_fil[0] > x_fil[-1]:
		x_fil = x_fil[::-1]
		y_fil = y_fil[::-1]
	ind = np.argmax(y_fil) 
	print("max ={} {}".format(x_fil[ind], y_fil[ind]))    
	#choose some points if they are:
	args = (y_fil >= 0.3) & (x_fil <= x_mean)
	x_ripple = x_fil[args]
	y_ripple = y_fil[args]

	print("sizes of ripple:{} {}".format(x_ripple.shape, y_ripple.shape))
	#return 0, 0, 0, 0, 0, 0, 0
	#find the first peak in a smoothed curve
	peaks, props = find_peaks(y_ripple, prominence=0.0001)
	print("peaks, props:", peaks, props)
	try:
		x_peak, y_peak = x_ripple[peaks[0]], y_ripple[peaks[0]]
	except:
		x_peak, y_peak = np.inf, np.inf
	length_x_peak_mean = x_mean - x_peak
	print('x_peak candidates=', x_ripple[peaks], 'y_peak candidates=', y_ripple[peaks])
	print('x_peak', x_peak, 'y_peak=', y_peak, "length_x_peak_mean=", length_x_peak_mean)
	
	return x_peak, y_peak, length_x_peak_mean

def pca_sort(x,y, M = 3, window = 100, is_cycle = False,  dcut = 0.05):
    x = np.asarray(x)
    y = np.asarray(y)
    window += 1 if window % 2 == 0 else 0

    xmin = x.min()
    xmax = x.max()
    xcm = 0.5*(xmax + xmin)
    inds = x < xcm
    ind_xy0 = y[inds].argmin()
    xy0 = x[inds][ind_xy0], y[inds][ind_xy0]
    #split points into 3 parts
    cut_x1, cut_x2 = xmin + .08*(xmax-xmin), xmin + .7*(xmax-xmin)


    inds = [(x < cut_x1 + dcut*(xmax-xmin)),
            (x > cut_x1 - dcut*(xmax-xmin)) & (x < cut_x2 + dcut*(xmax-xmin)),
            (x >= cut_x2 - dcut*(xmax-xmin))]

    ccenters = np.zeros((M, 2))
    color = plt.cm.rainbow(np.linspace(0, 1, M))
    for i, c in zip(range(M), color):
        x_, y_ = np.asarray(x[inds[i]]), np.asarray(y[inds[i]])
        ccenters[i,:] = x_.mean(), y_.mean()


    # 3. Hybrid solution: clusterization + 2-NN + PCA
    # Sort clusters: https://stackoverflow.com/a/37744549/2531400
    clf = NearestNeighbors(n_neighbors=2).fit(ccenters)
    G = clf.kneighbors_graph(mode='distance')

    T = nx.from_scipy_sparse_matrix(G)
    i0 = 0  # number of the first cluster
    if not is_cycle:
        min_dist = np.inf
        for i in range(M):
            dist = np.linalg.norm(ccenters[i] - xy0)
            if dist < min_dist:
                i0 = i; min_dist = dist

    order = list(nx.dfs_preorder_nodes(T, i0))

    X, Y = [], []
    for i in range(M):
        j = order[i]
        x_, y_ = np.asarray(x[inds[j]]), np.asarray(y[inds[j]])
        x_, y_, pca = PCAsort(x_, y_)
        x_, y_ = fromX(pca.inverse_transform(toX(x_, y_)))
        if i < M-1:
            x_, y_ = align_data(x_, y_, ccenters[order[i+1]] - ccenters[j])
        else:
            x_, y_ = align_data(x_, y_, ccenters[j] - ccenters[order[i-1]])
        X.append(list(x_)); Y.append(list(y_))
    return sum(X, []), sum(Y, [])
    return np.asarray(X), np.asarray(Y)

def Convolution(x, y, window=100, n_sigma=2):
    smootherx = ConvolutionSmoother(window_len=window, window_type='ones')
    smootherx.smooth(x)
    # generate intervals
    # lowx, upx = smootherx.get_intervals('sigma_interval', n_sigma=2)
    # operate smoothing
    smoothery = ConvolutionSmoother(window_len=window, window_type='ones')
    smoothery.smooth(y)
    # generate intervals
    # lowy, upy = smoothery.get_intervals('sigma_interval', n_sigma=2)
    return smootherx.smooth_data[0], smoothery.smooth_data[0];

def sort_and_smooth(x,y, x_mean, x_peak, M = 3, convolution = False, chebushev = True, Ninterp = 4000, Npolyfit = 50, window = 100, is_cycle = False, cut_1 = 0.08, cut_2 = 0.7, dcut = 0.05):
    window += 1 if window % 2 == 0 else 0
    xmin = x.min()
    xmax = x.max()
    inds = x < x_mean
    ind_xy0 = y[inds].argmin()
    xy0 = x[inds][ind_xy0], y[inds][ind_xy0]
    inds = x > x_mean
    ind_xyN = y[inds].argmin()
    xyN = x[inds][ind_xyN], y[inds][ind_xyN]
    length_x_clip_ends = xyN[0] - xy0[0]
    print('xy0=', xy0, 'xyN=', xyN, 'length=',length_x_clip_ends)
    wdw = max(window, Npolyfit)
    x = np.concatenate([xy0[0] + length_x_clip_ends*np.linspace(-0.01, 0, wdw), x, xyN[0] + length_x_clip_ends*np.linspace(0, 0.01, wdw)])
    y = np.concatenate([xy0[1] - np.linspace(0.1, 0, wdw), y, xyN[1] - np.linspace(0,0.1, wdw)])

    #split points into 3 parts
    dl1 = cut_1*(x_peak - xy0[0])
    dl2 = cut_2*(xmax - xmin) #cut_2*(xyN[0] - x_mean)
    cut_x1, cut_x2 = x_peak - dl1, x_peak + dl2 #    cut_x1, cut_x2 = xmin + cut_1*(xmax-xmin), xmin + cut_2*(xmax-xmin)

    xmm = {
        0: (-np.inf, cut_x1),
        1: (cut_x1, cut_x2),
        2: (cut_x2, np.inf),
    }
    ymm = np.min(y) + 0.01*(np.max(y) - np.min(y)) # to cut edges of the fitting curve
    inds = [(x < cut_x1 + dcut*(xmax-xmin)),
            (x > cut_x1 - dcut*(xmax-xmin)) & (x < cut_x2 + dcut*(xmax-xmin)),
            (x >= cut_x2 - dcut*(xmax-xmin))]
    print("cut_x1={} +/- {}, cut_x2={} +/- {}".format(cut_x1, dcut*(xmax-xmin), cut_x2, dcut*(xmax-xmin)))
    ccenters = np.zeros((M, 2))
    color = plt.cm.rainbow(np.linspace(0, 1, M))
    for i, c in zip(range(M), color):
        x_, y_ = np.asarray(x[inds[i]]), np.asarray(y[inds[i]])
        ccenters[i,:] = x_.mean(), y_.mean()

    # 3. Hybrid solution: clusterization + 2-NN + PCA
    # Sort clusters: https://stackoverflow.com/a/37744549/2531400
    clf = NearestNeighbors(n_neighbors=2).fit(ccenters)
    G = clf.kneighbors_graph(mode='distance')

    T = nx.from_scipy_sparse_matrix(G)
    i0 = 0  # number of the first cluster
    if not is_cycle:
        min_dist = np.inf
        for i in range(M):
            dist = np.linalg.norm(ccenters[i] - xy0)
            if dist < min_dist:
                i0 = i; min_dist = dist

    order = list(nx.dfs_preorder_nodes(T, i0))

    X, Y = [], []
    for i in range(M):
        j = order[i]
        x_, y_ = np.asarray(x[inds[j]]), np.asarray(y[inds[j]])
        # x_, y_ = x[kmeans.labels_ == j], y[kmeans.labels_ == j]
        x_, y_, pca = PCAsort(x_, y_)
        x_, y_ = fromX(pca.inverse_transform(toX(x_, y_)))
        if i < M-1:
            x_, y_ = align_data(x_, y_, ccenters[order[i+1]] - ccenters[j])
        else:
            x_, y_ = align_data(x_, y_, ccenters[j] - ccenters[order[i-1]])
            #print(ccenters[j], ccenters[order[i-1]])
        X.append(x_); Y.append(y_)

    XX, YY = np.array([]), np.array([])
    for i, c in zip(range(M), color):
        x_, y_ = X[i], Y[i]
        pca = PCA(2) #, whiten=True)
        x_, y_ = fromX(pca.fit_transform(toX(x_, y_)))

        xx = np.array([np.min(x_), np.max(x_), np.max(x_), np.min(x_), np.min(x_)])
        yy = np.array([np.min(y_), np.min(y_), np.max(y_), np.max(y_), np.min(y_)])

        print(f"Cluster {i} contains {X[i].size} points within \
            {np.min(x_):.2}<x<{np.max(x_):.2} & {np.min(y_):.2}<y<{np.max(y_):.2}, \
            ratio = {(np.max(x_)-np.min(x_))/(np.max(y_)-np.min(y_)):.3}")

        #'''Chebyshev polynomial fit
        if chebushev:
            p = np.polynomial.Chebyshev.fit(x_, y_, Npolyfit)
            x_ = np.linspace(x_[0], x_[-1], Ninterp//M); y_ = p(x_)
            print('Chebushev is done xminmax={} {}'.format(min(x_), max(x_)))

        if convolution: #Convolution Smoother fit
            N_x_ = len(x_)
            print("N_x_=",N_x_)
            x_, y_ = Convolution(x_, y_, window=N_x_//25, n_sigma=2)

        x_, y_ = fromX(pca.inverse_transform(toX(x_, y_)))

        mask = (x_ > xmm[i][0]) & (x_ < xmm[i][1]) & (y_ > ymm)
        x_, y_ = x_[mask], y_[mask]
        # Fill 1-D arrays
        XX = np.hstack((XX, x_))
        YY = np.hstack((YY, y_))
    return XX, YY, length_x_clip_ends, xy0, xyN, xmin, xmax

def find_smooth_curve_and_bounds(x, y, x_mean, alpha = 0.05, M = 3, convolution = False, chebushev = True, Ninterp = 4000, Npolyfit = 50, window = 100, is_cycle = False, cut_1 = 0.08, cut_2 = 0.7, dcut = 0.05):
    # calculate xmin, xmax
    xmin, xmax = min(x), max(x)
    # calculate xy0, xyN
    inds = x < x_mean
    ind_xy0 = y[inds].argmin()
    xy0 = x[inds][ind_xy0], y[inds][ind_xy0]
    inds = x > x_mean
    ind_xyN = y[inds].argmin()
    xyN = x[inds][ind_xyN], y[inds][ind_xyN]
    length_x_clip_ends = xyN[0] - xy0[0]
    # upper and lower lines
    lower_hull, upper_hull = find_min_max_curve(np.c_[x,y], alpha = alpha, p0=xy0, pN=xyN)
	
    print("lower_hull and upper_hull minmax<><><><>", min(lower_hull[:,0]), max(lower_hull[:,0]), min(upper_hull[:,0]), max(upper_hull[:,0]))

    x_peak, y_peak, length_x_peak_mean = find_first_peak(upper_hull[:,0], upper_hull[:,1], xy0[0], xmin, xmax, x_mean, 0)

    XX, YY, length_x_clip_ends, xy0, xyN, xmin, xmax = sort_and_smooth(x,y, x_mean, x_peak, M = M, convolution = convolution, chebushev = chebushev, Ninterp = Ninterp, Npolyfit = Npolyfit, window = window, is_cycle = is_cycle, cut_1 = cut_1, cut_2 = cut_2, dcut = dcut)



    delta_min_smooth, delta_max_smooth, delta_avg_fil, delta_avg_std_fil = calc_thickness(XX, YY, x_peak, x_mean, 'filtered')
    delta_min_lw, delta_max_lw, delta_avg_lw, delta_avg_std_lw = calc_thickness(lower_hull[:,0], lower_hull[:,1], x_peak, x_mean, 'lower_hull')
    delta_min_up, delta_max_up, delta_avg_up, delta_avg_std_up = calc_thickness(upper_hull[:,0], upper_hull[:,1], x_peak, x_mean, 'upper_hull')
    delta_min, delta_max, delta_avg, delta_avg_std = calc_thickness(x, y, x_peak, x_mean, 'sliced_x_y')
	
    lh_fil_x, lh_fil_y, lh_length, lh_xy0, lh_xyN, lh_xmin, lh_xmax = sort_and_smooth(lower_hull[:,0], lower_hull[:,1], x_mean, x_peak, M = M, convolution = convolution, chebushev = chebushev, Ninterp = Ninterp, Npolyfit = Npolyfit, window = window, is_cycle = is_cycle, cut_1 = cut_1, cut_2 = cut_2, dcut = dcut)
    lower_hull = np.c_[lh_fil_x,lh_fil_y]
    uh_fil_x, uh_fil_y, uh_length, uh_xy0, uh_xyN, uh_xmin, uh_xmax = sort_and_smooth(upper_hull[:,0], upper_hull[:,1], x_mean, x_peak, M = M, convolution = convolution, chebushev = chebushev, Ninterp = Ninterp, Npolyfit = Npolyfit, window = window, is_cycle = is_cycle, cut_1 = cut_1, cut_2 = cut_2, dcut = dcut)
    upper_hull = np.c_[uh_fil_x,uh_fil_y]
	
    return lower_hull, upper_hull, XX, YY, x_peak, y_peak, length_x_peak_mean, delta_min, delta_max, delta_min_smooth, delta_max_smooth, xy0, xyN, xmin, xmax

def Save1DArraysToFile(numpy_arrays, fn):
    lists = []
    with open(fn, 'w') as f:
        for n in range(len(numpy_arrays)):
            lists.append(numpy_arrays[n].tolist())
        f.write(json.dumps(lists))
        print('Successfully save file:', fn)

def passArray(input_data, fn, PointDataArrays=[ 'Points', 'Result'], CellDataArrays=['Volume']):
    # create a new 'Pass Arrays'
    passArrays1 = PassArrays(Input=input_data)
    passArrays1.PointDataArrays = PointDataArrays
    passArrays1.CellDataArrays = CellDataArrays

    # update the view to ensure updated data information
    spreadSheetView1.Update()

    ss_data = paraview.servermanager.Fetch(passArrays1)
    Np = ss_data.GetNumberOfPoints()
    print('Np=', Np)
    xr = []
    for ip in range(Np):
        xp = ss_data.GetPoint(ip)[0]
        rp = ss_data.GetPointData().GetArray('Result').GetValue(ip)
        xr.append((xp,rp))
    xr = np.array(xr)
    print('processing data size of x and y:', len(xr))
    
    Save1DArraysToFile([xr[:,0], xr[:,1]], fn)

    return xr[:,0], xr[:,1]

def gcd(a,b):
    while b: a, b = b, a%b
    return a
def simplify_fraction(x,y):
    return x // gcd(x,y), y // gcd(x,y)
# Read from arguments
#2 pvd file name (by defaults in the first file in a current directory)
#3 pvtu is swithed off by default
# ---------------------------------------------------------------------------------------------------------
start = timeit.default_timer()

parser = argparse.ArgumentParser()
parser._action_groups.pop()
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')

optional.add_argument("-infn", type=str, help="Provide the name of the input paraview files, please",
                    nargs='?', default='*.pvd')
optional.add_argument("-outfn", type=str, help="Provide the name of the output paraview files, please",
                    nargs='?', default='iso_surface')
required.add_argument("-maxlevel", type=int, help="Provide the maximum level of refinement",
                    nargs='?', default=10, required=True)
required.add_argument("-iter", type=int, help="Provide the iter argument level of refinement",
                    nargs='?', default=0, required=True)
required.add_argument("-volumetric_repr", type=int, help="Provide the volumetric_repr argument 1 if volumetric repr, 0 only surface",
                    nargs='?', default=0, required=True)
required.add_argument("-Nslice", type=int, help="Provide the Nslice argument",
                    nargs='?', default=6)


# parser.add_argument('--foo', action='store_const', const=2, default=42)
args = parser.parse_args()
print(args)
infn = args.infn
outfn = args.outfn
maxlevel = args.maxlevel
iter = args.iter
volumetric_repr = args.volumetric_repr
Nslice = args.Nslice

if volumetric_repr:
	outfn = 'iso_volume'
#Current PATH reading
if vtk_from_pvpython:
	path = os.path.abspath(os.getcwd())
else:
	path = '/Users/weugene/basilisk/work/tube/'

eprint("Current PATH=" + path)
if vtk_from_pvpython:
	if infn[-5::] == '.pvtu':
	# Find files with *.pvtu extension
		numbers = []
		file = ""
		for file in glob.glob(infn):
		    numbers.append(int(filter(lambda x: x.isdigit(), file)))
		file=file[0:-9]
		numbers.sort()
		N = len(numbers)
		infn = []
		for i in numbers:
		    infn.append('{}/{}{:04d}.pvtu'.format(path, file, i))
		print(infn)
		fn = path +  '/' + infn
		my_source = XMLPartitionedUnstructuredGridReader(FileName = fn)
	elif infn[-4::] == '.pvd':
	# create a new 'PVD Reader'
		fn = glob.glob(infn)
		print('Found files:',fn)
		id_fn = fn.index("dump2pvd_compressed.pvd")
		fn = fn[id_fn]
		print('Read the first one:',fn)
		my_source = PVDReader(FileName = path +  '/' + fn)
	else:
		eprint('Get Active Source: No pvd or pvtu files are provided')
		my_source = GetActiveSource()
else:
	eprint('Get Active Source')
	my_source = GetActiveSource()

#my_source.CellArrays = ['fs', 'f', 'l', 'l2', 'omega', 'u.x']
#     sys.exit()
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# get animation scene
animationScene1 = GetAnimationScene()

# get the time-keeper
timeKeeper1 = GetTimeKeeper()
timesteps = timeKeeper1.TimestepValues # 0, 0.1, 0.2 ...
print("timesteps=", timesteps)
try:
	NT = len(timesteps)
except:
    NT = 1
    timesteps = [timesteps]

print ("renderViews, axesGrid, SpreadSheetViews, layouts.. "),
# Create a new 'Render View'
renderView1 = CreateView('RenderView')

# Create a new 'SpreadSheet View'
spreadSheetView1 = CreateView('SpreadSheetView')

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------
# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.SplitHorizontal(0, 0.5)
layout1.AssignView(1, renderView1)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)

# defining of computational domains
Show(my_source, renderView1)
boundsDomain = my_source.GetDataInformation().GetBounds()

lDomain = boundsDomain[1] - boundsDomain[0]
print("boundsDomain of cube=", boundsDomain, " lDomain=", lDomain)


SetActiveView(renderView1)
# Properties modi bfied on animationScene1
animationScene1.AnimationTime = timesteps[0]
# Properties modified on timeKeeper1
timeKeeper1.Time = timesteps[0]
# ***************** CUT remain only a TUBE ****************************
# ***************** COARSE ESTIMATION of BORDERS OF BUBBLE ****************************

# create a new 'Iso Volume'
isoVolume1 = IsoVolume(Input=my_source)
isoVolume1.InputScalars = ['CELLS', 'f']
isoVolume1.ThresholdRange = [0.0, 0.5]

# create a new 'Connectivity'
connectivity1 = Connectivity(Input=isoVolume1)
connectivity1.ExtractionMode = 'Extract Largest Region'

Show(connectivity1, renderView1)

bounds = connectivity1.GetDataInformation().GetBounds()
print("bounds_coarse=",bounds)
center = [(bounds[0] + bounds[1])/2, (bounds[2] + bounds[3])/2,(bounds[4] + bounds[5])/2]
print('center_coarse=',center)

shift = 0.2
len_bub = bounds[1] - bounds[0]
len_min = max([bounds[0] - shift, 0])
len_max = min([bounds[1] + shift, lDomain])
length = len_max - len_min

s = "Coarse: len_min: {} len_max: {} len: {} len_bub: {} ".format(len_min, len_max, length, len_bub)
print (s)

# ***************** CLIP A BOX ****************************
# create a new 'Clip'
clip1 = Clip(Input=my_source)
clip1.ClipType = 'Box'
clip1.HyperTreeGridClipper = 'Plane'
clip1.Scalars = ['POINTS', 'f']
clip1.Value = 0.5
# init the 'Box' selected for 'ClipType'
clip1.ClipType.Position = [len_min, -0.6, -0.6]
clip1.ClipType.Length = [len_max, 1.2, 1.2]


# ***************** CELL DATA TO POINT DATA  from BOX ****************************
# create a new 'Cell Data to Point Data'
cellDatatoPointData1 = CellDatatoPointData(Input=clip1)

# ***************** CALCULATE OMEGA ****************************
# create a new 'Calculator' for OMEGA
calculator0 = Calculator(Input=cellDatatoPointData1)
calculator0.ResultArrayName = 'absOmega'
calculator0.Function = 'abs(omega)'

# ****************************************************************
# ***************** REFINED **************************************
# ***************** RESAMPLE TO IMAGE ****************************
# create a new 'Resample To Image'
resampleToImage1 = ResampleToImage(Input=calculator0)
resampleToImage1.UseInputBounds = 0
dx = lDomain/2.0**maxlevel
shift = 0.1
len_min -= shift
len_max -= shift
Nx = int((len_max - len_min)/dx)
Nyz = int(1.0/dx)
resampleDimension = [Nx, Nyz, Nyz]
print("resampleDimension:", resampleDimension)
resampleToImage1.SamplingDimensions = resampleDimension
resampleToImage1.SamplingBounds = [len_min, len_max, -0.5, 0.5, -0.5, 0.5]

# ***************** REFINED CONTOUR1 for VOLUME FRACTION f ****************************
# create a new 'Iso Volume'
isoVolume1 = IsoVolume(Input=resampleToImage1)
isoVolume1.InputScalars = ['POINTS', 'f']
isoVolume1.ThresholdRange = [0.0, 0.5]
print("Isovolume finished")

# create a new 'Connectivity'
connectivity1 = Connectivity(Input=isoVolume1)
connectivity1.ExtractionMode = 'Extract Largest Region'
connectivity1.ColorRegions = 0

Show(connectivity1, renderView1)

bounds = connectivity1.GetDataInformation().GetBounds()
print("bounds_refined=",bounds)
center = [(bounds[0] + bounds[1])/2, (bounds[2] + bounds[3])/2,(bounds[4] + bounds[5])/2]
print('center_refined=',center)

len_bub = bounds[1] - bounds[0]
len_min = max([bounds[0], 0])
len_max = min([bounds[1], lDomain])
length = len_max - len_min
s = "Refined: len_min: {} len_max: {} len: {} len_bub: {} ".format(len_min, len_max, length, len_bub)
print (s)
print("Connectivity1 finished")

# create a new 'Calculator'
calculator1 = Calculator(Input=connectivity1)
calculator1.Function = 'sqrt(coordsY^2+coordsZ^2)'

# ******************* X_mean calculation *****************
# create a new 'Integrate Variables'
integrateVariables1 = IntegrateVariables(Input=calculator1)

# show data from integrateVariables1
integrateVariables1Display = Show(integrateVariables1, renderView1, 'UnstructuredGridRepresentation')

# create a new 'Pass Arrays'
passArrays1 = PassArrays(Input=integrateVariables1)
passArrays1.PointDataArrays = [ 'Points', 'u.x']
passArrays1.CellDataArrays = ['Volume']

# update the view to ensure updated data information
spreadSheetView1.Update()

ss_data = paraview.servermanager.Fetch(passArrays1)
volume = ss_data.GetCellData().GetArray('Volume').GetValue(0)
x_mean = ss_data.GetPoint(0)[0]
u_mean = [ss_data.GetPointData().GetArray('u.x').GetValue(i)/volume for i in range(3)]
s = "refined_x_mean: {} refined_u_mean: {} volume: {}".format(x_mean, u_mean, volume)
print(s)

# create a new 'Extract Surface'
extractSurface1 = ExtractSurface(Input=calculator1)
#Show(extractSurface1, renderView1)
slices_x = []
slices_y = []
slices_x_raw = []
slices_y_raw = []
slices_name = []
slices_mode = []
slice_dash = []
xy0=[0,0]
xyN=[0,0]

#colorlist = ['rgb({:d},{:d}, {:d})'.format((37*i) %255, (31*i) % 255, (29*i) % 255) for i in range(Nslice)]
#

colorlist = """aliceblue, antiquewhite, aqua, aquamarine, azure,
                beige, bisque, black, blanchedalmond, blue,
                blueviolet, brown, burlywood, cadetblue,
                chartreuse, chocolate, coral, cornflowerblue,
                cornsilk, crimson, cyan, darkblue, darkcyan,
                darkgoldenrod, darkgray, darkgrey, darkgreen,
                darkkhaki, darkmagenta, darkolivegreen, darkorange,
                darkorchid, darkred, darksalmon, darkseagreen,
                darkslateblue, darkslategray, darkslategrey,
                darkturquoise, darkviolet, deeppink, deepskyblue,
                dimgray, dimgrey, dodgerblue, firebrick,
                floralwhite, forestgreen, fuchsia, gainsboro,
                ghostwhite, gold, goldenrod, gray, grey, green,
                greenyellow, honeydew, hotpink, indianred, indigo,
                ivory, khaki, lavender, lavenderblush, lawngreen,
                lemonchiffon, lightblue, lightcoral, lightcyan,
                lightgoldenrodyellow, lightgray, lightgrey,
                lightgreen, lightpink, lightsalmon, lightseagreen,
                lightskyblue, lightslategray, lightslategrey,
                lightsteelblue, lightyellow, lime, limegreen,
                linen, magenta, maroon, mediumaquamarine,
                mediumblue, mediumorchid, mediumpurple,
                mediumseagreen, mediumslateblue, mediumspringgreen,
                mediumturquoise, mediumvioletred, midnightblue,
                mintcream, mistyrose, moccasin, navajowhite, navy,
                oldlace, olive, olivedrab, orange, orangered,
                orchid, palegoldenrod, palegreen, paleturquoise,
                palevioletred, papayawhip, peachpuff, peru, pink,
                plum, powderblue, purple, red, rosybrown,
                royalblue, saddlebrown, salmon, sandybrown,
                seagreen, seashell, sienna, silver, skyblue,
                slateblue, slategray, slategrey, snow, springgreen,
                steelblue, tan, teal, thistle, tomato, turquoise,
                violet, wheat, white, whitesmoke, yellow,
                yellowgreen""".replace(",", "").split()
colorlist = shuffle(colorlist, random_state=0, n_samples=Nslice)
colorlist = colorlist + colorlist



# *********************** Find the x coordinate of the first peak x_peak ****************************
x,y = passArray(extractSurface1, fn="r_over_x_total_t={}.csv".format(timesteps[0]), PointDataArrays=[ 'Points', 'Result'], CellDataArrays=['Volume'])

lower_hull, upper_hull, XX, YY, x_peak, y_peak, length_x_peak_mean, delta_min, delta_max, delta_min_smooth, delta_max_smooth, xy0, xyN, xmin, xmax =  find_smooth_curve_and_bounds(x, y, x_mean, alpha = 0.05, M = 3, convolution = False, chebushev = True, Ninterp = 4000, Npolyfit = 50, window = 1000, is_cycle = False, cut_1 = 0.08, cut_2 = 0.7, dcut = 0.05)

ind = (YY > 0.05)
#XX = XX[ind]
#YY = YY[ind]
#XX = np.concatenate([ [xy0[0]], XX, [xyN[0]] ])
#YY = np.concatenate([ [xy0[1]], YY, [xyN[1]] ])
fn = "r_over_x_t={}.csv".format(timesteps[0])
Save1DArraysToFile([x, y, XX, YY, lower_hull, upper_hull], fn)
#read files as below:
#with open(fn, 'r') as f:
#	lists = json.load(f)
#    x, y, XX, YY, lower_hull, upper_hull = np.array(lists[0]), np.array(lists[1]), np.array(lists[2]), np.array(lists[3])

plot_graph([XX, [x_peak, x_peak], [x_mean, x_mean]], [YY, [0,0.5], [0,0.5]], \
		   ["smoothed", 'first peak', 'center of mass', "min edge", "max edge"], \
		   list_x_fill=[lower_hull[:,0], upper_hull[:,0]], list_y_fill=[lower_hull[:,1], upper_hull[:,1]], \
		   dash=['solid', 'dot', 'dot'], \
		   xtitle="x", ytitle="r", image_name=fn[:-3]+'pdf', mode=['lines', 'lines', 'lines'], \
		   colors=['red', 'black', 'black' ], yrange=[0,0.5], xrange=[xy0[0] - 0.1, xyN[0] + 0.1],\
		   marker_size=1, width=1000, height=500, path='./', yanchor='bottom', y0_anchor=0.01, xanchor='left', x0_anchor=0.3)

###
### ****************SLICES********************************
###
for n in range(Nslice):
	alpha = n*2*np.pi/Nslice
	#frac = Fraction(2 * n/Nslice); num = frac.numerator; den = frac.denominator
	num, den = simplify_fraction(2*n, Nslice)
	print('angle frac={}/{}'.format(num, den))
	# create a new 'Slice'
	slice1 = Slice(Input=extractSurface1)
	slice1.SliceType = 'Plane'
	slice1.HyperTreeGridSlicer = 'Plane'
	slice1.SliceOffsetValues = [0.0]

	# init the 'Plane' selected for 'SliceType'
	slice1.SliceType.Origin = [0, 0, 0]
	slice1.SliceType.Normal = [0.0, np.sin(alpha), np.cos(alpha)]

	# init the 'Plane' selected for 'HyperTreeGridSlicer'
	slice1.HyperTreeGridSlicer.Origin = [0, 0, 0]

	# ----------------------------------------------------------------
	# setup the visualization in view 'renderView1'
	# ----------------------------------------------------------------

	# show data from slice1
	slice1Display = Show(slice1, renderView1, 'GeometryRepresentation')

	# create a new 'Clip'
	clip2 = Clip(Input=slice1)
	clip2.ClipType = 'Plane'
	clip2.HyperTreeGridClipper = 'Plane'
	# clip2.Scalars = ['POINTS', 'Result']
	# clip2.Value = 0.23405700000000002
	# Properties modified on clip2.ClipType
	clip2.ClipType.Normal = [0.0, -np.cos(alpha), np.sin(alpha)]

	# init the 'Plane' selected for 'ClipType'
	clip2.ClipType.Origin = [0, 0, 0]

	# init the 'Plane' selected for 'HyperTreeGridClipper'
	clip2.HyperTreeGridClipper.Origin = [0, 0, 0]

	fn="slice_t={}_nalpha_{}.csv".format(timesteps[0], n)
	x, y = passArray(clip2, fn, PointDataArrays=[ 'Points', 'Result'], CellDataArrays=['Volume'])
	
	x_peak_, y_peak_, length_x_peak_mean_ = find_first_peak(x, y, xy0[0], xmin, xmax, x_mean, timesteps[0])
	print("slice n={} x_peak_= {} y_peak_= {} length_x_peak_mean_= {}".format(n, x_peak_, y_peak_, length_x_peak_mean_))
	XX, YY, length_x_ends, xy0, xyN, xmin, xmax = sort_and_smooth(x, y, x_mean, x_peak_, M = 3, convolution = True, chebushev = False, Ninterp = 4000, Npolyfit = 30, window = 50, is_cycle = False, cut_1 = 0.08, cut_2 = 0.7, dcut = 0.05)
	#XX = x; YY = y;
	slices_x_raw.append(x)
	slices_y_raw.append(y)
	slices_x.append(XX)
	slices_y.append(YY)
	
	if num == 0:
		slices_name.append("\u03B1=0")
	elif num == 1:
		slices_name.append("\u03B1=\u03C0/{}".format(den))
	else:
		slices_name.append("\u03B1={}\u03C0/{}".format(num,den))
	slices_mode.append('lines')
	slice_dash.append('solid')
	SaveData(path + '/slice_t={}_nalpha_{}.csv'.format(timesteps[0], n), proxy=clip2,\
		 CellDataArrays=['vtkGhostType'], FieldDataArrays=['TimeValue'], AddTime=1)
	fn = "sliceList_t={}_nalpha_{}.csv".format(timesteps[0], n)
	Save1DArraysToFile([ x, y, XX, YY ], fn)
	Delete(clip2)
	del clip2
	Delete(slice1)
	del slice1

for i in range(len(slices_x_raw)):
	#slices_x.append(slices_x_raw[i])
	#slices_y.append(slices_y_raw[i])
	slices_mode.append('markers')
	slice_dash.append(None)
slices_name[:] = slices_name[:] + slices_name[:]
fn = "slice_t={}.pdf".format(timesteps[0])
# mode=slices_mode, \ #xrange=[xy0[0] - 0.1, xyN[0] + 0.1], \ colors=['blue', 'red', 'green', 'black', 'aqua', 'indigo', 'pink', 'yellow' ]
plot_graph(slices_x, slices_y, slices_name, \
       #    mode =['markers', 'markers', 'markers', 'markers', 'markers'], \
	   mode = slices_mode, \
	   #list_x_fill=[upper_hull[:,0], lower_hull[:,0]], list_y_fill=[upper_hull[:,1], lower_hull[:,1]], \
	   dash=slice_dash, \
	   xtitle="x", ytitle="r", image_name=fn, yrange=[0.0,0.5], \
	   colors=[], \
	   #colors=['blue', 'red', 'green', 'black', 'blue', 'red', 'green', 'black' ], \
	   marker_size=0.5, width=1000, height=500, path='./', yanchor='bottom', y0_anchor=0.01, xanchor='left', x0_anchor=0.4)

# ***************** CLIP A BOX ****************************
# create a new 'Clip'
clip2 = Clip(Input=calculator1)
clip2.ClipType = 'Box'
clip2.HyperTreeGridClipper = 'Plane'
#clip2.Scalars = ['POINTS', 'f']
#clip2.Value = 0.5
# init the 'Box' selected for 'ClipType'
clip2.ClipType.Position = [x_peak, -0.6, -0.6]
clip2.ClipType.Length = [length_x_peak_mean, 1.2, 1.2]

# create a new 'Integrate Variables'
integrateVariables1 = IntegrateVariables(Input=clip2)

# show data from integrateVariables1
integrateVariables1Display = Show(integrateVariables1, renderView1, 'UnstructuredGridRepresentation')

# create a new 'Pass Arrays'
passArrays1 = PassArrays(Input=integrateVariables1)
passArrays1.PointDataArrays = ['Points']
passArrays1.CellDataArrays = ['Volume']

# update the view to ensure updated data information
spreadSheetView1.Update()

ss_data = paraview.servermanager.Fetch(passArrays1)
print('N_clip2=', ss_data.GetNumberOfPoints())
volumeB = ss_data.GetCellData().GetArray('Volume').GetValue(0)
rB = np.sqrt(volumeB/(np.pi*length_x_peak_mean))
delta_mean = 0.5 - rB

print('delta_min=', delta_min, 'delta_mean=', delta_mean, 'delta_max=',delta_max)

#save it in the file
fn = "delta_min_mean_max_over_t.txt"
with open(fn, 'a') as f:
	f.write("{} {} {} {} {} {} {} {} {} {}\n".format(timesteps[0], delta_min, delta_mean, delta_max, delta_min_smooth, delta_max_smooth, xy0[0], xy0[1], xyN[0], xyN[1]))
print('Successfully save file:', fn)

# cut tail
# create a new 'Clip'
clip3 = Clip(Input=extractSurface1)
clip3.ClipType = 'Plane'
clip3.HyperTreeGridClipper = 'Plane'
clip3.Scalars = ['POINTS', 'Result']
clip3.Value = 0.23469079123049602

# init the 'Plane' selected for 'ClipType'
clip3.ClipType.Origin = [x_peak, 0, 0]

# init the 'Plane' selected for 'HyperTreeGridClipper'
clip3.HyperTreeGridClipper.Origin = [x_peak, 0, 0]

# Properties modified on clip1.ClipType
clip3.ClipType.Normal = [1.0, 0.0, 0.0]

# create a new 'Transform'
transform1 = Transform(Input=clip3)
transform1.Transform = 'Transform'

# Properties modified on transform1.Transform
transform1.Transform.Translate = [-x_peak, 0.0, 0.0]

#Show the tail
Show(transform1, renderView1)

# create a new 'Extract Surface'
extractSurface2 = ExtractSurface(Input=transform1)

fn = "{}/res/{}_0_{:04d}.vtp".format(path, outfn + '_tail', iter)
SaveData(fn, proxy=extractSurface2)
print("Saved tail surface data of bubble:", fn)

if volumetric_repr == 1:
	#Show(calculator1, renderView1)
	fn = "{}/res/{}_0_{:04d}.pvtu".format(path, outfn, iter)
	SaveData(fn, proxy=calculator1, UseSubdirectory=0)
	print("Saved volumetric data of bubble:", fn)
else:
	fn = "{}/res/{}_0_{:04d}.vtp".format(path, outfn, iter)
	SaveData(fn, proxy=extractSurface1)
	print("Saved surface data of bubble:", fn)



# Freeing Memory
Delete(clip2)
del clip2
Delete(extractSurface1)
del extractSurface1
Delete(calculator1)
del calculator1
Delete(connectivity1)
del connectivity1
Delete(isoVolume1)
del isoVolume1
Delete(resampleToImage1)
del resampleToImage1
Delete(calculator0)
del calculator0
Delete(cellDatatoPointData1)
del cellDatatoPointData1
Delete(clip1)
del clip1



stop = timeit.default_timer()
total_time = stop - start

# output running time in a nice format.
mins, secs = divmod(total_time, 60)
hours, mins = divmod(mins, 60)


sys.stdout.write("Total running time: %d:%d:%d.\n" % (hours, mins, secs))


# tck, u = splprep([xx, yy], s=0.01)
# new_u = np.linspace(min(u), max(u), 500)
# x_new, y_new = splev(new_u, tck)
# x_new, y_new = xx, yy



# x_filtered = butter_lowpass_filter(new_u, x_new, cutoff=0.3, order=2)
# y_filtered = butter_lowpass_filter(new_u, y_new, cutoff=0.3, order=2)
#
# #add first and last points
# x_filtered = np.concatenate(([x_new[0]], x_filtered, [x_new[-1]]))
# y_filtered = np.concatenate(([y_new[0]], y_filtered, [y_new[-1]]))
# x_filtered = x_new
# y_filtered = y_new
