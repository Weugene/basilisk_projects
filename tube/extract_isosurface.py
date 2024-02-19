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
# from paraview.vtk.util import numpy_support # provides unique()
# from vtkmodules.numpy_interface.algorithms import * # provides volume()
import vtk
import numpy as np
import glob, os, sys
import logging
import timeit
import argparse
import re


logging.basicConfig(format='%(message)s')
log = logging.getLogger(__name__)
import matplotlib.pyplot as plt
# from scipy.signal import butter,filtfilt # for fft filtering
# from scipy.ndimage import gaussian_filter1d # gaussian filtering
# from scipy.signal import savgol_filter
import json
from scipy.spatial import Delaunay

import functools

# import __builtin__
# from inspect import getmodule
#
# print(getmodule(Show))
vtk_from_pvpython = True  # pvpython reads from file, otherwise from paraview GUI
# vtk_from_pvpython=False # pvpython reads from file, otherwise from paraview GUI

from matplotlib.pyplot import *
from scipy import interpolate

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate, ndimage

import networkx as nx
import statsmodels.api as sm
import random


import itertools
import operator
from scipy.signal import find_peaks
from tsmoothie.smoother import *
from fractions import Fraction  # convert a decimal number into fraction
import plotly.graph_objects as go
from pathlib import Path

from matplotlib.pyplot import *
from scipy import interpolate


def my_custom_timer(func):
    @functools.wraps(func)
    def wrapper_timer(*args, **kwargs):
        start_time = timeit.default_timer()
        print("Started {!r} {!r} ".format(func.__name__, ' line=' + str(sys._getframe().f_back.f_lineno))),
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
def HyperTreeGridToDualGrid(*args, **kwargs):
    return paraview.simple.HyperTreeGridToDualGrid(*args, **kwargs)


@my_custom_timer
def Connectivity(*args, **kwargs):
    return paraview.simple.Connectivity(*args, **kwargs)


@my_custom_timer
def Threshold(*args, **kwargs):
    return paraview.simple.Threshold(*args, **kwargs)


@my_custom_timer
def Cylinder(*args, **kwargs):
    return paraview.simple.Cylinder(*args, **kwargs)


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
def ExtractSurface(*args, **kwargs):
    return paraview.simple.ExtractSurface(*args, **kwargs)


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
def IntegrateVariables(*args, **kwargs):
    return paraview.simple.IntegrateVariables(*args, **kwargs)


@my_custom_timer
def Show(*args, **kwargs):
    return paraview.simple.Show(*args, **kwargs)


@my_custom_timer
def Hide(*args, **kwargs):
    return paraview.simple.Hide(*args, **kwargs)


@my_custom_timer
def GetDisplayProperties(*args, **kwargs):
    return paraview.simple.GetDisplayProperties(*args, **kwargs)


@my_custom_timer
def GetColorTransferFunction(*args, **kwargs):
    return paraview.simple.GetColorTransferFunction(*args, **kwargs)


@my_custom_timer
def GetOpacityTransferFunction(*args, **kwargs):
    return paraview.simple.GetOpacityTransferFunction(*args, **kwargs)


@my_custom_timer
def ColorBy(*args, **kwargs):
    return paraview.simple.ColorBy(*args, **kwargs)


@my_custom_timer
def GetScalarBar(*args, **kwargs):
    return paraview.simple.GetScalarBar(*args, **kwargs)


@my_custom_timer
def GetMaterialLibrary(*args, **kwargs):
    return paraview.simple.GetMaterialLibrary(*args, **kwargs)


@my_custom_timer
def CreateView(*args, **kwargs):
    return paraview.simple.CreateView(*args, **kwargs)


@my_custom_timer
def CreateLayout(*args, **kwargs):
    return paraview.simple.CreateLayout(*args, **kwargs)


@my_custom_timer
def GetAnimationScene(*args, **kwargs):
    return paraview.simple.GetAnimationScene(*args, **kwargs)


@my_custom_timer
def GetTimeKeeper(*args, **kwargs):
    return paraview.simple.GetTimeKeeper(*args, **kwargs)


@my_custom_timer
def SetActiveView(*args, **kwargs):
    return paraview.simple.SetActiveView(*args, **kwargs)


@my_custom_timer
def FindSource(*args, **kwargs):
    return paraview.simple.FindSource(*args, **kwargs)


@my_custom_timer
def SaveData(*args, **kwargs):
    return paraview.simple.SaveData(*args, **kwargs)


@my_custom_timer
def SaveScreenshot(*args, **kwargs):
    return paraview.simple.SaveScreenshot(*args, **kwargs)


@my_custom_timer
def Delete(*args, **kwargs):
    return paraview.simple.Delete(*args, **kwargs)


def plot_graph(list_x, list_y, names, xtitle, ytitle, image_name, list_x_fill=[], list_y_fill=[], mode=[], \
               dash=['solid', 'dot', 'dash', 'longdash'], \
               colors=['blue', 'red', 'hsv(120,100,100)', 'green', 'black'], \
               marker_size=15, xrange=[], yrange=[], \
               marker_style=['circle', 'triangle-up', 'triangle-down', 'square', 'diamond', 'cross', 'x-thin',
                             'cross-thin'], \
               width=1000, height=500, path='./', yanchor='center', y0_anchor=0.01, xanchor='left', x0_anchor=0.3):
    if mode == []:
        for i in range(len(list_x)):
            mode.append('lines+markers')

    while len(marker_style) < len(list_x):
        marker_style[:] = marker_style[:] + marker_style[:]
    figborderlinesize = 0.7
    legborderlinesize = 0.7
    yaxis = dict(
        tickfont=dict(
            family='Times New Roman',
            size=20,
            color='black'
        ),
        titlefont=dict(
            family='Times New Roman',
            size=25,
            color='black'
        ),
    )
    xaxis = dict(
        tickfont=dict(
            family='Times New Roman',
            size=20,
            color='black'
        ),
        titlefont=dict(
            family='Times New Roman',
            size=25,
            color='black'
        )
    )

    axis_style = dict(showline=True, gridwidth=1, gridcolor='lightgrey', linewidth=figborderlinesize, linecolor='black',
                      mirror=True, ticks='outside', tickfont=dict(family='Times New Roman', size=20, color='black'))
    bg_style = {'plot_bgcolor': 'rgba(255, 255, 255, 1)', 'paper_bgcolor': 'rgba(255, 255, 255, 1)', }

    fig = go.Figure()
    k = len(list_x)
    n_fill = len(list_x_fill)
    if len(list_x_fill) == 2 and len(list_y_fill) == 2:
        fig.add_trace(
            go.Scatter(x=list_x_fill[1], y=list_y_fill[1], name=names[k + 1], mode='lines', fillcolor='blueviolet',
                       line_color='blueviolet', fill='tozeroy'))  # fill to trace0 y
        fig.add_trace(
            go.Scatter(x=list_x_fill[0], y=list_y_fill[0], name=names[k], mode='lines', fillcolor='lightsteelblue',
                       line_color='indigo', fill='tozeroy'))  # fill down to xaxis
    for i, x in enumerate(list_x):
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
        width=width,
        height=height,
        xaxis_title=xtitle,
        yaxis_title=ytitle,
        yaxis=yaxis,
        xaxis=xaxis,
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
    # fig.show()
    fn = path + image_name
    print('Write image to file:', fn)
    fig.write_image(str(Path(fn)), engine="kaleido")
    print("Successfully generated:", fn)


def eprint(var):
    log.warning(var)


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
    for ia, ib, ic in tri.simplices:
        pa = points[ia]
        pb = points[ib]
        pc = points[ic]
        # Computing radius of triangle circumcircle
        # www.mathalino.com/reviewer/derivation-of-formulas/derivation-of-formula-for-radius-of-circumcircle
        a = np.sqrt((pa[0] - pb[0]) ** 2 + (pa[1] - pb[1]) ** 2)
        b = np.sqrt((pb[0] - pc[0]) ** 2 + (pb[1] - pc[1]) ** 2)
        c = np.sqrt((pc[0] - pa[0]) ** 2 + (pc[1] - pa[1]) ** 2)
        s = (a + b + c) / 2.0
        area = np.sqrt(s * np.abs(s - a) * np.abs(s - b) * np.abs(s - c))
        circum_r = a * b * c / (4.0 * (area + 1e-16))  # corrected by Weugene
        if circum_r < alpha:
            add_edge(edges, ia, ib)
            add_edge(edges, ib, ic)
            add_edge(edges, ic, ia)
    return edges


def find_min_max_curve(points, alpha, p0, pN):
    # Computing the alpha shape
    edges = alpha_shape(points, alpha=alpha, only_outer=True)
    # order edges
    edges = stitch_boundaries(edges)

    edge_points = np.zeros((len(edges), 2))
    k = 0
    for i, j in edges:
        edge_points[k, :] = points[[i, j], 0][0], points[[i, j], 1][0]
        k += 1
    inodes, jnodes = zip(*edges)
    min_x_ind = np.argmin(np.linalg.norm(edge_points - p0, axis=1))
    max_x_ind = np.argmin(np.linalg.norm(edge_points - pN, axis=1))
    print("min_x_ind={} max_x_ind={}".format(min_x_ind, max_x_ind))
    #     min_x_ind = np.argmin(edge_points[:, 0])
    #     max_x_ind = np.argmax(edge_points[:, 0])
    if min_x_ind < max_x_ind:
        lower_hull = edge_points[min_x_ind:max_x_ind + 1, :]
        upper_hull = np.concatenate([edge_points[max_x_ind:, :], edge_points[:min_x_ind + 1, :]])
    else:
        upper_hull = edge_points[max_x_ind:min_x_ind + 1, :]
        lower_hull = np.concatenate([edge_points[min_x_ind:, :], edge_points[:max_x_ind + 1, :]])
    return lower_hull, upper_hull


def find_edges_with(i, edge_set):
    i_first = [j for (x, j) in edge_set if x == i]
    i_second = [j for (j, x) in edge_set if x == i]
    return i_first, i_second


def stitch_boundaries(edges):
    edge_set = edges.copy()
    boundary_lst = []
    while len(edge_set) > 0:
        boundary = []
        edge0 = edge_set.pop()
        boundary.append(edge0)
        last_edge = edge0
        while len(edge_set) > 0:
            i, j = last_edge
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
    print("Estimated delta_{}".format(prefix), 'delta_min=', delta_min, 'delta_max=', delta_max, 'delta_avg=',
          delta_avg, 'delta_avg_std=', delta_avg_std, 'NOTE: here avg is calcculated differently')
    return delta_min, delta_max, delta_avg, delta_avg_std


def find_first_peak(x_fil, y_fil, x0, xmin, xmax, x_mean, time):
    if x_fil[0] > x_fil[-1]:
        x_fil = x_fil[::-1]
        y_fil = y_fil[::-1]
    ind = np.argmax(y_fil)
    print("max ={} {}".format(x_fil[ind], y_fil[ind]))
    # choose some points if they are:
    args = (y_fil >= 0.3) & (x_fil <= x_mean)
    x_ripple = x_fil[args]
    y_ripple = y_fil[args]

    print("sizes of ripple:{} {}".format(x_ripple.shape, y_ripple.shape))
    # return 0, 0, 0, 0, 0, 0, 0
    # find the first peak in a smoothed curve
    peaks, props = find_peaks(y_ripple, prominence=0.0001)
    print("peaks, props:", peaks, props)
    try:
        x_peak, y_peak = x_ripple[peaks[0]], y_ripple[peaks[0]]
    except:
        x_peak, y_peak = np.inf, np.inf
    length_x_peak_mean = x_mean - x_peak
    print('x_peak candidates=', x_ripple[peaks], 'y_peak candidates=', y_ripple[peaks])
    print('x_peak', x_peak, 'y_peak=', y_peak, "length_x_peak_mean=", length_x_peak_mean)
    if not x_peak:
        x_peak = x_mean
        y_peak = 0.5
        length_x_peak_mean = 0
    return x_peak, y_peak, length_x_peak_mean


def find_smooth_curve_and_bounds(x, y, x_mean, alpha=0.01):
    # calculate xmin, xmax, ymin, ymax
    xmin, xmax = min(x), max(x)
    ymin, ymax = min(y), max(y)
    # calculate xy0, xyN
    inds = x < x_mean
    print("xmin, x_mean, xmax", xmin, x_mean, xmax)
    print("ymin, ymax", ymin, ymax)
    ind_xy0 = y[inds].argmin()
    xy0 = x[inds][ind_xy0], y[inds][ind_xy0]
    inds = x > x_mean
    ind_xyN = y[inds].argmin()
    xyN = x[inds][ind_xyN], y[inds][ind_xyN]
    length_x_clip_ends = xyN[0] - xy0[0]
    # preprocess x, y arrays using histogram:
    N = 10000
    xedges = xmin + (xmax - xmin) * np.linspace(0, 1, N)
    yedges = ymin + (ymax - ymin) * np.linspace(0, 1, N)
    H, xedges, yedges = np.histogram2d(x, y, bins=(xedges, yedges))
    X, Y = np.meshgrid(xedges, yedges)
    # Histogram does not follow Cartesian convention (see Notes),
    # therefore transpose H for visualization purposes.
    H = H.T * (255.0 / H.max())
    N_non_zero = np.count_nonzero(H)
    coords = np.zeros((N_non_zero, 2))
    k = 0
    for i in range(len(xedges) - 1):
        for j in range(len(yedges) - 1):
            if H[i, j] > 0:
                coords[k, :] = X[i, j], Y[i, j]
                k += 1
    del X, Y, H
    # upper and lower lines
    lower_hull, upper_hull = find_min_max_curve(np.asarray(coords), alpha=alpha, p0=xy0, pN=xyN)

    print("lower_hull and upper_hull minmax<><><><>", min(lower_hull[:, 0]), max(lower_hull[:, 0]),
          min(upper_hull[:, 0]), max(upper_hull[:, 0]))

    x_peak, y_peak, length_x_peak_mean = find_first_peak(upper_hull[:, 0], upper_hull[:, 1], xy0[0], xmin, xmax, x_mean,
                                                         0)

    delta_min_lw, delta_max_lw, delta_avg_lw, delta_avg_std_lw = calc_thickness(lower_hull[:, 0], lower_hull[:, 1],
                                                                                x_peak, x_mean, 'lower_hull')
    delta_min_up, delta_max_up, delta_avg_up, delta_avg_std_up = calc_thickness(upper_hull[:, 0], upper_hull[:, 1],
                                                                                x_peak, x_mean, 'upper_hull')
    delta_min, delta_max, delta_avg, delta_avg_std = calc_thickness(x, y, x_peak, x_mean, 'sliced_x_y')

    return lower_hull, upper_hull, x_peak, y_peak, length_x_peak_mean, delta_min, delta_max, xy0, xyN, xmin, xmax


def Save1DArraysToFile(numpy_arrays, fn):
    lists = []
    with open(fn, 'w') as f:
        for n in range(len(numpy_arrays)):
            lists.append(numpy_arrays[n].tolist())
        f.write(json.dumps(lists))
        print('Successfully save file:', fn)


def get_x_over_R_array(input_data, fn, PointDataArrays, CellDataArrays):
    # create a new 'Pass Arrays'
    passArrays1 = PassArrays(Input=input_data)
    passArrays1.PointDataArrays = PointDataArrays
    passArrays1.CellDataArrays = CellDataArrays

    # update the view to ensure updated data information
    spreadSheetView1.Update()

    ss_data = Fetch(passArrays1)
    Np = ss_data.GetNumberOfPoints()
    print('Np=', Np)
    xr = []
    for ip in range(Np):
        regionId = ss_data.GetPointData().GetArray('RegionId').GetValue(ip)
        if regionId == 0:
            xp = ss_data.GetPoint(ip)[2]  # ONLY FOR HTG format, channel along Z axis
            rp = ss_data.GetPointData().GetArray('Result').GetValue(ip)
            xr.append((xp, rp))
    xr = np.array(xr)
    print('processing data size of x and y:', len(xr))

    Save1DArraysToFile([xr[:, 0], xr[:, 1]], fn)

    return xr[:, 0], xr[:, 1]


def gcd(a, b):
    while b: a, b = b, a % b
    return a


def simplify_fraction(x, y):
    return x // gcd(x, y), y // gcd(x, y)


def SavePvdFile(fn, source, print_text, times, three_dimension=0):
    # Split the extension from the path and normalise it to lowercase.
    file_extension = os.path.splitext(fn)[1].lower()
    # Regular expression to match the pattern
    # The structure is <directory>/<name>_<number>_<number>.<extension>
    # Where <directory> and <extension> can be of any length
    match = re.search(r"/([^/]+)_\d+_\d+\.[^.]+$", fn)
    if match:
        subn = match.group(1)
    else:
        subn = ""
    print('e:', subn)
    text1 = '''<VTKFile type="Collection" version="1.0" byte_order="LittleEndian" header_type="UInt64">\n\t<Collection>\n'''
    text3 = '''\t</Collection>
</VTKFile>'''
    text2 = ''
    for i, t in enumerate(times):
        text2 += "\t\t<DataSet timestep=\"{}\" part=\"0\" file=\"res/{}_0_{:04d}{}\"/>\n".format(t, subn, i, file_extension)
    with open(subn + '.pvd', 'w') as f:
        f.write(text1)
        f.write(text2)
        f.write(text3)
    if three_dimension:
        SaveData(fn, proxy=source)
    else:
        SaveData(fn, proxy=source)
    print("Saved", print_text, ':', fn)


def get_bounds(input) -> dict:
    info = input.GetDataInformation()
    Npoints = info.GetNumberOfPoints()
    Ncells = info.GetNumberOfCells()
    bounds = info.GetBounds()
    center = [(bounds[0] + bounds[1]) / 2, (bounds[2] + bounds[3]) / 2, (bounds[4] + bounds[5]) / 2]
    len_x = bounds[1] - bounds[0]
    len_y = bounds[3] - bounds[2]
    len_z = bounds[5] - bounds[4]

    return {
        "bounds": bounds,
        "center": center,
        "len_x": len_x,
        "len_y": len_y,
        "len_z": len_z,
        "n_points": Npoints,
        "n_cells": Ncells,
    }


def compute_volume_averaged_vars(integrateVariables) -> dict:
    passArrays1 = PassArrays(Input=integrateVariables)
    passArrays1.PointDataArrays = ['Points', 'u', 'f']
    passArrays1.CellDataArrays = ['Volume']

    ss_data = Fetch(passArrays1)
    volume = ss_data.GetCellData().GetArray('Volume').GetValue(0)
    x_mean = ss_data.GetPoint(0)
    u_mean = [ss_data.GetPointData().GetArray('u').GetValue(i) / volume for i in range(3)]
    f_mean = ss_data.GetPointData().GetArray('f').GetValue(0) / volume
    return {
        "Volume": volume,
        "x_mean": x_mean,
        "u_mean": u_mean,
        "f_mean": f_mean
    }

def single_compute_area_volume(connectivity, threshold_value: float, time: float):
    # create a new 'Threshold'
    threshold = Threshold(Input=connectivity, registrationName=f"IsoVolume_{threshold_value}")
    threshold.Scalars = ['CELLS', 'RegionId']
    threshold.ComponentMode = 'Selected'
    threshold.SelectedComponent = ''
    threshold.LowerThreshold = threshold_value
    threshold.UpperThreshold = threshold_value
    threshold.ThresholdMethod = 'Between'
    threshold.AllScalars = 1
    threshold.UseContinuousCellRange = 0
    threshold.Invert = 0
    threshold.MemoryStrategy = 'Mask Input'

    # UpdatePipeline(time=time, proxy=threshold)
    threshold.UpdatePipeline()

    # create a new 'Integrate Variables'
    integrateVolumetricVariables = IntegrateVariables(Input=threshold)
    integrateVolumetricVariables.DivideCellDataByVolume = 0

    # UpdatePipeline(time=time, proxy=integrateVolumetricVariables)
    integrateVolumetricVariables.UpdatePipeline()

    mean_vars = compute_volume_averaged_vars(integrateVolumetricVariables)

    # create a new 'Extract Surface'
    extractSurface = ExtractSurface(Input=threshold)
    extractSurface.PieceInvariant = 1
    extractSurface.NonlinearSubdivisionLevel = 1
    extractSurface.FastMode = 0
    extractSurface.RemoveGhostInterfaces = 1

    # UpdatePipeline(time=time, proxy=extractSurface)
    extractSurface.UpdatePipeline()

    # create a new 'Integrate Variables'
    integrateSurfaceVariables = IntegrateVariables(Input=extractSurface)
    integrateSurfaceVariables.DivideCellDataByVolume = 0

    # UpdatePipeline(time=time, proxy=integrateSurfaceVariables)
    integrateSurfaceVariables.UpdatePipeline()

    passArrays1 = PassArrays(Input=integrateSurfaceVariables)
    passArrays1.PointDataArrays = []
    passArrays1.CellDataArrays = ['Area']

    ss_data = Fetch(passArrays1)
    mean_vars["Area"] = ss_data.GetCellData().GetArray('Area').GetValue(0)

    # Extract bounds of a bubble
    mean_vars.update(get_bounds(threshold))

    # Do not delete volumetric representation for the biggest volume
    if threshold_value != 0:
        Delete(threshold)
        del threshold
    Delete(integrateVolumetricVariables)
    del integrateVolumetricVariables
    Delete(extractSurface)
    del extractSurface
    Delete(integrateSurfaceVariables)
    del integrateSurfaceVariables
    return mean_vars


def compute_area_volume(input, time) -> dict[int, dict]:
    hyperTreeGridToDualGrid = HyperTreeGridToDualGrid(Input=input)

    # create a new 'Iso Volume'
    isoVolume = IsoVolume(Input=hyperTreeGridToDualGrid)
    isoVolume.InputScalars = ['POINTS', 'f']
    isoVolume.ThresholdRange = [0.0, 0.5]

    # UpdatePipeline(time=time, proxy=isoVolume)
    isoVolume.UpdatePipeline()

    # create a new 'Connectivity'
    connectivity = Connectivity(Input=isoVolume)
    connectivity.ExtractionMode = 'Extract All Regions'
    connectivity.ColorRegions = 1
    connectivity.RegionIdAssignmentMode = 'Cell Count Descending'

    # UpdatePipeline(time=time, proxy=connectivity)
    connectivity.UpdatePipeline()

    info = connectivity.GetDataInformation().GetPointDataInformation()
    arrayInfo = info.GetArrayInformation("RegionId")
    print("arrayInfo of connectivity in compute_area_volume:", arrayInfo)
    region_id_range = arrayInfo.GetComponentRange(0)
    region_id_range = int(region_id_range[0]), int(region_id_range[1]),
    threshold_result = dict()
    for threshold_value in range(*region_id_range):
        print("threshold by RegionId", threshold_value)
        threshold_result[threshold_value] = single_compute_area_volume(
            connectivity=connectivity, threshold_value=threshold_value, time=time
        )

    # Sort result by volume of regions
    # Sort the dictionary by volume in descending order
    sorted_data = sorted(threshold_result.items(), key=lambda x: x[1]["Volume"], reverse=True)

    # Convert sorted list back to dictionary if needed
    threshold_result = {k: v for k, v in sorted_data}

    Delete(hyperTreeGridToDualGrid)
    del hyperTreeGridToDualGrid
    Delete(isoVolume)
    del isoVolume
    Delete(connectivity)
    del connectivity

    return threshold_result

# Read from arguments
# 2 pvd file name (by defaults in the first file in a current directory)
# 3 pvtu is swithed off by default
# ---------------------------------------------------------------------------------------------------------
start = timeit.default_timer()

parser = argparse.ArgumentParser()
parser._action_groups.pop()
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')

optional.add_argument("--infn", type=str, help="Provide the name of the input paraview files, please",
                      nargs='?', default='*.pvd')
optional.add_argument("--outPrefix", type=str, help="Provide the name of the output paraview files, please",
                      nargs='?', default='')
required.add_argument("--maxlevel", type=int, help="Provide the maximum level of refinement",
                      default=10, required=True)
required.add_argument("--iter", type=int, help="Provide the iter argument level of refinement",
                      default=0, required=True)
required.add_argument("--rangeColorbar", type=float, help="Provide max |u| for colorbar, please. 0 means automatic range",
                      default=0, required=True)
parser.add_argument("--picName", type=str, help="Provide the name for the outputed pictures, please",
                    nargs='?', default='pic')

args = parser.parse_args()
print(f"args: {args}")
infn = args.infn
out_prefix = args.outPrefix
maxlevel = args.maxlevel
iter = args.iter
range_colorbar = args.rangeColorbar
picName = args.picName

# Current PATH reading
if vtk_from_pvpython:
    path = os.path.abspath(os.getcwd())
else:
    path = os.path.join(os.getenv("HOME"), "basilisk/work/tube")

eprint("Current PATH=" + path)
if vtk_from_pvpython:
    if infn[-5::] == ".pvtu":
        my_source = XMLPartitionedUnstructuredGridReader(FileName=os.path.join(path, infn))
    elif infn[-4::] == '.pvd':
        my_source = PVDReader(FileName=os.path.join(path, infn))
    else:
        eprint('Get Active Source: No pvd or pvtu files are provided')
        my_source = GetActiveSource()
else:
    eprint('Get Active Source')
    my_source = GetActiveSource()

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# get animation scene
animationScene1 = GetAnimationScene()
print("AnimationTime=", animationScene1.AnimationTime)

# get the time-keeper
timeKeeper1 = GetTimeKeeper()
timesteps = timeKeeper1.TimestepValues  # 0, 0.1, 0.2 ...
print("timesteps=", timesteps)
try:
    NT = len(timesteps)
except:
    NT = 1
    timesteps = [timesteps]

print("renderViews, axesGrid, SpreadSheetViews, layouts.. ")
# Create a new 'Render View'
# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
# renderView1.ViewSize = [1840, 1156]
renderView1.InteractionMode = '2D'
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.OrientationAxesLabelColor = [1, 1, 1]
renderView1.OrientationAxesOutlineColor = [1, 1, 1]
renderView1.CenterOfRotation = [0.00012353062629699707, 0.0003523975610733032, 3.9594372510910034]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [1.290158870335943, 5.248223688485809, -0.040512771139578935]
renderView1.CameraFocalPoint = [1.284848907700422, 5.227816101479906, -0.023208475832503923]
renderView1.CameraViewUp = [0.9788702853248832, -0.10688095869808088, 0.17432562971565893]
renderView1.CameraViewAngle = 15.42391304347826
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 1
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.OSPRayMaterialLibrary = materialLibrary1

# AxesGrid property provides access to the AxesGrid object.
axesGrid = renderView1.AxesGrid
axesGrid.Visibility = 1
axesGrid.XTitle = 'Z'
axesGrid.YTitle = 'Y'
axesGrid.ZTitle = 'X'
axesGrid.XTitleFontSize = 20
axesGrid.YTitleFontSize = 20
axesGrid.ZTitleFontSize = 20
axesGrid.XLabelFontSize = 20
axesGrid.YLabelFontSize = 20
axesGrid.ZLabelFontSize = 20

# Edit the Properties of the AxesGrid
axesGrid.XAxisUseCustomLabels = 1  # 1 means true
axesGrid.YAxisUseCustomLabels = 1  # 1 means true
axesGrid.ZAxisUseCustomLabels = 1  # 1 means true

axesGrid.XTitleColor = [0.9, 0.9, 0.9]
axesGrid.YTitleColor = [0.9, 0.9, 0.9]
axesGrid.ZTitleColor = [0.9, 0.9, 0.9]

axesGrid.XLabelColor = [0.9, 0.9, 0.9]
axesGrid.YLabelColor = [0.9, 0.9, 0.9]
axesGrid.ZLabelColor = [0.9, 0.9, 0.9]

axesGrid.XAxisLabels = [-0.5, -0.25, 0.25, 0.5]  # np.around(np.linspace(-0.5,0.5,5),2)
axesGrid.YAxisLabels = [-0.5, -0.25, 0.25, 0.5]  # np.around(np.linspace(-0.5,0.5,5),2)
axesGrid.ZAxisLabels = np.around(np.arange(0, 30.2, 0.25), 2)

# Create a new 'SpreadSheet View'
spreadSheetView1 = CreateView('SpreadSheetView')

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------
# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
# layout1.SplitHorizontal(0, 0.5)
layout1.AssignView(0, renderView1)
layout1.SetSize(1840, 1156)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)

# defining of computational domains
my_source.UpdatePipeline()
Show(my_source, renderView1)
Hide(my_source, renderView1)

bounds_domain = get_bounds(my_source)
print({"info": "boundsDomain of cube", "bounds_domain": bounds_domain})

lDomain = bounds_domain["len_x"]

if lDomain < 1e-10:
    print("Error: lDomain is too small:", lDomain)
    sys.exit()
# Hide(my_source, renderView1)

###-----------------GENERATION of Cylinder-----------------------------------
### it is timeless therefore it is outside of the loop
if "create a new 'Cylinder'":
    print("Creating a cylinder ")
    cylinder1 = Cylinder()
    cylinder1.Resolution = 100
    cylinder1.Height = 30
    cylinder1.Capping = 0

    # create a new 'Clip'
    print("Clipping a created cylinder ")
    clip3 = Clip(Input=cylinder1)
    clip3.ClipType = 'Plane'
    clip3.HyperTreeGridClipper = 'Plane'
    clip3.Scalars = ['POINTS', 'Normals_Magnitude']
    clip3.Value = 1

    # init the 'Plane' selected for 'ClipType'
    clip3.ClipType.Normal = [0.0, 0.0, -1.0]
    clip3.ClipType.Origin = [0, 0, 0]


    # create a new 'Transform'
    print("Rotation of a created cylinder ")
    transform1 = Transform(Input=clip3)
    transform1.Transform = 'Transform'

    # init the 'Transform' selected for 'Transform'
    transform1.Transform.Translate = [0.0, 0.0, 0.5 * lDomain]  # ONLY FOR HTG format, channel along Z axis
    transform1.Transform.Rotate = [90.0, 0.0, 0.0]  # ONLY FOR HTG format, channel along Z axis

    # trace defaults for the display properties.
    transform1Display = Show(transform1, renderView1, 'UnstructuredGridRepresentation')
    transform1Display.Representation = 'Surface'
    transform1Display.AmbientColor = [1.0, 0.7843137254901961, 0.7529411764705882]
    transform1Display.ColorArrayName = [None, '']
    transform1Display.DiffuseColor = [1.0, 0.7843137254901961, 0.7529411764705882]
    transform1Display.Opacity = 0.85
    transform1Display.Specular = 1.0
    transform1Display.Luminosity = 35.0
    transform1Display.OSPRayUseScaleArray = 1
    transform1Display.OSPRayScaleArray = 'Normals'
    transform1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    transform1Display.OSPRayMaterial = 'copper'
    transform1Display.SelectOrientationVectors = 'None'
    transform1Display.ScaleFactor = 3.0
    transform1Display.SelectScaleArray = 'None'
    transform1Display.GlyphType = 'Arrow'
    transform1Display.GlyphTableIndexArray = 'None'
    transform1Display.GaussianRadius = 0.15
    transform1Display.SetScaleArray = ['POINTS', 'Normals']
    transform1Display.ScaleTransferFunction = 'PiecewiseFunction'
    transform1Display.OpacityArray = ['POINTS', 'Normals']
    transform1Display.OpacityTransferFunction = 'PiecewiseFunction'
    transform1Display.DataAxesGrid = 'GridAxesRepresentation'
    transform1Display.DataAxesGrid.GridColor = [0, 0, 0]

    transform1Display.PolarAxes = 'PolarAxesRepresentation'
    transform1Display.ScalarOpacityUnitDistance = 6.650076732513133

    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    transform1Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    transform1Display.ScaleTransferFunction.Points = [-2.220446049250313e-16, 0.0, 0.5, 0.0, 2.220446049250313e-16, 1.0,
                                                      0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    transform1Display.OpacityTransferFunction.Points = [-2.220446049250313e-16, 0.0, 0.5, 0.0, 2.220446049250313e-16, 1.0,
                                                        0.5, 0.0]

for timestep in timesteps:
    print("timestep:", timestep)
    # Properties modified on animationScene1
    print("animationScene1.AnimationTime", animationScene1.AnimationTime)
    animationScene1.AnimationTime = timestep
    # Properties modified on timeKeeper1
    timeKeeper1.Time = timestep

    ############################################################################################
    ### Compute averaged velocity, coordinate, area and volumes for each connectivity region ###
    ############################################################################################

    threshold_result = compute_area_volume(my_source, timestep)
    print("threshold_result", threshold_result)
    # get the largest bubble
    first_bubble = threshold_result[0]
    bounds = first_bubble["bounds"]

    # compare with the second bubble
    if threshold_result[1]["Volume"] / first_bubble["Volume"] > 0.1:
        bounds2 = threshold_result[1]["bounds"]
        bounds = min(bounds[0], bounds2[0]), max(bounds[1], bounds2[1]), \
            min(bounds[2], bounds2[2]), max(bounds[3], bounds2[3]), \
            min(bounds[4], bounds2[4]), max(bounds[5], bounds2[5])

    print("bounds", bounds)
    center = [(bounds[0] + bounds[1]) / 2, (bounds[2] + bounds[3]) / 2, (bounds[4] + bounds[5]) / 2]
    print("center", center)

    len_bub = threshold_result[0]["len_z"]
    len_min = max([bounds[0] - 3, 0])
    len_max = min([bounds[1] + 1, lDomain])
    length = len_max - len_min
    print(f"len_min: {len_min} len_max: {len_max} len: {length} len_bub: {len_bub} ")

    area = first_bubble["Area"]
    volume = first_bubble["Volume"]
    x_mean = first_bubble["x_mean"][2]  # ONLY FOR HTG format, channel along Z axis
    u_mean = first_bubble["u_mean"]
    print(f"First bubble: x_mean: {x_mean} u_mean: {u_mean} area: {area}")

    # ***************** SAVE ISOVOLUME ****************************
    hyperTreeGridToDualGrid1 = HyperTreeGridToDualGrid(Input=my_source)
    hyperTreeGridToDualGrid1.UpdatePipeline()

    # create a new 'Iso Volume'
    isoVolume1 = IsoVolume(Input=hyperTreeGridToDualGrid1)
    isoVolume1.InputScalars = ['POINTS', 'f']
    isoVolume1.ThresholdRange = [0.0, 0.5]
    isoVolume1.UpdatePipeline()

    fn = f"{path}/res/{out_prefix}iso_volume_0_{iter:04d}.vtu"
    SavePvdFile(fn, isoVolume1, 'volumetric data of bubble', timesteps, three_dimension=1)

    # ***************** SAVE ISOSURFACE ****************************
    # create a new 'Iso Volume'
    contour1 = Contour(Input=my_source)
    contour1.ContourBy = ['CELLS', 'f']
    contour1.ComputeNormals = 1
    contour1.ComputeGradients = 0
    contour1.ComputeScalars = 1
    contour1.OutputPointsPrecision = 'Same as input'
    contour1.GenerateTriangles = 1
    contour1.FastMode = 0
    contour1.Contourstrategy3D = 'Use Voxels'
    contour1.Isosurfaces = [0.5]
    contour1.PointMergeMethod = 'Uniform Binning'
    contour1.PointMergeMethod.Divisions = [50, 50, 50]
    contour1.PointMergeMethod.Numberofpointsperbucket = 8

    # create a new 'Connectivity'
    connectivity1 = Connectivity(Input=contour1)
    connectivity1.ColorRegions = 1
    connectivity1.RegionIdAssignmentMode = 'Cell Count Descending'
    connectivity1.ExtractionMode = 'Extract All Regions'

    # UpdatePipeline(time=timestep, proxy=connectivity1)
    connectivity1.UpdatePipeline()

    # create a new 'Calculator'
    calculator1 = Calculator(Input=connectivity1)
    calculator1.Function = 'sqrt(coordsX^2+coordsY^2)'  # ONLY FOR HTG format, channel along Z axis
    calculator1.ResultArrayName = 'Result'

    # UpdatePipeline(time=timestep, proxy=calculator1)
    calculator1.UpdatePipeline()

    fn = f"{path}/res/{out_prefix}iso_surface_0_{iter:04d}.vtp"
    SavePvdFile(fn, calculator1, 'surface data of bubble', timesteps)

    # color 'calculator1'
    calculator1Display = GetDisplayProperties(calculator1, renderView1)
    ColorBy(calculator1Display, ('POINTS', 'u', 'Magnitude'))

    # get color transfer function/color map for 'u'
    uLUT = GetColorTransferFunction('u')
    uLUT.RGBPoints = [0.003649621564209849, 0.0, 0.0, 0.5625, 0.24729264204713047, 0.0, 0.0, 1.0, 0.8041920709742089, 0.0,
                      1.0, 1.0, 1.082641237240404, 0.5, 1.0, 0.5, 1.3610904035065987, 1.0, 1.0, 0.0, 1.9179898324336777,
                      1.0, 0.0, 0.0, 2.1964389986998722, 0.5, 0.0, 0.0]
    uLUT.ColorSpace = 'RGB'
    uLUT.ScalarRangeInitialized = 1.0
    #         uxLUT.ApplyPreset('jet', True)
    len_bar = 0.5
    # get color legend/bar for uxLUT in view renderView1
    uLUTColorBar = GetScalarBar(uLUT, renderView1)
    uLUTColorBar.Orientation = 'Horizontal'
    uLUTColorBar.WindowLocation = 'Any Location'
    uLUTColorBar.ScalarBarLength = len_bar
    uLUTColorBar.Position = [0.7 - 0.5 * len_bar, 0.01]
    uLUTColorBar.Title = "mag(u)"
    uLUTColorBar.ComponentTitle = ""
    uLUTColorBar.TitleColor = [1, 1, 1]
    uLUTColorBar.LabelColor = [1, 1, 1]
    uLUTColorBar.LabelFormat = '%-#6.2g'
    uLUTColorBar.RangeLabelFormat = '%6.2g'
    # uxLUTColorBar.ScalarBarThickness = 16*2
    # uxLUTColorBar.TitleFontSize = 16*2
    # uxLUTColorBar.LabelFontSize = 16*2
    uLUTColorBar.ScalarBarThickness = 16
    uLUTColorBar.TitleFontSize = 16
    uLUTColorBar.LabelFontSize = 16

    # set color bar visibility
    uLUTColorBar.Visibility = 1

    # ----------------------------------------------------------------
    # ----------------------------------------------------------------
    info = connectivity1.GetDataInformation().GetPointDataInformation()
    arrayInfo = info.GetArrayInformation("u")
    range0 = arrayInfo.GetComponentRange(0)
    range1 = arrayInfo.GetComponentRange(1)
    range2 = arrayInfo.GetComponentRange(2)
    rgX, rgY, rgZ = np.max([abs(range0[0]), abs(range0[1])]), np.max([abs(range1[0]), abs(range1[1])]), np.max(
        [abs(range2[0]), abs(range2[1])])
    range_max = np.sqrt(rgX ** 2 + rgY ** 2 + rgZ ** 2)
    print(f"Colorbar: mag(u): [0, {range_max}]")
    print(f"Colorbar: Ux: {range0} Uy: {range1} Uz: {range2}")
    if range_colorbar:
        range_max = range_colorbar
        print(f"range_max is overrided by {range_max}")

    uLUT.RescaleTransferFunction(0, range_max)
    op = GetOpacityTransferFunction("u")
    op.RescaleTransferFunction(0, range_max)
    uLUTColorBar.UseCustomLabels = 1
    # labels from 100 to 200 in increments of 10
    uLUTColorBar.CustomLabels = np.around(np.linspace(0, range_max, 4), 1)
    print("CustomLabels=", uLUTColorBar.CustomLabels)

    # trace defaults for the display properties.
    calculator1Display.LookupTable = uLUT

    # *********************** Find the x coordinate of the first peak x_peak only for RegionId=0 ****************************
    xy0 = [0, 0]
    xyN = [0, 0]
    x, y = get_x_over_R_array(calculator1, fn="r_over_x_total_t={}.csv".format(timestep),
                              PointDataArrays=['Points', 'Result', "RegionId"], CellDataArrays=["Volume"])

    lower_hull, upper_hull, x_peak, y_peak, length_x_peak_mean, delta_min, delta_max, xy0, xyN, xmin, xmax = find_smooth_curve_and_bounds(
        x, y, x_mean, alpha=0.05)

    fn = "r_over_x_t={}.csv".format(timestep)
    Save1DArraysToFile([x, y, lower_hull, upper_hull], fn)
    # read files as below:
    # with open(fn, 'r') as f:
    #	lists = json.load(f)
    #    x, y, lower_hull, upper_hull = np.array(lists[0]), np.array(lists[1]), np.array(lists[2]), np.array(lists[3])

    # plot_graph([[x_peak, x_peak], [x_mean, x_mean]], [[0, 0.5], [0, 0.5]], \
    #            ['first peak', 'center of mass', "min edge", "max edge"], \
    #            list_x_fill=[lower_hull[:, 0], upper_hull[:, 0]], list_y_fill=[lower_hull[:, 1], upper_hull[:, 1]], \
    #            dash=['solid', 'dot', 'dot'], \
    #            xtitle="x", ytitle="r", image_name=fn[:-3] + 'pdf', mode=['lines', 'lines', 'lines'], \
    #            colors=['red', 'black', 'black'], yrange=[0, 0.5], xrange=[xmin - 0.1, xmax + 0.1], \
    #            marker_size=1, width=1000, height=500, path='./', yanchor='bottom', y0_anchor=0.01, xanchor='left',
    #            x0_anchor=0.3)

    # ***************** CLIP A BOX to calculate volume from x_peak to x_mean ****************************
    # create a new 'Clip'
    first_bubble_threshold = FindSource('IsoVolume_0')

    clip2 = Clip(Input=first_bubble_threshold)
    clip2.ClipType = 'Box'
    clip2.Scalars = ['POINTS', 'f']
    clip2.Value = 0.5
    clip2.ClipType.Position = [-0.6, -0.6, x_peak]  # ONLY FOR HTG format, channel along Z axis
    clip2.ClipType.Length = [1.2, 1.2, length_x_peak_mean]  # ONLY FOR HTG format, channel along Z axis

    # UpdatePipeline(time=timestep, proxy=clip2)
    clip2.UpdatePipeline()

    print("clip2:", get_bounds(clip2))

    # create a new 'Integrate Variables'
    integrateVariables1 = IntegrateVariables(Input=clip2)

    # UpdatePipeline(time=timestep, proxy=integrateVariables1)
    integrateVariables1.UpdatePipeline()

    # create a new 'Pass Arrays'
    passArrays1 = PassArrays(Input=integrateVariables1)
    passArrays1.PointDataArrays = ['Points']
    passArrays1.CellDataArrays = ['Volume']

    # update the view to ensure updated data information
    spreadSheetView1.Update()

    ss_data = Fetch(passArrays1)
    print('clip2 N_points=', ss_data.GetNumberOfPoints())
    print('clip2  ss_data.GetPointData()=', ss_data.GetPointData())
    print('clip2  ss_data.GetCellData()=', ss_data.GetCellData())
    volumeB = ss_data.GetCellData().GetArray('Volume').GetValue(0)
    rB = np.sqrt(volumeB / (np.pi * length_x_peak_mean))
    delta_mean = 0.5 - rB

    print('delta_min=', delta_min, 'delta_mean=', delta_mean, 'delta_max=', delta_max)

    fn = "for_excel_table.txt"
    if not Path(fn).exists():
        with open(fn, 'w') as f:
            f.write(
                "t	x_tail	x_peak	y_peak	x_mean	x_nose	x_nose_ISC	volume	UmeanV	delta_min	delta_mean	delta_max	delta_min_smooth	delta_max_smooth\n")

    with open(fn, 'a') as f:
        f.write(
            "{} {} {} {} {} {} {} {} {} {} {} {} {} {}\n".format(timestep, xy0[0], x_peak, y_peak, x_mean, xyN[0], '?',
                                                                 volume, u_mean[0], delta_min, delta_mean, delta_max, '0',
                                                                 '0'))
    print('Successfully save file:', fn)


    # ************************* SAVE CUT TAIL ***********************************************
    # create a new 'Clip'
    clip3 = Clip(Input=calculator1)
    clip3.ClipType = 'Plane'
    clip3.HyperTreeGridClipper = 'Plane'
    clip3.Scalars = ['POINTS', 'Result']
    clip3.Value = 0.23469079123049602
    x_cut = max(x_peak, xy0[0]) + 0.1
    clip3.ClipType.Origin = [0, 0, x_cut]  # ONLY FOR HTG format, channel along Z axis
    clip3.HyperTreeGridClipper.Origin = [0, 0, x_cut]  # ONLY FOR HTG format, channel along Z axis
    clip3.ClipType.Normal = [0.0, 0.0, 1.0]  # ONLY FOR HTG format, channel along Z axis

    # create a new 'Transform'
    transform2 = Transform(Input=clip3)
    transform2.Transform = 'Transform'

    # Properties modified on transform2.Transform
    transform2.Transform.Translate = [0.0, 0.0, -x_peak]  # ONLY FOR HTG format, channel along Z axis
    # create a new 'Extract Surface' convert vtu -> vtp
    extractSurface1 = ExtractSurface(Input=transform2)

    fn = f"{path}/res/{out_prefix}iso_surface_tail_0_{iter:04d}.vtp"
    SavePvdFile(fn, extractSurface1, 'tail surface data of bubble', timesteps)

    # ***************** SAVE ISOSURFACE PNG ****************************
    Show(calculator1, renderView1)
    renderView1.CameraViewUp = [1, 0, 0]  # ONLY FOR HTG format, channel along Z axis
    renderView1.CameraParallelScale = 1.3  # 1.4 0.5
    renderView1.CenterOfRotation = center
    renderView1.CameraFocalPoint = center
    renderView1.CameraPosition = [1.56, 4.5, center[2] - 4]  # ONLY FOR HTG format, channel along Z axis
    # update the view to ensure updated data information
    renderView1.Update()
    fn = path + "/" + picName + '_t=' + str(timestep) + '_ux.png'
    SaveScreenshot(fn, renderView1,
                   ImageResolution=[1900, 1077],
                   TransparentBackground=0,
                   CompressionLevel='2')
    Hide(calculator1)

    # ***************** SAVE HALF ISOSURFACE PNG ****************************
    clip_half = Clip(Input=calculator1)
    clip_half.ClipType = 'Box'
    clip_half.HyperTreeGridClipper = 'Plane'
    # clip_half.Scalars = ['POINTS', 'u']
    # clip_half.Value = 0.5
    clip_half.ClipType.Position = [-0.5, -0.5, 0]  # ONLY FOR HTG format, channel along Z axis
    clip_half.ClipType.Length = [1, 0.5, lDomain]  # ONLY FOR HTG format, channel along Z axis
    clip_half.UpdatePipeline()

    print("clip_half:", get_bounds(clip_half))

    clip_halfDisplay = Show(clip_half, renderView1, 'GeometryRepresentation')
    clip_halfDisplay.Representation = 'Surface'
    clip_halfDisplay.ColorArrayName = ['POINTS', 'u']
    clip_halfDisplay.LookupTable = uLUT
    clip_halfDisplay.Opacity = 1
    uLUTColorBar.Visibility = 1

    renderView1.CameraViewUp = [1, 0, 0]  # ONLY FOR HTG format, channel along Z axis
    renderView1.CameraParallelScale = 1.3  # 1.4 0.5
    renderView1.CenterOfRotation = center
    renderView1.CameraFocalPoint = center
    renderView1.CameraPosition = [1.56, 4.5, center[2] - 4]  # ONLY FOR HTG format, channel along Z axis
    # update the view to ensure updated data information
    renderView1.Update()

    fn = path + "/" + picName + '_t=' + str(timestep) + '_ux_half.png'
    SaveScreenshot(fn, renderView1,
                   ImageResolution=[1900, 1077],
                   TransparentBackground=0,
                   CompressionLevel='2')
    Hide(clip_half)
    # ***************** CONTOUR2 for LAMBDA2 l2 in whole domain ****************************

    # create a new 'Contour'
    contour2 = Contour(Input=my_source)
    contour2.ContourBy = ['POINTS', 'l2']
    contour2.Isosurfaces = [-4, -2.0, -1.0, -0.5, -0.25, -0.125]
    contour2.PointMergeMethod = 'Uniform Binning'
    contour2.PointMergeMethod.Divisions = [50, 50, 50]
    contour2.PointMergeMethod.Numberofpointsperbucket = 8

    # create a new 'Clip' to remove parasite lambda2 at the inlet
    clip4 = Clip(Input=contour2)
    clip4.ClipType = 'Plane'
    clip4.HyperTreeGridClipper = 'Plane'
    # init the 'Plane' selected for 'ClipType'
    clip4.ClipType.Origin = [0, 0, 1]  # ONLY FOR HTG format, channel along Z axis
    # init the 'Plane' selected for 'HyperTreeGridClipper'
    clip4.HyperTreeGridClipper.Origin = [0, 0, 1]  # ONLY FOR HTG format, channel along Z axis
    # Properties modified on clip1.ClipType
    clip4.ClipType.Normal = [0.0, 0.0, -1.0]  # ONLY FOR HTG format, channel along Z axis

    # create a new 'Extract Surface' convert vtu -> vtp
    extractSurface2 = ExtractSurface(Input=clip4)

    fn = f"{path}/res/{out_prefix}lambda2_0_{iter:04d}.vtp"
    SavePvdFile(fn, extractSurface2, 'lambda2 surface data', timesteps)

    # ***************** LAMBDA2 l2 inside bubble ****************************

    # create a new 'Threshold'
    threshold0 = Threshold(Input=clip4)
    threshold0.Scalars = ['POINTS', 'l2']
    # -2.0, -0.5
    threshold0.LowerThreshold = -2.0
    threshold0.UpperThreshold = -0.5
    threshold0.ThresholdMethod = 'Between'
    threshold0.AllScalars = 1

    # create a new 'Threshold'
    threshold1 = Threshold(Input=threshold0)
    threshold1.Scalars = ['POINTS', 'f']
    threshold1.LowerThreshold = 0
    threshold1.UpperThreshold = 0.5
    threshold1.ThresholdMethod = 'Between'
    threshold1.AllScalars = 1

    # create a new 'Extract Surface' convert vtu -> vtp
    extractSurface3 = ExtractSurface(Input=threshold1)

    fn = f"{path}/res/{out_prefix}lambda2_in_bubble_0_{iter:04d}.vtp"
    SavePvdFile(fn, extractSurface3, 'lambda2 in bubble', timesteps)

    # ****************** CONNECTIVITY1(f) as bubble contour AND CONTOUR2 (lambda2) ********************
    # get color transfer function/color map for 'l2'
    l2LUT = GetColorTransferFunction('l2')
    l2LUT.RGBPoints = [-2.0, 0.0, 1.0, 1.0, -1.55, 0.0, 0.0, 1.0, -1.5, 0.0, 0.0, 0.501960784314, -1.4499999999999997, 1.0,
                       0.0, 0.0, -1.0, 1.0, 1.0, 0.0]
    l2LUT.ColorSpace = 'RGB'
    l2LUT.ScalarRangeInitialized = 1.0
    # show data from connectivity1
    print("Showing transparent connectivity1.. ")
    connectivity1Display = Show(connectivity1, renderView1)
    # trace defaults for the display properties.
    connectivity1Display.Representation = 'Surface'
    # connectivity1Display.AmbientColor = [0.0392156862745098, 0.00784313725490196, 1.0]
    connectivity1Display.ColorArrayName = ['POINTS', '']
    connectivity1Display.DiffuseColor = [0.0392156862745098, 0.00784313725490196, 1.0]
    connectivity1Display.Opacity = 0.2
    connectivity1Display.Specular = 1.0
    connectivity1Display.SpecularPower = 1.0
    # connectivity1Display.Ambient = 0.21
    connectivity1Display.OSPRayScaleArray = 'f'
    connectivity1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    connectivity1Display.SelectOrientationVectors = 'None'
    connectivity1Display.SelectScaleArray = 'f'
    connectivity1Display.GlyphType = 'Arrow'
    connectivity1Display.GlyphTableIndexArray = 'f'
    connectivity1Display.SetScaleArray = ['POINTS', 'f']
    connectivity1Display.ScaleTransferFunction = 'PiecewiseFunction'
    connectivity1Display.OpacityArray = ['POINTS', 'f']
    connectivity1Display.OpacityTransferFunction = 'PiecewiseFunction'
    connectivity1Display.DataAxesGrid = 'GridAxesRepresentation'
    connectivity1Display.PolarAxes = 'PolarAxesRepresentation'
    connectivity1Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]
    connectivity1Display.ScaleTransferFunction.Points = [0.5, 0.0, 0.5, 0.0, 0.5001220703125, 1.0, 0.5, 0.0]
    connectivity1Display.OpacityTransferFunction.Points = [0.5, 0.0, 0.5, 0.0, 0.5001220703125, 1.0, 0.5, 0.0]

    # show data from contour2
    print("Showing contour2.. ")
    contour2Display = Show(threshold0, renderView1, 'GeometryRepresentation')
    # trace defaults for the display properties.
    contour2Display.Representation = 'Surface'
    contour2Display.ColorArrayName = ['POINTS', 'l2']
    contour2Display.LookupTable = l2LUT
    contour2Display.Opacity = 0.5
    contour2Display.Specular = 1.0
    contour2Display.OSPRayScaleArray = 'l2'
    contour2Display.OSPRayScaleFunction = 'PiecewiseFunction'
    contour2Display.SelectOrientationVectors = 'None'
    contour2Display.ScaleFactor = 0.2388025760650635
    contour2Display.SelectScaleArray = 'l2'
    contour2Display.GlyphType = 'Arrow'
    contour2Display.GlyphTableIndexArray = 'l2'
    # contour2Display.GaussianRadius = 0.011940128803253174
    contour2Display.SetScaleArray = ['POINTS', 'l2']
    contour2Display.ScaleTransferFunction = 'PiecewiseFunction'
    contour2Display.OpacityArray = ['POINTS', 'l2']
    contour2Display.OpacityTransferFunction = 'PiecewiseFunction'
    contour2Display.DataAxesGrid = 'GridAxesRepresentation'
    contour2Display.PolarAxes = 'PolarAxesRepresentation'
    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    contour2Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]
    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    contour2Display.ScaleTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, -0.9998779296875, 1.0, 0.5, 0.0]
    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    contour2Display.OpacityTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, -0.9998779296875, 1.0, 0.5, 0.0]
    connectivity1Display.SetScalarBarVisibility(renderView1, False)


    # update the view to ensure updated data information
    renderView1.CameraViewUp = [1, 0, 0]
    renderView1.CameraParallelScale = 1.8
    renderView1.CenterOfRotation = center
    renderView1.CameraFocalPoint = [0, 0, center[2] - 0.5]
    renderView1.CameraPosition = [1.56, 4.5, center[2] - 4]
    renderView1.Update()

    fn = path + "/" + picName + '_t=' + str(timestep) + '_Lambda2.png'
    SaveScreenshot(
        fn, renderView1,
        ImageResolution=[1900, 1077],
        TransparentBackground=0,
        CompressionLevel='2'
    )
    print('File=' + fn + ' generated successfully')

    renderView1.CameraFocalPoint = [0, 0, center[2] - 1]
    renderView1.CameraPosition = [1.56, 4.5, center[2] - 5]
    renderView1.Update()
    fn = path + "/" + picName + '_t=' + str(timestep) + '_Lambda2_small_shift.png'
    SaveScreenshot(
        fn, renderView1,
        ImageResolution=[1900, 1077],
        TransparentBackground=0,
        CompressionLevel='2'
    )
    print('File=' + fn + ' generated successfully')


    renderView1.CameraFocalPoint = [0, 0, center[2] - 2]
    renderView1.CameraPosition = [1.56, 4.5, center[2] - 6]
    renderView1.Update()
    fn = path + "/" + picName + '_t=' + str(timestep) + '_Lambda2_large_shift.png'
    SaveScreenshot(
        fn, renderView1,
        ImageResolution=[1900, 1077],
        TransparentBackground=0,
        CompressionLevel='2'
    )
    print('File=' + fn + ' generated successfully')

    print("Showing connectivity1 and hiding contour2.. ")
    Hide(threshold0, renderView1)  # hide lambda2 in whole domain
    # ****************** LAMBDA2 inside a bubble  ********************
    # trace defaults for the display properties.
    connectivity1Display.Representation = 'Surface'
    connectivity1Display.Opacity = 0.2

    # show data from threshold1
    print("Showing threshold1.. ")
    threshold1Display = Show(threshold1, renderView1, 'GeometryRepresentation')
    # get color transfer function/color map for 'l2'
    # l2LUT = GetColorTransferFunction('l2')
    # l2LUT.RGBPoints = [-1.0, 0.054901960784313725, 0.9411764705882353, 0.12941176470588237, -0.75, 0.865, 0.865, 0.865,
    #                    -0.5, 1.0, 1.0, 0.0]
    # l2LUT.ScalarRangeInitialized = 1.0
    # trace defaults for the display properties.
    threshold1Display.Representation = 'Surface'
    threshold1Display.ColorArrayName = ['POINTS', 'l2']
    threshold1Display.LookupTable = l2LUT
    threshold1Display.Opacity = 0.5
    threshold1Display.Specular = 1.0
    threshold1Display.OSPRayScaleArray = 'l2'
    threshold1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    threshold1Display.SelectOrientationVectors = 'None'
    threshold1Display.ScaleFactor = 0.2388025760650635
    threshold1Display.SelectScaleArray = 'l2'
    threshold1Display.GlyphType = 'Arrow'
    threshold1Display.GlyphTableIndexArray = 'l2'
    # threshold1Display.GaussianRadius = 0.011940128803253174
    threshold1Display.SetScaleArray = ['POINTS', 'l2']
    threshold1Display.ScaleTransferFunction = 'PiecewiseFunction'
    threshold1Display.OpacityArray = ['POINTS', 'l2']
    threshold1Display.OpacityTransferFunction = 'PiecewiseFunction'
    threshold1Display.DataAxesGrid = 'GridAxesRepresentation'
    threshold1Display.PolarAxes = 'PolarAxesRepresentation'
    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    threshold1Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]
    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    threshold1Display.ScaleTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, -0.9998779296875, 1.0, 0.5, 0.0]
    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    threshold1Display.OpacityTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, -0.9998779296875, 1.0, 0.5, 0.0]
    # show color legend
    threshold1Display.SetScalarBarVisibility(renderView1, False)
    connectivity1Display.SetScalarBarVisibility(renderView1, False)

    l2PWF = GetOpacityTransferFunction('l2')
    l2PWF.Points = [-1.0, 0.0, 0.5, 0.0, -0.5, 1.0, 0.5, 0.0]
    l2PWF.ScalarRangeInitialized = 1

    renderView1.CameraViewUp = [1, 0, 0]
    renderView1.CameraParallelScale = 1.3
    renderView1.CenterOfRotation = center
    renderView1.CameraFocalPoint = center
    renderView1.CameraPosition = [1.56, 4.5, center[2] - 4]
    renderView1.Update()

    fn = path + "/" + picName + '_t=' + str(timestep) + '_Lambda2_in_bubble.png'
    SaveScreenshot(fn, renderView1,
                   ImageResolution=[1900, 1077],
                   TransparentBackground=0,
                   CompressionLevel='2')
    print('File=' + fn + ' generated successfully')
    Hide(transform1, renderView1)  # hide cylinder
    Hide(extractSurface3, renderView1)  # hide lambda2 in bubble
    Hide(threshold1, renderView1)  # hide lambda2 in bubble
    Hide(connectivity1, renderView1)  # hide bubble contour

    # ******************************SLICE bubble*******************************************
    # *************************************************************************************
    # create a new 'Slice'
    slice1 = Slice(Input=hyperTreeGridToDualGrid1)
    slice1.SliceType = 'Plane'
    slice1.HyperTreeGridSlicer = 'Plane'
    # slice1.UseDual = 0
    # slice1.Crinkleslice = 1
    slice1.Triangulatetheslice = 1
    slice1.SliceOffsetValues = [0.0]
    slice1.PointMergeMethod = 'Uniform Binning'
    # init the 'Plane' selected for 'SliceType'
    slice1.SliceType.Origin = [0.0, 0.0, 0.5 * (len_min + len_max)]
    slice1.SliceType.Normal = [0.0, 1.0, 0.0]
    slice1.SliceType.Offset = 0.0
    # init the 'Plane' selected for 'HyperTreeGridSlicer'
    slice1.HyperTreeGridSlicer.Origin = [0.0, 0.0, 0.5 * (len_min + len_max)]
    slice1.HyperTreeGridSlicer.Normal = [0.0, 1.0, 0.0]
    slice1.HyperTreeGridSlicer.Offset = 0.0
    # init the 'Uniform Binning' selected for 'PointMergeMethod'
    slice1.PointMergeMethod.Divisions = [50, 50, 50]
    slice1.PointMergeMethod.Numberofpointsperbucket = 8
    slice1.UpdatePipeline()

    print("hyperTreeGridToDualGrid1:", get_bounds(hyperTreeGridToDualGrid1))
    print("slice1 before:", get_bounds(slice1))

    fn = f"{path}/res/{out_prefix}slice_0_{iter:04d}.vtp"
    SavePvdFile(fn, slice1, 'slice1 data', timesteps)

    # print("slice1 after SavePvdFile:", get_bounds(slice1))

    # pointSource = PointSource()
    # pointSource.Center = [0, 0, xmin]  # Center of the point cloud
    # pointSource.NumberOfPoints = 500  # Number of seed points
    # pointSource.Radius = 0.5  # Radius of the point cloud
    #
    # # create a new 'Stream Tracer'
    # streamTracer3 = StreamTracer(Input=slice1, SeedType=pointSource)
    # streamTracer3.Vectors = ['POINTS', 'u']
    # streamTracer3.InterpolatorType = 'Interpolator with Point Locator'
    # streamTracer3.SurfaceStreamlines = 1
    # streamTracer3.IntegrationDirection = 'BOTH'
    # streamTracer3.IntegratorType = 'Runge-Kutta 4-5'
    # streamTracer3.IntegrationStepUnit = 'Cell Length'
    # streamTracer3.InitialStepLength = 0.2
    # streamTracer3.MinimumStepLength = 0.01
    # streamTracer3.MaximumStepLength = 0.5
    # streamTracer3.MaximumSteps = 2000
    # streamTracer3.MaximumStreamlineLength = 20.0
    # streamTracer3.TerminalSpeed = 1e-12
    # streamTracer3.MaximumError = 1e-06
    # streamTracer3.ComputeVorticity = 1

    # print("streamTracer3:", get_bounds(streamTracer3))
    # print("slice1:", get_bounds(slice1))
    # # show data in view
    # streamTracer3Display = Show(streamTracer3, renderView1, 'GeometryRepresentation')
    # streamTracer3.UpdatePipeline()

    clip5 = Clip(Input=my_source)
    clip5.Scalars = ['CELLS', 'u']
    clip5.ClipType = 'Plane'
    clip5.ClipType.Origin = [0, 0, x_mean]
    clip5.ClipType.Normal = [0, 1, 0]
    clip5.HyperTreeGridClipper = 'Plane'
    clip5.HyperTreeGridClipper.Origin = [0, 0, x_mean]
    clip5.HyperTreeGridClipper.Normal = [0, 1, 0]
    clip5.UpdatePipeline()

    print("clip5:", get_bounds(clip5))

    threshold2 = Threshold(Input=clip5)
    threshold2.Scalars = ['POINTS', 'fs']
    threshold2.LowerThreshold = 0
    threshold2.UpperThreshold = 0.5
    threshold2.ThresholdMethod = 'Between'
    threshold2.AllScalars = 1

    threshold2Display = Show(threshold2, renderView1, 'GeometryRepresentation')
    threshold2Display.Representation = 'Surface'
    threshold2Display.ColorArrayName = ['CELLS', 'u']
    threshold2Display.LookupTable = uLUT
    threshold2Display.Opacity = 0.8
    threshold2Display.Specular = 1

    threshold2Display.OSPRayScaleArray = 'u'
    threshold2Display.OSPRayScaleFunction = 'PiecewiseFunction'
    threshold2Display.SelectOrientationVectors = 'None'
    threshold2Display.ScaleFactor = 0.25
    threshold2Display.SelectScaleArray = 'u'
    threshold2Display.GlyphType = 'Arrow'
    threshold2Display.GlyphTableIndexArray = 'u'
    threshold2Display.SetScaleArray = ['CELLS', 'u']
    threshold2Display.ScaleTransferFunction = 'PiecewiseFunction'
    threshold2Display.OpacityArray = ['CELLS', 'u']
    threshold2Display.OpacityTransferFunction = 'PiecewiseFunction'
    threshold2Display.DataAxesGrid = 'GridAxesRepresentation'
    threshold2Display.PolarAxes = 'PolarAxesRepresentation'
    threshold2Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]
    threshold2Display.ScaleTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, -0.9998779296875, 1.0, 0.5, 0.0]
    threshold2Display.OpacityTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, -0.9998779296875, 1.0, 0.5, 0.0]
    threshold2Display.SetScalarBarVisibility(renderView1, True)

    # uLUTColorBar.Visibility = 1

    print("clip5 after:", get_bounds(threshold2))

    # create a new 'Slice'
    slice1 = Slice(Input=contour1)
    slice1.SliceType = 'Plane'
    slice1.HyperTreeGridSlicer = 'Plane'
    slice1.Triangulatetheslice = 1
    slice1.SliceOffsetValues = [0.0]
    slice1.PointMergeMethod = 'Uniform Binning'
    # init the 'Plane' selected for 'SliceType'
    slice1.SliceType.Origin = [0.0, 0, 0.5 * (len_min + len_max)]
    slice1.SliceType.Normal = [0.0, 1.0, 0.0]
    slice1.SliceType.Offset = 0.1
    # init the 'Plane' selected for 'HyperTreeGridSlicer'
    slice1.HyperTreeGridSlicer.Origin = [0.0, 0, 0.5 * (len_min + len_max)]
    slice1.HyperTreeGridSlicer.Normal = [0.0, 1.0, 0.0]
    slice1.HyperTreeGridSlicer.Offset = 0.1
    # init the 'Uniform Binning' selected for 'PointMergeMethod'
    slice1.PointMergeMethod.Divisions = [50, 50, 50]
    slice1.PointMergeMethod.Numberofpointsperbucket = 32
    slice1.UpdatePipeline()

    print("Showing slice1Display of bubble.. ")
    slice1Display = Show(slice1, renderView1, 'GeometryRepresentation')
    slice1Display.Representation = 'Surface'
    slice1Display.ColorArrayName = ['POINTS', '']
    slice1Display.Opacity = 1
    slice1Display.Specular = 1.0
    slice1Display.LineWidth = 3.0
    slice1Display.SpecularPower = 1.0
    slice1Display.AmbientColor = [1.0, 1.0, 1.0]  # RGB for white
    slice1Display.DiffuseColor = [1.0, 1.0, 1.0]  # RGB for white

    print("slice1 after:", get_bounds(slice1))

    renderView1.InteractionMode = '2D'
    renderView1.CameraViewUp = [1, 0, 0]
    renderView1.CameraParallelScale = 1.2
    renderView1.CenterOfRotation = center
    renderView1.CameraFocalPoint = [-0.2, 0, center[2] - 1.8]
    renderView1.CameraPosition = [-0.2, 4.5, center[2] - 1.8]
    # renderView1.OrientationAxesVisibility = 1
    # Hide(transform1, renderView1)  # turn off cylinder
    renderView1.Update()

    fn = path + "/" + picName + '_t=' + str(timestep) + '_ux_slice.png'
    SaveScreenshot(
        fn, renderView1,
        ImageResolution=[1900, 500],
        TransparentBackground=0,
        CompressionLevel='2'
    )
    print('File=' + fn + ' generated successfully')
    Show(transform1, renderView1)  # turn on cylinder

    # create a new 'Pass Arrays'
    passArrays1 = PassArrays(Input=slice1)
    passArrays1.PointDataArrays = ['Points']
    passArrays1.CellDataArrays = []

    # update the view to ensure updated data information
    spreadSheetView1.Update()

    ss_data = Fetch(passArrays1)
    Np = ss_data.GetNumberOfPoints()
    print('Np=', Np)
    xr = []
    for ip in range(Np):
        xp = ss_data.GetPoint(ip)[0]
        yp = ss_data.GetPoint(ip)[1]
        xr.append((xp, yp))
    xr = np.array(xr)
    print('processing data size of x and y:', xr.shape[0])

    if xr.size:
        fn = "slice_t={}.csv".format(timestep)
        Save1DArraysToFile([xr[:, 0], xr[:, 1]], fn)

    # Freeing Memory
    Delete(slice1)
    del slice1
    Delete(threshold2)
    del threshold2
    Delete(clip5)
    del clip5
    Delete(extractSurface3)
    del extractSurface3
    Delete(threshold1)
    del threshold1
    Delete(threshold0)
    del threshold0
    Delete(extractSurface2)
    del extractSurface2
    Delete(clip4)
    del clip4
    Delete(contour2)
    del contour2
    Delete(clip_half)
    del clip_half
    Delete(extractSurface1)
    del extractSurface1
    Delete(transform2)
    del transform2
    Delete(clip3)
    del clip3
    Delete(clip2)
    del clip2
    Delete(first_bubble_threshold)
    del first_bubble_threshold
    Delete(calculator1)
    del calculator1
    Delete(connectivity1)
    del connectivity1
    Delete(contour1)
    del contour1
    Delete(isoVolume1)
    del isoVolume1
    Delete(hyperTreeGridToDualGrid1)
    del hyperTreeGridToDualGrid1

Delete(my_source)
del my_source

stop = timeit.default_timer()
total_time = stop - start

# output running time in a nice format.
mins, secs = divmod(total_time, 60)
hours, mins = divmod(mins, 60)

sys.stdout.write("Total running time: %d:%d:%d.\n" % (hours, mins, secs))
