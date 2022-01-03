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
import vtk
import numpy as np
import glob, os, sys
import logging
import timeit
import argparse
logging.basicConfig(format='%(message)s')
log = logging.getLogger(__name__)
import matplotlib.pyplot as plt
# from scipy.signal import butter,filtfilt # for fft filtering
# from scipy.ndimage import gaussian_filter1d # gaussian filtering
# from scipy.signal import savgol_filter
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
from scipy import interpolate, ndimage

import networkx as nx
import statsmodels.api as sm
import random

import itertools
import operator
from scipy.signal import find_peaks
from tsmoothie.smoother import *
from fractions import Fraction #convert a decimal number into fraction
import plotly.graph_objects as go
from pathlib import Path

from matplotlib.pyplot import *
from scipy import interpolate

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


def plot_graph(list_x, list_y, names, xtitle, ytitle, image_name, list_x_fill=[], list_y_fill=[], mode=[], \
			   dash=['solid', 'dot', 'dash', 'longdash'], \
			   colors=['blue', 'red', 'hsv(120,100,100)', 'green', 'black' ], \
			   marker_size=15, xrange =[], yrange = [], \
			   marker_style = ['circle', 'triangle-up', 'triangle-down','square', 'diamond', 'cross',  'x-thin', 'cross-thin' ], \
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
		area = np.sqrt(s * np.abs(s - a) * np.abs(s - b) * np.abs(s - c))
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
	print("Estimated delta_{}".format(prefix), 'delta_min=', delta_min, 'delta_max=', delta_max, 'delta_avg=', delta_avg, 'delta_avg_std=', delta_avg_std, 'NOTE: here avg is calcculated differently')
	return delta_min, delta_max, delta_avg, delta_avg_std


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
	if not x_peak:
		x_peak = x_mean
		y_peak = 0.5
		length_x_peak_mean = 0
	return x_peak, y_peak, length_x_peak_mean

def find_smooth_curve_and_bounds(x, y, x_mean, alpha = 0.01):
	# calculate xmin, xmax, ymin, ymax
	xmin, xmax = min(x), max(x)
	ymin, ymax = min(y), max(y)
	# calculate xy0, xyN
	inds = x < x_mean
	ind_xy0 = y[inds].argmin()
	xy0 = x[inds][ind_xy0], y[inds][ind_xy0]
	inds = x > x_mean
	ind_xyN = y[inds].argmin()
	xyN = x[inds][ind_xyN], y[inds][ind_xyN]
	length_x_clip_ends = xyN[0] - xy0[0]
	#preprocess x, y arrays using histogram:
	N = 10000
	xedges = xmin + (xmax - xmin)*np.linspace(0,1,N)
	yedges = ymin + (ymax - ymin)*np.linspace(0,1,N)
	H, xedges, yedges = np.histogram2d(x, y, bins=(xedges, yedges))
	X, Y = np.meshgrid(xedges, yedges)
	# Histogram does not follow Cartesian convention (see Notes),
	# therefore transpose H for visualization purposes.
	H = H.T *(255.0/H.max())
	N_non_zero = np.count_nonzero(H)
	coords = np.zeros((N_non_zero, 2))
	k = 0
	for i in range(len(xedges)-1):
		for j in range(len(yedges)-1):
			if H[i,j]>0:
				coords[k,:] = X[i,j], Y[i,j]
				k += 1
	del X, Y, H
	# upper and lower lines
	lower_hull, upper_hull = find_min_max_curve(np.asarray(coords), alpha = alpha, p0=xy0, pN=xyN)

	print("lower_hull and upper_hull minmax<><><><>", min(lower_hull[:,0]), max(lower_hull[:,0]), min(upper_hull[:,0]), max(upper_hull[:,0]))

	x_peak, y_peak, length_x_peak_mean = find_first_peak(upper_hull[:,0], upper_hull[:,1], xy0[0], xmin, xmax, x_mean, 0)


	delta_min_lw, delta_max_lw, delta_avg_lw, delta_avg_std_lw = calc_thickness(lower_hull[:,0], lower_hull[:,1], x_peak, x_mean, 'lower_hull')
	delta_min_up, delta_max_up, delta_avg_up, delta_avg_std_up = calc_thickness(upper_hull[:,0], upper_hull[:,1], x_peak, x_mean, 'upper_hull')
	delta_min, delta_max, delta_avg, delta_avg_std = calc_thickness(x, y, x_peak, x_mean, 'sliced_x_y')


	return lower_hull, upper_hull, x_peak, y_peak, length_x_peak_mean, delta_min, delta_max, xy0, xyN, xmin, xmax

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

def SavePvdFile(fn, source, print_text, three_dimension=0):

	basename = os.path.basename(fn)
	subn = basename[:-11]
	print('e:', subn)
	times = []
	for file in glob.glob('dump-*'):
		times.append(float(os.path.basename(file).split('-')[-1]) )
	times = sorted(times)
	text1='''<VTKFile type="Collection" version="1.0" byte_order="LittleEndian" header_type="UInt64">
    <Collection>'''
	text3='''    </Collection>
</VTKFile>'''
	text2 =''
	for i,t in enumerate(times):
		text2 += "        <DataSet timestep=\"{}\" part=\"0\" file=\"res/{}_0_{:04d}.vtp\"/>\n".format(t, subn, i)
	with open(subn +'.pvd' , 'w') as f:
		f.write(text1)
		f.write(text2)
		f.write(text3)
	if three_dimension:
		SaveData(fn, proxy=source, UseSubdirectory=0)
	else:
		SaveData(fn, proxy=source)
	print("Saved", print_text, ':', fn)

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
parser.add_argument("-picName", type=str, help="Provide the name for the outputed pictures, please",
					nargs='?', default='pic')

# parser.add_argument('--foo', action='store_const', const=2, default=42)
args = parser.parse_args()
print(args)
infn = args.infn
outfn = args.outfn
maxlevel = args.maxlevel
iter = args.iter
volumetric_repr = args.volumetric_repr
picName = args.picName

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
# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
# renderView1.ViewSize = [1840, 1156]
renderView1.InteractionMode = '2D'
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.OrientationAxesLabelColor = [1,1,1]
renderView1.OrientationAxesOutlineColor = [1,1,1]
renderView1.CenterOfRotation = [3.9594372510910034, 0.00012353062629699707, 0.0003523975610733032]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [-0.040512771139578935, 1.290158870335943, 5.248223688485809]
renderView1.CameraFocalPoint = [-0.023208475832503923, 1.284848907700422, 5.227816101479906]
renderView1.CameraViewUp = [0.17432562971565893, 0.9788702853248832, -0.10688095869808088]
renderView1.CameraViewAngle = 15.42391304347826
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 1.1256051366537485
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.OSPRayMaterialLibrary = materialLibrary1

# AxesGrid property provides access to the AxesGrid object.
axesGrid = renderView1.AxesGrid
axesGrid.Visibility = 1
axesGrid.XTitle = 'X'
axesGrid.YTitle = 'Y'
axesGrid.ZTitle = 'Z'
axesGrid.XTitleFontSize = 20
axesGrid.YTitleFontSize = 20
axesGrid.ZTitleFontSize = 20
axesGrid.XLabelFontSize = 20
axesGrid.YLabelFontSize = 20
axesGrid.ZLabelFontSize = 20

# Edit the Properties of the AxesGrid
axesGrid.XAxisUseCustomLabels = 1    # 1 means true
axesGrid.YAxisUseCustomLabels = 1    # 1 means true
axesGrid.ZAxisUseCustomLabels = 1    # 1 means true


axesGrid.XTitleColor = [0.9, 0.9, 0.9]
axesGrid.YTitleColor = [0.9, 0.9, 0.9]
axesGrid.ZTitleColor = [0.9, 0.9, 0.9]

axesGrid.XLabelColor = [0.9, 0.9, 0.9]
axesGrid.YLabelColor = [0.9, 0.9, 0.9]
axesGrid.ZLabelColor = [0.9, 0.9, 0.9]

axesGrid.XAxisLabels = np.around(np.arange(0,30.2,0.25),2)
axesGrid.YAxisLabels = [-0.5, -0.25, 0.25, 0.5]#np.around(np.linspace(-0.5,0.5,5),2)
axesGrid.ZAxisLabels = [-0.5, -0.25, 0.25, 0.5]#np.around(np.linspace(-0.5,0.5,5),2)

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
Show(my_source, renderView1)
renderView1.Update()

info = my_source.GetDataInformation()
boundsDomain = info.GetBounds()
Npoints = info.GetNumberOfPoints()
Ncells  = info.GetNumberOfCells()
print('Npoints=', Npoints, '  Ncells=', Ncells)

lDomain = boundsDomain[1] - boundsDomain[0]
print("boundsDomain of cube=", boundsDomain, " lDomain=", lDomain)

if lDomain<1e-10:
	print("Error: lDomain is too small:", lDomain)
	sys.exit()
Hide(my_source, renderView1)
###-----------------GENERATION of Cylinder-----------------------------------
### it is timeless therefore it is outside of the loop
# create a new 'Cylinder'
print ("creating a cylinder.. ")
cylinder1 = Cylinder()
cylinder1.Resolution = 100
cylinder1.Height = 30
cylinder1.Capping = 0

# create a new 'Clip'
clip3 = Clip(Input=cylinder1)
clip3.ClipType = 'Plane'
clip3.HyperTreeGridClipper = 'Plane'
clip3.Scalars = ['POINTS', 'Normals_Magnitude']
clip3.Value = 1

# init the 'Plane' selected for 'ClipType'
clip3.ClipType.Normal = [0.0, 0.0, 1.0]

# create a new 'Transform'
transform1 = Transform(Input=clip3)
transform1.Transform = 'Transform'

# init the 'Transform' selected for 'Transform'
transform1.Transform.Translate = [0.5*lDomain, 0.0, 0.0]
transform1.Transform.Rotate = [0.0, 0.0, 90.0]

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
transform1Display.DataAxesGrid.GridColor = [0,0,0]

transform1Display.PolarAxes = 'PolarAxesRepresentation'
transform1Display.ScalarOpacityUnitDistance = 6.650076732513133

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
transform1Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
transform1Display.ScaleTransferFunction.Points = [-2.220446049250313e-16, 0.0, 0.5, 0.0, 2.220446049250313e-16, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
transform1Display.OpacityTransferFunction.Points = [-2.220446049250313e-16, 0.0, 0.5, 0.0, 2.220446049250313e-16, 1.0, 0.5, 0.0]

# Show(transform1, renderView1)




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
Hide(connectivity1, renderView1)

bounds = connectivity1.GetDataInformation().GetBounds()
print("bounds_coarse=",bounds)
center = [(bounds[0] + bounds[1])/2, (bounds[2] + bounds[3])/2,(bounds[4] + bounds[5])/2]
print('center_coarse=',center)

len_bub = bounds[1] - bounds[0]
len_min = max([bounds[0] - 3, 0])
len_max = min([bounds[1] + 1, lDomain])
length = len_max - len_min

s = "Coarse: len_min: {} len_max: {} len: {} len_bub: {} ".format(len_min, len_max, length, len_bub)
print (s)

Delete(connectivity1)
del connectivity1
Delete(isoVolume1)
del isoVolume1

# ***************** CLIP A BOX ****************************
# create a new 'Clip'
clip1 = Clip(Input=my_source)
clip1.ClipType = 'Box'
clip1.HyperTreeGridClipper = 'Plane'
clip1.Scalars = ['POINTS', 'f']
clip1.Value = 0.5
# init the 'Box' selected for 'ClipType'
clip1.ClipType.Position = [len_min - 0.1, -0.6, -0.6]
clip1.ClipType.Length = [len_max + 0.2, 1.2, 1.2]

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
dy = lDomain/2.0**maxlevel
dx = 1.5*dy
Nx = int((len_max - len_min)/dx)
Nyz = int(1.0/dy)
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

# create a new 'Connectivity' to extract the biggest isovolume
connectivity1 = Connectivity(Input=isoVolume1)
# connectivity1.ExtractionMode = 'Extract Largest Region'
connectivity1.ColorRegions = 0
# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------
Show(connectivity1, renderView1)
Hide(connectivity1, renderView1)

# get color transfer function/color map for 'ux'
uxLUT = GetColorTransferFunction('ux')
uxLUT.RGBPoints = [0.003649621564209849, 0.0, 0.0, 0.5625, 0.24729264204713047, 0.0, 0.0, 1.0, 0.8041920709742089, 0.0, 1.0, 1.0, 1.082641237240404, 0.5, 1.0, 0.5, 1.3610904035065987, 1.0, 1.0, 0.0, 1.9179898324336777, 1.0, 0.0, 0.0, 2.1964389986998722, 0.5, 0.0, 0.0]
uxLUT.ColorSpace = 'RGB'
uxLUT.ScalarRangeInitialized = 1.0
#         uxLUT.ApplyPreset('jet', True)
len_bar = 0.5
# get color legend/bar for uxLUT in view renderView1
uxLUTColorBar = GetScalarBar(uxLUT, renderView1)
uxLUTColorBar.Orientation = 'Horizontal'
uxLUTColorBar.WindowLocation = 'AnyLocation'
uxLUTColorBar.ScalarBarLength = len_bar
uxLUTColorBar.Position = [0.7 - 0.5*len_bar, 0.01]
uxLUTColorBar.Title = '|u|'
uxLUTColorBar.ComponentTitle = ''
uxLUTColorBar.TitleColor = [1, 1, 1]
uxLUTColorBar.LabelColor = [1, 1, 1]
uxLUTColorBar.LabelFormat = '%-#6.2g'
uxLUTColorBar.RangeLabelFormat = '%6.2g'
# uxLUTColorBar.ScalarBarThickness = 16*2
# uxLUTColorBar.TitleFontSize = 16*2
# uxLUTColorBar.LabelFontSize = 16*2
uxLUTColorBar.ScalarBarThickness = 16
uxLUTColorBar.TitleFontSize = 16
uxLUTColorBar.LabelFontSize = 16


# set color bar visibility
uxLUTColorBar.Visibility = 1

# ----------------------------------------------------------------
# ----------------------------------------------------------------
info = connectivity1.GetDataInformation().DataInformation
arrayInfo = info.GetArrayInformation("u.x", vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS)
range0 = arrayInfo.GetComponentRange(0)
range1 = arrayInfo.GetComponentRange(1)
range2 = arrayInfo.GetComponentRange(2)
rgX, rgY, rgZ = np.max([abs(range0[0]),abs(range0[1])]), np.max([abs(range1[0]),abs(range1[1])]), np.max([abs(range2[0]),abs(range2[1])])
range_max = np.sqrt(rgX**2 + rgY**2 + rgZ**2)
print("velocity range_max= ", range_max)
print("velocity range_xyz= ", range0, range1, range2)
uxLUT.RescaleTransferFunction(0, range_max)
op = GetOpacityTransferFunction("ux")
op.RescaleTransferFunction(0, range_max)
uxLUTColorBar.UseCustomLabels = 1
# labels from 100 to 200 in increments of 10
uxLUTColorBar.CustomLabels = np.around(np.linspace(0, range_max, 4), 1)
#         uxLUTColorBar.CustomLabels = np.arange(0, range_max, 0.5)
print("CustomLabels=",uxLUTColorBar.CustomLabels )


bounds = connectivity1.GetDataInformation().GetBounds()
print("bounds_refined=",bounds)
center = [(bounds[0] + bounds[1])/2, (bounds[2] + bounds[3])/2, (bounds[4] + bounds[5])/2]
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
Hide(integrateVariables1, renderView1)
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
extractSurface1Display = Show(extractSurface1, renderView1)
# Hide(extractSurface1, renderView1)
# trace defaults for the display properties.
extractSurface1Display.Representation = 'Surface'
extractSurface1Display.ColorArrayName = ['POINTS', 'u.x']
extractSurface1Display.LookupTable = uxLUT
extractSurface1Display.Specular = 1.0
extractSurface1Display.SelectTCoordArray = 'None'
extractSurface1Display.SelectNormalArray = 'None'
extractSurface1Display.SelectTangentArray = 'None'
extractSurface1Display.OSPRayScaleArray = 'Result'
extractSurface1Display.OSPRayScaleFunction = 'PiecewiseFunction'
extractSurface1Display.SelectOrientationVectors = 'None'
extractSurface1Display.ScaleFactor = 0.31534409523010254
extractSurface1Display.SelectScaleArray = 'ux'
extractSurface1Display.GlyphType = 'Arrow'
extractSurface1Display.GlyphTableIndexArray = 'Result'
extractSurface1Display.GaussianRadius = 0.015767204761505126
extractSurface1Display.SetScaleArray = ['POINTS', 'ux']
extractSurface1Display.ScaleTransferFunction = 'PiecewiseFunction'
extractSurface1Display.OpacityArray = ['POINTS', 'ux']
extractSurface1Display.OpacityTransferFunction = 'PiecewiseFunction'
extractSurface1Display.DataAxesGrid = 'GridAxesRepresentation'
extractSurface1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
extractSurface1Display.ScaleTransferFunction.Points = [0.005237827916105199, 0.0, 0.5, 0.0, 0.46907547386229526, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
extractSurface1Display.OpacityTransferFunction.Points = [0.005237827916105199, 0.0, 0.5, 0.0, 0.46907547386229526, 1.0, 0.5, 0.0]




# *********************** Find the x coordinate of the first peak x_peak ****************************
xy0=[0,0]
xyN=[0,0]
x,y = passArray(extractSurface1, fn="r_over_x_total_t={}.csv".format(timesteps[0]), PointDataArrays=[ 'Points', 'Result'], CellDataArrays=['Volume'])


lower_hull, upper_hull, x_peak, y_peak, length_x_peak_mean, delta_min, delta_max, xy0, xyN, xmin, xmax =  find_smooth_curve_and_bounds(x, y, x_mean, alpha = 0.05)

fn = "r_over_x_t={}.csv".format(timesteps[0])
Save1DArraysToFile([x, y, lower_hull, upper_hull], fn)
#read files as below:
#with open(fn, 'r') as f:
#	lists = json.load(f)
#    x, y, lower_hull, upper_hull = np.array(lists[0]), np.array(lists[1]), np.array(lists[2]), np.array(lists[3])

plot_graph([[x_peak, x_peak], [x_mean, x_mean]], [[0,0.5], [0,0.5]], \
		   ['first peak', 'center of mass', "min edge", "max edge"], \
		   list_x_fill=[lower_hull[:,0], upper_hull[:,0]], list_y_fill=[lower_hull[:,1], upper_hull[:,1]], \
		   dash=['solid', 'dot', 'dot'], \
		   xtitle="x", ytitle="r", image_name=fn[:-3]+'pdf', mode=['lines', 'lines', 'lines'], \
		   colors=['red', 'black', 'black' ], yrange=[0,0.5], xrange=[xmin - 0.1, xmax + 0.1], \
		   marker_size=1, width=1000, height=500, path='./', yanchor='bottom', y0_anchor=0.01, xanchor='left', x0_anchor=0.3)


# ***************** CLIP A BOX to calculate volume from x_peak to x_mean ****************************
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


fn = "for_excel_table.txt"
if not Path(fn).exists():
	with open(fn, 'w') as f:
		f.write("t	x_tail	x_peak	y_peak	x_mean	x_nose	x_nose_ISC	volume	UmeanV	delta_min	delta_mean	delta_max	delta_min_smooth	delta_max_smooth\n")

with open(fn, 'a') as f:
	f.write("{} {} {} {} {} {} {} {} {} {} {} {} {} {}\n".format(timesteps[0], xy0[0], x_peak, y_peak, x_mean, xyN[0], '?', volume, u_mean[0], delta_min, delta_mean, delta_max, '0', '0' ))
print('Successfully save file:', fn)
# ************************* CUT TAIL ***********************************************
# create a new 'Clip'
clip3 = Clip(Input=extractSurface1)
clip3.ClipType = 'Plane'
clip3.HyperTreeGridClipper = 'Plane'
clip3.Scalars = ['POINTS', 'Result']
clip3.Value = 0.23469079123049602
x_cut = max(x_peak, xy0[0])
# init the 'Plane' selected for 'ClipType'
clip3.ClipType.Origin = [x_cut, 0, 0]

# init the 'Plane' selected for 'HyperTreeGridClipper'
clip3.HyperTreeGridClipper.Origin = [x_cut, 0, 0]

# Properties modified on clip1.ClipType
clip3.ClipType.Normal = [1.0, 0.0, 0.0]

# create a new 'Transform'
transform2 = Transform(Input=clip3)
transform2.Transform = 'Transform'

# Properties modified on transform2.Transform
transform2.Transform.Translate = [-x_peak, 0.0, 0.0]

#Show the tail
Show(transform2, renderView1)
Hide(transform2, renderView1)
# create a new 'Extract Surface'
extractSurface2 = ExtractSurface(Input=transform2)

fn = "{}/res/{}_0_{:04d}.vtp".format(path, outfn + '_tail', iter)
SavePvdFile(fn, extractSurface2, 'tail surface data of bubble')

if volumetric_repr == 1:
	#Show(calculator1, renderView1)
	fn = "{}/res/{}_0_{:04d}.pvtu".format(path, outfn, iter)
	SavePvdFile(fn, calculator1, 'volumetric data of bubble', three_dimension=1)
else:
	fn = "{}/res/{}_0_{:04d}.vtp".format(path, outfn, iter)
	SavePvdFile(fn, extractSurface1, 'surface data of bubble')

print("RenderView update..")
renderView1.CameraViewUp = [0.2, 1, 0]
renderView1.CameraParallelScale = 1.3 #1.4 0.5
renderView1.CenterOfRotation = center
renderView1.CameraFocalPoint = center
renderView1.CameraPosition = [center[0] - 4, 0.6, 4.5]
# update the view to ensure updated data information
renderView1.Update()
print ("end")
fn = path + "/" + picName+  '_t=' + str(timesteps[0]) +'_ux.png'
SaveScreenshot( fn, renderView1,
				ImageResolution=[1900, 1077],
				TransparentBackground=0,
				CompressionLevel='2' )
Hide(extractSurface1)
# ***************** CONTOUR2 for LAMBDA2 l2 ****************************
# create a new 'Contour'


contour2 = Contour(Input=resampleToImage1)
contour2.ContourBy = ['POINTS', 'l2']
contour2.Isosurfaces = [-4, -2.0, -1.0, -0.5, -0.25, -0.125]
contour2.PointMergeMethod = 'Uniform Binning'

# create a new 'Threshold'
threshold0 = Threshold(Input=contour2)
threshold0.Scalars = ['POINTS', 'l2']
threshold0.ThresholdRange = [-1.09, -0.49]

# ***************** LAMBDA2 inside a bubble: ISOVOLUME1 for VOLUME FRACTION f ****************************

# create a new 'Iso Volume'
isoVolume2 = IsoVolume(Input=connectivity1)
isoVolume2.InputScalars = ['POINTS', 'l2']
isoVolume2.ThresholdRange = [-2.0, -0.5]

# create a new 'Connectivity'
connectivity2 = Connectivity(Input=isoVolume2)
connectivity2.RegionIdAssignmentMode = 'Cell Count Descending'

# create a new 'Threshold'
threshold1 = Threshold(Input=connectivity2)
threshold1.Scalars = ['POINTS', 'RegionId']
threshold1.ThresholdRange = [0.0, 10.0]

# create a new 'Programmable Filter'
programmableFilter1 = ProgrammableFilter(Input=connectivity2)
programmableFilter1.Script = """
	regionArray = inputs[0].CellData['RegionId']
	regions = unique(regionArray)
	v = volume(inputs[0])
	a=[]
	for regionValue in regions:
	  regionVolume = sum(v[regionArray == regionValue])
	  a.append((regionValue,regionVolume))
	q=sorted(a, key=lambda x: (x[1],x[0]), reverse=True)
	print(q)
	Nbig=10
	ids=[]
	for n in range(Nbig):
	  ids.append(q[n][0])
	print(ids)
	"""
programmableFilter1.RequestInformationScript = ''
programmableFilter1.RequestUpdateExtentScript = ''
programmableFilter1.PythonPath = ''



# get color transfer function/color map for 'l2'
l2LUT = GetColorTransferFunction('l2')
l2LUT.RGBPoints = [-2.0, 0.0, 1.0, 1.0, -1.55, 0.0, 0.0, 1.0, -1.5, 0.0, 0.0, 0.501960784314, -1.4499999999999997, 1.0, 0.0, 0.0, -1.0, 1.0, 1.0, 0.0]
l2LUT.ColorSpace = 'RGB'
l2LUT.ScalarRangeInitialized = 1.0
#****************** CONNECTIVITY(f) AND CONTOUR2 (lambda2) ********************
# show data from connectivity1
print ("Showing connectivity1.. "),
connectivity1Display = Show(connectivity1, renderView1)

# trace defaults for the display properties.
connectivity1Display.Representation = 'Surface'
connectivity1Display.AmbientColor = [0.0392156862745098, 0.00784313725490196, 1.0]
connectivity1Display.ColorArrayName = ['POINTS', '']
connectivity1Display.DiffuseColor = [0.0392156862745098, 0.00784313725490196, 1.0]
connectivity1Display.Opacity = 0.53
connectivity1Display.Specular = 1.0
connectivity1Display.SpecularPower = 1.0
connectivity1Display.Ambient = 0.21
connectivity1Display.OSPRayScaleArray = 'f'
connectivity1Display.OSPRayScaleFunction = 'PiecewiseFunction'
connectivity1Display.SelectOrientationVectors = 'None'
connectivity1Display.ScaleFactor = 0.8
connectivity1Display.SelectScaleArray = 'f'
connectivity1Display.GlyphType = 'Arrow'
connectivity1Display.GlyphTableIndexArray = 'f'
connectivity1Display.GaussianRadius = 0.04
connectivity1Display.SetScaleArray = ['POINTS', 'f']
connectivity1Display.ScaleTransferFunction = 'PiecewiseFunction'
connectivity1Display.OpacityArray = ['POINTS', 'f']
connectivity1Display.OpacityTransferFunction = 'PiecewiseFunction'
connectivity1Display.DataAxesGrid = 'GridAxesRepresentation'
connectivity1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
connectivity1Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
connectivity1Display.ScaleTransferFunction.Points = [0.5, 0.0, 0.5, 0.0, 0.5001220703125, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
connectivity1Display.OpacityTransferFunction.Points = [0.5, 0.0, 0.5, 0.0, 0.5001220703125, 1.0, 0.5, 0.0]
print ("end")
# show data from contour2
print ("Showing contour2.. "),
contour2Display = Show(threshold0, renderView1, 'GeometryRepresentation')
# trace defaults for the display properties.
contour2Display.Representation = 'Surface'
contour2Display.ColorArrayName = ['POINTS', 'l2']
contour2Display.LookupTable = l2LUT
contour2Display.Opacity = 0.62
contour2Display.Specular = 1.0
contour2Display.OSPRayScaleArray = 'l2'
contour2Display.OSPRayScaleFunction = 'PiecewiseFunction'
contour2Display.SelectOrientationVectors = 'None'
contour2Display.ScaleFactor = 0.2388025760650635
contour2Display.SelectScaleArray = 'l2'
contour2Display.GlyphType = 'Arrow'
contour2Display.GlyphTableIndexArray = 'l2'
contour2Display.GaussianRadius = 0.011940128803253174
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
print ("end")

renderView1.CameraViewUp = [0.2, 1, 0]
renderView1.CameraParallelScale = 1.5 #1.4 0.5
renderView1.CenterOfRotation = center
renderView1.CameraFocalPoint = [center[0]-0.5, 0, 0]
renderView1.CameraPosition = [center[0] - 4, 0.6, 4.5]

fn = path + "/" + picName+  '_t=' + str(timesteps[0]) +'_noLambda2.png'
SaveScreenshot( fn, renderView1,
				ImageResolution=[1900, 1077],
				TransparentBackground=0,
				CompressionLevel='2' )
print('File=' + fn + ' generated succesfully')

fn = "{}/res/{}_0_{:04d}.vtp".format(path, 'lambda2', iter)
SavePvdFile(fn, contour2, 'lambda2 surface data')

#****************** LAMBDA2 inside a bubble  ********************

print ("Showing connectivity1 and hiding contour2.. "),
Hide(threshold0, renderView1) # hide lambda2 all
print ("end")
# trace defaults for the display properties.
connectivity1Display.Representation = 'Surface'
connectivity1Display.Opacity = 0.2

# show data from threshold1
print ("Showing threshold1.. "),
threshold1Display = Show(threshold1, renderView1, 'GeometryRepresentation')
# get color transfer function/color map for 'l2'
l2LUT = GetColorTransferFunction('l2')
l2LUT.RGBPoints = [-1.0, 0.054901960784313725, 0.9411764705882353, 0.12941176470588237, -0.75, 0.865, 0.865, 0.865, -0.5, 1.0, 1.0, 0.0]
l2LUT.ScalarRangeInitialized = 1.0
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
threshold1Display.GaussianRadius = 0.011940128803253174
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
print ("end")
l2PWF = GetOpacityTransferFunction('l2')
l2PWF.Points = [-1.0, 0.0, 0.5, 0.0, -0.5, 1.0, 0.5, 0.0]
l2PWF.ScalarRangeInitialized = 1

renderView1.CameraViewUp = [0.2, 1, 0]
renderView1.CameraParallelScale = 1.2
renderView1.CenterOfRotation = center
renderView1.CameraFocalPoint = center
renderView1.CameraPosition = [center[0] - 4, 0.6, 4.5]
print("RenderView update..")
renderView1.Update()

fn = path + "/" + picName+  '_t=' + str(timesteps[0]) +'_Lambda2_in_bubble.png'
SaveScreenshot( fn, renderView1,
				ImageResolution=[1900, 1077],
				TransparentBackground=0,
				CompressionLevel='2' )
print('File=' + fn + ' generated succesfully')

# create a new 'Extract Surface'
extractSurface2 = ExtractSurface(Input=threshold1)
Show(extractSurface2, renderView1, 'GeometryRepresentation')
Hide(extractSurface2, renderView1)

fn = "{}/res/{}_0_{:04d}.vtp".format(path, 'lambda2_in_bubble', iter)
SavePvdFile(fn, extractSurface2, 'lambda2 in bubble surface data')

# ******************************SLICE bubble*******************************************
# *************************************************************************************
# create a new 'Slice'
slice1 = Slice(Input=resampleToImage1)
slice1.SliceType = 'Plane'
slice1.HyperTreeGridSlicer = 'Plane'
slice1.SliceOffsetValues = [0]

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Origin = [0.5*(len_min + len_max), 0.0, 0.0]
slice1.SliceType.Normal = [0.0, 0.0, 1.0]
# init the 'Plane' selected for 'HyperTreeGridSlicer'
# slice1.HyperTreeGridSlicer.Origin = [10.999999761581421, 0.0, 0.0] #???
print ("Showing slice1Display of bubble.. ")
slice1Display = Show(slice1, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
slice1Display.Representation = 'Surface'
slice1Display.ColorArrayName = ['POINTS', 'u.x']
slice1Display.LookupTable = uxLUT
slice1Display.OSPRayScaleArray = 'f'
slice1Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice1Display.SelectOrientationVectors = 'None'
slice1Display.ScaleFactor = 0.3379646740552927
slice1Display.SelectScaleArray = 'f'
slice1Display.GlyphType = 'Arrow'
slice1Display.GlyphTableIndexArray = 'f'
slice1Display.GaussianRadius = 0.016898233702764633
slice1Display.SetScaleArray = ['POINTS', 'f']
slice1Display.ScaleTransferFunction = 'PiecewiseFunction'
slice1Display.OpacityArray = ['POINTS', 'f']
slice1Display.OpacityTransferFunction = 'PiecewiseFunction'
slice1Display.DataAxesGrid = 'GridAxesRepresentation'
slice1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
slice1Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
slice1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0000000000000002, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
slice1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0000000000000002, 1.0, 0.5, 0.0]

fn = "{}/res/{}_0_{:04d}.vtp".format(path, 'slice', iter)
SavePvdFile(fn, slice1, 'resampleToImage1 slice data')

# ******************************SLICE bubble*******************************************
# *************************************************************************************

# create a new 'Contour'
contour3 = Contour(Input=slice1)
contour3.ContourBy = ['POINTS', 'f']
contour3.Isosurfaces = [0.5]
contour3.PointMergeMethod = 'Uniform Binning'
Show(contour3, renderView1)

# create a new 'Pass Arrays'
passArrays1 = PassArrays(Input=contour3)
passArrays1.PointDataArrays = [ 'Points']
passArrays1.CellDataArrays = ['Volume']

# update the view to ensure updated data information
spreadSheetView1.Update()

ss_data = paraview.servermanager.Fetch(passArrays1)
Np = ss_data.GetNumberOfPoints()
print('Np=', Np)
xr = []
for ip in range(Np):
	xp = ss_data.GetPoint(ip)[0]
	yp = ss_data.GetPoint(ip)[1]
	xr.append((xp,yp))
xr = np.array(xr)
print('processing data size of x and y:', xr.shape[0])

fn = "slice_t={}.csv".format(timesteps[0])
Save1DArraysToFile([xr[:,0], xr[:,1]], fn)

# Freeing Memory
Delete(slice1)
del slice1
Delete(extractSurface2)
del extractSurface2
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
