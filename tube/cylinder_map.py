# state file generated using paraview version 5.8.0
#### import the simple module from the paraview
from paraview.simple import *
from paraview.vtk.util import numpy_support # provides unique()
from vtkmodules.numpy_interface.algorithms import * # provides volume()
from vtk.numpy_interface import dataset_adapter as dsa
# import vtk.numpy_interface.algorithms as algs
import vtk
import numpy as np
import glob, os, sys
import logging
from sys import argv
import timeit
import argparse
import pandas as pd
from vtk.numpy_interface import dataset_adapter as dsa
logging.basicConfig(format='%(message)s')
log = logging.getLogger(__name__)
from pathlib import Path

import plotly.graph_objects as go
import functools
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib import ticker, cm
from mpl_toolkits.axes_grid1 import make_axes_locatable

from scipy.interpolate import griddata

# plt.rcParams["font.family"] = "Times New Roman"
#import __builtin__
# from inspect import getmodule
#
# print(getmodule(Show))
vtk_from_pvpython=True # pvpython reads from file, otherwise from paraview GUI
delta_sym = '\u03B4'
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
def passArray(input_data, fn, PointDataArrays=[ 'Points'], CellDataArrays=['Volume']):
	# create a new 'Pass Arrays'
	passArrays1 = PassArrays(Input=input_data)
	passArrays1.PointDataArrays = PointDataArrays
	passArrays1.CellDataArrays = CellDataArrays

	ss_data = paraview.servermanager.Fetch(passArrays1)
	Np = ss_data.GetNumberOfPoints()
	print('Np=', Np)
	xp = np.zeros(Np)
	yp = np.zeros(Np)
	zp = np.zeros(Np)
	print("gathering points is started...")

	for ip in range(Np):
		xp[ip] = ss_data.GetPoint(ip)[0]
		yp[ip] = ss_data.GetPoint(ip)[1]
		zp[ip] = ss_data.GetPoint(ip)[2]
	print("gathering points is finished...")

	tab = pd.DataFrame({'x':xp, 'y':yp, 'z':zp})
	tab['r'] = np.sqrt(tab['y']**2 + tab['z']**2)

	tab['theta'] = np.arctan2(zp, yp)
	tab['theta'][tab['theta']<0] += 2*np.pi # it ranges [0, 2*pi]
	tab.to_csv(fn, index=False)
	return tab

def my_ceil(a, precision=0):
	return np.round(a + 0.5 * 10**(-precision), precision)

def my_floor(a, precision=0):
	return np.round(a - 0.5 * 10**(-precision), precision)
# Create and show figure
'''
x,y are 1D array
z is 2D array
'''
def plot_2D_color(xx, yy, zz, colorbartitle, xtitle, ytitle, title, path, fn, Nround=2,
				  colorscale='Jet', width=500, height=200, inds=None, symmetricalcolorscale=False,
				  zrange=None):
	if inds is None or len(inds) != 4:
		print("inds is all")
		inds = [0, len(xx), 0, len(yy)]

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
	x = np.copy(xx)
	y = np.copy(yy)
	z = np.copy(zz)
	x = x[inds[0]:inds[1]]
	y = y[inds[2]:inds[3]]
	z = z[inds[2]:inds[3], inds[0]:inds[1]]
	print(len(x), z.shape)
	fig = go.Figure()

	fig.add_trace(go.Heatmap(x=x,
							 y=y,
							 z=z,
							 colorscale=colorscale,
#							  colorbar=dict(
#										 title=delta_sym+"x",
#										 titlefont=dict(
#											 size=14,
#											 family='Times New Roman')
#									  )
							 colorbar=dict(
								title=colorbartitle,
								tickvals=np.linspace(my_floor(z.min(), Nround), my_ceil(z.max(), Nround), 5),
								thickness=25,
								titleside="top",
								tickmode="array",
								ticks="outside",
								titlefont=dict(
									size=25,
									family='Times New Roman')
							 )
	))

	layout = go.Layout(xaxis=go.layout.XAxis(
		title=go.layout.xaxis.Title(
			text=xtitle,
			font=dict(size=25,
					  family='Times New Roman')
		)),
	yaxis=go.layout.YAxis(
		title=go.layout.yaxis.Title(
			text=ytitle,
			font=dict(size=25,
					  family='Times New Roman')
		)
	))

	fig.update_layout(layout)

	fig.update_layout(font_family="Times New Roman",
					  width=width, height=height,
					  margin=dict(r=0, b=0, l=0, t=0, pad=0.0),
					  autosize=False)
	zmin = my_floor(z.min(), Nround)
	zmax = my_ceil(z.max(), Nround)
	print('zminmax?:', zmin, zmax)
	if symmetricalcolorscale:
		tmp = max(abs(zmin), abs(zmax))
		zmin = -tmp
		zmax = tmp
	if zrange!=None:
		zmin = zrange[0]
		zmax = zrange[1]
	print('zminmax exact:', zmin, zmax)
	fig.data[0].update(zmin=zmin, zmax=zmax)
	fig.update_layout(
		title=title,
		autosize=False,
		margin=dict(
				l=0,
				r=0,
				b=0,
				t=0,
				pad=0
			),
		#	 paper_bgcolor="LightSteelBlue"
	)
	fig.update_layout(
		yaxis = yaxis,
		xaxis = xaxis,
#		 showlegend=True
	)
	fig.show()
	fn = path + fn
	print("Successfully saved files:", fn)
	fig.write_image(str(Path(fn)), engine="kaleido")

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


def eprint(var):
	log.warning(var)

def sort_tab(tab):
	x = tab['Points:0']
	y = tab['Points:1']
	z = tab['Points:2']
	theta = np.arctan2(z, y)
	theta[theta<0] += 2*np.pi # it ranges [0, 2*pi]
	tab['theta'] = theta
	tab.sort_values(by=['theta'], inplace=True)

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [914, 491]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.CenterOfRotation = [13.507097721099854, 1.2822449207305908e-05, 2.4259090423583984e-05]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [8.4737468746048, 2.203557788303081, 3.1096608416862654]
renderView1.CameraFocalPoint = [13.507097721099857, 1.2822449207046548e-05, 2.4259090423096053e-05]
renderView1.CameraViewUp = [0.11333514905529504, 0.8879455127959932, -0.44576665453359654]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 1.6340497078444391
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.OSPRayMaterialLibrary = materialLibrary1

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.AssignView(0, renderView1)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------


start = timeit.default_timer()

parser = argparse.ArgumentParser()
parser._action_groups.pop()
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')

optional.add_argument("-infn", type=str, help="Provide the name of the input paraview files, please",
					nargs='?', default='save_isosurface.pvd')
parser.add_argument("-timeList", type=int, help="Provide the list of step, please",
					nargs='+', default=[0])
parser.add_argument("-nt", type=int, help="Provide the list of time step, please",
					nargs='?', default=0)

# parser.add_argument('--foo', action='store_const', const=2, default=42)
args = parser.parse_args()
print('arguments are:', args)
infn = args.infn
timeList = args.timeList
nt = args.nt
# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

#Current PATH reading
if vtk_from_pvpython:
	path = os.path.abspath(os.getcwd())
else:
	path = '/home/e.sharaborin/basilisk/work/tube/'


# create a new 'PVD Reader'
print("Reading source file:", path + '/' + infn)
save_isosurfacepvd = PVDReader(FileName=path + '/' + infn)
#save_isosurfacepvd.CellArrays = ['vtkGhostType']
#save_isosurfacepvd.PointArrays = ['absOmega', 'f', 'fs', 'l', 'l2', 'omega', 'p', 'u.x', 'vtkValidPointMask', 'vtkGhostType', 'Result']

# Create a new 'SpreadSheet View'
spreadSheetView1 = CreateView('SpreadSheetView')


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

if len(timeList) == 1 and timeList[0] == 0:
	timeList = range(0,NT)
if nt != 0:
	timeList = range(0,NT,nt)
if NT != 0:
	print("NT=", NT, " timeList=", timeList)
else:
	print("ERROR: can't read pvd file NT=0")

# Compute areas and colors
N = 150
ngridx = 1000
ngridy = 200
vmin=0
vmax=0.15
Nc = 200
v1 = np.linspace(vmin, vmax, 10)

for i in timeList:
	t = timesteps[i]
	t_round = np.round(t, 3)
	print("iteration:{} time={}".format(i, t))
	SetActiveView(renderView1)
	# Properties modified on animationScene1
	animationScene1.AnimationTime = t
	# Properties modified on timeKeeper1
	timeKeeper1.Time = t
	# update the view to ensure updated data information
	spreadSheetView1.Update()
	Show(save_isosurfacepvd)

	boundsDomain = save_isosurfacepvd.GetDataInformation().GetBounds()
	boundsDomain = save_isosurfacepvd.GetDataInformation().GetBounds()

	lDomain = boundsDomain[1] - boundsDomain[0]
	x_mean = 0.5*(boundsDomain[1] + boundsDomain[0])
	x_tail = boundsDomain[0]
	x_nose = boundsDomain[1]
	print("boundsDomain of my source=", boundsDomain, "lDomain=", lDomain, "x_mean=", x_mean, "x_tail=", x_tail, "x_nose=", x_nose)

	# ----------------------------------------------------------------
	# setup the visualization in view 'renderView1'
	# ----------------------------------------------------------------
	fn = path + '/xyz_total_t={}.csv'.format(t)
	tab = passArray(save_isosurfacepvd, fn, PointDataArrays=[ 'Points'], CellDataArrays=['Volume'])


	delta = 0.5 - tab['r'].values
	args = delta < vmax
	delta = delta[args]

	x = tab['x'].values[args]
	xmin, xmax = x.min(), x.max()
	print("xmin={} xmax={}".format(xmin, xmax))
	r = tab['r'].values[args]
	theta = tab['theta'].values[args]
	area = 20 * r**2

	deltamin, deltamax = delta.min(), delta.max()
	print("deltamin={} deltamax={}".format(deltamin, deltamax))
	colors = delta

	fig = plt.figure(figsize=(6,5))
	ax = fig.add_subplot(111, projection='polar')
	c = ax.scatter(theta, r, c=colors, s=area, cmap='hsv', alpha=0.75)
	fig.suptitle('t={}'.format(t_round), fontsize=16)
	fn = path + '/xyz_total_t_polar={}.png'.format(t)
	fig.savefig(fn, rasterized=True, dpi=300)

	limits = [[x_tail, x_nose], [0, x_nose - x_tail]]
	names = ['xyz_total_tric', 'xyz_total_rel_tric']
	ylabel = ['$x$', '$|x - x_{nose}|$']
	for k in range(2):
		xmin_, xmax_ = limits[k][0], limits[k][1]
		fig, ax = plt.subplots()
		ax.set_title('t={}'.format(t_round), fontsize=18)
		ax.set_xlabel('\u03B8', fontsize=16)
		ax.set_ylabel(ylabel[k], fontsize=16)
		levels = np.linspace(vmin, vmax, Nc + 1)
		if k==0:
			xx = x
		else:
			xx = x_nose - x
	# 	ax.tricontour(theta, xx, delta, levels=levels, linewidths=0.5, colors='k', locator=ticker.LogLocator())
		cntr = ax.tricontourf(theta, xx, delta, levels=levels, cmap="hsv", vmin=vmin, vmax=vmax)
		divider = make_axes_locatable(ax)
		cax = divider.append_axes("right", size="5%", pad=0.05)
		cbar = fig.colorbar(cntr, cax=cax, ticks=v1)
		cbar.ax.set_yticklabels(["{:4.2f}".format(i) for i in v1]) # add the labels
		cbar.ax.set_title("$\delta$", ha='left', x=0)
		ax.set(xlim=(0, 2*np.pi), ylim=(xmin_, xmax_))

		fn = path + '/{}_t={}.png'.format(names[k], t)
		plt.savefig(fn, rasterized=True, dpi=300)


	points = np.c_[theta, x]
	values = delta
	grid_x, grid_y = np.mgrid[0:2*np.pi:Nc*1j, xmin:xmax:Nc*1j]
	grid_z2 = griddata(points, values, (grid_x, grid_y), method='cubic')

	limits = [[xmin, xmax], [0, xmax - xmin]]
	names = ['xyz_total_imshow', 'xyz_total_rel_imshow']
	ylabel = ['$x$', '$|x - x_{nose}|$']
	for k in range(2):
		xmin_, xmax_ = limits[k][0], limits[k][1]
		title = 't={}'.format(t_round)
		fig, ax = plt.subplots()
		ax.set_title(title, fontsize=18)
		ax.set_xlabel('\u03B8', fontsize=16)
		ax.set_ylabel(ylabel[k], fontsize=16)
		print('xmin={}, xmax={}'.format(xmin_, xmax_))
		pos = ax.imshow(grid_z2.T, extent=(0, 2*np.pi, xmin_, xmax_), origin='lower', vmin=vmin, vmax=vmax)
		# create an axes on the right side of ax. The width of cax will be 5%
		# of ax and the padding between cax and ax will be fixed at 0.05 inch.
		divider = make_axes_locatable(ax)
		cax = divider.append_axes("right", size="5%", pad=0.05)


		cbar=fig.colorbar(pos, cax=cax, ticks=v1)
		cbar.ax.set_yticklabels(["{:4.2f}".format(i) for i in v1]) # add the labels
		cbar.ax.set_title("$\delta$", fontsize=16)
		ax.set(xlim=(0, 2*np.pi), ylim=(xmin_, xmax_))
		fn = path + '/{}_t={}.png'.format(names[k], t)
		fig.savefig(fn, rasterized=True, dpi=300)

		# fn = '/xyz_plotly_{}_t={}.png'.format(names[k], t)
		# xx, yy = np.linspace(0, 2*np.pi, Nc), np.linspace(xmin_, xmax_, Nc)
		# plot_2D_color(xx, yy, grid_z2.T, '\u03B4', '\u03B8', ylabel[k], title, path, fn, Nround=2,
		# 				  colorscale='Jet', width=500, height=600, inds=None, symmetricalcolorscale=False,
		# 				  zrange=None)

