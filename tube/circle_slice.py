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
#import __builtin__
# from inspect import getmodule
#
# print(getmodule(Show))
vtk_from_pvpython=True # pvpython reads from file, otherwise from paraview GUI

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
save_isosurfacepvd = PVDReader(FileName=path + '/' + infn)
#save_isosurfacepvd.CellArrays = ['vtkGhostType']
#save_isosurfacepvd.PointArrays = ['absOmega', 'f', 'fs', 'l', 'l2', 'omega', 'p', 'u.x', 'vtkValidPointMask', 'vtkGhostType', 'Result']

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

# create a new 'Slice'
slice1 = Slice(Input=save_isosurfacepvd)
slice1.SliceType = 'Cylinder'
slice1.HyperTreeGridSlicer = 'Plane'
slice1.SliceOffsetValues = [0.0]

# init the 'Cylinder' selected for 'SliceType'
slice1.SliceType.Axis = [1.0, 0.0, 0.0]
slice1.SliceType.Radius = 0.25

# init the 'Plane' selected for 'HyperTreeGridSlicer'
slice1.HyperTreeGridSlicer.Origin = [0, 0, 0]

# show data from slice1
slice1Display = Show(slice1, renderView1, 'GeometryRepresentation')

x = []
y = []
z = []
for i in timeList:
	print("iteration:{} time={}".format(i, timesteps[i]))
	SetActiveView(renderView1)
	# Properties modified on animationScene1
	animationScene1.AnimationTime = timesteps[i]
	# Properties modified on timeKeeper1
	timeKeeper1.Time = timesteps[i]
	
	boundsDomain = save_isosurfacepvd.GetDataInformation().GetBounds()

	lDomain = boundsDomain[1] - boundsDomain[0]
	x_mean = 0.5*(boundsDomain[1] + boundsDomain[0])
	x_tail = boundsDomain[0]
	x_nose = boundsDomain[1]
	print("boundsDomain of my source=", boundsDomain, "lDomain=", lDomain, "x_mean=", x_mean, "x_tail=", x_tail, "x_nose=", x_nose)

	# create a new 'Clip'
	clip1 = Clip(Input=slice1)
	clip1.ClipType = 'Plane'
	clip1.HyperTreeGridClipper = 'Plane'
	clip1.Scalars = ['POINTS', 'Result']
	clip1.Value = 0

	# init the 'Plane' selected for 'ClipType'
	clip1.ClipType.Origin = [x_mean, 0, 0]

	# init the 'Plane' selected for 'HyperTreeGridClipper'
	clip1.HyperTreeGridClipper.Origin = [x_mean, 0, 0]
	
	# show data from save_isosurfacepvd
	clip1Display = Show(clip1, renderView1, 'GeometryRepresentation')

	# ----------------------------------------------------------------
	# setup the visualization in view 'renderView1'
	# ----------------------------------------------------------------
	fn = path + '/circle_slice.csv'
	SaveData(fn, proxy=clip1,\
		 CellDataArrays=[], FieldDataArrays=[], AddTime=0)
	tab = pd.read_csv(fn, usecols=["Points:0","Points:1","Points:2", "u.x:0","u.x:1","u.x:2"])
	sort_tab(tab) # calculate theta
	x_tail_ref = tab['Points:0'].max()
	tab['x_tail_rel_nose'] = (tab['Points:0'] - x_nose).abs()
	tab['x_tail_rel'] = (tab['Points:0'] - x_tail_ref).abs()
	tab['t'] = timesteps[i]
	print(tab.head())
#'''
#	if i == 0:
#		with pd.ExcelWriter(fn) as writer:  
#			tab.to_excel(writer, index=False, sheet_name='{}'.format(timesteps[i]))
#	else:
#		with pd.ExcelWriter(fn, engine='openpyxl', mode='a') as writer:  
#			tab.to_excel(writer, index=False, sheet_name='{}'.format(timesteps[i]))
#'''

	xi = tab['t'].values.tolist()
	yi = tab['theta'].values.tolist()
	zi = tab['x_tail_rel'].values.tolist()
#	x += xi
#	y += yi
#	z += zi

	fn = path + '/piece_table_t={}.csv'.format(timesteps[i])
	tab.to_csv(fn, index=False)

#df = pd.DataFrame({'t':x, 'theta':y, 'x_tail_rel':z})
#fn = path + '/table_x_tail_rel_over_t_theta.csv'
#df.to_csv(fn, index=False)

#fig = go.Figure(data=[go.Mesh3d(x=x, y=y, z=z, color='yellowgreen', opacity=0.9)])
#fig.update_layout(scene = dict(
#                    xaxis_title='t',
#                    yaxis_title='\u03B8',
#                    zaxis_title='|x<sub>tail</sub> - x<sub>nose</sub>|'),
#                    width=600, height=550,
#                    margin=dict(r=0, b=0, l=0, t=0, pad=0.0),
#					autosize=False)
# Default parameters which are used when `layout.scene.camera` is not provided
#camera = dict(
#    up=dict(x=0, y=0, z=1),
#    center=dict(x=0, y=0, z=0),
#    eye=dict(x=0, y=1.25, z=1.25)
#)

#fig.update_layout(scene_camera=camera)

#fig.show()
#fn = path + "/surface_x_tail_over_t_l.pdf"
#fig.write_image(str(Path(fn)), engine="kaleido")	

