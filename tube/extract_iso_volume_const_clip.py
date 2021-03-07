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
#from vtk.numpy_interface import dataset_adapter as dsa
# import vtk.numpy_interface.algorithms as algs
import vtk
#import numpy as np
import glob, os, sys
import logging
from sys import argv
import timeit
import argparse
#from vtk.numpy_interface import dataset_adapter as dsa
logging.basicConfig(format='%(message)s')
log = logging.getLogger(__name__)

import functools
import __builtin__
# from inspect import getmodule
#
# print(getmodule(Show))
vtk_from_pvpython=True # pvpython reads from file, otherwise from paraview GUI
# vtk_from_pvpython=False # pvpython reads from file, otherwise from paraview GUI

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
                    nargs='?', default='iso_volume')
required.add_argument("-maxlevel", type=int, help="Provide the maximum level of refinement",
                    nargs='?', default=10, required=True)
required.add_argument("-iter", type=int, help="Provide the iter argument level of refinement",
                    nargs='?', default=0, required=True)
required.add_argument("-volumetric_repr", type=int, help="Provide the volumetric_repr argument 1 if volumetric repr, 0 only surface",
                    nargs='?', default=0, required=True)

# parser.add_argument('--foo', action='store_const', const=2, default=42)
args = parser.parse_args()
print(args)
infn = args.infn
outfn = args.outfn
maxlevel = args.maxlevel
iter = args.iter
volumetric_repr = args.volumetric_repr
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
spreadSheetView1.ColumnToSort = ''
spreadSheetView1.BlockSize = 1024L
# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------
# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.SplitHorizontal(0, 0.5)
layout1.AssignView(1, renderView1)

print ("end")
# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)

# defining of computational domains
Show(my_source, renderView1)
boundsDomain = my_source.GetDataInformation().GetBounds()

lDomain = boundsDomain[1] - boundsDomain[0]
print("boundsDomain of my source=", boundsDomain, " lDomain=", lDomain)


SetActiveView(renderView1)
# Properties modi bfied on animationScene1
animationScene1.AnimationTime = timesteps[0]
# Properties modified on timeKeeper1
timeKeeper1.Time = timesteps[0]
# ***************** CUT remain only a TUBE ****************************


# ***************** CLIP A BOX of a BOX ****************************
# create a new 'Clip'
clip1 = Clip(Input=my_source)
clip1.ClipType = 'Box'
clip1.HyperTreeGridClipper = 'Plane'
clip1.Scalars = ['POINTS', 'f']
clip1.Value = 0.5
# init the 'Box' selected for 'ClipType'
clip1.ClipType.Position = [0, -0.6, -0.6]
clip1.ClipType.Length = [lDomain, 1.2, 1.2]

# ***************** CELL DATA TO POINT DATA  from BOX of BOX ****************************
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
Nx = int(lDomain/dx)
Nyz = int(1.0/dx)
resampleDimension = [Nx, Nyz, Nyz]
print("resampleDimension:", resampleDimension)
resampleToImage1.SamplingDimensions = resampleDimension
resampleToImage1.SamplingBounds = [0, lDomain, -0.5, 0.5, -0.5, 0.5]

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

print("Connectivity1 finished")
#fn = path + '/res/save_isovolume_0_' + iter + '.pvtu'
#SaveData(fn, proxy=connectivity1, 
#        PointDataArrays=['f', 'fs', 'l', 'l2', 'omega', 'u.x', 'vtkGhostType', 'vtkValidPointMask'],
#	CellDataArrays=['vtkGhostType'],
#	WriteTimeSteps=1
#)

if volumetric_repr == 1:
	Show(connectivity1, renderView1)
	fn = "{}/res/{}_0_{:04d}.pvtu".format(path, outfn, iter)
	SaveData(fn, proxy=connectivity1, 
		#PointDataArrays=['f', 'fs', 'l', 'l2', 'omega', 'u.x'],
		UseSubdirectory=0
	)
	print("Saved volumetric data of bubble:", fn)
else:
	# create a new 'Extract Surface'
	extractSurface1 = ExtractSurface(Input=connectivity1)
	Show(extractSurface1, renderView1)
	fn = "{}/res/{}_0_{:04d}.vtp".format(path, outfn, iter)
	SaveData(fn, proxy=extractSurface1
		#PointDataArrays=['f', 'fs', 'l', 'l2', 'omega', 'u.x']
	)
	print("Saved surface data of bubble:", fn)



# Freeing Memory
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
