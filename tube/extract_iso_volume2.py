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
from contextlib import redirect_stdout

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
parser.add_argument("-infn", type=str, help="Provide the name of the input paraview files, please",
                    nargs='?', default='*.pvd')
parser.add_argument("-outfn", type=str, help="Provide the name of the output paraview files, please",
                    nargs='?', default='iso_volume.pvd')
parser.add_argument("-maxlevel", type=int, help="Provide the maximum level of refinement",
                    nargs='?', default=10)

# parser.add_argument('--foo', action='store_const', const=2, default=42)
args = parser.parse_args()
print(args)
infn = args.infn
maxlevel = args.maxlevel

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
for i in timesteps:
	# Properties modified on animationScene1
	animationScene1.AnimationTime = timesteps[i]
	# Properties modified on timeKeeper1
	timeKeeper1.Time = timesteps[i]
	# ***************** CUT remain only a TUBE ****************************
	# create a new 'Clip' a long box
	clip0 = Clip(Input=my_source)
	clip0.ClipType = 'Box'
	clip0.HyperTreeGridClipper = 'Plane'
	clip0.Scalars = ['POINTS', '']
	clip0.Value = 0.5
	# init the 'Box' selected for 'ClipType'
	clip0.ClipType.Position = [0, -0.6, -0.6]
	clip0.ClipType.Length = [lDomain, 1.2, 1.2]

	# ***************** COARSE ESTIMATION of CENTER Xcg Umean ****************************
	# create a new 'Iso Volume'
	isoVolume0 = IsoVolume(Input=clip0)
	isoVolume0.InputScalars = ['POINTS', 'f']
	isoVolume0.ThresholdRange = [0.0, 0.5]

	 # create a new 'Cell Data to Point Data'
	cellDatatoPointData0 = CellDatatoPointData(Input=isoVolume0)

	# create a new 'Integrate Variables'
	integrateVariables1 = IntegrateVariables(Input=cellDatatoPointData0)

	# show data from integrateVariables1
	integrateVariables1Display = Show(integrateVariables1, spreadSheetView1, 'SpreadSheetRepresentation')

	# create a new 'Pass Arrays'
	passArrays1 = PassArrays(Input=integrateVariables1)
	passArrays1.PointDataArrays = [ 'Points']
	passArrays1.CellDataArrays = []

	# update the view to ensure updated data information
	spreadSheetView1.Update()

	ss_data = Fetch(passArrays1)
	print('N=', ss_data.GetNumberOfPoints())
	#         x_sum = ss_data.GetPointData().GetArray('Point').GetValue(0)
	x_mean = ss_data.GetPoint(0)[0]
	#         u_mean = ss_data.GetPointData().GetArray('u.x').GetValue(0)
	#         print('u_mean:',u_mean)
	s = "time: {} i: {} coarse_x_mean: {}  ".format(timesteps[i], 0, x_mean)
	print (s)


	del ss_data
	Delete (passArrays1)
	del passArrays1
	Delete (integrateVariables1)
	del integrateVariables1
	Delete (cellDatatoPointData0)
	del cellDatatoPointData0
	Delete(isoVolume0)
	del isoVolume0

	bounds = [x_mean - 1.7, x_mean + 1.7, -0.5, 0.5, -0.5, 0.5  ]
	print("bounds_rude=", bounds)
	center = [(bounds[0] + bounds[1])/2, (bounds[2] + bounds[3])/2,(bounds[4] + bounds[5])/2]
	print('center_rude=', center)
	len_bub = bounds[1] - bounds[0]
	len_min = max([bounds[0] - 2, 0.5])
	len_max = min([bounds[1] + 2, lDomain - 0.5])

	length = len_max - len_min
	print('COARSE ESTIMATION: len_min=',len_min,' len_max=',len_max,' len=',length, ' len_bub=', len_bub)
	if abs(length) > 1e30:
	    print ("ERROR in reading: out of box... go to the next time step...")
	# *****************END OF A COARSE ESTIMATION ********

	# ***************** CLIP A BOX of a BOX ****************************
	# create a new 'Clip'
	clip1 = Clip(Input=clip0)
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
	Nx = int(lDomain/dx)//10
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
	Show(connectivity1, renderView1)
	print("Connectivity1 finished")
	fn = path + '/save_isovolume.pvd'
	#SaveData(fn, proxy=connectivity1, 
	#	PointDataArrays=['f', 'fs', 'l', 'l2', 'omega', 'u.x', 'vtkGhostType', 'vtkValidPointMask'],
	#	CellDataArrays=['vtkGhostType'],
	#	WriteTimeSteps=1
	#)
	# save data
	
	#fn = path + '/save_isovolume.pvd'
	#with open(fn, 'w') as f:
	#    with redirect_stdout(f):
	#	print("<VTKFile type=\"Collection\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n<Collection>")
	#	for t in enumerate(timesteps):
 	#	    print("<DataSet timestep=\"" + t + "\" part=\"0\" file=\"save_isovolume/save_isovolume_0_" + i + ".vtu"/>")
	#	print(" </Collection>\n</VTKFile>")

	#fn = path + '/save_isovolume.pvtu'
	#SaveData(fn, 
	#	proxy=connectivity1,
	#	#PointDataArrays=['f', 'fs', 'l', 'l2', 'omega', 'u.x'],
	#   	Writealltimestepsasfileseries=1,
	#    	Lasttimestep=9
	#)
	# save data
	fn = path + '/save_isovolume.pvtu'
	SaveData('/home/e.sharaborin/basilisk/work/tube/res/as2.pvtu', proxy=cellDatatoPointData1, PointDataArrays=['f', 'fs', 'l', 'l2', 'omega', 'u.x'],
    		UseSubdirectory=0
	)
	print("Saved data of bubble:", fn)

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
	Delete(clip0)
	del clip0



stop = timeit.default_timer()
total_time = stop - start

# output running time in a nice format.
mins, secs = divmod(total_time, 60)
hours, mins = divmod(mins, 60)


sys.stdout.write("Total running time: %d:%d:%d.\n" % (hours, mins, secs))
