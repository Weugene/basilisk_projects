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
import vtk
import numpy as np
import glob, os, sys
import logging
from sys import argv
import timeit
import argparse
from vtk.numpy_interface import dataset_adapter as dsa

logging.basicConfig(format='%(message)s')
log = logging.getLogger(__name__)

import functools
from math import sin
import __builtin__
# from inspect import getmodule
#
# print(getmodule(Show))

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
#1 video Name
#2 pvd file name (by defaults in the first file in a current directory)
#3 pvtu is swithed off by default
# ---------------------------------------------------------------------------------------------------------
start = timeit.default_timer()

parser = argparse.ArgumentParser()
parser.add_argument("-vidName", type=str, help="Provide the name for the outputed video, please",
                    nargs='?', default='lambda2')
parser.add_argument("-picName", type=str, help="Provide the name for the outputed pictures, please",
                    nargs='?', default='lambda2')
parser.add_argument("-filenames", type=str, help="Provide the name of the outputing paraview files, please",
                    nargs='?', default='*.pvd')
parser.add_argument("-frameRate", type=int, help="Provide the frame rate for the video, please",
                    nargs='?', default=10)
parser.add_argument("-frameWindow", type=int, help="Provide the frame window (can be overwritten further), please",
                    nargs='+', default=[0,0])
parser.add_argument("-timeList", type=int, help="Provide the list of step, please",
                    nargs='+', default=[0])
parser.add_argument("-nt", type=int, help="Provide the list of time step, please",
                    nargs='?', default=0)
parser.add_argument("-maxlevel", type=int, help="Provide the maximum level of refinement",
                    nargs='?', default=12)
parser.add_argument("-nslices", type=int, help="Provide the number of slices",
                    nargs='?', default=10)

parser.add_argument("-viewSize", type=int, help="Provide the view size of the output pictures, please",
                    nargs='+', default=[3108, 1168])
parser.add_argument("-noVideo", type=bool, help="Provide the no video Mode",
                    nargs='?', default=False)
parser.add_argument("-noPic", type=bool, help="Provide the no Picture Mode",
                    nargs='?', default=False)
parser.add_argument("-noData", type=bool, help="Provide the no Data Exporting Mode",
                    nargs='?', default=False)
parser.add_argument("-noLambda2", type=bool, help="Provide the no Lambda2 Mode",
                    nargs='?', default=False)
# parser.add_argument('--foo', action='store_const', const=2, default=42)
args = parser.parse_args()
print(args)
filenames = args.filenames
vidName = args.vidName
picName = args.picName
frameRate = args.frameRate
frameWindow = args.frameWindow
maxlevel = args.maxlevel
viewSize = args.viewSize
noVideo = args.noVideo
noPic = args.noPic
noData = args.noData
noLambda2 = args.noLambda2
timeList = args.timeList
nt = args.nt
nslices = args.nslices

if len(frameWindow) != 2 or frameWindow[0]<0 or frameWindow[1] < 0:
    eprint('Error in frameWindow' + frameWindow)
    sys.exit()

#Current PATH reading
path = os.path.abspath(os.getcwd())
eprint("Current PATH=" + path)

if filenames[-5::] == '.pvtu':
# Find files with *.pvtu extension
    numbers = []
    file = ""
    for file in glob.glob(filenames):
        numbers.append(int(filter(lambda x: x.isdigit(), file)))
    file=file[0:-9]
    numbers.sort()
    N = len(numbers)
    filenames = []
    for i in numbers:
        filenames.append('{}/{}{:04d}.pvtu'.format(path, file, i))
    print(filenames)
    fn = path +  '/' + filenames
    my_source = XMLPartitionedUnstructuredGridReader(FileName = fn)
elif filenames[-4::] == '.pvd':
# create a new 'PVD Reader'
    fn = glob.glob(filenames)
    print('Found files:',fn)
    fn = fn[0]
    print('Read the first one:',fn)
    my_source = PVDReader(FileName = path +  '/' + fn)
else:
    eprint('No pvd or pvtu files are provided')
    my_source = GetActiveSource()
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

if len(timeList) == 1 and timeList[0] == 0:
    timeList = range(0,NT)
if nt != 0:
    timeList = range(0,NT,nt)

print("NT=", NT, " timeList=", timeList)

if frameWindow[0] == 0 and frameWindow[1] == 0:
    frameWindow[0] = 0
    frameWindow[1] = NT-1
    print('frameWindow is updated:',frameWindow)

print ("renderViews, axesGrid, SpreadSheetViews, layouts.. "),
# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.InteractionMode = '2D'
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.StereoType = 'Crystal Eyes'
renderView1.BackEnd = 'OSPRay raycaster'
# get the material library
materialLibrary1 = GetMaterialLibrary()
renderView1.OSPRayMaterialLibrary = materialLibrary1
renderView1.CameraFocalDisk = 1.0
# AxesGrid property provides access to the AxesGrid object.
axesGrid = renderView1.AxesGrid
axesGrid.Visibility = 1
axesGrid.XTitle = 'X'
axesGrid.YTitle = 'Y'
axesGrid.ZTitle = 'Z'
axesGrid.XTitleFontSize = 20
axesGrid.YTitleFontSize = 20
axesGrid.ZTitleFontSize = 20
axesGrid.XLabelFontSize = 18
axesGrid.YLabelFontSize = 18
axesGrid.ZLabelFontSize = 18

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

# Create a new 'Render View'
renderView2 = CreateView('RenderView')
renderView2.ViewSize = [1168, 1168]
renderView2.InteractionMode = '2D'
renderView2.AxesGrid = 'GridAxes3DActor'
renderView2.OrientationAxesVisibility = 0
renderView2.StereoType = 'Crystal Eyes'
renderView2.BackEnd = 'OSPRay raycaster'
# copy grid settings from renderView1
renderView2.AxesGrid = axesGrid
# get the material library
materialLibrary2 = GetMaterialLibrary()
renderView2.OSPRayMaterialLibrary = materialLibrary2
renderView2.CameraFocalDisk = 1.0
SetActiveView(None)

# Create a new 'SpreadSheet View'
spreadSheetView1 = CreateView('SpreadSheetView')
spreadSheetView1.ColumnToSort = ''
spreadSheetView1.BlockSize = 1024L
# uncomment following to set a specific view size
# spreadSheetView1.ViewSize = [400, 400]

# Create a new 'SpreadSheet View'
spreadSheetView2 = CreateView('SpreadSheetView')
spreadSheetView2.ColumnToSort = ''
spreadSheetView2.BlockSize = 1024L
# uncomment following to set a specific view size
# spreadSheetView2.ViewSize = [400, 400]

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------
# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.SplitHorizontal(0, 0.500000)
layout1.AssignView(1, renderView1)
layout1.AssignView(2, renderView2)
# create new layout object 'Layout #2'
layout2 = CreateLayout(name='Layout #2')
layout2.AssignView(0, spreadSheetView1)

# create new layout object 'Layout #3'
layout3 = CreateLayout(name='Layout #3')
layout3.AssignView(0, spreadSheetView2)

print ("end")
# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)

# defining of computational domains
Show(my_source, renderView1)
Hide(my_source, renderView1)
boundsDomain = my_source.GetDataInformation().GetBounds()

lDomain = boundsDomain[1] - boundsDomain[0]
print("boundsDomain of my source=", boundsDomain, " lDomain=", lDomain)

print ("creating a cylinder.. "),
###-----------------GENERATION of Cylinder-----------------------------------
### it is timeless therefore it is outside of the loop
# create a new 'Cylinder'
cylinder1 = Cylinder()
cylinder1.Resolution = 100
cylinder1.Height = lDomain
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
for rv in [renderView1, renderView2]:
    transform1Display = Show(transform1, rv, 'UnstructuredGridRepresentation')
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
    transform1Display.PolarAxes = 'PolarAxesRepresentation'
    transform1Display.ScalarOpacityUnitDistance = 6.650076732513133

    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    transform1Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    transform1Display.ScaleTransferFunction.Points = [-2.220446049250313e-16, 0.0, 0.5, 0.0, 2.220446049250313e-16, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    transform1Display.OpacityTransferFunction.Points = [-2.220446049250313e-16, 0.0, 0.5, 0.0, 2.220446049250313e-16, 1.0, 0.5, 0.0]

    Show(transform1, rv)

# create a new 'Slice'
slice2 = Slice(Input=transform1)
slice2.SliceType = 'Plane'
slice2.HyperTreeGridSlicer = 'Plane'
slice2.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice2.SliceType.Origin = [0.5*lDomain, 0.0, -1e-5]
slice2.SliceType.Normal = [0.0, 0.0, 1.0]

# init the 'Plane' selected for 'HyperTreeGridSlicer'
slice2.HyperTreeGridSlicer.Origin = [0.5*lDomain, 0.0, -0.21650634706020355]

# show data from slice2
slice2Display = Show(slice2, renderView1, 'GeometryRepresentation')
Show(slice2)

# trace defaults for the display properties.
slice2Display.Representation = 'Surface'
slice2Display.AmbientColor = [0.0, 0.0, 0.0]
slice2Display.ColorArrayName = [None, '']
slice2Display.DiffuseColor = [0.0, 0.0, 0.0]
slice2Display.LineWidth = 4.0
slice2Display.OSPRayScaleArray = 'Normals'
slice2Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice2Display.SelectOrientationVectors = 'None'
slice2Display.ScaleFactor = 3.0
slice2Display.SelectScaleArray = 'None'
slice2Display.GlyphType = 'Arrow'
slice2Display.GlyphTableIndexArray = 'None'
slice2Display.GaussianRadius = 0.15
slice2Display.SetScaleArray = ['POINTS', 'Normals']
slice2Display.ScaleTransferFunction = 'PiecewiseFunction'
slice2Display.OpacityArray = ['POINTS', 'Normals']
slice2Display.OpacityTransferFunction = 'PiecewiseFunction'
slice2Display.DataAxesGrid = 'GridAxesRepresentation'
slice2Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
slice2Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
slice2Display.ScaleTransferFunction.Points = [-1.6653345369377348e-16, 0.0, 0.5, 0.0, 1.6653345369377348e-16, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
slice2Display.OpacityTransferFunction.Points = [-1.6653345369377348e-16, 0.0, 0.5, 0.0, 1.6653345369377348e-16, 1.0, 0.5, 0.0]
print ("end")
###-----------------GENERATION of a CLIP and Convert to Point data-----------------------------------
global my_stats # output

if not noPic:
    for i in timeList:
        print("in loop iteration:" + str(i))
        renderView1.ViewSize = [3108, 1168]
        SetActiveView(renderView1)
        # Properties modified on animationScene1
        animationScene1.AnimationTime = timesteps[i]
        # Properties modified on timeKeeper1
        timeKeeper1.Time = timesteps[i]
        if i==0:
            fn = str(timesteps[i]) + 'stats.txt'
            my_stats = open(fn, 'w')
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

        # update the view to ensure updated data information
        spreadSheetView1.Update()

        ss_data = Fetch(passArrays1)
        print('N=', ss_data.GetNumberOfPoints())
#         x_sum = ss_data.GetPointData().GetArray('Point').GetValue(0)
        x_mean = ss_data.GetPoint(0)[0]
#         u_mean = ss_data.GetPointData().GetArray('u.x').GetValue(0)
#         print('u_mean:',u_mean)
        s = "time: {} i: {} coarse_x_mean: {}  ".format(timesteps[i], i, x_mean)
        print (s)
        my_stats.write(s)

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
        len_min = max(bounds[0] - 2, 0.5)
        len_max = min(bounds[1] + 2, 29.5)

        length = len_max - len_min
        print('COARSE ESTIMATION: len_min=',len_min,' len_max=',len_max,' len=',length, ' len_bub=', len_bub)
        if abs(length) > 1e30:
            print ("ERROR in reading: out of box... go to the next time step...")
            continue


# ***************** CLIP A BOX of a BOX ****************************
        # create a new 'Clip'
        clip1 = Clip(Input=clip0)
        clip1.ClipType = 'Box'
        clip1.HyperTreeGridClipper = 'Plane'
        clip1.Scalars = ['POINTS', 'f']
        clip1.Value = 0.5
        # init the 'Box' selected for 'ClipType'
        clip1.ClipType.Position = [len_min - 0.1, -0.6, -0.6]
        clip1.ClipType.Length = [length + 0.2, 1.2, 1.2]

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
        Nx = int(length/(lDomain/2.0**maxlevel))
        Nyz = int(1.0/(lDomain/2.0**maxlevel))
        resampleDimension = [Nx, Nyz, Nyz]
        print("resampleDimension:", resampleDimension)
        resampleToImage1.SamplingDimensions = resampleDimension
        resampleToImage1.SamplingBounds = [len_min, len_max, -0.5, 0.5, -0.5, 0.5]

# ***************** REFINED CONTOUR1 for VOLUME FRACTION f ****************************
        # create a new 'Iso Volume'
        isoVolume1 = IsoVolume(Input=resampleToImage1)
        isoVolume1.InputScalars = ['POINTS', 'f']
        isoVolume1.ThresholdRange = [0.0, 0.5]

        # create a new 'Connectivity'
        connectivity1 = Connectivity(Input=isoVolume1)
        connectivity1.ExtractionMode = 'Extract Largest Region'
        connectivity1.ColorRegions = 0
        Show(connectivity1, renderView1)
        Hide(connectivity1, renderView1)

# ***************** REFINED BOUNDS bounds, center, len_min, len_max ****************************
        bounds = connectivity1.GetDataInformation().GetBounds()
        print("bounds_refined=",bounds)
        center = [(bounds[0] + bounds[1])/2, (bounds[2] + bounds[3])/2,(bounds[4] + bounds[5])/2]
        print('center_refined=',center)

        len_bub = bounds[1] - bounds[0]
        len_min = max(bounds[0] - 2, 0.5)
        len_max = min(bounds[1] + 2, 29.5)
        length = len_max - len_min

        # create a new 'Integrate Variables'
        integrateVariables1 = IntegrateVariables(Input=connectivity1)

        # show data from integrateVariables1
        integrateVariables1Display = Show(integrateVariables1, spreadSheetView1, 'SpreadSheetRepresentation')

        # create a new 'Pass Arrays'
        passArrays1 = PassArrays(Input=integrateVariables1)
        passArrays1.PointDataArrays = [ 'Points', 'u.x']

        # update the view to ensure updated data information
        spreadSheetView1.Update()

        ss_data = Fetch(passArrays1)
        print('N=', ss_data.GetNumberOfPoints())
        x_mean = ss_data.GetPoint(0)[0]
        u_mean = ss_data.GetPointData().GetArray('u.x').GetValue(0)
        s = "refined_x_mean: {} refined_u_mean: {} len_min: {} len_max: {} len: {} len_bub: {} \n".format(x_mean, u_mean, len_min, len_max, length, len_bub)
        print (s)
        my_stats.write(s)

# ***************** CONTOUR2 for LAMBDA2 l2 ****************************
        # create a new 'Contour'
        contour2 = Contour(Input=resampleToImage1)
        contour2.ContourBy = ['POINTS', 'l2']
        contour2.Isosurfaces = [-2.0, -1.0]
        contour2.PointMergeMethod = 'Uniform Binning'

        # create a new 'Slice'
        slice1 = Slice(Input=connectivity1)
        slice1.SliceType = 'Plane'
        slice1.HyperTreeGridSlicer = 'Plane'
        slice1.SliceOffsetValues = [0]

        # init the 'Plane' selected for 'SliceType'
        slice1.SliceType.Origin = [0.5*(len_min + len_max), 0.0, 0.0]
        slice1.SliceType.Normal = [0.0, 0.0, 1.0]
# ***************** LAMBDA2 inside a bubble: ISOVOLUME1 for VOLUME FRACTION f ****************************
#         # create a new 'Iso Volume'
#         isoVolume1 = IsoVolume(Input=connectivity1)
#         isoVolume1.InputScalars = ['POINTS', 'f']
#         isoVolume1.ThresholdRange = [0.0, 0.5]
#
#         # Find unique region IDs
#         regionArray = isoVolume1.CellData['RegionId']
#         regions = unique(regionArray)
#         v = volume(inp)
#         a=[]
#         for regionValue in regions:
#           regionVolume = sum(v[regionArray == regionValue])
#           a.append((regionValue,regionVolume))
#         #  print(regionValue, regionVolume)
#         xx = sorted(a,reverse=True,key=lambda x: (x[1]))
#         print(xx[0:10])
#
#         # create a new 'Integrate Variables'
#         integrateVariables1 = IntegrateVariables(Input=isoVolume1)
#
#         # show data from integrateVariables1
#         integrateVariables1Display = Show(integrateVariables1, spreadSheetView1, 'SpreadSheetRepresentation')
#
#         # create a new 'Pass Arrays'
#         passArrays1 = PassArrays(Input=integrateVariables1)
#         passArrays1.PointDataArrays = [ 'Points', 'u.x']
#
#         # update the view to ensure updated data information
#         spreadSheetView1.Update()
#
#         ss_data = Fetch(passArrays1)
#         print('N=', ss_data.GetNumberOfPoints())
#         x_mean = ss_data.GetPoint(0)[0]
#         u_mean = ss_data.GetPointData().GetArray('u.x').GetValue(0)
#         s = "refined_x_mean: {}  u_mean: {} \n".format(x_mean, u_mean)
#         print (s)
#         my_stats.write(s)

        # create a new 'Contour'
        contour5 = Contour(Input=connectivity1)
        contour5.ContourBy = ['POINTS', 'l2']
        contour5.Isosurfaces = [-1.5]
        contour5.PointMergeMethod = 'Uniform Binning'

# ----------------------------------------------------------------
        # setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

        # get color transfer function/color map for 'ux'
        uxLUT = GetColorTransferFunction('ux')
        uxLUT.RGBPoints = [0.3730543491279064, 0.0, 0.0, 0.5625, 0.5106463634876334, 0.0, 0.0, 1.0, 0.8251430154745502, 0.0, 1.0, 1.0, 0.9823910318856668, 0.5, 1.0, 0.5, 1.139639048296783, 1.0, 1.0, 0.0, 1.4541357002836997, 1.0, 0.0, 0.0, 1.6113837166948162, 0.5, 0.0, 0.0]
        uxLUT.ColorSpace = 'RGB'
        uxLUT.ScalarRangeInitialized = 1.0
#         uxLUT.ApplyPreset('jet', True)
        len_bar = 0.5
        # get color legend/bar for uxLUT in view renderView1
        uxLUTColorBar = GetScalarBar(uxLUT, renderView1)
        uxLUTColorBar.Orientation = 'Horizontal'
        uxLUTColorBar.WindowLocation = 'AnyLocation'
        uxLUTColorBar.ScalarBarLength = len_bar
        uxLUTColorBar.Position = [0.7 - 0.5*len_bar, 0.02]
        uxLUTColorBar.Title = '|u|'
        uxLUTColorBar.ComponentTitle = ''
        uxLUTColorBar.TitleColor = [1, 1, 1]
        uxLUTColorBar.LabelColor = [1, 1, 1]
        uxLUTColorBar.LabelFormat = '%-#6.2g'
        uxLUTColorBar.RangeLabelFormat = '%6.2g'
        uxLUTColorBar.ScalarBarThickness = 16*2
        uxLUTColorBar.TitleFontSize = 16*2
        uxLUTColorBar.LabelFontSize = 16*2


        # set color bar visibility
        uxLUTColorBar.Visibility = 1


        len_bar = 0.5
        uxLUTColorBar2 = GetScalarBar(uxLUT, renderView2)
        uxLUTColorBar2.Orientation = 'Horizontal'
        uxLUTColorBar2.WindowLocation = 'AnyLocation'
        uxLUTColorBar2.Position = [0.5 - 0.5*len_bar, 0.02]
        uxLUTColorBar2.Title = '|u|'
        uxLUTColorBar2.ComponentTitle = ''
        uxLUTColorBar2.TitleColor = [1, 1, 1]
        uxLUTColorBar2.LabelColor = [1, 1, 1]
        uxLUTColorBar2.LabelFormat = '%-#6.2g'
        uxLUTColorBar2.RangeLabelFormat = '%6.2g'
        uxLUTColorBar2.ScalarBarLength = len_bar

        # set color bar visibility
        uxLUTColorBar2.Visibility = 0



        # get color transfer function/color map for 'l2'
        l2LUT = GetColorTransferFunction('l2')
        l2LUT.RGBPoints = [-2.0, 0.0, 1.0, 1.0, -1.55, 0.0, 0.0, 1.0, -1.5, 0.0, 0.0, 0.501960784314, -1.4499999999999997, 1.0, 0.0, 0.0, -1.0, 1.0, 1.0, 0.0]
        l2LUT.ColorSpace = 'RGB'
        l2LUT.ScalarRangeInitialized = 1.0

        for rv in [renderView1, renderView2]:
            # show data from connectivity1
            print ("showing connectivity1Display.. "),
            connectivity1Display = Show(connectivity1, rv, 'GeometryRepresentation')
            # trace defaults for the display properties.
            connectivity1Display.Representation = 'Surface'
            connectivity1Display.ColorArrayName = ['POINTS', 'u.x']
            connectivity1Display.LookupTable = uxLUT
            connectivity1Display.Specular = 1.0
            connectivity1Display.Luminosity = 49.0
            connectivity1Display.Ambient = 0.13
            connectivity1Display.OSPRayScaleArray = 'f'
            connectivity1Display.OSPRayScaleFunction = 'PiecewiseFunction'
            connectivity1Display.SelectOrientationVectors = 'None'
            connectivity1Display.ScaleFactor = 0.3379646740552927
            connectivity1Display.SelectScaleArray = 'f'
            connectivity1Display.GlyphType = 'Arrow'
            connectivity1Display.GlyphTableIndexArray = 'f'
            connectivity1Display.GaussianRadius = 0.016898233702764633
            connectivity1Display.SetScaleArray = ['POINTS', 'f']
            connectivity1Display.ScaleTransferFunction = 'PiecewiseFunction'
            connectivity1Display.OpacityArray = ['POINTS', 'f']
            connectivity1Display.OpacityTransferFunction = 'PiecewiseFunction'
            connectivity1Display.DataAxesGrid = 'GridAxesRepresentation'
            connectivity1Display.PolarAxes = 'PolarAxesRepresentation'
            print ("end")
        # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
#         connectivity1Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

        # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
#         connectivity1Display.ScaleTransferFunction.Points = [0.5, 0.0, 0.5, 0.0, 0.5001220703125, 1.0, 0.5, 0.0]

        # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
#         connectivity1Display.OpacityTransferFunction.Points = [0.5, 0.0, 0.5, 0.0, 0.5001220703125, 1.0, 0.5, 0.0]

        # setup the color legend parameters for each legend in this view

        # show color legend
        connectivity1Display.SetScalarBarVisibility(renderView1, True)

        # ----------------------------------------------------------------
        # setup color maps and opacity mapes used in the visualization
        # note: the Get..() functions create a new object, if needed
        # ----------------------------------------------------------------
        passArrays0 = PassArrays(Input=connectivity1)
        passArrays0.PointDataArrays = [ 'Points']

        ss_data = Fetch(passArrays0)
        numPoints = ss_data.GetNumberOfPoints()
        print ("Number of points on contour:", numPoints)
        data=np.zeros((numPoints,3))
        for x in range(numPoints):
            data[x] = ss_data.GetPoint(x)
#             data[i] = connectivity1.GetPointData().GetArray('f').GetValue(x)
#         rad = np.linalg.norm(data[],axis=1)
        rad = np.sqrt(data[:,1]**2 + data[:,2]**2)
        i_max = np.argmax(rad)
        x_maxy = data[i_max,0]
#         data = data[(data[:,0]>=x_maxy) & (data[:,0]<=x_mean)]
        if x_maxy < x_mean:
            print('x_maxy < x_mean')
            rad = rad[(data[:,0]>=x_maxy) & (data[:,0]<=x_mean)]
        else:
            print('x_maxy >= x_mean')
            rad = rad[(data[:,0]>=x_mean) & (data[:,0]<=x_maxy)]
        print('     x_maxy:',x_maxy, " x_mean:", x_mean)
        rad_mean = rad.mean()
        s = "RAD_mean: " + str(round(rad_mean,6)) + " delta: " + str(round(0.5 - rad_mean,6)) + "\n"
        print (s)
        my_stats.write(s)

        print("resampleToImage1.GetDataInformation.."),
        info = resampleToImage1.GetDataInformation().DataInformation
        arrayInfo = info.GetArrayInformation("u.x", vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS)
        print('end')
        range0 = arrayInfo.GetComponentRange(0)
        range1 = arrayInfo.GetComponentRange(1)
        range2 = arrayInfo.GetComponentRange(2)
        rgX, rgY, rgZ = max(abs(range0[0]),abs(range0[1])), max(abs(range1[0]),abs(range1[1])), max(abs(range2[0]),abs(range2[1]))
        range_max = np.sqrt(rgX**2 + rgY**2 + rgZ**2)
        print("velocity range_max= ", range_max)
        print("velocity range_xyz= ", range0, range1, range2)
        uxLUT.RescaleTransferFunction(0, range_max)
        op = GetOpacityTransferFunction("ux")
        op.RescaleTransferFunction(0, range_max)
        uxLUTColorBar.UseCustomLabels = 1
        # labels from 100 to 200 in increments of 10
        uxLUTColorBar.CustomLabels = np.around(np.linspace(0, range_max, 4),1)
#         uxLUTColorBar.CustomLabels = np.arange(0, range_max, 0.5)
        print("CustomLabels=",uxLUTColorBar.CustomLabels )

        uxLUTColorBar2.UseCustomLabels = 1
        # labels from 100 to 200 in increments of 10
        uxLUTColorBar2.CustomLabels = np.around(np.linspace(0, range_max, 4),1)
#         uxLUTColorBar2.CustomLabels = np.arange(0, range_max, 0.5)
        print("CustomLabels=",uxLUTColorBar2.CustomLabels )

        # ----------------------------------------------------------------
        # ----------------------------------------------------------------

#         renderView1.CameraPosition = [center[0], center[1], 6]
#         renderView1.CameraFocalPoint = center
        renderView1.CameraViewUp = [0.2, 1, 0]
        renderView1.CameraParallelScale = 1
        renderView1.CenterOfRotation = [0.5*(len_min + len_max), 0.0, 0.0]
        renderView1.CameraFocalPoint = [0.5*(bounds[0]+center[0]), 0, 0]
        renderView1.CameraPosition = [0.5*(len_min + len_max) - 4, 0.6, 4.5]

        # update the view to ensure updated data information
        renderView1.Update()

#****************** CONNECTIVITY(f) AND U MAGNITUDE ********************
        # show data from connectivity1
        fn = path + "/" + picName +  '_t=' + str(timesteps[i]) +'ux.png'
        SaveScreenshot( fn, renderView1,
        #      ImageResolution=[2316, 2204],
            TransparentBackground=0,
            CompressionLevel='2' )
        print('File=' + fn + ' generated succesfully')


#****************** CONNECTIVITY(f) AND CONTOUR2 (lambda2) ********************
        # set color bar visibility
        uxLUTColorBar.Visibility = 0

        # show data from connectivity1
        print ("Showing connectivity1.. "),
        connectivity1Display = Show(connectivity1, renderView1, 'GeometryRepresentation')

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
        contour2Display = Show(contour2, renderView1, 'GeometryRepresentation')
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
        print ("end")
        print("RenderView update..")
        renderView1.Update()
        print ("end")
        fn = path + "/" + picName+  '_t=' + str(timesteps[i]) +'_noLambda2.png'
        SaveScreenshot( fn, renderView1,
                    #      ImageResolution=[2316, 2204],
                        TransparentBackground=0,
                        CompressionLevel='2' )
        print('File=' + fn + ' generated succesfully')


#****************** LAMBDA2 inside a bubble  ********************

        print ("Showing connectivity1 and hiding contour2.. "),
        Show(connectivity1, renderView1) # bubble shape with opacity
        Hide(contour2, renderView1) # hide lambda2 all
        print ("end")
        # trace defaults for the display properties.
        connectivity1Display.Representation = 'Surface'
        connectivity1Display.Opacity = 0.1

        # show data from contour5
        print ("Showing contour5.. "),
        contour5Display = Show(contour5, renderView1, 'GeometryRepresentation')
        # get color transfer function/color map for 'l2'
        l2LUT = GetColorTransferFunction('l2')
        l2LUT.RGBPoints = [-1.0, 0.054901960784313725, 0.9411764705882353, 0.12941176470588237, -0.75, 0.865, 0.865, 0.865, -0.5, 1.0, 1.0, 0.0]
        l2LUT.ScalarRangeInitialized = 1.0
        # trace defaults for the display properties.
        contour5Display.Representation = 'Surface'
        contour5Display.ColorArrayName = ['POINTS', 'l2']
        contour5Display.LookupTable = l2LUT
        contour5Display.Opacity = 0.62
        contour5Display.Specular = 1.0
        contour5Display.OSPRayScaleArray = 'l2'
        contour5Display.OSPRayScaleFunction = 'PiecewiseFunction'
        contour5Display.SelectOrientationVectors = 'None'
        contour5Display.ScaleFactor = 0.2388025760650635
        contour5Display.SelectScaleArray = 'l2'
        contour5Display.GlyphType = 'Arrow'
        contour5Display.GlyphTableIndexArray = 'l2'
        contour5Display.GaussianRadius = 0.011940128803253174
        contour5Display.SetScaleArray = ['POINTS', 'l2']
        contour5Display.ScaleTransferFunction = 'PiecewiseFunction'
        contour5Display.OpacityArray = ['POINTS', 'l2']
        contour5Display.OpacityTransferFunction = 'PiecewiseFunction'
        contour5Display.DataAxesGrid = 'GridAxesRepresentation'
        contour5Display.PolarAxes = 'PolarAxesRepresentation'

        # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
        contour5Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

        # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
        contour5Display.ScaleTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, -0.9998779296875, 1.0, 0.5, 0.0]

        # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
        contour5Display.OpacityTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, -0.9998779296875, 1.0, 0.5, 0.0]

        # show color legend
        contour5Display.SetScalarBarVisibility(renderView1, False)
        print ("end")
        l2PWF = GetOpacityTransferFunction('l2')
        l2PWF.Points = [-1.0, 0.0, 0.5, 0.0, -0.5, 1.0, 0.5, 0.0]
        l2PWF.ScalarRangeInitialized = 1

        renderView1.CameraViewUp = [0.2, 1, 0]
        renderView1.CameraParallelScale = 0.8
        renderView1.CenterOfRotation = [center[0], 0.0, 0.0]
        renderView1.CameraFocalPoint = [center[0], 0, 0]
#         renderView1.CameraFocalPoint = [0.5*(len_min + len_max) , -0.25, -1.3]
        renderView1.CameraPosition = [center[0] - 2, 0.6, 4.5]
        print("RenderView update..")
        renderView1.Update()
        print ("end")
        fn = path + "/" + picName+  '_t=' + str(timesteps[i]) +'_Lambda2_in_bubble.png'
        SaveScreenshot( fn, renderView1,
                    #      ImageResolution=[2316, 2204],
                        TransparentBackground=0,
                        CompressionLevel='2' )
        print('File=' + fn + ' generated succesfully')
#****************** CONNECTIVITY1(f) AND TRACERS ********************
        print('hiding connectivity1, contour2, contour5..')
        Hide(connectivity1, renderView1)
        Hide(contour2, renderView1)
        Hide(contour5, renderView1)
        print ("end")

        selection=SelectPoints()
        selection.QueryString='(pointIsNear([(' + str(lDomain) + ',0,0),], ' + str(lDomain) + ', inputs))'
        selection.FieldType = 'POINT'
        selection.UpdatePipelineInformation()

        extractSelection1 = ExtractSelection(Input=connectivity1, Selection=selection)
        extractSelection1.UpdatePipeline()



        # create a query selection
#         sel = QuerySelect(InsideOut=0, QueryString='(pointIsNear([(' + str(lDomain) + ',0,0),], ' + str(lDomain) + ', inputs))', FieldType='POINT')
#         # create a new 'Extract Selection'
#         extractSelection1 = ExtractSelection(Input=connectivity1, Selection=sel)
#         extractSelection1Display = Show(extractSelection1, renderView1, 'UnstructuredGridRepresentation')
#         Show(extractSelection1, spreadSheetView1, 'SpreadSheetRepresentation')


#         Show(extractSelection1, spreadSheetView1, 'SpreadSheetRepresentation')
#         Hide(extractSelection1, spreadSheetView1)
#         u_tip = 1.29#delete me
#         info = resampleToImage1.GetDataInformation().DataInformation
#         arrayInfo = info.GetArrayInformation("u.x", vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS)
#         print('end')
#         range0 = arrayInfo.GetComponentRange(0)
#         range1 = arrayInfo.GetComponentRange(1)
#         range2 = arrayInfo.GetComponentRange(2)

        info = extractSelection1.GetDataInformation().DataInformation
        print('info Number of Points=',info.GetNumberOfPoints ())
        arrayInfo = info.GetArrayInformation("u.x", vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS)
        u_tip = arrayInfo.GetComponentRange(0)[0]
        print("u_tip = ", u_tip)


        # create a new 'Calculator'
        calculator1 = Calculator(Input=resampleToImage1)
        calculator1.ResultArrayName = 'deltaU'
        calculator1.Function = 'u.x - iHat*' + str(u_tip)
        print('calculator: ', calculator1.Function)

        deltaULUT = GetColorTransferFunction('deltaU')
        deltaULUT.RGBPoints = [0.3730543491279064, 0.0, 0.0, 0.5625, 0.5106463634876334, 0.0, 0.0, 1.0, 0.8251430154745502, 0.0, 1.0, 1.0, 0.9823910318856668, 0.5, 1.0, 0.5, 1.139639048296783, 1.0, 1.0, 0.0, 1.4541357002836997, 1.0, 0.0, 0.0, 1.6113837166948162, 0.5, 0.0, 0.0]
        deltaULUT.ColorSpace = 'RGB'
        deltaULUT.ScalarRangeInitialized = 1.0

#         deltaULUT.ApplyPreset('jet', True)
        len_bar = 0.5
        # get color legend/bar for deltaULUT in view renderView1
        deltaULUTColorBar = GetScalarBar(deltaULUT, renderView1)
        deltaULUTColorBar.Orientation = 'Horizontal'
        deltaULUTColorBar.WindowLocation = 'AnyLocation'
        deltaULUTColorBar.ScalarBarLength = len_bar
        deltaULUTColorBar.Position = [0.5 - 0.5*len_bar, 0.02]
        deltaULUTColorBar.Title = '|u - u_tip|'
        deltaULUTColorBar.ComponentTitle = ''
        deltaULUTColorBar.TitleColor = [1, 1, 1]
        deltaULUTColorBar.LabelColor = [1, 1, 1]
        deltaULUTColorBar.LabelFormat = '%-#6.1g'
        deltaULUTColorBar.RangeLabelFormat = '%6.1g'
        deltaULUTColorBar.ScalarBarThickness = 16*2
        deltaULUTColorBar.TitleFontSize = 16*2
        deltaULUTColorBar.LabelFontSize = 16*2
        # set color bar visibility
        deltaULUTColorBar.Visibility = 1


        connectivity1Display.Opacity = 0.23
        # create a new 'Stream Tracer'
        streamTracer1 = StreamTracer(Input=calculator1, SeedType='High Resolution Line Source')
        streamTracer1.Vectors = ['POINTS', 'deltaU']
        streamTracer1.SurfaceStreamlines = 1
        streamTracer1.InitialStepLength = 0.05
        streamTracer1.MaximumStreamlineLength = length
        streamTracer1.MaximumSteps = 20000
        streamTracer1.TerminalSpeed = 0.0001

        # init the 'High Resolution Line Source' selected for 'SeedType'
        streamTracer1.SeedType.Point1 = [bounds[1] + 1, -0.5, 0.0]
        streamTracer1.SeedType.Point2 = [bounds[1] + 1, 0.5, 0.0]
        streamTracer1.SeedType.Resolution = 50

        # show data from streamTracer1
        print ("Showing streamTracer1Display.. "),
        streamTracer1Display = Show(streamTracer1, renderView1, 'GeometryRepresentation')

        # trace defaults for the display properties.
        streamTracer1Display.Representation = 'Surface'
        streamTracer1Display.ColorArrayName = ['POINTS', 'deltaU']
        streamTracer1Display.LookupTable = deltaULUT
        streamTracer1Display.OSPRayScaleArray = 'AngularVelocity'
        streamTracer1Display.OSPRayScaleFunction = 'PiecewiseFunction'
        streamTracer1Display.SelectOrientationVectors = 'Normals'
        streamTracer1Display.ScaleFactor = 0.3379646740552927
        streamTracer1Display.SelectScaleArray = 'AngularVelocity'
        streamTracer1Display.GlyphType = 'Arrow'
        streamTracer1Display.GlyphTableIndexArray = 'AngularVelocity'
        streamTracer1Display.GaussianRadius = 0.03949999809265137
        streamTracer1Display.SetScaleArray = ['POINTS', 'AngularVelocity']
        streamTracer1Display.ScaleTransferFunction = 'PiecewiseFunction'
        streamTracer1Display.OpacityArray = ['POINTS', 'AngularVelocity']
        streamTracer1Display.OpacityTransferFunction = 'PiecewiseFunction'
        streamTracer1Display.DataAxesGrid = 'GridAxesRepresentation'
        streamTracer1Display.PolarAxes = 'PolarAxesRepresentation'

        streamTracer2 = StreamTracer(Input=calculator1, SeedType='High Resolution Line Source')
        streamTracer2.Vectors = ['POINTS', 'deltaU']
        streamTracer2.SurfaceStreamlines = 1
        streamTracer2.InitialStepLength = 0.01
        streamTracer2.MaximumStreamlineLength = length
        streamTracer2.MaximumSteps = 20000
        streamTracer2.TerminalSpeed = 0.0001

        # init the 'High Resolution Line Source' selected for 'SeedType'
        streamTracer2.SeedType.Point1 = [bounds[0] - 1, -0.5, 0.0]
        streamTracer2.SeedType.Point2 = [bounds[0] - 1, 0.5, 0.0]
        streamTracer2.SeedType.Resolution = 50
        # show data from streamTracer2
        print ("Showing streamTracer2Display.. "),
        streamTracer2Display = Show(streamTracer2, renderView1, 'GeometryRepresentation')

        # trace defaults for the display properties.
        streamTracer2Display.Representation = 'Surface'
        streamTracer2Display.ColorArrayName = ['POINTS', 'deltaU']
        streamTracer2Display.LookupTable = deltaULUT
        streamTracer2Display.OSPRayScaleArray = 'AngularVelocity'
        streamTracer2Display.OSPRayScaleFunction = 'PiecewiseFunction'
        streamTracer2Display.SelectOrientationVectors = 'Normals'
        streamTracer2Display.ScaleFactor = 0.3379646740552927
        streamTracer2Display.SelectScaleArray = 'AngularVelocity'
        streamTracer2Display.GlyphType = 'Arrow'
        streamTracer2Display.GlyphTableIndexArray = 'AngularVelocity'
        streamTracer2Display.GaussianRadius = 0.03949999809265137
        streamTracer2Display.SetScaleArray = ['POINTS', 'AngularVelocity']
        streamTracer2Display.ScaleTransferFunction = 'PiecewiseFunction'
        streamTracer2Display.OpacityArray = ['POINTS', 'AngularVelocity']
        streamTracer2Display.OpacityTransferFunction = 'PiecewiseFunction'
        streamTracer2Display.DataAxesGrid = 'GridAxesRepresentation'
        streamTracer2Display.PolarAxes = 'PolarAxesRepresentation'

        streamTracer3 = StreamTracer(Input=calculator1, SeedType='High Resolution Line Source')
        streamTracer3.Vectors = ['POINTS', 'deltaU']
        streamTracer3.SurfaceStreamlines = 1
        streamTracer3.InitialStepLength = 0.01
        streamTracer3.MaximumStreamlineLength = length
        streamTracer3.MaximumSteps = 20000
        streamTracer3.TerminalSpeed = 0.0001

        # init the 'High Resolution Line Source' selected for 'SeedType'
        streamTracer3.SeedType.Point1 = [center[0], -0.4, 0.0]
        streamTracer3.SeedType.Point2 = [center[0], 0.4, 0.0]
        streamTracer3.SeedType.Resolution = 50

        # show data from streamTracer3
        print ("Showing streamTracer3Display.. "),
        streamTracer3Display = Show(streamTracer3, renderView1, 'GeometryRepresentation')

        # trace defaults for the display properties.
        streamTracer3Display.Representation = 'Surface'
        streamTracer3Display.ColorArrayName = ['POINTS', 'deltaU']
        streamTracer3Display.LookupTable = deltaULUT
        streamTracer3Display.OSPRayScaleArray = 'AngularVelocity'
        streamTracer3Display.OSPRayScaleFunction = 'PiecewiseFunction'
        streamTracer3Display.SelectOrientationVectors = 'Normals'
        streamTracer3Display.ScaleFactor = 0.3379646740552927
        streamTracer3Display.SelectScaleArray = 'AngularVelocity'
        streamTracer3Display.GlyphType = 'Arrow'
        streamTracer3Display.GlyphTableIndexArray = 'AngularVelocity'
        streamTracer3Display.GaussianRadius = 0.03949999809265137
        streamTracer3Display.SetScaleArray = ['POINTS', 'AngularVelocity']
        streamTracer3Display.ScaleTransferFunction = 'PiecewiseFunction'
        streamTracer3Display.OpacityArray = ['POINTS', 'AngularVelocity']
        streamTracer3Display.OpacityTransferFunction = 'PiecewiseFunction'
        streamTracer3Display.DataAxesGrid = 'GridAxesRepresentation'
        streamTracer3Display.PolarAxes = 'PolarAxesRepresentation'

        print ("Showing slice1Display.. "),
        slice1Display = Show(slice1, renderView1, 'GeometryRepresentation')

        # trace defaults for the display properties.
        slice1Display.Representation = 'Surface'
        slice1Display.AmbientColor = [0.0, 0.0, 0.0]
        slice1Display.ColorArrayName = ['POINTS', '']
        slice1Display.DiffuseColor = [0.0, 0.0, 0.0]
        slice1Display.LineWidth = 4.0
        slice1Display.OSPRayScaleArray = 'f'
        slice1Display.OSPRayScaleFunction = 'PiecewiseFunction'
        slice1Display.SelectOrientationVectors = 'None'
        slice1Display.ScaleFactor = 0.3379646740552927
        slice1Display.SelectScaleArray = 'f'
        slice1Display.GlyphType = 'Arrow'
        slice1Display.GlyphTableIndexArray = 'f'
        slice1Display.GaussianRadius = 0.003989625424146652
        slice1Display.SetScaleArray = ['POINTS', 'f']
        slice1Display.ScaleTransferFunction = 'PiecewiseFunction'
        slice1Display.OpacityArray = ['POINTS', 'f']
        slice1Display.OpacityTransferFunction = 'PiecewiseFunction'
        slice1Display.DataAxesGrid = 'GridAxesRepresentation'
        slice1Display.PolarAxes = 'PolarAxesRepresentation'

        print ("calculator1.GetDataInformation.. ")
        info = calculator1.GetDataInformation().DataInformation
        arrayInfo = info.GetArrayInformation("deltaU", vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS)
        range0 = arrayInfo.GetComponentRange(0)
        range1 = arrayInfo.GetComponentRange(1)
        range2 = arrayInfo.GetComponentRange(2)
        print ("end")
        rgX, rgY, rgZ = max(abs(range0[0]),abs(range0[1])), max(abs(range1[0]),abs(range1[1])), max(abs(range2[0]),abs(range2[1]))
        range_max = np.sqrt(rgX**2 + rgY**2 + rgZ**2)
        print("deltaU range_max= ", range_max)
        print("deltaU range_xyz= ", range0, range1, range2)
        deltaULUT.RescaleTransferFunction(0, range_max)
        op = GetOpacityTransferFunction("deltaU")
        op.RescaleTransferFunction(0, range_max)
        deltaULUTColorBar.UseCustomLabels = 1
        # labels from 100 to 200 in increments of 10
        deltaULUTColorBar.CustomLabels = np.around(np.linspace(0, range_max, 4),1)
#         uxLUTColorBar.CustomLabels = np.arange(0, range_max, 0.5)
        print("CustomLabels deltaU=",uxLUTColorBar.CustomLabels )

        renderView1.CameraPosition = [center[0], center[1], 6]
        renderView1.CameraFocalPoint = center
        renderView1.CameraParallelScale = 1
        renderView1.CameraViewUp = [0, 1, 0]

        print("RenderView update..")
        renderView1.Update()
        print ("end")
        fn = path + "/" + picName +  '_t=' + str(timesteps[i]) +'_tracer.png'
        SaveScreenshot( fn, renderView1,
        #      ImageResolution=[2316, 2204],
            TransparentBackground=0,
            CompressionLevel='2' )
        Hide(streamTracer1, renderView1)
        Hide(streamTracer2, renderView1)
        Hide(streamTracer3, renderView1)
        # set color bar visibility
        uxLUTColorBar.Visibility = 0
        deltaULUTColorBar.Visibility = 0

        print('File=' + fn + ' generated succesfully')

#***************** Like in the article SIDE SLICE and contour(f)*********************
        # create a new 'Slice'
        slice3 = Slice(Input=resampleToImage1)
        slice3.SliceType = 'Plane'
        slice3.HyperTreeGridSlicer = 'Plane'
        slice3.SliceOffsetValues = [0]

        # init the 'Plane' selected for 'SliceType'
        slice3.SliceType.Origin = [0.5*(len_min + len_max), 0.0, 0.0]
        slice3.SliceType.Normal = [0.0, 0.0, 1.0]
        # init the 'Plane' selected for 'HyperTreeGridSlicer'
        # slice3.HyperTreeGridSlicer.Origin = [10.999999761581421, 0.0, 0.0] #???
        print ("Showing slice3Display.. "),
        slice3Display = Show(slice3, renderView1, 'GeometryRepresentation')

        # trace defaults for the display properties.
        slice3Display.Representation = 'Surface'
        slice3Display.ColorArrayName = ['POINTS', 'u.x']
        slice3Display.LookupTable = uxLUT
        slice3Display.OSPRayScaleArray = 'f'
        slice3Display.OSPRayScaleFunction = 'PiecewiseFunction'
        slice3Display.SelectOrientationVectors = 'None'
        slice3Display.ScaleFactor = 0.3379646740552927
        slice3Display.SelectScaleArray = 'f'
        slice3Display.GlyphType = 'Arrow'
        slice3Display.GlyphTableIndexArray = 'f'
        slice3Display.GaussianRadius = 0.016898233702764633
        slice3Display.SetScaleArray = ['POINTS', 'f']
        slice3Display.ScaleTransferFunction = 'PiecewiseFunction'
        slice3Display.OpacityArray = ['POINTS', 'f']
        slice3Display.OpacityTransferFunction = 'PiecewiseFunction'
        slice3Display.DataAxesGrid = 'GridAxesRepresentation'
        slice3Display.PolarAxes = 'PolarAxesRepresentation'

        # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
        slice3Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

        # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
        slice3Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0000000000000002, 1.0, 0.5, 0.0]

        # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
        slice3Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0000000000000002, 1.0, 0.5, 0.0]

        # show data from connectivity1
        print ("Showing connectivity1Display.. "),
        connectivity1Display = Show(connectivity1, renderView1, 'GeometryRepresentation')
        # trace defaults for the display properties.
        connectivity1Display.Representation = 'Surface'
        connectivity1Display.ColorArrayName = ['POINTS', 'u.x']
        connectivity1Display.LookupTable = uxLUT
        connectivity1Display.Opacity = 1
        connectivity1Display.Specular = 0.8
        connectivity1Display.SpecularPower = 100.0
        connectivity1Display.Ambient = 0.1
        connectivity1Display.OSPRayScaleArray = 'f'
        connectivity1Display.OSPRayScaleFunction = 'PiecewiseFunction'
        connectivity1Display.SelectOrientationVectors = 'None'
        connectivity1Display.ScaleFactor = 0.3379646740552927
        connectivity1Display.SelectScaleArray = 'f'
        connectivity1Display.GlyphType = 'Arrow'
        connectivity1Display.GlyphTableIndexArray = 'f'
        connectivity1Display.GaussianRadius = 0.016898233702764633
        connectivity1Display.SetScaleArray = ['POINTS', 'f']
        connectivity1Display.ScaleTransferFunction = 'PiecewiseFunction'
        connectivity1Display.OpacityArray = ['POINTS', 'f']
        connectivity1Display.OpacityTransferFunction = 'PiecewiseFunction'
        connectivity1Display.DataAxesGrid = 'GridAxesRepresentation'
        connectivity1Display.PolarAxes = 'PolarAxesRepresentation'
        uxLUTColorBar.Position = [0.5 - 0.5*len_bar, 0.02]
        uxLUTColorBar.Visibility = 1
        print("RenderView update..")
        renderView1.Update()
        print ("end")
        fn = path + "/" + picName +  '_t=' + str(timesteps[i]) +'_uxSide.png'
        SaveScreenshot( fn, renderView1,
        #      ImageResolution=[2316, 2204],
            TransparentBackground=0,
            CompressionLevel='2' )
#         Hide(streamTracer1, renderView1)
        Hide(slice3, renderView1)
        Hide(slice1, renderView1)

        print('File=' + fn + ' generated succesfully')
#**************************************************************
#**************** SLICES of U along a bubble ******************
#**************************************************************
        ss = GetSources()
        for so in ss:
          Hide(ss[so])

         # create a new 'Slice'
        slice4 = Slice(Input=resampleToImage1)
        slice4.SliceType = 'Plane'
        slice4.HyperTreeGridSlicer = 'Plane'
        slice4.SliceOffsetValues = [0]

         # create a new 'Slice' contour(f)
        slice5 = Slice(Input=connectivity1)
        slice5.SliceType = 'Plane'
        slice5.HyperTreeGridSlicer = 'Plane'
        slice5.SliceOffsetValues = [0]

        uxLUTColorBar.ScalarBarThickness = 16
        uxLUTColorBar.TitleFontSize = 16
        uxLUTColorBar.LabelFontSize = 16
        for k,sl in enumerate(np.linspace(bounds[0]+0.1, bounds[1]-0.1, nslices)):
            print("slice x=", sl)

            # create a new 'Text'
            text1 = Text()
            text1.Text = 'x=' + str(round(sl,2)) + " l/l_b=" + str(round((sl-bounds[0])/len_bub,2))

            renderView1.CameraPosition = [sl-3, 0.0, 0.0]
            renderView1.CameraFocalPoint = [sl, 0.0, 0.0]
            renderView1.CameraViewUp = [0.0, 0.0, 1.0]
            renderView1.CameraParallelScale = 14.017845768876187

            renderView1.ViewSize = [2048, 2048]
            renderView1.CameraViewUp = [0.0, 0.0, 1.0]
            renderView1.ResetCamera(sl, sl, -0.5, 0.5, -0.5, 0.5)

            renderView2.CameraViewUp = [0, 1, 0]
            renderView2.ResetCamera(bounds[0], bounds[1], -0.5, 0.5, -0.5, 0.5)
            SetActiveView(renderView2)
            camera = GetActiveCamera()
            camera.Azimuth(10)

            print("renderView2:",bounds[0], bounds[1], -0.5, 0.5, -0.5, 0.5)

            # init the 'Plane' selected for 'SliceType'
            slice4.SliceType.Origin = [sl, 0.0, 0.0]
            slice4.SliceType.Normal = [1.0, 0.0, 0.0]

            # show data from text1
            text1Display = Show(text1, renderView2, 'TextSourceRepresentation')

            # trace defaults for the display properties.
            text1Display.WindowLocation = 'AnyLocation'
            text1Display.Position = [0.3, 0.2]


            for rv in [renderView1,renderView2]:
                print ("Showing slice4Display.. "),
                slice4Display = Show(slice4, rv, 'GeometryRepresentation')

                # trace defaults for the display properties.
                slice4Display.Representation = 'Surface'
                slice4Display.ColorArrayName = ['POINTS', 'u.x']
                slice4Display.LookupTable = uxLUT
                slice4Display.OSPRayScaleArray = 'f'
                slice4Display.OSPRayScaleFunction = 'PiecewiseFunction'
                slice4Display.SelectOrientationVectors = 'None'
                slice4Display.ScaleFactor = 0.3379646740552927
                slice4Display.SelectScaleArray = 'f'
                slice4Display.GlyphType = 'Arrow'
                slice4Display.GlyphTableIndexArray = 'f'
                slice4Display.GaussianRadius = 0.016898233702764633
                slice4Display.SetScaleArray = ['POINTS', 'f']
                slice4Display.ScaleTransferFunction = 'PiecewiseFunction'
                slice4Display.OpacityArray = ['POINTS', 'f']
                slice4Display.OpacityTransferFunction = 'PiecewiseFunction'
                slice4Display.DataAxesGrid = 'GridAxesRepresentation'
                slice4Display.PolarAxes = 'PolarAxesRepresentation'

                # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
                slice4Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

                # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
                slice4Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0000000000000002, 1.0, 0.5, 0.0]

                # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
                slice4Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0000000000000002, 1.0, 0.5, 0.0]


            # init the 'Plane' selected for 'SliceType'
            slice5.SliceType.Origin = [sl, 0.0, 0.0]
            slice5.SliceType.Normal = [1.0, 0.0, 0.0]
            print ("Showing slice5Display.. "),
            slice5Display = Show(slice5, renderView1, 'GeometryRepresentation')

            # trace defaults for the display properties.
            slice5Display.Representation = 'Surface'
            slice5Display.AmbientColor = [0.0, 0.0, 0.0]
            slice5Display.ColorArrayName = ['POINTS', '']
            slice5Display.DiffuseColor = [0.0, 0.0, 0.0]
            slice5Display.LineWidth = 4.0
            slice5Display.OSPRayScaleArray = 'f'
            slice5Display.OSPRayScaleFunction = 'PiecewiseFunction'
            slice5Display.SelectOrientationVectors = 'None'
            slice5Display.ScaleFactor = 0.3379646740552927
            slice5Display.SelectScaleArray = 'f'
            slice5Display.GlyphType = 'Arrow'
            slice5Display.GlyphTableIndexArray = 'f'
            slice5Display.GaussianRadius = 0.003989625424146652
            slice5Display.SetScaleArray = ['POINTS', 'f']
            slice5Display.ScaleTransferFunction = 'PiecewiseFunction'
            slice5Display.OpacityArray = ['POINTS', 'f']
            slice5Display.OpacityTransferFunction = 'PiecewiseFunction'
            slice5Display.DataAxesGrid = 'GridAxesRepresentation'
            slice5Display.PolarAxes = 'PolarAxesRepresentation'

            print ("Showing transform1.. "),
            Show(transform1, renderView1)
            uxLUTColorBar.Visibility = 1

            # show color bar/color legend
            connectivity1Display.SetScalarBarVisibility(renderView2, True)
            # get color legend/bar for uxLUT in view renderView2
            uxLUTColorBar2 = GetScalarBar(uxLUT, renderView2)

            uxLUTColorBar2.ComponentTitle = ''
            uxLUTColorBar2.Title = '|u|'
            uxLUTColorBar2.ComponentTitle = ''
            uxLUTColorBar2.TitleColor = [1, 1, 1]
            uxLUTColorBar2.LabelColor = [1, 1, 1]
            uxLUTColorBar2.LabelFormat = '%-#6.2g'
            uxLUTColorBar2.RangeLabelFormat = '%6.2g'

            len_bar = 0.5
            # change scalar bar placement
            uxLUTColorBar2.Orientation = 'Horizontal'
            uxLUTColorBar2.WindowLocation = 'AnyLocation'
            uxLUTColorBar2.Position = [0.5 - 0.5*len_bar, 0.02]
            uxLUTColorBar2.ScalarBarLength = len_bar

            print ("Showing connectivity1.. "),
            Show(connectivity1, renderView2, 'GeometryRepresentation')
            print ("Showing slice4.. "),
            Show(slice4, renderView2, 'GeometryRepresentation')
            print ("Showing transform1.. "),
            Show(transform1, renderView2, 'GeometryRepresentation')
            renderView2.CameraParallelScale = 1.8
            print("RenderView 1,2 update..")
            renderView1.Update()
            renderView2.Update()
            print ("end")
            fn = path + "/" + picName +  '_t=' + str(timesteps[i]) +'_n=' + str(k) + '_uxSlice.png'

            SaveScreenshot( fn, layout1, SaveAllViews=1,
                 ImageResolution=[4096, 2048],
                TransparentBackground=0,
                CompressionLevel='2' )

            SetActiveView(renderView2)
            camera = GetActiveCamera()
            camera.Azimuth(-10)
            Delete(text1)
            del text1
            print('File=' + fn + ' generated succesfully')
        Hide(slice4, renderView1)
        Hide(slice5, renderView1)
        Hide(connectivity1, renderView2)
        Hide(slice5, renderView2)


#**************************************************************
#**************** SLICES of OMEGA along a bubble **************
#**************************************************************
        ss = GetSources()
        for so in ss:
          Hide(ss[so])

        # get color transfer function/color map for 'omega'
        omegaLUT = GetColorTransferFunction('absOmega')
        omegaLUT.RGBPoints = [0.3730543491279064, 0.0, 0.0, 0.5625, 0.5106463634876334, 0.0, 0.0, 1.0, 0.8251430154745502, 0.0, 1.0, 1.0, 0.9823910318856668, 0.5, 1.0, 0.5, 1.139639048296783, 1.0, 1.0, 0.0, 1.4541357002836997, 1.0, 0.0, 0.0, 1.6113837166948162, 0.5, 0.0, 0.0]
        omegaLUT.ColorSpace = 'RGB'
        omegaLUT.ScalarRangeInitialized = 1.0
        # omegaLUT.ApplyPreset('jet', True)
        len_bar = 0.5

        # get color legend/bar for omegaLUT in view renderView1
        omegaLUTColorBar = GetScalarBar(omegaLUT, renderView2)
        omegaLUTColorBar.Orientation = 'Horizontal'
        omegaLUTColorBar.WindowLocation = 'AnyLocation'
        omegaLUTColorBar.ScalarBarLength = len_bar
        omegaLUTColorBar.Position = [0.5 - 0.5*len_bar, 0.02]
        omegaLUTColorBar.Title = '|Omega|'
        omegaLUTColorBar.ComponentTitle = ''
        omegaLUTColorBar.TitleColor = [1, 1, 1]
        omegaLUTColorBar.LabelColor = [1, 1, 1]
        omegaLUTColorBar.LabelFormat = '%-#6.2g'
        omegaLUTColorBar.RangeLabelFormat = '%6.2g'
        omegaLUTColorBar.ScalarBarThickness = 16
        omegaLUTColorBar.TitleFontSize = 16
        omegaLUTColorBar.LabelFontSize = 16

        print ("Showing resampleToImage1.GetDataInformation.. "),
        info = resampleToImage1.GetDataInformation().DataInformation
        arrayInfo = info.GetArrayInformation("absOmega", vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS)
        range0 = arrayInfo.GetComponentRange(0)
        range_max = np.max(np.abs(range0))
        print("|Omega| range= ", range0)
        omegaLUT.RescaleTransferFunction(0, range_max)
        op = GetOpacityTransferFunction("absOmega")
        op.RescaleTransferFunction(0, range_max)
        omegaLUTColorBar.CustomLabels = np.ceil(np.linspace(0, range_max, 5))
        print("CustomLabels=",omegaLUTColorBar.CustomLabels )


        for k,sl in enumerate(np.linspace(bounds[0]+0.1, bounds[1]-0.1, nslices)):
            print("slice x=", sl)

            # create a new 'Text'
            text1 = Text()
            text1.Text = 'x=' + str(round(sl,2)) + " l/l_b=" + str(round((sl-bounds[0])/len_bub,2))

            renderView1.ViewSize = [2048, 2048]
            renderView1.CameraViewUp = [0, 0, 1]
            renderView1.ResetCamera(sl, sl, -0.5, 0.5, -0.5, 0.5)

            renderView2.CameraViewUp = [0, 1, 0]
            renderView2.ResetCamera(bounds[0], bounds[1], -0.5, 0.5, -0.5, 0.5)
            SetActiveView(renderView2)
            camera = GetActiveCamera()
            camera.Azimuth(10)

            print("renderView2:",bounds[0], bounds[1], -0.5, 0.5, -0.5, 0.5)

            # init the 'Plane' selected for 'SliceType'
            slice4.SliceType.Origin = [sl, 0.0, 0.0]
            slice4.SliceType.Normal = [1.0, 0.0, 0.0]

            # show data from text1
            text1Display = Show(text1, renderView2, 'TextSourceRepresentation')

            # trace defaults for the display properties.
            text1Display.WindowLocation = 'AnyLocation'
            text1Display.Position = [0.3, 0.2]




            slice4Display = [0,0]
            for id,rv in enumerate([renderView1,renderView2]):
                print ("Showing slice4Display.. "),
                slice4Display[id] = Show(slice4, rv, 'GeometryRepresentation')

                # trace defaults for the display properties.
                slice4Display[id].Representation = 'Surface'
                slice4Display[id].ColorArrayName = ['POINTS', 'absOmega']
                slice4Display[id].LookupTable = omegaLUT
                slice4Display[id].OSPRayScaleArray = 'f'
                slice4Display[id].OSPRayScaleFunction = 'PiecewiseFunction'
                slice4Display[id].SelectOrientationVectors = 'None'
                slice4Display[id].ScaleFactor = 0.3379646740552927
                slice4Display[id].SelectScaleArray = 'f'
                slice4Display[id].GlyphType = 'Arrow'
                slice4Display[id].GlyphTableIndexArray = 'f'
                slice4Display[id].GaussianRadius = 0.016898233702764633
                slice4Display[id].SetScaleArray = ['POINTS', 'f']
                slice4Display[id].ScaleTransferFunction = 'PiecewiseFunction'
                slice4Display[id].OpacityArray = ['POINTS', 'f']
                slice4Display[id].OpacityTransferFunction = 'PiecewiseFunction'
                slice4Display[id].DataAxesGrid = 'GridAxesRepresentation'
                slice4Display[id].PolarAxes = 'PolarAxesRepresentation'

                # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
                slice4Display[id].OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

                # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
                slice4Display[id].ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0000000000000002, 1.0, 0.5, 0.0]

                # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
                slice4Display[id].OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0000000000000002, 1.0, 0.5, 0.0]


            # init the 'Plane' selected for 'SliceType'
            print ("Showing slice5.. "),
            slice5.SliceType.Origin = [sl, 0.0, 0.0]
            slice5.SliceType.Normal = [1.0, 0.0, 0.0]

            slice5Display = Show(slice5, renderView1, 'GeometryRepresentation')

            # trace defaults for the display properties.
            slice5Display.Representation = 'Surface'
            slice5Display.AmbientColor = [0.0, 0.0, 0.0]
            slice5Display.ColorArrayName = ['POINTS', '']
            slice5Display.DiffuseColor = [0.0, 0.0, 0.0]
            slice5Display.LineWidth = 4.0
            slice5Display.OSPRayScaleArray = 'f'
            slice5Display.OSPRayScaleFunction = 'PiecewiseFunction'
            slice5Display.SelectOrientationVectors = 'None'
            slice5Display.ScaleFactor = 0.3379646740552927
            slice5Display.SelectScaleArray = 'f'
            slice5Display.GlyphType = 'Arrow'
            slice5Display.GlyphTableIndexArray = 'f'
            slice5Display.GaussianRadius = 0.003989625424146652
            slice5Display.SetScaleArray = ['POINTS', 'f']
            slice5Display.ScaleTransferFunction = 'PiecewiseFunction'
            slice5Display.OpacityArray = ['POINTS', 'f']
            slice5Display.OpacityTransferFunction = 'PiecewiseFunction'
            slice5Display.DataAxesGrid = 'GridAxesRepresentation'
            slice5Display.PolarAxes = 'PolarAxesRepresentation'
            print ("Showing transform1.. "),
            Show(transform1, renderView1)
            info = slice4.GetDataInformation().DataInformation
            arrayInfo = info.GetArrayInformation("absOmega", vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS)
            range0 = arrayInfo.GetComponentRange(0)
            range_max_slice = np.max(np.abs(range0))
            print("|Omega| range= ", range0, " range_max_slice=", range_max_slice)
#             slice4Display[0].SetScalarBarVisibility(renderView1, True)

            omegaLUT1 = GetColorTransferFunction('absOmega', slice4Display[0], separate=True)
            omegaLUT1.RGBPoints = [0.3730543491279064, 0.0, 0.0, 0.5625, 0.5106463634876334, 0.0, 0.0, 1.0, 0.8251430154745502, 0.0, 1.0, 1.0, 0.9823910318856668, 0.5, 1.0, 0.5, 1.139639048296783, 1.0, 1.0, 0.0, 1.4541357002836997, 1.0, 0.0, 0.0, 1.6113837166948162, 0.5, 0.0, 0.0]
            omegaLUT1.ColorSpace = 'RGB'
            omegaLUT1.ScalarRangeInitialized = 1.0
            # get separate opacity transfer function/opacity map for 'ux'
            omegaPWF1 = GetOpacityTransferFunction('absOmega', slice4Display[0], separate=True)
            omegaPWF1.Points = [2.2554501288418073e-07, 0.0, 0.5, 0.0, 2.1026964870830196, 1.0, 0.5, 0.0]
            omegaPWF1.ScalarRangeInitialized = 1

#             # get color legend/bar for omegaLUT in view renderView1
            omegaLUTColorBar1 = GetScalarBar(omegaLUT1, renderView1)
            omegaLUT1.RescaleTransferFunction(0, range_max_slice)
            op = GetOpacityTransferFunction("absOmega")
            op.RescaleTransferFunction(0, range_max_slice)
            omegaLUTColorBar1.UseCustomLabels = 1
            omegaLUTColorBar1.CustomLabels = np.ceil(np.linspace(0, range_max_slice, 5))
            print("CustomLabels |Omega|=",omegaLUTColorBar1.CustomLabels )
            len_bar = 0.5
            omegaLUTColorBar1.Orientation = 'Horizontal'
            omegaLUTColorBar1.WindowLocation = 'AnyLocation'
            omegaLUTColorBar1.Position = [0.5 - 0.5*len_bar, 0.02]
            omegaLUTColorBar1.Title = '|Omega|'
            omegaLUTColorBar1.ComponentTitle = ''
            omegaLUTColorBar1.LabelFormat = '%-#6.2g'
            omegaLUTColorBar1.RangeLabelFormat = '%6.2g'
            omegaLUTColorBar1.TitleColor = [1, 1, 1]
            omegaLUTColorBar1.LabelColor = [1, 1, 1]
            omegaLUTColorBar1.ScalarBarLength = len_bar
            omegaLUTColorBar1.ScalarBarThickness = 16
            omegaLUTColorBar1.TitleFontSize = 16
            omegaLUTColorBar1.LabelFontSize = 16

            # set color bar visibility
            omegaLUTColorBar1.Visibility = 1
            omegaLUTColorBar.Visibility = 1
            uxLUTColorBar.Visibility = 0
#             omegaLUTColorBar.CustomLabels = np.ceil(np.linspace(0, range_max, 5))


            slice4Display[0].LookupTable = omegaLUT1
#             print ("Showing connectivity1 renderView2.. "),
#             Show(connectivity1, renderView2, 'GeometryRepresentation') # Weugene: ??? was commented deliberately
#             print ("end")
            # show data from connectivity1
            print ("Showing connectivity1.. "),
            connectivity1Display = Show(connectivity1, renderView2, 'GeometryRepresentation')
            # trace defaults for the display properties.

            connectivity1Display.ColorArrayName = ['POINTS', 'absOmega']
            connectivity1Display.LookupTable = omegaLUT
#             connectivity1Display.Opacity = 1
#             connectivity1Display.Specular = 0.8
#             connectivity1Display.SpecularPower = 100.0
#             connectivity1Display.Ambient = 0.1
#             connectivity1Display.OSPRayScaleArray = 'f'
#             connectivity1Display.OSPRayScaleFunction = 'PiecewiseFunction'
#             connectivity1Display.SelectOrientationVectors = 'None'
#             connectivity1Display.ScaleFactor = 0.3379646740552927
#             connectivity1Display.SelectScaleArray = 'f'
#             connectivity1Display.GlyphType = 'Arrow'
#             connectivity1Display.GlyphTableIndexArray = 'f'
#             connectivity1Display.GaussianRadius = 0.016898233702764633
#             connectivity1Display.SetScaleArray = ['POINTS', 'f']
#             connectivity1Display.ScaleTransferFunction = 'PiecewiseFunction'
#             connectivity1Display.OpacityArray = ['POINTS', 'f']
#             connectivity1Display.OpacityTransferFunction = 'PiecewiseFunction'
#             connectivity1Display.DataAxesGrid = 'GridAxesRepresentation'
#             connectivity1Display.PolarAxes = 'PolarAxesRepresentation'
            # show color bar/color legend
            connectivity1Display.SetScalarBarVisibility(renderView2, True)
            print ("Showing slice4, transform1 renderView2.. "),
            Show(slice4, renderView2, 'GeometryRepresentation')
            Show(transform1, renderView2, 'GeometryRepresentation')
            renderView2.CameraParallelScale = 1.8
            renderView1.Update()
            renderView2.Update()
            fn = path + "/" + picName +  '_t=' + str(timesteps[i]) +'_n=' + str(k) + '_omegaSlice.png'


            SaveScreenshot( fn,  layout1, SaveAllViews=1,
                 ImageResolution=[4096, 2048],
                TransparentBackground=0,
                CompressionLevel='2' )

#             Hide(connectivity1, renderView2)
#             Hide(slice4, renderView2)
#             Hide(slice5, renderView2)
            SetActiveView(renderView2)
            camera = GetActiveCamera()
            camera.Azimuth(-10)
            # set color bar visibility
            omegaLUTColorBar1.Visibility = 0
            omegaLUTColorBar.Visibility = 0
            uxLUTColorBar.Visibility = 0

            Delete(text1)
            del text1
#             Hide(linearExtrusion1, renderView2)
            print('File=' + fn + ' generated succesfully')



        print("Contours")
        # show data in view
        contour3Display = Show(slice1, spreadSheetView1, 'SpreadSheetRepresentation')
        fn = path + '/' + 'contour_t=' + str(timesteps[i]) + '.csv'
        ExportView(fn, view=spreadSheetView1)
        print('File=' + fn + ' generated succesfully')

        #print("Tracers")
        # show data in view
        #contour3Display = Show(streamTracer1, spreadSheetView2, 'SpreadSheetRepresentation')
        #fn = path + '/' + 'tracer_t=' + str(timesteps[i]) + '.csv'
        #ExportView(fn, view=spreadSheetView2)
        #print('File=' + fn + ' generated succesfully')


        # Freeing Memory
        Delete(slice5)
        del slice5
        Delete(slice4)
        del slice4
        Delete(slice3)
        del slice3
        Delete(streamTracer3)
        del streamTracer3
        Delete(streamTracer2)
        del streamTracer2
        Delete(streamTracer1)
        del streamTracer1
        Delete(calculator1)
        del calculator1
        Delete(extractSelection1)
        del extractSelection1
        Delete(contour5)
        del contour5
        Delete(isoVolume1)
        del isoVolume1
        Delete(slice1)
        del slice1
        Delete(contour2)
        del contour2
        Delete(connectivity1)
        del connectivity1
        Delete(contour1)
        del contour1
        Delete(resampleToImage1)
        del resampleToImage1
        Delete(calculator0)
        del calculator0
        Delete(cellDatatoPointData1)
        del cellDatatoPointData1
        Delete(clip1)
        del clip1
#         Delete(threshold0)
#         del threshold0
#         Delete(resampleToImage0)
#         del resampleToImage0
        Delete(clip0)
        del clip0



stop = timeit.default_timer()
total_time = stop - start

# output running time in a nice format.
mins, secs = divmod(total_time, 60)
hours, mins = divmod(mins, 60)
my_stats.close()

sys.stdout.write("Total running time: %d:%d:%d.\n" % (hours, mins, secs))
