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
from vtk.numpy_interface import dataset_adapter as dsa
logging.basicConfig(format='%(message)s')
log = logging.getLogger(__name__)

import functools
# import __builtin__
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
#1 video Name
#2 pvd file name (by defaults in the first file in a current directory)
#3 pvtu is swithed off by default
# ---------------------------------------------------------------------------------------------------------
start = timeit.default_timer()

parser = argparse.ArgumentParser()
parser.add_argument("-vidName", type=str, help="Provide the name for the outputed video, please",
                    nargs='?', default='video')
parser.add_argument("-picName", type=str, help="Provide the name for the outputed pictures, please",
                    nargs='?', default='pic')
parser.add_argument("-filenames", type=str, help="Provide the name of the outputing paraview files, please",
                    nargs='?', default='save*.pvd')
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
parser.add_argument("-lDomain", type=float, help="Provide the length of the channel",
                    nargs='?', default=30)


parser.add_argument("-viewSize", type=int, help="Provide the view size of the output pictures, please",
                    nargs='+', default=[1900, 1078])
parser.add_argument("-noVideo", type=bool, help="Provide the no video Mode",
                    nargs='?', default=False)
parser.add_argument("-noPic", type=bool, help="Provide the no Picture Mode",
                    nargs='?', default=False)
parser.add_argument("-noData", type=bool, help="Provide the no Data Exporting Mode",
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
lDomain = args.lDomain
viewSize = args.viewSize
noVideo = args.noVideo
noPic = args.noPic
noData = args.noData
timeList = args.timeList
nt = args.nt

scale_size = 1.38#2048./1078.
if len(frameWindow) != 2 or frameWindow[0]<0 or frameWindow[1] < 0:
    eprint('Error in frameWindow' + frameWindow)
    sys.exit()

#Current PATH reading
if vtk_from_pvpython:
    path = os.path.abspath(os.getcwd())

eprint("Current PATH=" + path)
if vtk_from_pvpython:
    if filenames[-4::] == '.pvd':
    # create a new 'PVD Reader'
        fn = glob.glob(filenames)
        print('Found files:',fn)
        fn = fn[0]
        print('Read the first one:',fn)
        my_source = PVDReader(FileName = path +  '/' + fn)
    else:
        eprint('Get Active Source: No pvd or pvtu files are provided')
        my_source = GetActiveSource()
else:
    eprint('Get Active Source')
    my_source = GetActiveSource()

my_source.CellArrays = [ 'f', 'omega', 'u.x', 'p']
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
if NT != 0:
    print("NT=", NT, " timeList=", timeList)
else:
    print("ERROR: can't read pvd file NT=0")
if frameWindow[0] == 0 and frameWindow[1] == 0:
    frameWindow[0] = 0
    frameWindow[1] = NT-1
    print('frameWindow is updated:',frameWindow)

print ("renderViews, axesGrid, SpreadSheetViews, layouts.. ")


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


# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------
# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.AssignView(0, renderView1)
layout1.SetSize(1840, 1156)
# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)


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

Show(transform1, renderView1)


print ("end")
###-----------------GENERATION of a CLIP and Convert to Point data-----------------------------------
global my_stats # output

if not noPic:
    for i in timeList:
        try:
            print("in loop iteration:" + str(i))
            SetActiveView(renderView1)
            # Properties modified on animationScene1
            animationScene1.AnimationTime = timesteps[i]
            # Properties modified on timeKeeper1
            timeKeeper1.Time = timesteps[i]
            if i==0 and vtk_from_pvpython:
                fn = str(timesteps[i]) + 'stats.txt'
                my_stats = open(fn, 'w')

# ***************** CELL DATA TO POINT DATA  from BOX of BOX ****************************
            # create a new 'Cell Data to Point Data'
            cellDatatoPointData1 = CellDatatoPointData(Input=my_source)

# ***************** CALCULATE OMEGA ****************************
            # create a new 'Calculator' for OMEGA
            calculator0 = Calculator(Input=cellDatatoPointData1)
            calculator0.ResultArrayName = 'absOmega'
            calculator0.Function = 'abs(omega)'

            Show(transform1, renderView1)
            calculator0Display = Show(calculator0, renderView1)
            # defining of computational domains
            bounds = calculator0.GetDataInformation().GetBounds()
            print("bounds_refined=",bounds)
            center = [(bounds[0] + bounds[1])/2, (bounds[2] + bounds[3])/2,(bounds[4] + bounds[5])/2]
            print('center_refined=',center)

            len_bub = bounds[1] - bounds[0]
            len_min = bounds[0] - 2
            len_max = bounds[1] + 2
            length = len_max - len_min

            print("len_bub=", len_bub, "len_min=", len_min, "len_max=", len_max, "length of the domain=", length)

# ----------------------------------------------------------------
        # setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

            # get color transfer function/color map for 'ux'
            # uxLUT = GetColorTransferFunction('ux')
            # uxLUT.RGBPoints = [0.3730543491279064, 0.0, 0.0, 0.5625, 0.5106463634876334, 0.0, 0.0, 1.0, 0.8251430154745502, 0.0, 1.0, 1.0, 0.9823910318856668, 0.5, 1.0, 0.5, 1.139639048296783, 1.0, 1.0, 0.0, 1.4541357002836997, 1.0, 0.0, 0.0, 1.6113837166948162, 0.5, 0.0, 0.0]
            # uxLUT.ColorSpace = 'RGB'
            # uxLUT.ScalarRangeInitialized = 1.0
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
            info = calculator0.GetDataInformation().DataInformation
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
            uxLUTColorBar.CustomLabels = np.around(np.linspace(0, range_max, 4),1)
            #         uxLUTColorBar.CustomLabels = np.arange(0, range_max, 0.5)
            print("CustomLabels=",uxLUTColorBar.CustomLabels )


            # trace defaults for the display properties.
            calculator0Display.Representation = 'Surface'
            calculator0Display.ColorArrayName = ['POINTS', 'u.x']
            calculator0Display.LookupTable = uxLUT
            calculator0Display.Specular = 1.0
            calculator0Display.SelectTCoordArray = 'None'
            calculator0Display.SelectNormalArray = 'None'
            calculator0Display.SelectTangentArray = 'None'
            calculator0Display.OSPRayScaleArray = 'Result'
            calculator0Display.OSPRayScaleFunction = 'PiecewiseFunction'
            calculator0Display.SelectOrientationVectors = 'None'
            calculator0Display.ScaleFactor = 0.31534409523010254
            calculator0Display.SelectScaleArray = 'ux'
            calculator0Display.GlyphType = 'Arrow'
            calculator0Display.GlyphTableIndexArray = 'Result'
            calculator0Display.GaussianRadius = 0.015767204761505126
            calculator0Display.SetScaleArray = ['POINTS', 'ux']
            calculator0Display.ScaleTransferFunction = 'PiecewiseFunction'
            calculator0Display.OpacityArray = ['POINTS', 'ux']
            calculator0Display.OpacityTransferFunction = 'PiecewiseFunction'
            calculator0Display.DataAxesGrid = 'GridAxesRepresentation'
            calculator0Display.PolarAxes = 'PolarAxesRepresentation'

            # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
            calculator0Display.ScaleTransferFunction.Points = [0.005237827916105199, 0.0, 0.5, 0.0, 0.46907547386229526, 1.0, 0.5, 0.0]

            # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
            calculator0Display.OpacityTransferFunction.Points = [0.005237827916105199, 0.0, 0.5, 0.0, 0.46907547386229526, 1.0, 0.5, 0.0]

            axesGrid.XTitleFontSize = 20
            axesGrid.YTitleFontSize = 20
            axesGrid.ZTitleFontSize = 20
            axesGrid.XLabelFontSize = 20
            axesGrid.YLabelFontSize = 20
            axesGrid.ZLabelFontSize = 20
            renderView1.CameraViewUp = [0.2, 1, 0]
            renderView1.CameraParallelScale = 1. #1.4 0.5
            renderView1.CenterOfRotation = center
            renderView1.CameraFocalPoint = center
            renderView1.CameraPosition = [center[0] - 4, 0.6, 4.5]
            # update the view to ensure updated data information
            renderView1.Update()


#****************** CONNECTIVITY(f) AND U MAGNITUDE ********************
            # show data from connectivity1
            fn = path + "/" + picName +  '_t=' + str(timesteps[i]) +'ux.png'
            SaveScreenshot( fn, renderView1,
                ImageResolution=[1900, 1078],
                TransparentBackground=0,
                CompressionLevel='2' )
            print('File=' + fn + ' generated succesfully')
            # break



#***************** Like in the article SIDE SLICE and contour(f)*********************
        # create a new 'Slice'
        # slice3 = Slice(Input=calculator0)
        # slice3.SliceType = 'Plane'
        # slice3.HyperTreeGridSlicer = 'Plane'
        # slice3.SliceOffsetValues = [0]
        #
        # # init the 'Plane' selected for 'SliceType'
        # slice3.SliceType.Origin = [0.5*(len_min + len_max), 0.0, 0.0]
        # slice3.SliceType.Normal = [0.0, 0.0, 1.0]
        # # init the 'Plane' selected for 'HyperTreeGridSlicer'
        # # slice3.HyperTreeGridSlicer.Origin = [10.999999761581421, 0.0, 0.0] #???
        # print ("Showing slice3Display.. "),
        # slice3Display = Show(slice3, renderView1, 'GeometryRepresentation')
        #
        # # trace defaults for the display properties.
        # slice3Display.Representation = 'Surface'
        # slice3Display.ColorArrayName = ['POINTS', 'u.x']
        # slice3Display.LookupTable = uxLUT
        # slice3Display.OSPRayScaleArray = 'f'
        # slice3Display.OSPRayScaleFunction = 'PiecewiseFunction'
        # slice3Display.SelectOrientationVectors = 'None'
        # slice3Display.ScaleFactor = 0.3379646740552927
        # slice3Display.SelectScaleArray = 'f'
        # slice3Display.GlyphType = 'Arrow'
        # slice3Display.GlyphTableIndexArray = 'f'
        # slice3Display.GaussianRadius = 0.016898233702764633
        # slice3Display.SetScaleArray = ['POINTS', 'f']
        # slice3Display.ScaleTransferFunction = 'PiecewiseFunction'
        # slice3Display.OpacityArray = ['POINTS', 'f']
        # slice3Display.OpacityTransferFunction = 'PiecewiseFunction'
        # slice3Display.DataAxesGrid = 'GridAxesRepresentation'
        # slice3Display.PolarAxes = 'PolarAxesRepresentation'
        #
        # # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
        # slice3Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]
        #
        # # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
        # slice3Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0000000000000002, 1.0, 0.5, 0.0]
        #
        # # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
        # slice3Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0000000000000002, 1.0, 0.5, 0.0]


            uxLUTColorBar.Position = [0.5 - 0.5*len_bar, 0.02]
            uxLUTColorBar.Visibility = 1
            uxLUTColorBar.ScalarBarThickness = int(16*scale_size)
            uxLUTColorBar.TitleFontSize = int(16*scale_size)
            uxLUTColorBar.LabelFontSize = int(16*scale_size)
            print("RenderView update..")
            axesGrid.XTitleFontSize = int(20*scale_size)
            axesGrid.YTitleFontSize = int(20*scale_size)
            axesGrid.ZTitleFontSize = int(20*scale_size)
            axesGrid.XLabelFontSize = int(20*scale_size)
            axesGrid.YLabelFontSize = int(20*scale_size)
            axesGrid.ZLabelFontSize = int(20*scale_size)
            renderView1.CameraPosition = [center[0]-3, 0.0, 0.0]
            renderView1.CameraFocalPoint = [center[0], 0.0, 0.0]
            renderView1.CameraViewUp = [0.0, 0.0, 1.0]
            renderView1.CameraParallelScale = 1.9
            # renderView1.ViewSize = [2048, 2048]
            renderView1.ResetCamera(center[0], center[0], -0.5, 0.5, -0.5, 0.5)
            # update the view to ensure updated data information
            renderView1.Update()
            print ("end")
            fn = path + "/" + picName +  '_t=' + str(timesteps[i]) +'_ux_tail.png'
            SaveScreenshot( fn, renderView1,
                ImageResolution=[2048, 2048],
                TransparentBackground=0,
                CompressionLevel='2' )
#         Hide(streamTracer1, renderView1)
#         Hide(slice3, renderView1)

            print('File=' + fn + ' generated succesfully')


            # Freeing Memory


            # Delete(slice3)
            # del slice3
            Delete(calculator0)
            del calculator0
            Delete(cellDatatoPointData1)
            del cellDatatoPointData1
        except:
            print('ERROR: the time step does not exist')



stop = timeit.default_timer()
total_time = stop - start

# output running time in a nice format.
mins, secs = divmod(total_time, 60)
hours, mins = divmod(mins, 60)
if vtk_from_pvpython:
    my_stats.close()

sys.stdout.write("Total running time: %d:%d:%d.\n" % (hours, mins, secs))
