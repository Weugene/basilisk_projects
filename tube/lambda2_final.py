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
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()
import glob, os, sys
import logging
from sys import argv
import timeit
logging.basicConfig(format='%(message)s')
log = logging.getLogger(__name__)

def eprint(var):
    log.warning(var)
# Read from arguments
#1 video Name
#2 pvd file name (by defaults in the first file in a current directory)
#3 pvtu is swithed off by default

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-vidName", type=str, help="Provide the name for the outputed video, please",
                    nargs='?', default='lambda2_')
parser.add_argument("-picName", type=str, help="Provide the name for the outputed pictures, please",
                    nargs='?', default='lambda2_')
parser.add_argument("-filenames", type=str, help="Provide the name of the outputing paraview files, please",
                    nargs='?', default='*.pvd')
parser.add_argument("-frameRate", type=int, help="Provide the frame rate for the video, please",
                    nargs='?', default=10)
parser.add_argument("-frameWindow", type=int, help="Provide the frame window (can be overwritten further), please",
                    nargs='+', default=[0,0])
parser.add_argument("-viewSize", type=int, help="Provide the view size of the output pictures, please",
                    nargs='+', default=[2152, 862])
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
viewSize = args.viewSize
noVideo = args.noVideo
noPic = args.noPic
noData = args.noData

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
    sys.exit()
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# get animation scene
animationScene1 = GetAnimationScene()

# get the time-keeper
timeKeeper1 = GetTimeKeeper()
timesteps = timeKeeper1.TimestepValues # 0, 0.1, 0.2 ...
NT = len(timesteps)
print("NT=", NT, " timesteps=", timesteps)

if frameWindow[0] == 0 and frameWindow[1] == 0:
    frameWindow[0] = 0
    frameWindow[1] = NT-1
    print('frameWindow is updated:',frameWindow)

# ---------------------------------------------------------------------------------------------------------
start = timeit.default_timer()

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [2150, 1354]
renderView1.InteractionMode = '2D'
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.CenterOfRotation = [15.0, 0.0, 0.0]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [10.15847112215707, 0.5977000151160986, 4.5392388055413]
renderView1.CameraFocalPoint = [16.142004301693653, -0.27836824652388686, -1.3233473689572015]
renderView1.CameraViewUp = [0.2035575809744966, 0.9771141634582637, 0.06174320041161353]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 1.8015930673699396
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

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# # create a new 'PVD Reader'
# my_source = PVDReader(FileName='/Users/weugene/basilisk/work/tube/dump2pvd_compressed.pvd')
# my_source.CellArrays = ['p', 'fs', 'f', 'l', 'l2', 'omega', 'u.x']

# create a new 'Cylinder'
cylinder1 = Cylinder()
cylinder1.Resolution = 140
cylinder1.Height = 30.0
cylinder1.Capping = 0

# create a new 'Clip'
clip2 = Clip(Input=cylinder1)
clip2.ClipType = 'Plane'
clip2.HyperTreeGridClipper = 'Plane'
clip2.Scalars = ['POINTS', 'Normals_Magnitude']
clip2.Value = 1

# init the 'Plane' selected for 'ClipType'
clip2.ClipType.Normal = [0.0, 0.0, 1.0]

# create a new 'Transform'
transform1 = Transform(Input=clip2)
transform1.Transform = 'Transform'

# init the 'Transform' selected for 'Transform'
transform1.Transform.Translate = [15.0, 0.0, 0.0]
transform1.Transform.Rotate = [0.0, 0.0, 90.0]

# create a new 'Contour' for estimation of boundaries
contour0 = Contour(Input=my_source)
contour0.ContourBy = ['POINTS', 'f']
contour0.Isosurfaces = [0.5]
contour0.PointMergeMethod = 'Uniform Binning'

# show data in view
contour0Display = Show(contour0, renderView1)
Hide(contour0, renderView1)
boundsContour0 = contour0.GetDataInformation().GetBounds()
print("bounds contour0=", boundsContour0)

# create a new 'Clip'
clip1 = Clip(Input=my_source)
clip1.ClipType = 'Box'
clip1.HyperTreeGridClipper = 'Plane'
clip1.Scalars = ['POINTS', 'f']
clip1.Value = 0.5

# init the 'Box' selected for 'ClipType'
len_min = boundsContour0[0] - 3.5
len_max = boundsContour0[1] + 1.2
length = len_max -len_min
clip1.ClipType.Position = [len_min, -0.5, -0.5]

clip1.ClipType.Length = [length, 1.0, 1.0]

# init the 'Plane' selected for 'HyperTreeGridClipper'
# clip1.HyperTreeGridClipper.Origin = [15.0, 0.0, 0.0]

# create a new 'Cell Data to Point Data'
cellDatatoPointData1 = CellDatatoPointData(Input=clip1)
cellDatatoPointData1.CellDataArraytoprocess = ['f', 'fs', 'l', 'l2', 'omega', 'p', 'u.x']

# create a new 'Delaunay 3D'
delaunay3D1 = Delaunay3D(Input=cellDatatoPointData1)
delaunay3D1.AlphaVerts = 1

# create a new 'Contour'
contour1 = Contour(Input=delaunay3D1)
contour1.ContourBy = ['POINTS', 'f']
contour1.Isosurfaces = [0.5]
contour1.PointMergeMethod = 'Uniform Binning'

# create a new 'Contour'
contour2 = Contour(Input=delaunay3D1)
contour2.ContourBy = ['POINTS', 'l2']
contour2.Isosurfaces = [-1.0]
contour2.PointMergeMethod = 'Uniform Binning'

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from contour1
contour1Display = Show(contour1, renderView1, 'GeometryRepresentation')
contour1Display.RescaleTransferFunctionToDataRange(False, True)

# get color transfer function/color map for 'ux'
uxLUT = GetColorTransferFunction('ux')

info = contour1.GetDataInformation().DataInformation
arrayInfo = info.GetArrayInformation("u.x", vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS)
range0 = arrayInfo.GetComponentRange(0)
range1 = arrayInfo.GetComponentRange(1)
range2 = arrayInfo.GetComponentRange(2)
range = [np.sqrt(range0[0]**2+range1[0]**2+range2[0]**2), np.sqrt(range0[1]**2+range1[1]**2+range2[1]**2)]
print("velocity range= ", range)
print("velocity range1= ", range1)
print("velocity range2= ", range2)
print("velocity range3= ", range3)
tf = GetColorTransferFunction("ux")
tf.RescaleTransferFunction(0, range[1])
op = GetOpacityTransferFunction("ux")
op.RescaleTransferFunction(0, range[1])

# uxLUT.RGBPoints = [0, 0.0, 0.0, 0.5625, 0.23363290986076246, 0.0, 0.0, 1.0, 0.7676515253591777, 0.0, 1.0, 1.0, 1.0346603074343204, 0.5, 1.0, 0.5, 1.301669089509462, 1.0, 1.0, 0.0, 1.8356877050078775, 1.0, 0.0, 0.0, 2.1026964870830196, 0.5, 0.0, 0.0]
uxLUT.ColorSpace = 'RGB'
uxLUT.ScalarRangeInitialized = 1.0
# Rescale transfer function set minimal velocity
uxLUT.RescaleTransferFunction(0.0, 1.0)

# trace defaults for the display properties.
contour1Display.Representation = 'Surface'
contour1Display.ColorArrayName = ['POINTS', 'u.x']
contour1Display.LookupTable = uxLUT
contour1Display.Specular = 1.0
contour1Display.Luminosity = 49.0
contour1Display.Ambient = 0.13
contour1Display.OSPRayScaleArray = 'f'
contour1Display.OSPRayScaleFunction = 'PiecewiseFunction'
contour1Display.SelectOrientationVectors = 'None'
contour1Display.ScaleFactor = 0.3379646740552927
contour1Display.SelectScaleArray = 'f'
contour1Display.GlyphType = 'Arrow'
contour1Display.GlyphTableIndexArray = 'f'
contour1Display.GaussianRadius = 0.016898233702764633
contour1Display.SetScaleArray = ['POINTS', 'f']
contour1Display.ScaleTransferFunction = 'PiecewiseFunction'
contour1Display.OpacityArray = ['POINTS', 'f']
contour1Display.OpacityTransferFunction = 'PiecewiseFunction'
contour1Display.DataAxesGrid = 'GridAxesRepresentation'
contour1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
contour1Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
contour1Display.ScaleTransferFunction.Points = [0.5, 0.0, 0.5, 0.0, 0.5001220703125, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
contour1Display.OpacityTransferFunction.Points = [0.5, 0.0, 0.5, 0.0, 0.5001220703125, 1.0, 0.5, 0.0]

# show data from contour2
contour2Display = Show(contour2, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
contour2Display.Representation = 'Surface'
contour2Display.AmbientColor = [0.5764705882352941, 0.7607843137254902, 1.0]
contour2Display.ColorArrayName = ['POINTS', '']
contour2Display.DiffuseColor = [0.5764705882352941, 0.7607843137254902, 1.0]
contour2Display.Specular = 1.0
contour2Display.Luminosity = 46.0
contour2Display.Ambient = 0.15
contour2Display.OSPRayScaleArray = 'l2'
contour2Display.OSPRayScaleFunction = 'PiecewiseFunction'
contour2Display.SelectOrientationVectors = 'None'
contour2Display.ScaleFactor = 0.6500720447488617
contour2Display.SelectScaleArray = 'l2'
contour2Display.GlyphType = 'Arrow'
contour2Display.GlyphTableIndexArray = 'l2'
contour2Display.GaussianRadius = 0.03250360223744308
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

# show data from transform1
transform1Display = Show(transform1, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
transform1Display.Representation = 'Surface'
transform1Display.ColorArrayName = [None, '']
transform1Display.OSPRayScaleArray = 'Normals'
transform1Display.OSPRayScaleFunction = 'PiecewiseFunction'
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

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
transform1Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
transform1Display.ScaleTransferFunction.Points = [-2.220446049250313e-16, 0.0, 0.5, 0.0, 2.220446049250313e-16, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
transform1Display.OpacityTransferFunction.Points = [-2.220446049250313e-16, 0.0, 0.5, 0.0, 2.220446049250313e-16, 1.0, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for uxLUT in view renderView1
uxLUTColorBar = GetScalarBar(uxLUT, renderView1)
uxLUTColorBar.Orientation = 'Horizontal'
uxLUTColorBar.WindowLocation = 'AnyLocation'
uxLUTColorBar.Position = [0.5355813953488372, 0.10487444608567209]
uxLUTColorBar.Title = 'u'
uxLUTColorBar.ComponentTitle = 'Magnitude'
uxLUTColorBar.LabelFormat = '%-#6.2g'
uxLUTColorBar.RangeLabelFormat = '%6.2g'
uxLUTColorBar.ScalarBarLength = 0.3300000000000003

# set color bar visibility
uxLUTColorBar.Visibility = 1

# show color legend
contour1Display.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get opacity transfer function/opacity map for 'ux'
# uxPWF = GetOpacityTransferFunction('ux')
# uxPWF.Points = [2.2554501288418073e-07, 0.0, 0.5, 0.0, 2.1026964870830196, 1.0, 0.5, 0.0]
# uxPWF.ScalarRangeInitialized = 1

# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(my_source)
# ----------------------------------------------------------------
Hide(cylinder1,renderView1)
Render()
SaveScreenshot( "lambda2_final.png", renderView1,
#      ImageResolution=[2316, 2204],
    TransparentBackground=0,
#     magnification=1,
    # PNG options
    CompressionLevel='2' )
print('File=' + "lambda2_final.png" + ' generated succesfully')