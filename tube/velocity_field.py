# trace generated using paraview version 5.8.0
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

import glob, os, sys
import logging
from sys import argv
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
                    nargs='?', default='vid_')
parser.add_argument("-picName", type=str, help="Provide the name for the outputed pictures, please",
                    nargs='?', default='velocity_field_')
parser.add_argument("-filenames", type=str, help="Provide the name of the outputing paraview files, please",
                    nargs='?', default='*.pvd')
parser.add_argument("-frameRate", type=int, help="Provide the frame rate for the video, please",
                    nargs='?', default=10)
parser.add_argument("-frameWindow", type=int, help="Provide the frame window (can be overwritten further), please",
                    nargs='+', default=[0,0])
parser.add_argument("-contourList", type=int, help="Provide the list of contour (can be overwritten further), please. Has higher priority",
                    nargs='+', default=[0,0,0])
parser.add_argument("-contourStartEndStep", type=int, help="Provide the start, end and step of output (can be overwritten further), please",
                    nargs='+', default=[0,0,0])
parser.add_argument("-lDomain", type=float, help="Provide the domain size, please",
                    nargs='?', default=20.0)
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
contourList = args.contourList
contourStartEndStep = args.contourStartEndStep
lDomain = args.lDomain
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
print("NT=", NT)

if frameWindow[0] == 0 and frameWindow[1] == 0:
    frameWindow[0] = 0
    frameWindow[1] = NT-1
    print('frameWindow is updated:',frameWindow)

if sum(n < 0 for n in contourList):
    eprint('Error in contourList:' + contourList)
    sys.exit()

if len(contourStartEndStep) != 3 or sum(n < 0 for n in contourStartEndStep):
    eprint('Error in contourStartEndStep' + str(contourStartEndStep))
    sys.exit()

notDefContourStartEndStep = ( contourStartEndStep[0] == 0 and contourStartEndStep[1] == 0 )
notDefContourList = ( contourList[0] == 0 and contourList[1] == 0 )

if notDefContourStartEndStep and notDefContourList:
    contourStartEndStep[0] = 0
    contourStartEndStep[1] = NT - 1
    contourStartEndStep[2] = 1
    contourList = range(contourStartEndStep[0], contourStartEndStep[1], contourStartEndStep[2])
    print('contourStartEndStep is updated:',contourStartEndStep)
    print('contourList is updated:',contourList)
elif notDefContourList:
    contourList = range(contourStartEndStep[0], contourStartEndStep[1], contourStartEndStep[2])
    print('contourList is updated:',contourList)
# go to first step
animationScene1.GoToFirst()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
renderView1.ViewSize = [1200, 1024]

contour0 = Contour(Input=my_source)

# Properties modified on contour0 for estimation of clip
contour0.ContourBy = ['POINTS', 'f']
contour0.Isosurfaces = [0.5]
contour0.PointMergeMethod = 'Uniform Binning'
contour0Display = Show(contour0, renderView1, 'GeometryRepresentation')
# hide data in view
Hide(contour0, renderView1)
bounds = contour0.GetDataInformation().GetBounds()
print("bounds=", bounds)

# # calculation of bounds
boundsDomain = my_source.GetDataInformation().GetBounds()
print("bounds of the whole domain=",boundsDomain)
centerDomain = [(boundsDomain[0] + boundsDomain[1])/2, (boundsDomain[2] + boundsDomain[3])/2,(boundsDomain[4] + boundsDomain[5])/2]
print('center of the whole domain=',centerDomain)

# get layout
layout1 = GetLayout()

# show data in view
my_sourceDisplay = Show(my_source, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
my_sourceDisplay.Representation = 'Surface'
my_sourceDisplay.ColorArrayName = [None, '']
my_sourceDisplay.Opacity = 0.85
my_sourceDisplay.OSPRayScaleArray = 'f'
my_sourceDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
my_sourceDisplay.SelectOrientationVectors = 'None'
my_sourceDisplay.ScaleFactor = 3.0
my_sourceDisplay.SelectScaleArray = 'None'
my_sourceDisplay.GlyphType = 'Arrow'
my_sourceDisplay.GlyphTableIndexArray = 'None'
my_sourceDisplay.GaussianRadius = 0.15
my_sourceDisplay.SetScaleArray = ['POINTS', 'f']
my_sourceDisplay.ScaleTransferFunction = 'PiecewiseFunction'
my_sourceDisplay.OpacityArray = ['POINTS', 'f']
my_sourceDisplay.OpacityTransferFunction = 'PiecewiseFunction'
my_sourceDisplay.DataAxesGrid = 'GridAxesRepresentation'
my_sourceDisplay.PolarAxes = 'PolarAxesRepresentation'
my_sourceDisplay.ScalarOpacityUnitDistance = 0.48577906707444607

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
my_sourceDisplay.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# reset view to fit data
renderView1.ResetCamera()

# get the material library
materialLibrary1 = GetMaterialLibrary()

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Cell Data to Point Data'
cellDatatoPointData1 = CellDatatoPointData(Input=my_source)
cellDatatoPointData1.CellDataArraytoprocess = ['f', 'fs', 'l', 'p', 'residual_of_p', 'u.x']

# show data in view
cellDatatoPointData1Display = Show(cellDatatoPointData1, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
cellDatatoPointData1Display.Representation = 'Surface'
cellDatatoPointData1Display.ColorArrayName = [None, '']
cellDatatoPointData1Display.Opacity = 0.85
cellDatatoPointData1Display.OSPRayScaleArray = 'f'
cellDatatoPointData1Display.OSPRayScaleFunction = 'PiecewiseFunction'
cellDatatoPointData1Display.SelectOrientationVectors = 'None'
cellDatatoPointData1Display.ScaleFactor = 3.0
cellDatatoPointData1Display.SelectScaleArray = 'None'
cellDatatoPointData1Display.GlyphType = 'Arrow'
cellDatatoPointData1Display.GlyphTableIndexArray = 'None'
cellDatatoPointData1Display.GaussianRadius = 0.15
cellDatatoPointData1Display.SetScaleArray = ['POINTS', 'f']
cellDatatoPointData1Display.ScaleTransferFunction = 'PiecewiseFunction'
cellDatatoPointData1Display.OpacityArray = ['POINTS', 'f']
cellDatatoPointData1Display.OpacityTransferFunction = 'PiecewiseFunction'
cellDatatoPointData1Display.DataAxesGrid = 'GridAxesRepresentation'
cellDatatoPointData1Display.PolarAxes = 'PolarAxesRepresentation'
cellDatatoPointData1Display.ScalarOpacityUnitDistance = 0.48577906707444607

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
cellDatatoPointData1Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# hide data in view
Hide(my_source, renderView1)

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Threshold'
threshold1 = Threshold(Input=cellDatatoPointData1)
# Properties modified on threshold1
threshold1.Scalars = ['POINTS', 'fs']
threshold1.ThresholdRange = [0.0, 0.95]

# show data in view
threshold1Display = Show(threshold1, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
threshold1Display.Representation = 'Surface'
threshold1Display.ColorArrayName = [None, '']
threshold1Display.Opacity = 0.85
threshold1Display.OSPRayScaleArray = 'f'
threshold1Display.OSPRayScaleFunction = 'PiecewiseFunction'
threshold1Display.SelectOrientationVectors = 'None'
threshold1Display.ScaleFactor = 3.0
threshold1Display.SelectScaleArray = 'None'
threshold1Display.GlyphType = 'Arrow'
threshold1Display.GlyphTableIndexArray = 'None'
threshold1Display.GaussianRadius = 0.15
threshold1Display.SetScaleArray = ['POINTS', 'f']
threshold1Display.ScaleTransferFunction = 'PiecewiseFunction'
threshold1Display.OpacityArray = ['POINTS', 'f']
threshold1Display.OpacityTransferFunction = 'PiecewiseFunction'
threshold1Display.DataAxesGrid = 'GridAxesRepresentation'
threshold1Display.PolarAxes = 'PolarAxesRepresentation'
threshold1Display.ScalarOpacityUnitDistance = 0.30986217978916686

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
threshold1Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# hide data in view
Hide(cellDatatoPointData1, renderView1)

# create a new 'Clip'
clip1 = Clip(Input=threshold1)
clip1.ClipType = 'Plane'
clip1.HyperTreeGridClipper = 'Plane'
clip1.Scalars = ['POINTS', 'f']
clip1.Value = 0.5

# Properties modified on clip2.ClipType
clip1.ClipType.Origin = [bounds[0] - 0.1, 0.0, 0.0]
clip1.ClipType.Normal = [-1.0, 0.0, 0.0]

# init the 'Plane' selected for 'HyperTreeGridClipper'
clip1.HyperTreeGridClipper.Origin = [bounds[0] - 0.1, 0.0, 0.0]

# show data in view
clip1Display = Show(clip1, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
clip1Display.Representation = 'Surface'
clip1Display.ColorArrayName = [None, '']
clip1Display.Opacity = 0.85
clip1Display.OSPRayScaleArray = 'f'
clip1Display.OSPRayScaleFunction = 'PiecewiseFunction'
clip1Display.SelectOrientationVectors = 'None'
clip1Display.ScaleFactor = 1.5
clip1Display.SelectScaleArray = 'None'
clip1Display.GlyphType = 'Arrow'
clip1Display.GlyphTableIndexArray = 'None'
clip1Display.GaussianRadius = 0.075
clip1Display.SetScaleArray = ['POINTS', 'f']
clip1Display.ScaleTransferFunction = 'PiecewiseFunction'
clip1Display.OpacityArray = ['POINTS', 'f']
clip1Display.OpacityTransferFunction = 'PiecewiseFunction'
clip1Display.DataAxesGrid = 'GridAxesRepresentation'
clip1Display.PolarAxes = 'PolarAxesRepresentation'
clip1Display.ScalarOpacityUnitDistance = 0.19649440830229037

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
clip1Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# hide data in view
Hide(threshold1, renderView1)

# hide data in view
Hide(clip1, renderView1)

# update the view to ensure updated data information
renderView1.Update()


# create a new 'Clip'
clip2 = Clip(Input=clip1)
clip2.ClipType = 'Plane'
clip2.HyperTreeGridClipper = 'Plane'
clip2.Scalars = ['POINTS', 'f']
clip2.Value = 0.5

# Properties modified on clip2.ClipType
clip2.ClipType.Origin = [bounds[1] + 0.1, 0.0, 0.0]
# init the 'Plane' selected for 'HyperTreeGridClipper'
clip2.HyperTreeGridClipper.Origin = [bounds[1] + 0.1, 0.0, 0.0]
clip2.ClipType.Normal = [1.0, 0.0, 0.0]

# show data in view
clip2Display = Show(clip2, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
clip2Display.Representation = 'Surface'
clip2Display.ColorArrayName = [None, '']
clip2Display.Opacity = 0.85
clip2Display.OSPRayScaleArray = 'f'
clip2Display.OSPRayScaleFunction = 'PiecewiseFunction'
clip2Display.SelectOrientationVectors = 'None'
clip2Display.ScaleFactor = 1.3
clip2Display.SelectScaleArray = 'None'
clip2Display.GlyphType = 'Arrow'
clip2Display.GlyphTableIndexArray = 'None'
clip2Display.GaussianRadius = 0.065
clip2Display.SetScaleArray = ['POINTS', 'f']
clip2Display.ScaleTransferFunction = 'PiecewiseFunction'
clip2Display.OpacityArray = ['POINTS', 'f']
clip2Display.OpacityTransferFunction = 'PiecewiseFunction'
clip2Display.DataAxesGrid = 'GridAxesRepresentation'
clip2Display.PolarAxes = 'PolarAxesRepresentation'
clip2Display.ScalarOpacityUnitDistance = 0.17615700667391

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
clip2Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(clip1)

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=clip2.ClipType)

# toggle 3D widget visibility (only when running from the GUI)
Show3DWidgets(proxy=clip1.ClipType)

# set active source
SetActiveSource(clip2)

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=clip1.ClipType)

# toggle 3D widget visibility (only when running from the GUI)
Show3DWidgets(proxy=clip2.ClipType)

# Properties modified on clip2.ClipType
clip2.ClipType.Normal = [1.0, 0.0, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Delaunay 3D'
delaunay3D1 = Delaunay3D(Input=clip2)

# show data in view
delaunay3D1Display = Show(delaunay3D1, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
delaunay3D1Display.Representation = 'Surface'
delaunay3D1Display.ColorArrayName = [None, '']
delaunay3D1Display.Opacity = 0.85
delaunay3D1Display.OSPRayScaleArray = 'f'
delaunay3D1Display.OSPRayScaleFunction = 'PiecewiseFunction'
delaunay3D1Display.SelectOrientationVectors = 'None'
delaunay3D1Display.ScaleFactor = 0.2
delaunay3D1Display.SelectScaleArray = 'None'
delaunay3D1Display.GlyphType = 'Arrow'
delaunay3D1Display.GlyphTableIndexArray = 'None'
delaunay3D1Display.GaussianRadius = 0.01
delaunay3D1Display.SetScaleArray = ['POINTS', 'f']
delaunay3D1Display.ScaleTransferFunction = 'PiecewiseFunction'
delaunay3D1Display.OpacityArray = ['POINTS', 'f']
delaunay3D1Display.OpacityTransferFunction = 'PiecewiseFunction'
delaunay3D1Display.DataAxesGrid = 'GridAxesRepresentation'
delaunay3D1Display.PolarAxes = 'PolarAxesRepresentation'
delaunay3D1Display.ScalarOpacityUnitDistance = 0.0351274886866612

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
delaunay3D1Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# hide data in view
Hide(clip2, renderView1)

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Contour'
contour1 = Contour(Input=delaunay3D1)
contour1.ContourBy = ['POINTS', 'f']
contour1.Isosurfaces = [0.5]
contour1.PointMergeMethod = 'Uniform Binning'

# show data in view
contour1Display = Show(contour1, renderView1, 'GeometryRepresentation')

# get color transfer function/color map for 'f'
fLUT = GetColorTransferFunction('f')
fLUT.RGBPoints = [0.5, 0.23137254902, 0.298039215686, 0.752941176471, 0.50006103515625, 0.865, 0.865, 0.865, 0.5001220703125, 0.705882352941, 0.0156862745098, 0.149019607843]
fLUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
contour1Display.Representation = 'Surface'
contour1Display.ColorArrayName = ['POINTS', 'f']
contour1Display.LookupTable = fLUT
contour1Display.OSPRayScaleArray = 'f'
contour1Display.OSPRayScaleFunction = 'PiecewiseFunction'
contour1Display.SelectOrientationVectors = 'None'
contour1Display.ScaleFactor = 0.2
contour1Display.SelectScaleArray = 'f'
contour1Display.GlyphType = 'Arrow'
contour1Display.GlyphTableIndexArray = 'f'
contour1Display.GaussianRadius = 0.01
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

# show color bar/color legend
contour1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# get opacity transfer function/opacity map for 'f'
fPWF = GetOpacityTransferFunction('f')
fPWF.Points = [0.5, 0.0, 0.5, 0.0, 0.5001220703125, 1.0, 0.5, 0.0]
fPWF.ScalarRangeInitialized = 1

# hide data in view
Hide(delaunay3D1, renderView1)

# Hide orientation axes
renderView1.OrientationAxesVisibility = 0

# set scalar coloring
ColorBy(contour1Display, ('POINTS', 'u.x', 'Magnitude'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(fLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
contour1Display.RescaleTransferFunctionToDataRange(True, False)

# hide color bar/color legend
contour1Display.SetScalarBarVisibility(renderView1, False)

# get color transfer function/color map for 'ux'
uxLUT = GetColorTransferFunction('ux')
uxLUT.RGBPoints = [2.2554501288418073e-07, 0.23137254902, 0.298039215686, 0.752941176471, 1.0513483563140162, 0.865, 0.865, 0.865, 2.1026964870830196, 0.705882352941, 0.0156862745098, 0.149019607843]
uxLUT.ScalarRangeInitialized = 1.0

# get color legend/bar for uxLUT in view renderView1
uxLUTColorBar = GetScalarBar(uxLUT, renderView1)
uxLUTColorBar.WindowLocation = 'AnyLocation'
uxLUTColorBar.Title = 'u'
uxLUTColorBar.ComponentTitle = 'Magnitude'
uxLUTColorBar.AutomaticLabelFormat = 0
uxLUTColorBar.LabelFormat = '%-#6.5g'
uxLUTColorBar.RangeLabelFormat = '%6.5g'
uxLUTColorBar.ScalarBarLength = 0.3300000000000002

# get opacity transfer function/opacity map for 'ux'
uxPWF = GetOpacityTransferFunction('ux')
uxPWF.Points = [2.2554501288418073e-07, 0.0, 0.5, 0.0, 2.1026964870830196, 1.0, 0.5, 0.0]
uxPWF.ScalarRangeInitialized = 1
# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
uxLUT.ApplyPreset('jet', True)

# create a new 'Slice'
slice1 = Slice(Input=threshold1)
slice1.SliceType = 'Plane'
slice1.HyperTreeGridSlicer = 'Plane'
slice1.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Origin = [0, 0.0, 0.0]

# init the 'Plane' selected for 'HyperTreeGridSlicer'
slice1.HyperTreeGridSlicer.Origin = [0, 0.0, 0.0]

# Properties modified on slice1.SliceType
slice1.SliceType.Normal = [0.0, 0.0, 1.0]

# get display properties
slice1Display = GetDisplayProperties(slice1, view=renderView1)
# set scalar coloring
ColorBy(slice1Display, ('POINTS', 'u.x', 'Magnitude'))
# hide color bar/color legend
slice1Display.SetScalarBarVisibility(renderView1, False)

if not noPic:
    for i in contourList:
        print("in loop iteration:" + str(i) + " t:" + str(timesteps[i]))
        # Properties modified on animationScene1
        animationScene1.AnimationTime = timesteps[i]
        # Properties modified on timeKeeper1
        timeKeeper1.Time = timesteps[i]
        # create a new 'Contour'
        contour0 = Contour(Input=my_source)

        # Properties modified on contour0
        contour0.ContourBy = ['POINTS', 'f']
        contour0.Isosurfaces = [0.5]
        contour0.PointMergeMethod = 'Uniform Binning'
        contour0Display = Show(contour0, renderView1, 'GeometryRepresentation')
        # hide data in view
        Hide(contour0, renderView1)
    #     cylinder1 = FindSource('contour3')
        bounds = contour0.GetDataInformation().GetBounds()
        print("bounds=", bounds)
#         center = [(bounds[0] + bounds[1])/2, (bounds[2] + bounds[3])/2,(bounds[4] + bounds[5])/2]
        center = [bounds[0] + 0.0, (bounds[2] + bounds[3])/2,(bounds[4] + bounds[5])/2]
        print('center=', center)

        # Properties modified on clip2.ClipType
        clip1.ClipType.Origin = [bounds[0] - 0.1, 0.0, 0.0]
        clip1.ClipType.Normal = [-1.0, 0.0, 0.0]
        # init the 'Plane' selected for 'HyperTreeGridClipper'
        clip1.HyperTreeGridClipper.Origin = [bounds[0] - 0.1, 0.0, 0.0]

        # Properties modified on clip2.ClipType
        clip2.ClipType.Origin = [bounds[1] + 0.1, 0.0, 0.0]
        clip2.ClipType.Normal = [1.0, 0.0, 0.0]
        # init the 'Plane' selected for 'HyperTreeGridClipper'
        clip2.HyperTreeGridClipper.Origin = [bounds[1] + 0.1, 0.0, 0.0]


        renderView1.Update()
        renderView1.InteractionMode = '2D'
        renderView1.CameraPosition = [center[0] + 0.2, center[1], 6]
        renderView1.CameraFocalPoint = center
        renderView1.CameraParallelScale = 1
#         renderView1.CameraPosition = [14.85133751567109, -0.009534123838279692, 58.013718000153766]
#         renderView1.CameraFocalPoint = [14.85133751567109, -0.014148889638709456, 0.0012655406649076384]
#         renderView1.CameraParallelScale = 0.8897494071437293

#         renderView1.InteractionMode = '2D'
#         renderView1.CameraPosition = [14.85133751567109, -0.014148889638709456, 58.013718000153766]
#         renderView1.CameraFocalPoint = [14.85133751567109, -0.014148889638709456, 0.0012655406649076384]
#         renderView1.CameraParallelScale = 0.8897494071437293
        # hide color bar/color legend
        contour1Display.SetScalarBarVisibility(renderView1, False)

        # update the view to ensure updated data information
        renderView1.Update()
        fn = path + "/" + picName + str(i) +  '_t=' + str(timesteps[i]) +'.png'
        SaveScreenshot( fn, renderView1,
    #      ImageResolution=[2316, 2204],
    #     TransparentBackground=1,
    #     magnification=1,
        # PNG options
        CompressionLevel='2' )
        print('File=' + fn + ' generated succesfully')







































# # trace generated using paraview version 5.8.0
# #
# # To ensure correct image size when batch processing, please search
# # for and uncomment the line `# renderView*.ViewSize = [*,*]`
#
# #### import the simple module from the paraview
# from paraview.simple import *
# #### disable automatic camera reset on 'Show'
# paraview.simple._DisableFirstRenderCameraReset()
#
# import glob, os, sys
# import logging
# from sys import argv
# logging.basicConfig(format='%(message)s')
# log = logging.getLogger(__name__)
#
# def eprint(var):
#     log.warning(var)
#
# # Read from arguments
# #1 video Name
# #2 pvd file name (by defaults in the first file in a current directory)
# #3 pvtu is swithed off by default
#
# import argparse
# parser = argparse.ArgumentParser()
# parser.add_argument("-vidName", type=str, help="Provide the name for the outputed video, please",
#                     nargs='?', default='vid_')
# parser.add_argument("-picName", type=str, help="Provide the name for the outputed pictures, please",
#                     nargs='?', default='velocity_field_')
# parser.add_argument("-filenames", type=str, help="Provide the name of the outputing paraview files, please",
#                     nargs='?', default='*.pvd')
# parser.add_argument("-frameRate", type=int, help="Provide the frame rate for the video, please",
#                     nargs='?', default=10)
# parser.add_argument("-frameWindow", type=int, help="Provide the frame window (can be overwritten further), please",
#                     nargs='+', default=[0,0])
# parser.add_argument("-contourList", type=int, help="Provide the list of contour (can be overwritten further), please. Has higher priority",
#                     nargs='+', default=[0,0,0])
# parser.add_argument("-contourStartEndStep", type=int, help="Provide the start, end and step of output (can be overwritten further), please",
#                     nargs='+', default=[0,0,0])
# parser.add_argument("-lDomain", type=float, help="Provide the domain size, please",
#                     nargs='?', default=20.0)
# parser.add_argument("-viewSize", type=int, help="Provide the view size of the output pictures, please",
#                     nargs='+', default=[2152, 862])
# parser.add_argument("-noVideo", type=bool, help="Provide the no video Mode",
#                     nargs='?', default=False)
# parser.add_argument("-noPic", type=bool, help="Provide the no Picture Mode",
#                     nargs='?', default=False)
# parser.add_argument("-noData", type=bool, help="Provide the no Data Exporting Mode",
#                     nargs='?', default=False)
# # parser.add_argument('--foo', action='store_const', const=2, default=42)
# args = parser.parse_args()
# print(args)
# filenames = args.filenames
# vidName = args.vidName
# picName = args.picName
# frameRate = args.frameRate
# frameWindow = args.frameWindow
# contourList = args.contourList
# contourStartEndStep = args.contourStartEndStep
# lDomain = args.lDomain
# viewSize = args.viewSize
# noVideo = args.noVideo
# noPic = args.noPic
# noData = args.noData
#
# if len(frameWindow) != 2 or frameWindow[0]<0 or frameWindow[1] < 0:
#     eprint('Error in frameWindow' + frameWindow)
#     sys.exit()
#
# #Current PATH reading
# path = os.path.abspath(os.getcwd())
# eprint("Current PATH=" + path)
#
# if filenames[-5::] == '.pvtu':
# # Find files with *.pvtu extension
#     numbers = []
#     file = ""
#     for file in glob.glob(filenames):
#         numbers.append(int(filter(lambda x: x.isdigit(), file)))
#     file=file[0:-9]
#     numbers.sort()
#     N = len(numbers)
#     filenames = []
#     for i in numbers:
#         filenames.append('{}/{}{:04d}.pvtu'.format(path, file, i))
#     print(filenames)
#     fn = path +  '/' + filenames
#     my_source = XMLPartitionedUnstructuredGridReader(FileName = fn)
# elif filenames[-4::] == '.pvd':
# # create a new 'PVD Reader'
#     fn = glob.glob(filenames)
#     print('Found files:',fn)
#     fn = fn[0]
#     print('Read the first one:',fn)
#     my_source = PVDReader(FileName = path +  '/' + fn)
# else:
#     eprint('No pvd or pvtu files are provided')
#     sys.exit()
# #### disable automatic camera reset on 'Show'
# paraview.simple._DisableFirstRenderCameraReset()
#
# # get animation scene
# animationScene1 = GetAnimationScene()
#
# # get the time-keeper
# timeKeeper1 = GetTimeKeeper()
# timesteps = timeKeeper1.TimestepValues # 0, 0.1, 0.2 ...
# NT = len(timesteps)
# print("NT=", NT)
#
# if frameWindow[0] == 0 and frameWindow[1] == 0:
#     frameWindow[0] = 0
#     frameWindow[1] = NT-1
#     print('frameWindow is updated:',frameWindow)
#
# if sum(n < 0 for n in contourList):
#     eprint('Error in contourList:' + contourList)
#     sys.exit()
#
# if len(contourStartEndStep) != 3 or sum(n < 0 for n in contourStartEndStep):
#     eprint('Error in contourStartEndStep' + str(contourStartEndStep))
#     sys.exit()
#
# notDefContourStartEndStep = ( contourStartEndStep[0] == 0 and contourStartEndStep[1] == 0 )
# notDefContourList = ( contourList[0] == 0 and contourList[1] == 0 )
#
# if notDefContourStartEndStep and notDefContourList:
#     contourStartEndStep[0] = 0
#     contourStartEndStep[1] = NT - 1
#     contourStartEndStep[2] = 1
#     contourList = range(contourStartEndStep[0], contourStartEndStep[1], contourStartEndStep[2])
#     print('contourStartEndStep is updated:',contourStartEndStep)
#     print('contourList is updated:',contourList)
# elif notDefContourList:
#     contourList = range(contourStartEndStep[0], contourStartEndStep[1], contourStartEndStep[2])
#     print('contourList is updated:',contourList)
#
#
#
# # calculation of bounds
# boundsDomain = [0,30, -15, 15, -15, 15]  #my_source.GetDataInformation().GetBounds()
# print("bounds of the whole domain=",boundsDomain)
# centerDomain = [(boundsDomain[0] + boundsDomain[1])/2, (boundsDomain[2] + boundsDomain[3])/2,(boundsDomain[4] + boundsDomain[5])/2]
# print('center of the whole domain=',centerDomain)
#
#
# # get active view
# renderView1 = GetActiveViewOrCreate('RenderView')
# # uncomment following to set a specific view size
# renderView1.ViewSize = [1200, 1024]
#
# # get layout
# layout1 = GetLayout()
#
# # hide data in view
# Hide(my_source, renderView1)
#
# # create a new 'Threshold'
# threshold1 = Threshold(Input=my_source)
#
# # Properties modified on threshold1
# threshold1.Scalars = ['POINTS', 'fs']
# threshold1.ThresholdRange = [0.0, 0.99]
#
# # show data in view
# threshold1Display = Show(threshold1, renderView1, 'UnstructuredGridRepresentation')
#
# # trace defaults for the display properties.
# threshold1Display.Representation = 'Surface'
# threshold1Display.ColorArrayName = [None, '']
# threshold1Display.Opacity = 0.85
# threshold1Display.OSPRayScaleArray = 'fs'
# threshold1Display.OSPRayScaleFunction = 'PiecewiseFunction'
# threshold1Display.SelectOrientationVectors = 'None'
# threshold1Display.ScaleFactor = 3.0
# threshold1Display.SelectScaleArray = 'None'
# threshold1Display.GlyphType = 'Arrow'
# threshold1Display.GlyphTableIndexArray = 'None'
# threshold1Display.GaussianRadius = 0.15
# threshold1Display.SetScaleArray = ['POINTS', 'fs']
# threshold1Display.ScaleTransferFunction = 'PiecewiseFunction'
# threshold1Display.OpacityArray = ['POINTS', 'fs']
# threshold1Display.OpacityTransferFunction = 'PiecewiseFunction'
# threshold1Display.DataAxesGrid = 'GridAxesRepresentation'
# threshold1Display.PolarAxes = 'PolarAxesRepresentation'
# threshold1Display.ScalarOpacityUnitDistance = 0.3071941883753823
#
# # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
# threshold1Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]
#
# # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
# threshold1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 0.9855154608384954, 1.0, 0.5, 0.0]
#
# # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
# threshold1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 0.9855154608384954, 1.0, 0.5, 0.0]
#
# # reset view to fit data
# renderView1.ResetCamera()
#
# # get the material library
# materialLibrary1 = GetMaterialLibrary()
#
# # update the view to ensure updated data information
# renderView1.Update()
#
# # create a new 'Slice'
# slice1 = Slice(Input=threshold1)
# slice1.SliceType = 'Plane'
# slice1.HyperTreeGridSlicer = 'Plane'
# slice1.SliceOffsetValues = [0.0]
#
# # init the 'Plane' selected for 'SliceType'
# slice1.SliceType.Origin = centerDomain
#
# # init the 'Plane' selected for 'HyperTreeGridSlicer'
# slice1.HyperTreeGridSlicer.Origin = centerDomain
#
# # Properties modified on slice1.SliceType
# slice1.SliceType.Normal = [0.0, 0.0, 1.0]
#
# # show data in view
# slice1Display = Show(slice1, renderView1, 'GeometryRepresentation')
#
# # trace defaults for the display properties.
# slice1Display.Representation = 'Surface'
# slice1Display.ColorArrayName = [None, '']
# slice1Display.OSPRayScaleArray = 'fs'
# slice1Display.OSPRayScaleFunction = 'PiecewiseFunction'
# slice1Display.SelectOrientationVectors = 'None'
# slice1Display.ScaleFactor = 3.0
# slice1Display.SelectScaleArray = 'None'
# slice1Display.GlyphType = 'Arrow'
# slice1Display.GlyphTableIndexArray = 'None'
# slice1Display.GaussianRadius = 0.15
# slice1Display.SetScaleArray = ['POINTS', 'fs']
# slice1Display.ScaleTransferFunction = 'PiecewiseFunction'
# slice1Display.OpacityArray = ['POINTS', 'fs']
# slice1Display.OpacityTransferFunction = 'PiecewiseFunction'
# slice1Display.DataAxesGrid = 'GridAxesRepresentation'
# slice1Display.PolarAxes = 'PolarAxesRepresentation'
#
# # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
# slice1Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]
#
# # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
# slice1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 0.9499012345679013, 1.0, 0.5, 0.0]
#
# # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
# slice1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 0.9499012345679013, 1.0, 0.5, 0.0]
#
# # update the view to ensure updated data information
# renderView1.Update()
#
# # hide data in view
# Hide(threshold1, renderView1)
#
# # set scalar coloring
# ColorBy(slice1Display, ('CELLS', 'u.x', 'Magnitude'))
#
# # rescale color and/or opacity maps used to include current data range
# slice1Display.RescaleTransferFunctionToDataRange(True, False)
#
# # show color bar/color legend
# slice1Display.SetScalarBarVisibility(renderView1, True)
#
# # get color transfer function/color map for 'ux'
# uxLUT = GetColorTransferFunction('ux')
# uxLUT.RGBPoints = [2.2554501288418073e-07, 0.23137254902, 0.298039215686, 0.752941176471, 1.0513483563140162, 0.865, 0.865, 0.865, 2.1026964870830196, 0.705882352941, 0.0156862745098, 0.149019607843]
# uxLUT.ScalarRangeInitialized = 1.0
#
# # get opacity transfer function/opacity map for 'ux'
# uxPWF = GetOpacityTransferFunction('ux')
# uxPWF.Points = [2.2554501288418073e-07, 0.0, 0.5, 0.0, 2.1026964870830196, 1.0, 0.5, 0.0]
# uxPWF.ScalarRangeInitialized = 1
#
# # set active source
# SetActiveSource(threshold1)
#
# # toggle 3D widget visibility (only when running from the GUI)
# Hide3DWidgets(proxy=slice1.SliceType)
#
# # create a new 'Contour'
# contour1 = Contour(Input=threshold1)
#
# # Properties modified on contour1
# contour1.ContourBy = ['POINTS', 'f']
# contour1.Isosurfaces = [0.5]
# contour1.PointMergeMethod = 'Uniform Binning'
#
# # show data in view
# contour1Display = Show(contour1, renderView1, 'GeometryRepresentation')
#
# # get color transfer function/color map for 'f'
# fLUT = GetColorTransferFunction('f')
# fLUT.RGBPoints = [0.5, 0.23137254902, 0.298039215686, 0.752941176471, 0.50006103515625, 0.865, 0.865, 0.865, 0.5001220703125, 0.705882352941, 0.0156862745098, 0.149019607843]
# fLUT.ScalarRangeInitialized = 1.0
#
# # trace defaults for the display properties.
# contour1Display.Representation = 'Surface'
# contour1Display.ColorArrayName = ['POINTS', 'f']
# contour1Display.LookupTable = fLUT
# contour1Display.OSPRayScaleArray = 'f'
# contour1Display.OSPRayScaleFunction = 'PiecewiseFunction'
# contour1Display.SelectOrientationVectors = 'None'
# contour1Display.ScaleFactor = 0.3379618165202066
# contour1Display.SelectScaleArray = 'f'
# contour1Display.GlyphType = 'Arrow'
# contour1Display.GlyphTableIndexArray = 'f'
# contour1Display.GaussianRadius = 0.01689809082601033
# contour1Display.SetScaleArray = ['POINTS', 'f']
# contour1Display.ScaleTransferFunction = 'PiecewiseFunction'
# contour1Display.OpacityArray = ['POINTS', 'f']
# contour1Display.OpacityTransferFunction = 'PiecewiseFunction'
# contour1Display.DataAxesGrid = 'GridAxesRepresentation'
# contour1Display.PolarAxes = 'PolarAxesRepresentation'
#
# # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
# contour1Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]
#
# # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
# contour1Display.ScaleTransferFunction.Points = [0.5, 0.0, 0.5, 0.0, 0.5001220703125, 1.0, 0.5, 0.0]
#
# # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
# contour1Display.OpacityTransferFunction.Points = [0.5, 0.0, 0.5, 0.0, 0.5001220703125, 1.0, 0.5, 0.0]
#
# # show color bar/color legend
# contour1Display.SetScalarBarVisibility(renderView1, True)
#
# # update the view to ensure updated data information
# renderView1.Update()
#
# # get opacity transfer function/opacity map for 'f'
# fPWF = GetOpacityTransferFunction('f')
# fPWF.Points = [0.5, 0.0, 0.5, 0.0, 0.5001220703125, 1.0, 0.5, 0.0]
# fPWF.ScalarRangeInitialized = 1
#
# # set scalar coloring
# ColorBy(contour1Display, ('CELLS', 'u.x', 'Magnitude'))
#
# # Hide the scalar bar for this color map if no visible data is colored by it.
# HideScalarBarIfNotNeeded(fLUT, renderView1)
#
# # rescale color and/or opacity maps used to include current data range
# contour1Display.RescaleTransferFunctionToDataRange(True, False)
#
# # show color bar/color legend
# contour1Display.SetScalarBarVisibility(renderView1, True)
#
# # Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
# uxLUT.ApplyPreset('jet', True)
#
# # get color legend/bar for uxLUT in view renderView1
# uxLUTColorBar = GetScalarBar(uxLUT, renderView1)
# uxLUTColorBar.WindowLocation = 'AnyLocation'
# uxLUTColorBar.Title = 'u'
# uxLUTColorBar.ComponentTitle = 'Magnitude'
# uxLUTColorBar.AutomaticLabelFormat = 0
# uxLUTColorBar.ScalarBarLength = 0.3300000000000002
#
# # Properties modified on uxLUTColorBar
# uxLUTColorBar.LabelFormat = '%-#6.2f'
# uxLUTColorBar.RangeLabelFormat = '%6.2f'
#
# # Rescale transfer function
# uxLUT.RescaleTransferFunction(0.0, 3.2)
#
# # Rescale transfer function
# uxPWF.RescaleTransferFunction(0.0, 3.2)
#
# # Hide orientation axes
# renderView1.OrientationAxesVisibility = 0
#
# #change interaction mode for render view
# renderView1.InteractionMode = '2D'
#
# # create a new 'Plane'
# plane1 = Plane()
#
# # Properties modified on plane1
# plane1.Origin = [0.0, 0.0, 0.0]
# plane1.Point1 = [boundsDomain[1], 0.0, 0.0]
# plane1.Point2 = [0.0, 0.0, 1.0]
#
# # show data in view
# plane1Display = Show(plane1, renderView1, 'GeometryRepresentation')
#
# # trace defaults for the display properties.
# plane1Display.Representation = 'Surface'
# plane1Display.ColorArrayName = [None, '']
# plane1Display.OSPRayScaleArray = 'Normals'
# plane1Display.OSPRayScaleFunction = 'PiecewiseFunction'
# plane1Display.SelectOrientationVectors = 'None'
# plane1Display.ScaleFactor = 1.5
# plane1Display.SelectScaleArray = 'None'
# plane1Display.GlyphType = 'Arrow'
# plane1Display.GlyphTableIndexArray = 'None'
# plane1Display.GaussianRadius = 0.075
# plane1Display.SetScaleArray = ['POINTS', 'Normals']
# plane1Display.ScaleTransferFunction = 'PiecewiseFunction'
# plane1Display.OpacityArray = ['POINTS', 'Normals']
# plane1Display.OpacityTransferFunction = 'PiecewiseFunction'
# plane1Display.DataAxesGrid = 'GridAxesRepresentation'
# plane1Display.PolarAxes = 'PolarAxesRepresentation'
#
# # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
# plane1Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]
#
# # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
# plane1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 0, 1.0, 0.5, 0.0]
#
# # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
# plane1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 0, 1.0, 0.5, 0.0]
#
# # update the view to ensure updated data information
# renderView1.Update()
#
# # Properties modified on plane1
# plane1.XResolution = 10
# plane1.YResolution = 10
#
# # update the view to ensure updated data information
# renderView1.Update()
#
# #### saving camera placements for all active views
#
# # current camera placement for renderView1
# renderView1.InteractionMode = '2D'
# renderView1.CameraPosition = [16.12945558130605, -0.05378359910981157, 7.128400042325374]
# renderView1.CameraFocalPoint = [16.12945558130605, -0.05378359910981157, 0.0]
# renderView1.CameraParallelScale = 1.8449656920634208
#
# #### uncomment the following to render all views
# # RenderAllViews()
# # alternatively, if you want to write images, you can use SaveScreenshot(...).rite images, you can use SaveScreenshot(...).
#
#
# renderView1.Update()
# if not noPic:
#     for i in contourList:
#         print("in loop iteration:" + str(i))
#         # Properties modified on animationScene1
#         animationScene1.AnimationTime = timesteps[i]
#         # Properties modified on timeKeeper1
#         timeKeeper1.Time = timesteps[i]
#     #     Hide(line1, renderView1)
#     #     Hide(line2, renderView1)
#
#     #     cylinder1 = FindSource('contour3')
#         bounds = contour1.GetDataInformation().GetBounds()
#         print("bounds=", bounds)
# #         center = [(bounds[0] + bounds[1])/2, (bounds[2] + bounds[3])/2,(bounds[4] + bounds[5])/2]
#         center = [bounds[0] + 0.4, (bounds[2] + bounds[3])/2,(bounds[4] + bounds[5])/2]
#         print('center=', center)
#
#         renderView1.InteractionMode = '2D'
#         renderView1.CameraPosition = [center[0], center[1], 6]
#         renderView1.CameraFocalPoint = center
#         renderView1.CameraParallelScale = 0.7
# #         renderView1.CameraPosition = [14.85133751567109, -0.009534123838279692, 58.013718000153766]
# #         renderView1.CameraFocalPoint = [14.85133751567109, -0.014148889638709456, 0.0012655406649076384]
# #         renderView1.CameraParallelScale = 0.8897494071437293
#
# #         renderView1.InteractionMode = '2D'
# #         renderView1.CameraPosition = [14.85133751567109, -0.014148889638709456, 58.013718000153766]
# #         renderView1.CameraFocalPoint = [14.85133751567109, -0.014148889638709456, 0.0012655406649076384]
# #         renderView1.CameraParallelScale = 0.8897494071437293
#         # hide color bar/color legend
#         contour1Display.SetScalarBarVisibility(renderView1, False)
#
#         # update the view to ensure updated data information
#         renderView1.Update()
#         fn = path + "/" + picName + str(i) +  '_t=' + str(timesteps[i]) +'.png'
#         SaveScreenshot( fn, renderView1,
#     #      ImageResolution=[2316, 2204],
#     #     TransparentBackground=1,
#     #     magnification=1,
#         # PNG options
#         CompressionLevel='2' )
#         print('File=' + fn + ' generated succesfully')
#
# CreateLayout('Layout #2')