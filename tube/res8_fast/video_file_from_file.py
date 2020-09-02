# trace generated using paraview version 5.7.0
#
# To ensure correct image size when batch processing, please search
# for and uncomment the line `# renderView*.ViewSize = [*,*]`
# specify M step
#### import the simple module from the paraview
from paraview.simple import *
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
                    nargs='?', default='pic_')
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
                    nargs='?', default=30.0)
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

renderView1 = GetActiveViewOrCreate('RenderView')
print('contour1...')
# create a new 'Contour'
contour1 = Contour(Input=my_source)
contour1.ContourBy = ['POINTS', 'f']
contour1.Isosurfaces = [0.5]
contour1.PointMergeMethod = 'Uniform Binning'

# show data in view
contour1Display = Show(contour1, renderView1)
# contour1Display = Show(contour1, renderView1, 'GeometryRepresentation')

# get color transfer function/color map for 'f'
fLUT = GetColorTransferFunction('f')
fLUT.RGBPoints = [0.0, 0.23137254902, 0.298039215686, 0.752941176471, 0.5, 0.865, 0.865, 0.865, 1.0, 0.705882352941, 0.0156862745098, 0.149019607843]
fLUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
contour1Display.Representation = 'Surface'
contour1Display.ColorArrayName = ['POINTS', 'f']
contour1Display.LookupTable = fLUT
contour1Display.OSPRayScaleArray = 'f'
contour1Display.OSPRayScaleFunction = 'PiecewiseFunction'
contour1Display.SelectOrientationVectors = 'None'
contour1Display.ScaleFactor = 0.336715683006052
contour1Display.SelectScaleArray = 'f'
contour1Display.GlyphType = 'Arrow'
contour1Display.GlyphTableIndexArray = 'f'
contour1Display.GaussianRadius = 0.016835784150302596
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


# get opacity transfer function/opacity map for 'f'
fPWF = GetOpacityTransferFunction('f')
fPWF.ScalarRangeInitialized = 1

# turn off scalar coloring
ColorBy(contour1Display, None)

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(fLUT, renderView1)

# change solid color
contour1Display.AmbientColor = [0.00784313725490196, 0.00784313725490196, 1.0]
contour1Display.DiffuseColor = [0.00784313725490196, 0.00784313725490196, 1.0]

# Properties modified on contour1Display
contour1Display.Specular = 0.8

# Properties modified on contour1Display
contour1Display.Specular = 0.87

# Properties modified on contour1Display
contour1Display.Specular = 0.89

# Properties modified on contour1Display
contour1Display.Specular = 1.0

print('contour1Display parameters are set')
# update the view to ensure updated data information
renderView1.Update()

print('contour1 created succesfully')
print('cylinder fs...')

# create a new 'Contour for solid'
contour2 = Contour(Input=my_source)
contour2.ContourBy = ['POINTS', 'fs']
contour2.Isosurfaces = [0.5]
contour2.PointMergeMethod = 'Uniform Binning'

# show data in view
contour2Display = Show(contour2, renderView1)
# contour2Display = Show(contour2, renderView1, 'GeometryRepresentation')

# get color transfer function/color map for 'fs'
fsLUT = GetColorTransferFunction('fs')
fsLUT.RGBPoints = [0.5, 0.23137254902, 0.298039215686, 0.752941176471, 0.50006103515625, 0.865, 0.865, 0.865, 0.5001220703125, 0.705882352941, 0.0156862745098, 0.149019607843]
fsLUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
contour2Display.Representation = 'Surface'
contour2Display.ColorArrayName = ['POINTS', 'fs']
contour2Display.LookupTable = fsLUT
contour2Display.OSPRayScaleArray = 'fs'
contour2Display.OSPRayScaleFunction = 'PiecewiseFunction'
contour2Display.SelectOrientationVectors = 'None'
contour2Display.ScaleFactor = 3.0
contour2Display.SelectScaleArray = 'fs'
contour2Display.GlyphType = 'Arrow'
contour2Display.GlyphTableIndexArray = 'fs'
contour2Display.GaussianRadius = 0.15
contour2Display.SetScaleArray = ['POINTS', 'fs']
contour2Display.ScaleTransferFunction = 'PiecewiseFunction'
contour2Display.OpacityArray = ['POINTS', 'fs']
contour2Display.OpacityTransferFunction = 'PiecewiseFunction'
contour2Display.DataAxesGrid = 'GridAxesRepresentation'
contour2Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
contour2Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
contour2Display.ScaleTransferFunction.Points = [0.5, 0.0, 0.5, 0.0, 0.5001220703125, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
contour2Display.OpacityTransferFunction.Points = [0.5, 0.0, 0.5, 0.0, 0.5001220703125, 1.0, 0.5, 0.0]

# show color bar/color legend
contour2Display.SetScalarBarVisibility(renderView1, False)

# update the view to ensure updated data information
renderView1.Update()

# get opacity transfer function/opacity map for 'fs'
fsPWF = GetOpacityTransferFunction('fs')
fsPWF.Points = [0.5, 0.0, 0.5, 0.0, 0.5001220703125, 1.0, 0.5, 0.0]
fsPWF.ScalarRangeInitialized = 1
# Properties modified on contour2Display
contour2Display.Opacity = 0.15

# Properties modified on contour2Display
contour2Display.Opacity = 0.21

# Properties modified on contour2Display
contour2Display.Luminosity = 44.0

# Properties modified on contour2Display
contour2Display.Specular = 0.32
print('contour2 created succesfully')

# # create a new 'Cylinder'
# cylinder1 = Cylinder()
# print('cylinder created...')
# # Properties modified on cylinder1
# cylinder1.Resolution = 60
# cylinder1.Height = 30.0
# print('cylinder Resolution and Height...')
# # show data in view
# cylinder1Display = Show(cylinder1, renderView1)
# # cylinder1Display = Show(cylinder1, renderView1, 'GeometryRepresentation')
# print('cylinder show...')
# # trace defaults for the display properties.
# cylinder1Display.Representation = 'Surface'
# cylinder1Display.ColorArrayName = [None, '']
# cylinder1Display.OSPRayScaleArray = 'Normals'
# cylinder1Display.OSPRayScaleFunction = 'PiecewiseFunction'
# cylinder1Display.SelectOrientationVectors = 'None'
# cylinder1Display.ScaleFactor = 3.0
# cylinder1Display.SelectScaleArray = 'None'
# cylinder1Display.GlyphType = 'Arrow'
# cylinder1Display.GlyphTableIndexArray = 'None'
# cylinder1Display.GaussianRadius = 0.15
# cylinder1Display.SetScaleArray = ['POINTS', 'Normals']
# cylinder1Display.ScaleTransferFunction = 'PiecewiseFunction'
# cylinder1Display.OpacityArray = ['POINTS', 'Normals']
# cylinder1Display.OpacityTransferFunction = 'PiecewiseFunction'
# cylinder1Display.DataAxesGrid = 'GridAxesRepresentation'
# cylinder1Display.PolarAxes = 'PolarAxesRepresentation'
#
# # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
# cylinder1Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]
#
# # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
# cylinder1Display.ScaleTransferFunction.Points = [-0.9983690977096558, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]
#
# # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
# cylinder1Display.OpacityTransferFunction.Points = [-0.9983690977096558, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]
# print('cylinder Display parameters are set...')
# # update the view to ensure updated data information
# renderView1.Update()
# print('cylinder Display parameters are updated...')
#
# # Properties modified on cylinder1Display
# cylinder1Display.Position = [15.0, 0.0, 0.0]
#
# # Properties modified on cylinder1Display.DataAxesGrid
# cylinder1Display.DataAxesGrid.Position = [15.0, 0.0, 0.0]
#
# # Properties modified on cylinder1Display.PolarAxes
# cylinder1Display.PolarAxes.Translation = [15.0, 0.0, 0.0]
#
# # Properties modified on cylinder1Display
# cylinder1Display.Orientation = [0.0, 0.0, 90.0]
#
# # Properties modified on cylinder1Display.PolarAxes
# cylinder1Display.PolarAxes.Orientation = [0.0, 0.0, 90.0]
#
# # Properties modified on cylinder1Display
# cylinder1Display.Specular = 0.82
#
# # Properties modified on cylinder1Display
# cylinder1Display.Specular = 0.49
#
# # Properties modified on cylinder1Display
# cylinder1Display.Ambient = 0.14
#
# # Properties modified on cylinder1Display
# cylinder1Display.Opacity = 0.19
# print('cylinder display, location, bright etc. parameters are updated...')
# set active source
SetActiveSource(contour1)
SetActiveSource(contour2)
# SetActiveSource(cylinder1)
print('contour and cylinder are active')

# get animation scene
animationScene1 = GetAnimationScene()
print('Get animationScene1')
# # get the time-keeper
# timeKeeper1 = GetTimeKeeper()


animationScene1.GoToFirst()
timeKeeper1.Time = timesteps[0]
bounds = contour1.GetDataInformation().GetBounds()
print("boundsF=",bounds)
centerF = [(bounds[0] + bounds[1])/2, (bounds[2] + bounds[3])/2,(bounds[4] + bounds[5])/2]
print('centerFirst=',centerF)

animationScene1.GoToLast()
timeKeeper1.Time = timesteps[-1]
bounds = contour1.GetDataInformation().GetBounds()
print("boundsLast=",bounds)
centerL = [(bounds[0] + bounds[1])/2, (bounds[2] + bounds[3])/2,(bounds[4] + bounds[5])/2]
print('centerLast=',centerL)

center = [0.5*(centerF[i] + centerL[i]) for i in range(3)]

renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [center[0], center[1], 6]
renderView1.CameraFocalPoint = center
renderView1.CameraParallelScale = 0.7

# current camera placement for renderView1
# renderView1.InteractionMode = '2D'
# renderView1.CameraPosition = [lDomain/2., 0, 6]
# renderView1.CameraFocalPoint = [lDomain/2., 0, 0]
# renderView1.CameraParallelScale = 2

# update the view to ensure updated data information
renderView1.Update()
# save animation
if vidName[-4::] != '.avi':
    fn = path + '/' + vidName + '.avi'
else:
     fn = path + '/' + vidName
if not noVideo:
    SaveAnimation(fn, renderView1, ImageResolution=[2048, 512], FrameRate=frameRate, FrameWindow=frameWindow)
    print("File:" + fn + " successfully saved")
# hide data in view
# Hide(cylinder1, renderView1)
Hide(contour1, renderView1)
################################ PICTURE ################################################################

paraview.simple._DisableFirstRenderCameraReset()
print("slice1...")
# create a new 'Slice'
slice1 = Slice(Input=my_source)
slice1.SliceType = 'Plane'
# slice1.HyperTreeGridSlicer = 'Plane'
slice1.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Origin = [15.0, 0.0, 0.0]

# init the 'Plane' selected for 'HyperTreeGridSlicer'
# slice1.HyperTreeGridSlicer.Origin = [15.0, 0.0, 0.0]

# Properties modified on slice1.SliceType
slice1.SliceType.Normal = [0.0, 0.0, 1.0]

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
renderView1.ViewSize = viewSize

# get layout
layout1 = GetLayout()

# show data in view
slice1Display = Show(slice1, renderView1)

# trace defaults for the display properties.
slice1Display.Representation = 'Surface'
slice1Display.ColorArrayName = [None, '']
slice1Display.OSPRayScaleArray = 'f'
slice1Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice1Display.SelectOrientationVectors = 'None'
slice1Display.ScaleFactor = 3.0
slice1Display.SelectScaleArray = 'None'
slice1Display.GlyphType = 'Arrow'
slice1Display.GlyphTableIndexArray = 'None'
slice1Display.GaussianRadius = 0.15
slice1Display.SetScaleArray = ['POINTS', 'f']
slice1Display.ScaleTransferFunction = 'PiecewiseFunction'
slice1Display.OpacityArray = ['POINTS', 'f']
slice1Display.OpacityTransferFunction = 'PiecewiseFunction'
slice1Display.DataAxesGrid = 'GridAxesRepresentation'
slice1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
slice1Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]
print("slice1 finished")
# # reset view to fit data
# renderView1.ResetCamera()
#
# #changing interaction mode based on data extents
# renderView1.CameraPosition = [15.0, 0.0, 10000.0]
# renderView1.CameraFocalPoint = [15.0, 0.0, 0.0]

# get the material library
materialLibrary1 = GetMaterialLibrary()
print("contour3...")
# create a new 'Contour'
contour3 = Contour(Input=slice1)
contour3.ContourBy = ['POINTS', 'f']
contour3.Isosurfaces = [0.5]
contour3.PointMergeMethod = 'Uniform Binning'

# show data in view
contour3Display = Show(contour3, renderView1)

# get color transfer function/color map for 'f'
fLUT = GetColorTransferFunction('f')
fLUT.RGBPoints = [0.0, 0.23137254902, 0.298039215686, 0.752941176471, 0.5, 0.865, 0.865, 0.865, 1.0, 0.705882352941, 0.0156862745098, 0.149019607843]
fLUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
contour3Display.Representation = 'Surface'
contour3Display.ColorArrayName = ['POINTS', 'f']
contour3Display.LookupTable = fLUT
contour3Display.OSPRayScaleArray = 'f'
contour3Display.OSPRayScaleFunction = 'PiecewiseFunction'
contour3Display.SelectOrientationVectors = 'None'
contour3Display.ScaleFactor = 0.3329599454085324
contour3Display.SelectScaleArray = 'f'
contour3Display.GlyphType = 'Arrow'
contour3Display.GlyphTableIndexArray = 'f'
contour3Display.GaussianRadius = 0.01664799727042662
contour3Display.SetScaleArray = ['POINTS', 'f']
contour3Display.ScaleTransferFunction = 'PiecewiseFunction'
contour3Display.OpacityArray = ['POINTS', 'f']
contour3Display.OpacityTransferFunction = 'PiecewiseFunction'
contour3Display.DataAxesGrid = 'GridAxesRepresentation'
contour3Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
contour3Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
contour3Display.ScaleTransferFunction.Points = [0.5, 0.0, 0.5, 0.0, 0.5001220703125, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
contour3Display.OpacityTransferFunction.Points = [0.5, 0.0, 0.5, 0.0, 0.5001220703125, 1.0, 0.5, 0.0]

# hide data in view
Hide(slice1, renderView1)

# show color bar/color legend
contour3Display.SetScalarBarVisibility(renderView1, False)

# update the view to ensure updated data information
renderView1.Update()

# get opacity transfer function/opacity map for 'f'
fPWF = GetOpacityTransferFunction('f')
fPWF.ScalarRangeInitialized = 1

# turn off scalar coloring
ColorBy(contour3Display, None)

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(fLUT, renderView1)

# change solid color
contour3Display.AmbientColor = [0.0, 0.0, 0.0]
contour3Display.DiffuseColor = [0.0, 0.0, 0.0]
print("contoured")
# Hide orientation axes
renderView1.OrientationAxesVisibility = 0

# reset view to fit data bounds
# renderView1.ResetCamera(8.06708301404, 11.3966824681, -0.447185845316, 0.445647527971, 0.0, 0.0)

# set active source
SetActiveSource(my_source)
print("SetActiveSource(my_source)")
print("line1...")
# create a new 'Line'
line1 = Line()

# Properties modified on line1
line1.Point1 = [0.0, 0.5, 0.0]
line1.Point2 = [lDomain, 0.5, 0.0]

# show data in view
line1Display = Show(line1, renderView1)
print("line1 show")
# trace defaults for the display properties.
line1Display.Representation = 'Surface'
line1Display.ColorArrayName = [None, '']
line1Display.OSPRayScaleArray = 'Texture Coordinates'
line1Display.OSPRayScaleFunction = 'PiecewiseFunction'
line1Display.SelectOrientationVectors = 'None'
line1Display.ScaleFactor = 1.5
line1Display.SelectScaleArray = 'None'
line1Display.GlyphType = 'Arrow'
line1Display.GlyphTableIndexArray = 'None'
line1Display.GaussianRadius = 0.075
line1Display.SetScaleArray = ['POINTS', 'Texture Coordinates']
line1Display.ScaleTransferFunction = 'PiecewiseFunction'
line1Display.OpacityArray = ['POINTS', 'Texture Coordinates']
line1Display.OpacityTransferFunction = 'PiecewiseFunction'
line1Display.DataAxesGrid = 'GridAxesRepresentation'
line1Display.PolarAxes = 'PolarAxesRepresentation'

# change solid color
line1Display.AmbientColor = [0.0, 0.0, 0.0]
line1Display.DiffuseColor = [0.0, 0.0, 0.0]

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
line1Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]
print("line1 is finished")
# update the view to ensure updated data information
# renderView1.Update()
print("line2...")
# create a new 'Line'
line2 = Line()

# set active source
SetActiveSource(line1)
print("line1 is set active")
# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=line2)
print("Hide3DWidgets for line2")
# set active source
SetActiveSource(line2)
print("line2 is set active")
# Properties modified on line2
line2.Point1 = [0.0, -0.5, 0.0]
line2.Point2 = [lDomain, -0.5, 0.0]

# show data in view
line2Display = Show(line2, renderView1)

# trace defaults for the display properties.
line2Display.Representation = 'Surface'
line2Display.ColorArrayName = [None, '']
line2Display.OSPRayScaleArray = 'Texture Coordinates'
line2Display.OSPRayScaleFunction = 'PiecewiseFunction'
line2Display.SelectOrientationVectors = 'None'
line2Display.ScaleFactor = 1.5
line2Display.SelectScaleArray = 'None'
line2Display.GlyphType = 'Arrow'
line2Display.GlyphTableIndexArray = 'None'
line2Display.GaussianRadius = 0.075
line2Display.SetScaleArray = ['POINTS', 'Texture Coordinates']
line2Display.ScaleTransferFunction = 'PiecewiseFunction'
line2Display.OpacityArray = ['POINTS', 'Texture Coordinates']
line2Display.OpacityTransferFunction = 'PiecewiseFunction'
line2Display.DataAxesGrid = 'GridAxesRepresentation'
line2Display.PolarAxes = 'PolarAxesRepresentation'

# change solid color
line2Display.AmbientColor = [0.0, 0.0, 0.0]
line2Display.DiffuseColor = [0.0, 0.0, 0.0]

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
line2Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# current camera placement for renderView1
# renderView1.InteractionMode = '2D'
# renderView1.CameraPosition = [10, 0, 6]
# renderView1.CameraFocalPoint = [10, 0, 0]
# renderView1.CameraParallelScale = 0.12

# update the view to ensure updated data information
renderView1.Update()
print("line2 is finished")
if not noPic:
    for i in contourList:
        print("in loop iteration:" + str(i))
        # Properties modified on animationScene1
        animationScene1.AnimationTime = timesteps[i]
        # Properties modified on timeKeeper1
        timeKeeper1.Time = timesteps[i]
    #     Hide(line1, renderView1)
    #     Hide(line2, renderView1)

    #     cylinder1 = FindSource('contour3')
        bounds = contour3.GetDataInformation().GetBounds()
        print("bounds=",bounds)
        center = [(bounds[0] + bounds[1])/2, (bounds[2] + bounds[3])/2,(bounds[4] + bounds[5])/2]
        print('center=',center)

        renderView1.InteractionMode = '2D'
        renderView1.CameraPosition = [center[0], center[1], 6]
        renderView1.CameraFocalPoint = center
        renderView1.CameraParallelScale = 0.7

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

CreateLayout('Layout #2')

# set active view
SetActiveView(None)

# Create a new 'SpreadSheet View'
spreadSheetView1 = CreateView('SpreadSheetView')
spreadSheetView1.ColumnToSort = ''
spreadSheetView1.BlockSize = 1024L
# uncomment following to set a specific view size
# spreadSheetView1.ViewSize = [400, 400]

# show data in view
contour3Display_1 = Show(contour3, spreadSheetView1)
# contour3Display_1 = Show(contour3, spreadSheetView1, 'SpreadSheetRepresentation')

# get layout
layout2 = GetLayoutByName("Layout #2")

# assign view to a particular cell in the layout
# AssignViewToLayout(view=spreadSheetView1, layout=layout2, hint=0)
if not noData:
    for i in contourList:
        # Properties modified on animationScene1
        animationScene1.AnimationTime = timesteps[i]
        # Properties modified on timeKeeper1
        timeKeeper1.Time = timesteps[i]
        # export view
        fn = path + '/' + 'data_i=' + str(i) + '_t=' + str(timesteps[i]) + '.csv'
        ExportView(fn, view=spreadSheetView1)
        print('File=' + fn + ' generated succesfully')


