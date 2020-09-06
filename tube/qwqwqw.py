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
logging.basicConfig(format='%(message)s')
log = logging.getLogger(__name__)

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
                    nargs='?', default='lambda2_')
parser.add_argument("-picName", type=str, help="Provide the name for the outputed pictures, please",
                    nargs='?', default='lambda2_')
parser.add_argument("-filenames", type=str, help="Provide the name of the outputing paraview files, please",
                    nargs='?', default='*.pvd')
parser.add_argument("-frameRate", type=int, help="Provide the frame rate for the video, please",
                    nargs='?', default=10)
parser.add_argument("-frameWindow", type=int, help="Provide the frame window (can be overwritten further), please",
                    nargs='+', default=[0,0])
parser.add_argument("-timeList", type=int, help="Provide the list of step, please",
                    nargs='+', default=[0,-1])

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
viewSize = args.viewSize
noVideo = args.noVideo
noPic = args.noPic
noData = args.noData
noLambda2 = args.noLambda2
timeList = args.timeList

if len(frameWindow) != 2 or frameWindow[0]<0 or frameWindow[1] < 0:
    eprint('Error in frameWindow' + frameWindow)
    sys.exit()
# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [3108, 1168]
renderView1.InteractionMode = '2D'
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.StereoType = 'Crystal Eyes'
renderView1.BackEnd = 'OSPRay raycaster'
# get the material library
materialLibrary1 = GetMaterialLibrary()
renderView1.OSPRayMaterialLibrary = materialLibrary1
renderView1.CameraViewUp = [0.2, 1, 0]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 1.13

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.AssignView(0, renderView1)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
#Current PATH reading
path = "/Users/weugene/basilisk/work/tube/res22_adaptffsu/"
# path = os.path.abspath(os.getcwd())
eprint("Current PATH=" + path)
# create a new 'PVD Reader'
# my_source = PVDReader(FileName=path+'dump2pvd_compressed.pvd')
my_source = GetActiveSource()
my_source.CellArrays = ['p', 'fs', 'f', 'l', 'residual_of_p', 'l2', 'omega', 'u.x']

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
timeList = range(0,NT)
print("NT=", NT, " timeList=", timeList)

if frameWindow[0] == 0 and frameWindow[1] == 0:
    frameWindow[0] = 0
    frameWindow[1] = NT-1
    print('frameWindow is updated:',frameWindow)




###-----------------GENERATION of Cylinder-----------------------------------
### it is timeless therefore it is outside of the loop
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
transform1Display.PolarAxes = 'PolarAxesRepresentation'
transform1Display.ScalarOpacityUnitDistance = 6.650076732513133

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
transform1Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
transform1Display.ScaleTransferFunction.Points = [-2.220446049250313e-16, 0.0, 0.5, 0.0, 2.220446049250313e-16, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
transform1Display.OpacityTransferFunction.Points = [-2.220446049250313e-16, 0.0, 0.5, 0.0, 2.220446049250313e-16, 1.0, 0.5, 0.0]

SetActiveSource(transform1)
# Show()
Render()
# SetActiveSource(transform1)
# renderView1.Update()

###-----------------GENERATION of a CLIP and Convert to Point data-----------------------------------
if not noPic:
    for i in timeList:
        print("in loop iteration:" + str(i))
        # Properties modified on animationScene1
        animationScene1.AnimationTime = timesteps[i]
        # Properties modified on timeKeeper1
        timeKeeper1.Time = timesteps[i]
        contour0 = Contour(Input=my_source)

        # Properties modified on contour0
        contour0.ContourBy = ['POINTS', 'f']
        contour0.Isosurfaces = [0.5]
        contour0.PointMergeMethod = 'Uniform Binning'
        contour0Display = Show(contour0, renderView1, 'GeometryRepresentation')
        # hide data in view
        Hide(contour0, renderView1)
        bounds = contour0.GetDataInformation().GetBounds()
        print("bounds=",bounds)
        center = [(bounds[0] + bounds[1])/2, (bounds[2] + bounds[3])/2,(bounds[4] + bounds[5])/2]
        print('center=',center)
        len_min = max(bounds[0] - 3.5,0.5)
        len_max = min(bounds[1] + 1.2,29)
#         len_min = bounds[0] - 0.5   # correct it!!!
#         len_max = bounds[1] + 0.5   # correct it!!!
        length = len_max - len_min
        print('len_min=',len_min,' len_max=',len_max,' len=',length)

# ***************** CLIP A BOX ****************************
        # create a new 'Clip'
        clip1 = Clip(Input=my_source)
        clip1.ClipType = 'Box'
        clip1.HyperTreeGridClipper = 'Plane'
        clip1.Scalars = ['POINTS', 'f']
        clip1.Value = 0.5

        # init the 'Box' selected for 'ClipType'
        clip1.ClipType.Position = [len_min - 0.2, -0.55, -0.55]
        clip1.ClipType.Length = [length + 0.4, 1.1, 1.1]

        # init the 'Plane' selected for 'HyperTreeGridClipper'
        clip1.HyperTreeGridClipper.Origin = [15.0, 0.0, 0.0] #???

# ***************** CELL DATA TO POINT DATA ****************************
        # create a new 'Cell Data to Point Data'
        cellDatatoPointData1 = CellDatatoPointData(Input=clip1)
#         cellDatatoPointData1.CellDataArraytoprocess = ['f', 'fs', 'l', 'l2', 'omega', 'p', 'u.x']

# ***************** RESAMPLE TO IMAGE ****************************
        # create a new 'Resample To Image'
        resampleToImage1 = ResampleToImage(Input=cellDatatoPointData1)
        resampleToImage1.UseInputBounds = 0
        resampleToImage1.SamplingDimensions = [800, 100, 100]
        resampleToImage1.SamplingBounds = [len_min, len_max, -0.5, 0.5, -0.5, 0.5]

# ***************** CONTOUR1 for VOLUME FRACTION f ****************************
        # create a new 'Contour'
        contour1 = Contour(Input=resampleToImage1)
        contour1.ContourBy = ['POINTS', 'f']
        contour1.Isosurfaces = [0.5]
        contour1.PointMergeMethod = 'Uniform Binning'

# ***************** CONTOUR2 for LAMBDA2 l2 ****************************
        # create a new 'Contour'
        contour2 = Contour(Input=resampleToImage1)
        contour2.ContourBy = ['POINTS', 'l2']
        contour2.Isosurfaces = [-1.0, -2.0]
        contour2.PointMergeMethod = 'Uniform Binning'

# ----------------------------------------------------------------
        # setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

        # get color transfer function/color map for 'ux'
        uxLUT = GetColorTransferFunction('ux')
#         uxLUT.RGBPoints = [0.3730543491279064, 0.0, 0.0, 0.5625, 0.5106463634876334, 0.0, 0.0, 1.0, 0.8251430154745502, 0.0, 1.0, 1.0, 0.9823910318856668, 0.5, 1.0, 0.5, 1.139639048296783, 1.0, 1.0, 0.0, 1.4541357002836997, 1.0, 0.0, 0.0, 1.6113837166948162, 0.5, 0.0, 0.0]
#         uxLUT.ColorSpace = 'RGB'
        uxLUT.ScalarRangeInitialized = 1.0
        uxLUT.ApplyPreset('jet', True)

        # get color legend/bar for uxLUT in view renderView1
        uxLUTColorBar = GetScalarBar(uxLUT, renderView1)
        uxLUTColorBar.Orientation = 'Horizontal'
        uxLUTColorBar.WindowLocation = 'AnyLocation'
        uxLUTColorBar.Position = [0.4, 0.1]
        uxLUTColorBar.Title = 'u'
        uxLUTColorBar.ComponentTitle = 'Magnitude'
        uxLUTColorBar.LabelFormat = '%-#6.2g'
        uxLUTColorBar.RangeLabelFormat = '%6.2g'
        uxLUTColorBar.ScalarBarLength = 0.3


        # set color bar visibility
        uxLUTColorBar.Visibility = 1

        # get color transfer function/color map for 'l2'
        l2LUT = GetColorTransferFunction('l2')
        l2LUT.RGBPoints = [-2.0, 0.0, 1.0, 1.0, -1.55, 0.0, 0.0, 1.0, -1.5, 0.0, 0.0, 0.501960784314, -1.4499999999999997, 1.0, 0.0, 0.0, -1.0, 1.0, 1.0, 0.0]
        l2LUT.ColorSpace = 'RGB'
        l2LUT.ScalarRangeInitialized = 1.0

        # show data from contour1
        contour1Display = Show(contour1, renderView1, 'GeometryRepresentation')
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

        # setup the color legend parameters for each legend in this view

        # show color legend
        contour1Display.SetScalarBarVisibility(renderView1, True)

        # ----------------------------------------------------------------
        # setup color maps and opacity mapes used in the visualization
        # note: the Get..() functions create a new object, if needed
        # ----------------------------------------------------------------
        info = contour1.GetDataInformation().DataInformation
        arrayInfo = info.GetArrayInformation("u.x", vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS)
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
#         uxLUTColorBar.CustomLabels = np.linspace(0, range_max, 5)
        uxLUTColorBar.CustomLabels = np.arange(0, range_max, 0.5)
        print("CustomLbels=",uxLUTColorBar.CustomLabels )

        # get opacity transfer function/opacity map for 'ux'
        uxPWF = GetOpacityTransferFunction('ux')
        uxPWF.Points = [0.3730543491279064, 0.0, 0.5, 0.0, 1.6113837166948162, 1.0, 0.5, 0.0]
        uxPWF.ScalarRangeInitialized = 1

        # get opacity transfer function/opacity map for 'l2'
        l2PWF = GetOpacityTransferFunction('l2')
        l2PWF.Points = [-2.0, 0.0, 0.5, 0.0, -1.0, 1.0, 0.5, 0.0]
        l2PWF.ScalarRangeInitialized = 1
        # ----------------------------------------------------------------
        # finally, restore active source
        SetActiveSource(contour1)
        # ----------------------------------------------------------------

        renderView1.InteractionMode = '2D'
#         renderView1.CameraPosition = [center[0], center[1], 6]
#         renderView1.CameraFocalPoint = center
        renderView1.CenterOfRotation = [0.5*(len_min + len_max), 0.0, 0.0]
        renderView1.CameraFocalPoint = [0.5*(len_min + len_max) , 0, 0]
#         renderView1.CameraFocalPoint = [0.5*(len_min + len_max) , -0.25, -1.3]
        renderView1.CameraPosition = [0.5*(len_min + len_max) - 2, 0.6, 4.5]
        renderView1.CameraViewUp = [0.2, 1, 0]
        renderView1.CameraFocalDisk = 1.0
#         renderView1.CameraParallelScale = 1.8

        # update the view to ensure updated data information
        renderView1.Update()

#****************** CONTOUR1(f) AND U MAGNITUDE ********************
        # show data from contour1
        fn = path + "/" + picName + str(i) +  '_t=' + str(timesteps[i]) +'ux.png'
        SaveScreenshot( fn, renderView1,
        #      ImageResolution=[2316, 2204],
            TransparentBackground=0,
            CompressionLevel='2' )
        print('File=' + fn + ' generated succesfully')
#****************** CONTOUR1(f) AND CONTOUR2 (lambda2) ********************
        # set color bar visibility
        uxLUTColorBar.Visibility = 0

        # show data from contour1
        contour1Display = Show(contour1, renderView1, 'GeometryRepresentation')

        # trace defaults for the display properties.
        contour1Display.Representation = 'Surface'
        contour1Display.AmbientColor = [0.0392156862745098, 0.00784313725490196, 1.0]
        contour1Display.ColorArrayName = ['POINTS', '']
        contour1Display.DiffuseColor = [0.0392156862745098, 0.00784313725490196, 1.0]
        contour1Display.Opacity = 0.53
        contour1Display.Specular = 1.0
        contour1Display.SpecularPower = 1.0
        contour1Display.Ambient = 0.21
        contour1Display.OSPRayScaleArray = 'f'
        contour1Display.OSPRayScaleFunction = 'PiecewiseFunction'
        contour1Display.SelectOrientationVectors = 'None'
        contour1Display.ScaleFactor = 0.8
        contour1Display.SelectScaleArray = 'f'
        contour1Display.GlyphType = 'Arrow'
        contour1Display.GlyphTableIndexArray = 'f'
        contour1Display.GaussianRadius = 0.04
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
        Show()
        Render()
#         renderView1.Update()
        fn = path + "/" + picName + str(i) +  '_t=' + str(timesteps[i]) +'_noLambda2.png'
        SaveScreenshot( fn, renderView1,
                    #      ImageResolution=[2316, 2204],
                        TransparentBackground=0,
                        CompressionLevel='2' )
        print('File=' + fn + ' generated succesfully')
#****************** CONTOUR1(f) AND TRACERS ********************
        Hide(contour2, renderView1)
        # set color bar visibility
        uxLUTColorBar.Visibility = 1
        contour1Display.Opacity = 0.23
        # create a new 'Stream Tracer'
        streamTracer1 = StreamTracer(Input=resampleToImage1,
            SeedType='High Resolution Line Source')
        streamTracer1.Vectors = ['POINTS', 'u.x']
        streamTracer1.SurfaceStreamlines = 1
        streamTracer1.MaximumStreamlineLength = 8

        # init the 'High Resolution Line Source' selected for 'SeedType'
        streamTracer1.SeedType.Point1 = [bounds[1] + 1, -0.5, 0.0]
        streamTracer1.SeedType.Point2 = [bounds[1] + 1, 0.5, 0.0]
        streamTracer1.SeedType.Resolution = 50
        streamTracer1Display = Show(streamTracer1, renderView1, 'GeometryRepresentation')
        # show data from streamTracer1
        streamTracer1Display = Show(streamTracer1, renderView1, 'GeometryRepresentation')

        # trace defaults for the display properties.
        streamTracer1Display.Representation = 'Surface'
        streamTracer1Display.ColorArrayName = ['POINTS', 'u.x']
        streamTracer1Display.LookupTable = uxLUT
        streamTracer1Display.OSPRayScaleArray = 'AngularVelocity'
        streamTracer1Display.OSPRayScaleFunction = 'PiecewiseFunction'
        streamTracer1Display.SelectOrientationVectors = 'Normals'
        streamTracer1Display.ScaleFactor = 0.7899999618530273
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

        # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
        streamTracer1Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

        # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
        streamTracer1Display.ScaleTransferFunction.Points = [-9.447468533227408, 0.0, 0.5, 0.0, 7.903744217703089, 1.0, 0.5, 0.0]

        # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
        streamTracer1Display.OpacityTransferFunction.Points = [-9.447468533227408, 0.0, 0.5, 0.0, 7.903744217703089, 1.0, 0.5, 0.0]
        Show()
        Render()
#         renderView1.Update()
        fn = path + "/" + picName + str(i) +  '_t=' + str(timesteps[i]) +'_tracer.png'
        SaveScreenshot( fn, renderView1,
        #      ImageResolution=[2316, 2204],
            TransparentBackground=0,
            CompressionLevel='2' )
        Hide(streamTracer1, renderView1)
        print('File=' + fn + ' generated succesfully')


stop = timeit.default_timer()
total_time = stop - start

# output running time in a nice format.
mins, secs = divmod(total_time, 60)
hours, mins = divmod(mins, 60)

sys.stdout.write("Total running time: %d:%d:%d.\n" % (hours, mins, secs))