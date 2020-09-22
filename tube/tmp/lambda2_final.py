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
parser.add_argument("-maxlevel", type=int, help="Provide the maximum level of refinement",
                    nargs='?', default=12)
parser.add_argument("-nslices", type=int, help="Provide the number of slices",
                    nargs='?', default=10)
parser.add_argument("-lDomain", type=float, help="Provide the Domain size",
                    nargs='?', default=20)
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
parser.add_argument("-isolambda2", type=int, help="Provide the isosurfaces for lambda2 criterium",
                    nargs='+', default=[-2])
# parser.add_argument('--foo', action='store_const', const=2, default=42)
args = parser.parse_args()
print(args)
filenames = args.filenames
vidName = args.vidName
picName = args.picName
frameRate = args.frameRate
frameWindow = args.frameWindow
lDomain = args.lDomain
maxlevel = args.maxlevel
viewSize = args.viewSize
noVideo = args.noVideo
noPic = args.noPic
noData = args.noData
noLambda2 = args.noLambda2
timeList = args.timeList
nslices = args.nslices
isolambda2 = args.isolambda2
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
print("timesteps=", timesteps)
try:
    NT = len(timesteps)
except:
    NT = 1
    timesteps = [timesteps]

if len(timeList) == 1 and timeList[0] == 0:
    timeList = range(0,NT)
print("NT=", NT, " timeList=", timeList)

if frameWindow[0] == 0 and frameWindow[1] == 0:
    frameWindow[0] = 0
    frameWindow[1] = NT-1
    print('frameWindow is updated:',frameWindow)

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


# Create a new 'Render View'
renderView2 = CreateView('RenderView')
renderView2.ViewSize = [1168, 1168]
renderView2.InteractionMode = '2D'
renderView2.AxesGrid = 'GridAxes3DActor'
renderView2.OrientationAxesVisibility = 0
renderView2.StereoType = 'Crystal Eyes'
renderView2.BackEnd = 'OSPRay raycaster'
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

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)

###-----------------GENERATION of Cylinder-----------------------------------
### it is timeless therefore it is outside of the loop
# create a new 'Cylinder'
cylinder1 = Cylinder()
cylinder1.Resolution = 140
cylinder1.Height = 30.0
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
transform1.Transform.Translate = [15.0, 0.0, 0.0]
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
slice2.SliceType.Origin = [15.0, 0.0, -1e-5]
slice2.SliceType.Normal = [0.0, 0.0, 1.0]

# init the 'Plane' selected for 'HyperTreeGridSlicer'
slice2.HyperTreeGridSlicer.Origin = [15.0, 0.0, -0.21650634706020355]

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

###-----------------GENERATION of a CLIP and Convert to Point data-----------------------------------
if not noPic:
    for i in timeList:
        print("in loop iteration:" + str(i))
        renderView1.ViewSize = [3108, 1168]
        SetActiveView(renderView1)
        # Properties modified on animationScene1
        animationScene1.AnimationTime = timesteps[i]
        # Properties modified on timeKeeper1
        timeKeeper1.Time = timesteps[i]

        # create a new 'Threshold'
        threshold1 = Threshold(Input=my_source)
        threshold1.Scalars = ['POINTS', 'f']
        threshold1.ThresholdRange = [0.0, 0.]

        contour0 = Contour(Input=threshold1)
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
        len_bub = bounds[1] - bounds[0]
        len_min = max(bounds[0] - 2, 0.5)
        len_max = min(bounds[1] + 2, 29)
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
        clip1.ClipType.Position = [len_min - 0.1, -0.6, -0.6]
        clip1.ClipType.Length = [length + 0.2, 1.2, 1.2]
#         clip1.ClipType.Position = [len_min - 0.2, -0.55, -0.55]
#         clip1.ClipType.Length = [length + 0.4, 1.1, 1.1]

        # init the 'Plane' selected for 'HyperTreeGridClipper'
#         clip1.HyperTreeGridClipper.Origin = [15.0, 0.0, 0.0] #???

# ***************** CELL DATA TO POINT DATA ****************************
        # create a new 'Cell Data to Point Data'
        cellDatatoPointData1 = CellDatatoPointData(Input=clip1)
#         cellDatatoPointData1.CellDataArraytoprocess = ['f', 'fs', 'l', 'l2', 'omega', 'p', 'u.x']

# ***************** RESAMPLE TO IMAGE ****************************
        # create a new 'Resample To Image'
        resampleToImage1 = ResampleToImage(Input=cellDatatoPointData1)
        resampleToImage1.UseInputBounds = 0
        Nx = int(length/(lDomain/2.0**maxlevel))
        Nyz = int(1.0/(lDomain/2.0**maxlevel))
        resampleDimension = [Nx, Nyz, Nyz]
        print("resampleDimension:", resampleDimension)
        resampleToImage1.SamplingDimensions = resampleDimension
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
        contour2.Isosurfaces = isolambda2
        contour2.PointMergeMethod = 'Uniform Binning'

        # create a new 'Slice'
        slice1 = Slice(Input=contour1)
        slice1.SliceType = 'Plane'
        slice1.HyperTreeGridSlicer = 'Plane'
        slice1.SliceOffsetValues = [0]

        # init the 'Plane' selected for 'SliceType'
        slice1.SliceType.Origin = [0.5*(len_min + len_max), 0.0, 0.0]
        slice1.SliceType.Normal = [0.0, 0.0, 1.0]

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
        uxLUTColorBar.Position = [0.5 - 0.5*len_bar, 0.02]
        uxLUTColorBar.Title = 'u'
        uxLUTColorBar.ComponentTitle = 'Magnitude'
        uxLUTColorBar.LabelFormat = '%-#6.2g'
        uxLUTColorBar.RangeLabelFormat = '%6.2g'
        uxLUTColorBar.ScalarBarThickness = 16*2
        uxLUTColorBar.TitleFontSize = 16*2
        uxLUTColorBar.LabelFontSize = 16*2


        # set color bar visibility
        uxLUTColorBar.Visibility = 1

        # get color transfer function/color map for 'omega'
        omegaLUT = GetColorTransferFunction('omega')
        omegaLUT.RGBPoints = [0.3730543491279064, 0.0, 0.0, 0.5625, 0.5106463634876334, 0.0, 0.0, 1.0, 0.8251430154745502, 0.0, 1.0, 1.0, 0.9823910318856668, 0.5, 1.0, 0.5, 1.139639048296783, 1.0, 1.0, 0.0, 1.4541357002836997, 1.0, 0.0, 0.0, 1.6113837166948162, 0.5, 0.0, 0.0]
        omegaLUT.ColorSpace = 'RGB'
        omegaLUT.ScalarRangeInitialized = 1.0
        # omegaLUT.ApplyPreset('jet', True)
        len_bar = 0.5
        # get color legend/bar for omegaLUT in view renderView1
        omegaLUTColorBar = GetScalarBar(omegaLUT, renderView1)
        omegaLUTColorBar.Orientation = 'Horizontal'
        omegaLUTColorBar.WindowLocation = 'AnyLocation'
        omegaLUTColorBar.ScalarBarLength = len_bar
        omegaLUTColorBar.Position = [0.5 - 0.5*len_bar, 0.02]
        omegaLUTColorBar.Title = 'omega'
        omegaLUTColorBar.ComponentTitle = ''
        omegaLUTColorBar.LabelFormat = '%-#6.2g'
        omegaLUTColorBar.RangeLabelFormat = '%6.2g'


        # set color bar visibility
        omegaLUTColorBar.Visibility = 0

        len_bar = 0.5
        uxLUTColorBar2 = GetScalarBar(omegaLUT, renderView2)
        uxLUTColorBar2.Orientation = 'Horizontal'
        uxLUTColorBar2.WindowLocation = 'AnyLocation'
        uxLUTColorBar2.Position = [0.5 - 0.5*len_bar, 0.02]
        uxLUTColorBar2.Title = 'u'
        uxLUTColorBar2.ComponentTitle = 'Magnitude'
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
            # show data from contour1
            contour1Display = Show(contour1, rv, 'GeometryRepresentation')
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
#         contour1Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

        # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
#         contour1Display.ScaleTransferFunction.Points = [0.5, 0.0, 0.5, 0.0, 0.5001220703125, 1.0, 0.5, 0.0]

        # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
#         contour1Display.OpacityTransferFunction.Points = [0.5, 0.0, 0.5, 0.0, 0.5001220703125, 1.0, 0.5, 0.0]

        # setup the color legend parameters for each legend in this view

        # show color legend
        contour1Display.SetScalarBarVisibility(renderView1, True)

        # ----------------------------------------------------------------
        # setup color maps and opacity mapes used in the visualization
        # note: the Get..() functions create a new object, if needed
        # ----------------------------------------------------------------
        info = resampleToImage1.GetDataInformation().DataInformation
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
        uxLUTColorBar.CustomLabels = np.around(np.linspace(0, range_max, 4),2)
#         uxLUTColorBar.CustomLabels = np.arange(0, range_max, 0.5)
        print("CustomLabels=",uxLUTColorBar.CustomLabels )

        uxLUTColorBar2.UseCustomLabels = 1
        # labels from 100 to 200 in increments of 10
        uxLUTColorBar2.CustomLabels = np.around(np.linspace(0, range_max, 4),2)
#         uxLUTColorBar2.CustomLabels = np.arange(0, range_max, 0.5)
        print("CustomLabels=",uxLUTColorBar2.CustomLabels )


        info = resampleToImage1.GetDataInformation().DataInformation
        arrayInfo = info.GetArrayInformation("omega", vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS)
        range0 = arrayInfo.GetComponentRange(0)
        range_max = np.max(np.abs(range0))
        print("omega range= ", range0)
        omegaLUT.RescaleTransferFunction(-range_max, range_max)
        op = GetOpacityTransferFunction("omega")
        op.RescaleTransferFunction(-range_max, range_max)
        omegaLUTColorBar.UseCustomLabels = 1
        # labels from 100 to 200 in increments of 10
#         omegaLUTColorBar.CustomLabels = np.linspace(0, range_max, 5)
        omegaLUTColorBar.CustomLabels = np.around(np.linspace(-range_max, range_max, 5),2)
        print("CustomLabels=",omegaLUTColorBar.CustomLabels )



        # get opacity transfer function/opacity map for 'ux'
#         uxPWF = GetOpacityTransferFunction('ux')
#         uxPWF.Points = [0.3730543491279064, 0.0, 0.5, 0.0, 1.6113837166948162, 1.0, 0.5, 0.0]
#         uxPWF.ScalarRangeInitialized = 1
#
#          # get opacity transfer function/opacity map for 'omega'
#         omegaPWF = GetOpacityTransferFunction('omega')
#         omegaPWF.Points = [0.3730543491279064, 0.0, 0.5, 0.0, 1.6113837166948162, 1.0, 0.5, 0.0]
#         omegaPWF.ScalarRangeInitialized = 1
#
#         # get opacity transfer function/opacity map for 'l2'
#         l2PWF = GetOpacityTransferFunction('l2')
#         l2PWF.Points = [-2.0, 0.0, 0.5, 0.0, -1.0, 1.0, 0.5, 0.0]
#         l2PWF.ScalarRangeInitialized = 1
        # ----------------------------------------------------------------
        # ----------------------------------------------------------------

#         renderView1.CameraPosition = [center[0], center[1], 6]
#         renderView1.CameraFocalPoint = center
        renderView1.CameraViewUp = [0.2, 1, 0]
        renderView1.CameraParallelScale = 1.5
        renderView1.CenterOfRotation = [0.5*(len_min + len_max), 0.0, 0.0]
        renderView1.CameraFocalPoint = [0.5*(len_min + len_max) , 0, 0]
#         renderView1.CameraFocalPoint = [0.5*(len_min + len_max) , -0.25, -1.3]
        renderView1.CameraPosition = [0.5*(len_min + len_max) - 2, 0.6, 4.5]
#         renderView1.CameraParallelScale = 1.8

        # update the view to ensure updated data information
        renderView1.Update()

#****************** CONTOUR1(f) AND U MAGNITUDE ********************
        # show data from contour1
        fn = path + "/" + picName +  '_t=' + str(timesteps[i]) +'ux.png'
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

        renderView1.Update()
        fn = path + "/" + picName+  '_t=' + str(timesteps[i]) +'_noLambda2.png'
        SaveScreenshot( fn, renderView1,
                    #      ImageResolution=[2316, 2204],
                        TransparentBackground=0,
                        CompressionLevel='2' )
        print('File=' + fn + ' generated succesfully')
#****************** CONTOUR1(f) AND TRACERS ********************
        Hide(contour1, renderView1)
        Hide(contour2, renderView1)
        # set color bar visibility
        uxLUTColorBar.Visibility = 1
        contour1Display.Opacity = 0.23
        # create a new 'Stream Tracer'
        streamTracer1 = StreamTracer(Input=resampleToImage1, SeedType='High Resolution Line Source')
        streamTracer1.Vectors = ['POINTS', 'u.x']
        streamTracer1.SurfaceStreamlines = 1
        streamTracer1.InitialStepLength = 0.05
        streamTracer1.MaximumStreamlineLength = 7.800000000000001
        streamTracer1.TerminalSpeed = 0.0001

        # init the 'High Resolution Line Source' selected for 'SeedType'
        streamTracer1.SeedType.Point1 = [bounds[1] + 1, -0.5, 0.0]
        streamTracer1.SeedType.Point2 = [bounds[1] + 1, 0.5, 0.0]
        streamTracer1.SeedType.Resolution = 50
        # show data from streamTracer1
        streamTracer1Display = Show(streamTracer1, renderView1, 'GeometryRepresentation')

        # trace defaults for the display properties.
        streamTracer1Display.Representation = 'Surface'
        streamTracer1Display.ColorArrayName = ['POINTS', 'u.x']
        streamTracer1Display.LookupTable = uxLUT
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


        # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
#         slice1Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

        # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
#         slice1Display.ScaleTransferFunction.Points = [0.5, 0.0, 0.5, 0.0, 0.5001220703125, 1.0, 0.5, 0.0]

        # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
#         slice1Display.OpacityTransferFunction.Points = [0.5, 0.0, 0.5, 0.0, 0.5001220703125, 1.0, 0.5, 0.0]

        # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
#         streamTracer1Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

        # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
#         streamTracer1Display.ScaleTransferFunction.Points = [-9.447468533227408, 0.0, 0.5, 0.0, 7.903744217703089, 1.0, 0.5, 0.0]

        # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
#         streamTracer1Display.OpacityTransferFunction.Points = [-9.447468533227408, 0.0, 0.5, 0.0, 7.903744217703089, 1.0, 0.5, 0.0]

        renderView1.CameraPosition = [center[0], center[1], 6]
        renderView1.CameraFocalPoint = center
        renderView1.CameraParallelScale = 1
        renderView1.CameraViewUp = [0, 1, 0]
#         len_min = max(bounds[0] - 1, 0.5)
#         len_max = min(bounds[1] + 1, 29)
#         renderView1.CameraPosition = [center[0] , center[1], 6]
#         renderView1.CameraFocalPoint = center
#         renderView1.CameraParallelScale = 1
        #         renderView1.CenterOfRotation = [0.5*(len_min + len_max), 0.0, 0.0]
        #         renderView1.CameraFocalPoint = [0.5*(len_min + len_max) , 0, 0]
        #         renderView1.CameraFocalPoint = [0.5*(len_min + len_max) , -0.25, -1.3]
        #         renderView1.CameraPosition = [0.5*(len_min + len_max) , 0, 4.5]


        renderView1.Update()
        fn = path + "/" + picName +  '_t=' + str(timesteps[i]) +'_tracer.png'
        SaveScreenshot( fn, renderView1,
        #      ImageResolution=[2316, 2204],
            TransparentBackground=0,
            CompressionLevel='2' )
        Hide(streamTracer1, renderView1)

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

        # show data from contour1
        contour1Display = Show(contour1, renderView1, 'GeometryRepresentation')
        # trace defaults for the display properties.
        contour1Display.Representation = 'Surface'
        contour1Display.ColorArrayName = ['POINTS', 'u.x']
        contour1Display.LookupTable = uxLUT
        contour1Display.Opacity = 1
        contour1Display.Specular = 0.8
        contour1Display.SpecularPower = 100.0
        contour1Display.Ambient = 0.1
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

        uxLUTColorBar.Visibility = 1
        renderView1.Update()
        fn = path + "/" + picName +  '_t=' + str(timesteps[i]) +'_uxSide.png'
        SaveScreenshot( fn, renderView1,
        #      ImageResolution=[2316, 2204],
            TransparentBackground=0,
            CompressionLevel='2' )
#         Hide(streamTracer1, renderView1)
        Hide(slice3, renderView1)
        Hide(slice1, renderView1)

        print('File=' + fn + ' generated succesfully')

         # create a new 'Slice'
        slice4 = Slice(Input=resampleToImage1)
        slice4.SliceType = 'Plane'
        slice4.HyperTreeGridSlicer = 'Plane'
        slice4.SliceOffsetValues = [0]
         # create a new 'Slice' contour(f)
        slice5 = Slice(Input=contour1)
        slice5.SliceType = 'Plane'
        slice5.HyperTreeGridSlicer = 'Plane'
        slice5.SliceOffsetValues = [0]
#**************** SLICES of U along a bubble ******************
        ss = GetSources()
        for s in ss:
          Hide(ss[s])
        uxLUTColorBar.ScalarBarThickness = 16
        uxLUTColorBar.TitleFontSize = 16
        uxLUTColorBar.LabelFontSize = 16
        for k,sl in enumerate(np.linspace(bounds[0]+0.1, bounds[1]-0.1, nslices)):
#         for k,sl in enumerate(np.linspace(len_min+1, len_max-1, nslices)):
            print("slice x=", sl)
            renderView1.CameraPosition = [sl - 3, 0, 0]
            renderView1.CameraFocalPoint = [sl, 0, 0]
            renderView1.CameraParallelScale = 0.7
            renderView1.CameraViewUp = [0, 0, 1]

            # init the 'Plane' selected for 'SliceType'
            slice4.SliceType.Origin = [sl, 0.0, 0.0]
            slice4.SliceType.Normal = [1.0, 0.0, 0.0]

            # init the 'Plane' selected for 'HyperTreeGridSlicer'
            # slice4.HyperTreeGridSlicer.Origin = [10.999999761581421, 0.0, 0.0] #???

            slice4Display = Show(slice4, renderView1, 'GeometryRepresentation')

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
            Show(transform1, renderView1)
            uxLUTColorBar.Visibility = 1
            omegaLUTColorBar.Visibility = 0
            renderView1.Update()
            fn = path + "/" + picName +  '_t=' + str(timesteps[i]) +'_n=' + str(k) + '_uxSlice.png'
            SaveScreenshot( fn, renderView1,
                 ImageResolution=[2316, 2204],
                TransparentBackground=0,
                CompressionLevel='2' )

            print('File=' + fn + ' generated succesfully')
        Hide(slice4,renderView1)
        Hide(slice5,renderView1)





#**************** SLICES of Omega along a bubble ******************
        ss = GetSources()
        for s in ss:
          Hide(ss[s])

#         slice4 = Slice(Input=resampleToImage1)
#         slice4.SliceType = 'Plane'
#         slice4.HyperTreeGridSlicer = 'Plane'
#         slice4.SliceOffsetValues = [0]
#
#         # create a new 'Slice' contour(f)
#         slice5 = Slice(Input=contour1)
#         slice5.SliceType = 'Plane'
#         slice5.HyperTreeGridSlicer = 'Plane'
#         slice5.SliceOffsetValues = [0]

        uxLUTColorBar.ScalarBarThickness = 16
        uxLUTColorBar.TitleFontSize = 16
        uxLUTColorBar.LabelFontSize = 16
        for k,sl in enumerate(np.linspace(bounds[0]+0.1, bounds[1]-0.1, nslices)):
            print("slice x=", sl)

            # create a new 'Text'
            text1 = Text()
            text1.Text = 'x=' + str(round(sl,2)) + " l/l_b=" + str(round((sl-bounds[0])/len_bub,2))

#             renderView1.CameraPosition = [sl - 3, 0, 0]
#             renderView1.CameraFocalPoint = [sl, 0, 0]
#             renderView1.CameraParallelScale = 0.7
            renderView1.CameraViewUp = [0, 0, 1]
            renderView1.ViewSize = [2048, 2048]
            renderView1.ResetCamera(sl, sl, -0.5, 0.5, -0.5, 0.5)

#             renderView2.CameraPosition = [center[0], center[1], 6]
#             renderView2.CameraFocalPoint = center
#             renderView2.CameraParallelScale = 1.7
            renderView2.CameraViewUp = [0, 1, 0]
            renderView2.ResetCamera(bounds[0], bounds[1], -0.5, 0.5, -0.5, 0.5)
            SetActiveView(renderView2)
            camera = GetActiveCamera()
            camera.Azimuth(10)

            print("renderView2:",bounds[0], bounds[1], -0.5, 0.5, -0.5, 0.5)

            # create a new 'Slice'


            # init the 'Plane' selected for 'SliceType'
            slice4.SliceType.Origin = [sl, 0.0, 0.0]
            slice4.SliceType.Normal = [1.0, 0.0, 0.0]

            # init the 'Plane' selected for 'HyperTreeGridSlicer'
            # slice4.HyperTreeGridSlicer.Origin = [10.999999761581421, 0.0, 0.0] #???

            # show data from text1
            text1Display = Show(text1, renderView2, 'TextSourceRepresentation')

            # trace defaults for the display properties.
            text1Display.WindowLocation = 'AnyLocation'
            text1Display.Position = [0.3, 0.2]
            renderView2.Update()
#             textBounds = text1.GetDataInformation().GetBounds() # doesnot work
#             print("textBounds=",textBounds)

            for rv in [renderView1,renderView2]:
                slice4Display = Show(slice4, rv, 'GeometryRepresentation')

                # trace defaults for the display properties.
                slice4Display.Representation = 'Surface'
                slice4Display.ColorArrayName = ['POINTS', 'u.x']
                slice4Display.LookupTable = omegaLUT
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
            Show(transform1, renderView1)
            # set color bar visibility
            uxLUTColorBar.Visibility = 0
            omegaLUTColorBar.Visibility = 1

#             uxLUTColorBar2.Visibility = 1
            # get display properties
#             contour1Display = GetDisplayProperties(contour2, view=renderView1)

            # show color bar/color legend
            contour1Display.SetScalarBarVisibility(renderView2, True)
            # get color legend/bar for uxLUT in view renderView2
            uxLUTColorBar2 = GetScalarBar(uxLUT, renderView2)
            uxLUTColorBar2.WindowLocation = 'UpperLeftCorner'

            uxLUTColorBar2.ComponentTitle = ''
            uxLUTColorBar2.Title = 'u'
            uxLUTColorBar2.ComponentTitle = 'Magnitude'
            uxLUTColorBar2.LabelFormat = '%-#6.2g'
            uxLUTColorBar2.RangeLabelFormat = '%6.2g'

            len_bar = 0.5
            # change scalar bar placement
            uxLUTColorBar2.Orientation = 'Horizontal'
            uxLUTColorBar2.WindowLocation = 'AnyLocation'
            uxLUTColorBar2.Position = [0.5 - 0.5*len_bar, 0.02]
            uxLUTColorBar2.ScalarBarLength = len_bar

            Show(contour1, renderView2, 'GeometryRepresentation')
            Show(slice4, renderView2, 'GeometryRepresentation')
            Show(transform1, renderView2, 'GeometryRepresentation')
            renderView1.Update()
            renderView2.Update()
            fn = path + "/" + picName +  '_t=' + str(timesteps[i]) +'_n=' + str(k) + '_uxSlice.png'


            SaveScreenshot( fn,  layout1, SaveAllViews=1,
            #renderView1,
                 ImageResolution=[4096, 2048],
                TransparentBackground=0,
                CompressionLevel='2' )
#             uxLUTColorBar2.Visibility = 0
            Hide(slice4, renderView1)
            Hide(slice5, renderView1)
            Hide(contour1, renderView2)
            Hide(slice5, renderView2)
            SetActiveView(renderView2)
            camera = GetActiveCamera()
            camera.Azimuth(-10)
#             Delete(slice4)
#             del slice4
#             Delete(slice5)
#             del slice5
            Delete(text1)
            del text1
#             Hide(linearExtrusion1, renderView2)
            print('File=' + fn + ' generated succesfully')





#**************** SLICES of Omega along a bubble ******************
        ss = GetSources()
        for s in ss:
          Hide(ss[s])

        slice4 = Slice(Input=resampleToImage1)
        slice4.SliceType = 'Plane'
        slice4.HyperTreeGridSlicer = 'Plane'
        slice4.SliceOffsetValues = [0]

        # create a new 'Slice' contour(f)
        slice5 = Slice(Input=contour1)
        slice5.SliceType = 'Plane'
        slice5.HyperTreeGridSlicer = 'Plane'
        slice5.SliceOffsetValues = [0]
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

            # create a new 'Slice'


            # init the 'Plane' selected for 'SliceType'
            slice4.SliceType.Origin = [sl, 0.0, 0.0]
            slice4.SliceType.Normal = [1.0, 0.0, 0.0]

            # init the 'Plane' selected for 'HyperTreeGridSlicer'
            # slice4.HyperTreeGridSlicer.Origin = [10.999999761581421, 0.0, 0.0] #???

            # show data from text1
            text1Display = Show(text1, renderView2, 'TextSourceRepresentation')

            # trace defaults for the display properties.
            text1Display.WindowLocation = 'AnyLocation'
            text1Display.Position = [0.3, 0.2]
            renderView2.Update()
            textBounds = text1.GetDataInformation().GetBounds()
            print("textBounds=",textBounds)

            for rv in [renderView1,renderView2]:
                slice4Display = Show(slice4, rv, 'GeometryRepresentation')

                # trace defaults for the display properties.
                slice4Display.Representation = 'Surface'
                slice4Display.ColorArrayName = ['POINTS', 'omega']
                slice4Display.LookupTable = omegaLUT
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
            Show(transform1, renderView1)
            # set color bar visibility
            uxLUTColorBar.Visibility = 0
            omegaLUTColorBar.Visibility = 1

#             uxLUTColorBar2.Visibility = 1
            # get display properties
#             contour1Display = GetDisplayProperties(contour2, view=renderView1)

            # show color bar/color legend
            contour1Display.SetScalarBarVisibility(renderView2, True)
            # get color legend/bar for uxLUT in view renderView2
            uxLUTColorBar2 = GetScalarBar(uxLUT, renderView2)
            uxLUTColorBar2.WindowLocation = 'UpperLeftCorner'

            uxLUTColorBar2.ComponentTitle = ''
            uxLUTColorBar2.Title = 'u'
            uxLUTColorBar2.ComponentTitle = 'Magnitude'
            uxLUTColorBar2.LabelFormat = '%-#6.2g'
            uxLUTColorBar2.RangeLabelFormat = '%6.2g'

            len_bar = 0.5
            # change scalar bar placement
            uxLUTColorBar2.Orientation = 'Horizontal'
            uxLUTColorBar2.WindowLocation = 'AnyLocation'
            uxLUTColorBar2.Position = [0.5 - 0.5*len_bar, 0.02]
            uxLUTColorBar2.ScalarBarLength = len_bar

            Show(contour1, renderView2, 'GeometryRepresentation')
            Show(slice4, renderView2, 'GeometryRepresentation')
            Show(transform1, renderView2, 'GeometryRepresentation')
            renderView1.Update()
            renderView2.Update()
            fn = path + "/" + picName +  '_t=' + str(timesteps[i]) +'_n=' + str(k) + '_omegaSlice.png'


            SaveScreenshot( fn,  layout1, SaveAllViews=1,
            #renderView1,
                 ImageResolution=[4096, 2048],
                TransparentBackground=0,
                CompressionLevel='2' )
#             uxLUTColorBar2.Visibility = 0
            Hide(slice4, renderView1)
            Hide(slice5, renderView1)
            Hide(contour1, renderView2)
            Hide(slice5, renderView2)
            SetActiveView(renderView2)
            camera = GetActiveCamera()
            camera.Azimuth(-10)
#             Delete(slice4)
#             del slice4
#             Delete(slice5)
#             del slice5
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

        print("Tracers")
        # show data in view
        contour3Display = Show(streamTracer1, spreadSheetView2, 'SpreadSheetRepresentation')
        fn = path + '/' + 'tracer_t=' + str(timesteps[i]) + '.csv'
        ExportView(fn, view=spreadSheetView2)
        print('File=' + fn + ' generated succesfully')
        # Freeing Memory
        Delete(streamTracer1)
        del streamTracer1
        Delete(slice1)
        del slice1
        Delete(contour2)
        del contour2
        Delete(contour1)
        del contour1
        Delete(resampleToImage1)
        del resampleToImage1
        Delete(cellDatatoPointData1)
        del cellDatatoPointData1
        Delete(clip1)
        del clip1
        Delete(contour0)
        del contour0




stop = timeit.default_timer()
total_time = stop - start

# output running time in a nice format.
mins, secs = divmod(total_time, 60)
hours, mins = divmod(mins, 60)

sys.stdout.write("Total running time: %d:%d:%d.\n" % (hours, mins, secs))