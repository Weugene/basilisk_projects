# state file generated using paraview version 5.9.0

#### import the simple module from the paraview
from paraview.simple import *
import glob, os, sys
import logging
import timeit
import argparse
logging.basicConfig(format='%(message)s')
log = logging.getLogger(__name__)
def eprint(var):
    log.warning(var)
import json
from scipy.spatial import Delaunay

import functools
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()
# --------------------------------User functions -----------------
def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')
# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

start = timeit.default_timer()

parser = argparse.ArgumentParser()
parser._action_groups.pop()
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
optional.add_argument("-filename", type=str, help="Provide the name of the input paraview files, please",
                      nargs='?', default='output_res.pvd')
optional.add_argument("-frameRate", type=int, help="Provide the frameRate, please",
                      nargs='?', default=10)
parser.add_argument("-picName", type=str, help="Provide the name for the output pictures, please",
                    nargs='?', default='pone')
parser.add_argument('-noPic', type=str2bool, nargs='?', const=True, default=False)
parser.add_argument('-noVideo', type=str2bool, nargs='?', const=True, default=False)
parser.add_argument("-timeList", type=int, help="Provide the list of step, please",
                    nargs='+', default=[0])
optional.add_argument("-nt", type=int, help="Provide the nt argument",
                      nargs='?', default=0)
optional.add_argument("-fontsize", type=int, help="Provide the font size argument",
                      nargs='?', default=32)
parser.add_argument("-viewSize", type=int, help="Provide the view size of the output pictures, please",
                    nargs='+', default=[1600, 1298])

args = parser.parse_args()
print(args)
filename = args.filename
frameRate = args.frameRate
viewSize = args.viewSize
picName = args.picName
noPic = args.noPic
noVideo = args.noVideo
timeList = args.timeList
nt = args.nt
fontsize = args.fontsize

#Current PATH reading
path = os.path.abspath(os.getcwd())

eprint("Current PATH=" + path)

if filename[-5::] == '.pvtu':
    # Find files with *.pvtu extension
    numbers = []
    file = ""
    for file in glob.glob(filename):
        numbers.append(int(filter(lambda x: x.isdigit(), file)))
    file=file[0:-9]
    numbers.sort()
    N = len(numbers)
    infn = []
    for i in numbers:
        infn.append('{}/{}{:04d}.pvtu'.format(path, file, i))
    print(infn)
    fn = path +  '/' + infn
    output_respvd = XMLPartitionedUnstructuredGridReader(FileName = fn)
elif filename[-4::] == '.pvd':
    # create a new 'PVD Reader'
    print('Read the first one:',filename)
    output_respvd = PVDReader(FileName = path +  '/' + filename)
else:
    eprint('Get Active Source: No pvd or pvtu files are provided')
    output_respvd = GetActiveSource()

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

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'

renderView1 = CreateView('RenderView')
renderView1.ViewSize = viewSize
# renderView1.ViewSize = [1644, 1298]
renderView1.InteractionMode = '2D'
# renderView1.AxesGrid = 'GridAxes3DActor'
# renderView1.OrientationAxesVisibility = 0
# renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [-0.8, 0, 10]
renderView1.CameraFocalPoint = [-0.8, 0, 0.0]
# renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 5
renderView1.OrientationAxesVisibility = 0
# renderView1.BackEnd = 'OSPRay raycaster'
renderView1.OSPRayMaterialLibrary = materialLibrary1


SetActiveView(None)

# defining of computational domains
Show(output_respvd, renderView1)
Hide(output_respvd, renderView1)
boundsDomain = output_respvd.GetDataInformation().GetBounds()

lDomain = boundsDomain[1] - boundsDomain[0]
print("boundsDomain of my source=", boundsDomain, " lDomain=", lDomain)


# renderView1.Update()
# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.AssignView(0, renderView1)
layout1.SetSize(viewSize[0], viewSize[1])

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------
renderView1.ResetCamera(1.3*boundsDomain[0] , boundsDomain[1], boundsDomain[2], boundsDomain[3], boundsDomain[4], boundsDomain[5])
renderView1.CameraParallelScale = lDomain/2.0
renderView1.Update()
# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------
latex_names = {'alpha_doc':'$\\alpha$', 'u':'$|\\mathbf{u}|$', 'mu_cell':'$\\mu$', 'T':'$T$'}

for i in timeList:
    print("in loop iteration:" + str(i))
    SetActiveView(renderView1)
    # Properties modified on animationScene1
    animationScene1.AnimationTime = timesteps[i]
    # Properties modified on timeKeeper1
    timeKeeper1.Time = timesteps[i]
    for feature in ['alpha_doc', 'u', 'mu_cell', 'T']:
        os.makedirs(f'{path}/{picName}/{feature}', exist_ok=True)
        # create a new 'HyperTreeGrid To Dual Grid'
        hyperTreeGridToDualGrid1 = HyperTreeGridToDualGrid(registrationName='HyperTreeGridToDualGrid1', Input=output_respvd)

        # create a new 'Contour'
        contour1 = Contour(registrationName='Contour1', Input=hyperTreeGridToDualGrid1)
        contour1.ContourBy = ['POINTS', 'fs']
        contour1.Isosurfaces = [0.5]
        contour1.PointMergeMethod = 'Uniform Binning'

        # create a new 'Contour'
        contour2 = Contour(registrationName='Contour2', Input=hyperTreeGridToDualGrid1)
        contour2.ContourBy = ['POINTS', 'f']
        contour2.Isosurfaces = [0.5]
        contour2.PointMergeMethod = 'Uniform Binning'

        # ----------------------------------------------------------------
        # setup the visualization in view 'renderView1'
        # ----------------------------------------------------------------

        # show data from output_respvd
        output_respvdDisplay = Show(output_respvd, renderView1, 'GeometryRepresentation')

        # get color transfer function/color map for feature
        featureLUT = GetColorTransferFunction(feature)
        featureLUT.RGBPoints = [0.0, 0.231373, 0.298039, 0.752941, 1.5, 0.865003, 0.865003, 0.865003, 3.0, 0.705882, 0.0156863, 0.14902]
        featureLUT.ScalarRangeInitialized = 1.0
        # Rescale transfer function
        if feature == 'alpha_doc':
            featureLUT.RescaleTransferFunction(0.0, 1.0)

        # trace defaults for the display properties.
        output_respvdDisplay.Representation = 'Surface'
        output_respvdDisplay.ColorArrayName = ['CELLS', feature]
        output_respvdDisplay.LookupTable = featureLUT
        output_respvdDisplay.SelectTCoordArray = 'None'
        output_respvdDisplay.SelectNormalArray = 'None'
        output_respvdDisplay.SelectTangentArray = 'None'
        output_respvdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
        output_respvdDisplay.SelectOrientationVectors = 'None'
        output_respvdDisplay.SelectScaleArray = 'fs_face'
        output_respvdDisplay.GlyphType = 'Arrow'
        output_respvdDisplay.GlyphTableIndexArray = 'fs_face'
        output_respvdDisplay.GaussianRadius = 0.05
        output_respvdDisplay.SetScaleArray = [None, '']
        output_respvdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
        output_respvdDisplay.OpacityArray = [None, '']
        output_respvdDisplay.OpacityTransferFunction = 'PiecewiseFunction'
        output_respvdDisplay.DataAxesGrid = 'GridAxesRepresentation'
        output_respvdDisplay.PolarAxes = 'PolarAxesRepresentation'

        # show data from hyperTreeGridToDualGrid1
        hyperTreeGridToDualGrid1Display = Show(hyperTreeGridToDualGrid1, renderView1, 'UnstructuredGridRepresentation')

        # get opacity transfer function/opacity map for feature
        featurePWF = GetOpacityTransferFunction(feature)
        featurePWF.Points = [0.0, 0.0, 0.5, 0.0, 3.0, 1.0, 0.5, 0.0]
        featurePWF.ScalarRangeInitialized = 1
        # Rescale transfer function
        if feature == 'alpha_doc':
            featurePWF.RescaleTransferFunction(0.0, 1.0)

        # trace defaults for the display properties.
        hyperTreeGridToDualGrid1Display.Representation = 'Surface'
        hyperTreeGridToDualGrid1Display.ColorArrayName = ['POINTS', feature]
        hyperTreeGridToDualGrid1Display.LookupTable = featureLUT
        hyperTreeGridToDualGrid1Display.SelectTCoordArray = 'None'
        hyperTreeGridToDualGrid1Display.SelectNormalArray = 'None'
        hyperTreeGridToDualGrid1Display.SelectTangentArray = 'None'
        hyperTreeGridToDualGrid1Display.OSPRayScaleArray = 'fs_face'
        hyperTreeGridToDualGrid1Display.OSPRayScaleFunction = 'PiecewiseFunction'
        hyperTreeGridToDualGrid1Display.SelectOrientationVectors = 'None'
        hyperTreeGridToDualGrid1Display.SelectScaleArray = 'fs_face'
        hyperTreeGridToDualGrid1Display.GlyphType = 'Arrow'
        hyperTreeGridToDualGrid1Display.GlyphTableIndexArray = 'fs_face'
        hyperTreeGridToDualGrid1Display.GaussianRadius = 0.05
        hyperTreeGridToDualGrid1Display.SetScaleArray = ['POINTS', 'fs_face']
        hyperTreeGridToDualGrid1Display.ScaleTransferFunction = 'PiecewiseFunction'
        hyperTreeGridToDualGrid1Display.OpacityArray = ['POINTS', 'fs_face']
        hyperTreeGridToDualGrid1Display.OpacityTransferFunction = 'PiecewiseFunction'
        hyperTreeGridToDualGrid1Display.DataAxesGrid = 'GridAxesRepresentation'
        hyperTreeGridToDualGrid1Display.PolarAxes = 'PolarAxesRepresentation'
        hyperTreeGridToDualGrid1Display.ScalarOpacityFunction = featurePWF
        hyperTreeGridToDualGrid1Display.ScalarOpacityUnitDistance = 0.6912841351258522
        hyperTreeGridToDualGrid1Display.OpacityArrayName = ['POINTS', 'fs_face']


        # get display properties
        hyperTreeGridToDualGrid2Display = GetDisplayProperties(hyperTreeGridToDualGrid1, view=renderView1)
        # rescale color and/or opacity maps used to exactly fit the current data range
        if feature != 'alpha_doc':
            hyperTreeGridToDualGrid2Display.RescaleTransferFunctionToDataRange(False, True)

        # show data from contour2
        contour2Display = Show(contour2, renderView1, 'GeometryRepresentation')

        # trace defaults for the display properties.
        contour2Display.Representation = 'Surface'
        contour2Display.AmbientColor = [0.0, 0.0, 1.0]
        contour2Display.ColorArrayName = ['POINTS', '']
        contour2Display.DiffuseColor = [0.0, 0.0, 1.0]
        contour2Display.LineWidth = 5.0
        contour2Display.SelectTCoordArray = 'None'
        contour2Display.SelectNormalArray = 'None'
        contour2Display.SelectTangentArray = 'None'
        contour2Display.OSPRayScaleArray = 'f'
        contour2Display.OSPRayScaleFunction = 'PiecewiseFunction'
        contour2Display.SelectOrientationVectors = 'None'
        contour2Display.ScaleFactor = 0.18000437021255494
        contour2Display.SelectScaleArray = 'f'
        contour2Display.GlyphType = 'Arrow'
        contour2Display.GlyphTableIndexArray = 'f'
        contour2Display.GaussianRadius = 0.009000218510627747
        contour2Display.SetScaleArray = ['POINTS', 'f']
        contour2Display.ScaleTransferFunction = 'PiecewiseFunction'
        contour2Display.OpacityArray = ['POINTS', 'f']
        contour2Display.OpacityTransferFunction = 'PiecewiseFunction'
        contour2Display.DataAxesGrid = 'GridAxesRepresentation'
        contour2Display.PolarAxes = 'PolarAxesRepresentation'

        # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
        contour2Display.ScaleTransferFunction.Points = [0.5, 0.0, 0.5, 0.0, 0.5001220703125, 1.0, 0.5, 0.0]

        # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
        contour2Display.OpacityTransferFunction.Points = [0.5, 0.0, 0.5, 0.0, 0.5001220703125, 1.0, 0.5, 0.0]

        # show data from contour1
        contour1Display = Show(contour1, renderView1, 'GeometryRepresentation')

        # trace defaults for the display properties.
        contour1Display.Representation = 'Surface'
        contour1Display.AmbientColor = [0.3333333333333333, 0.0, 0.0]
        contour1Display.ColorArrayName = ['POINTS', '']
        contour1Display.DiffuseColor = [0.3333333333333333, 0.0, 0.0]
        contour1Display.LineWidth = 5.0
        contour1Display.SelectTCoordArray = 'None'
        contour1Display.SelectNormalArray = 'None'
        contour1Display.SelectTangentArray = 'None'
        contour1Display.OSPRayScaleArray = 'fs'
        contour1Display.OSPRayScaleFunction = 'PiecewiseFunction'
        contour1Display.SelectOrientationVectors = 'None'
        contour1Display.ScaleFactor = 0.6996097326278687
        contour1Display.SelectScaleArray = 'fs'
        contour1Display.GlyphType = 'Arrow'
        contour1Display.GlyphTableIndexArray = 'fs'
        contour1Display.GaussianRadius = 0.034980486631393436
        contour1Display.SetScaleArray = ['POINTS', 'fs']
        contour1Display.ScaleTransferFunction = 'PiecewiseFunction'
        contour1Display.OpacityArray = ['POINTS', 'fs']
        contour1Display.OpacityTransferFunction = 'PiecewiseFunction'
        contour1Display.DataAxesGrid = 'GridAxesRepresentation'
        contour1Display.PolarAxes = 'PolarAxesRepresentation'

        # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
        contour1Display.ScaleTransferFunction.Points = [0.5, 0.0, 0.5, 0.0, 0.5001220703125, 1.0, 0.5, 0.0]

        # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
        contour1Display.OpacityTransferFunction.Points = [0.5, 0.0, 0.5, 0.0, 0.5001220703125, 1.0, 0.5, 0.0]

        # setup the color legend parameters for each legend in this view

        # get color legend/bar for featureLUT in view renderView1
        featureLUTColorBar = GetScalarBar(featureLUT, renderView1)
        featureLUTColorBar.WindowLocation = 'UpperLeftCorner'
        # featureLUTColorBar.Position = [0.0069051881724923395, 0.3364396537793629]
        featureLUTColorBar.Title = latex_names[feature]
        featureLUTColorBar.ComponentTitle = ''
        featureLUTColorBar.TitleFontFamily = 'Times'
        featureLUTColorBar.LabelFontFamily = 'Times'
        featureLUTColorBar.ScalarBarLength = 0.8
        featureLUTColorBar.LabelFormat = '%-#6.2g'
        featureLUTColorBar.RangeLabelFormat = '%6.2g'
        featureLUTColorBar.ScalarBarThickness = fontsize
        featureLUTColorBar.TitleFontSize = fontsize
        featureLUTColorBar.LabelFontSize = fontsize

        # set color bar visibility
        featureLUTColorBar.Visibility = 1

        # show color legend
        output_respvdDisplay.SetScalarBarVisibility(renderView1, True)

        # show color legend
        hyperTreeGridToDualGrid1Display.SetScalarBarVisibility(renderView1, True)

        renderView1.Update()
        fn = f'{path}/{picName}/{feature}/{picName}_{feature}_t={timesteps[i]}.png'
        SaveScreenshot( fn,  layout1, SaveAllViews=1,
                        ImageResolution=viewSize,
                        TransparentBackground=0,
                        CompressionLevel='2' )
        # turn off color bar visibility
        featureLUTColorBar.Visibility = 0
        Delete(contour2)
        del contour2
        Delete(contour1)
        del contour1
        Delete(hyperTreeGridToDualGrid1)
        del hyperTreeGridToDualGrid1
