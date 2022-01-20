# state file generated using paraview version 5.9.0

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()
import argparse
import glob, os, sys

import logging
# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument("-pvdfile", type=str, help="Provide the name for the pvd file, please",
                    nargs='?', default='')
parser.add_argument("-maxlevel", type=int, help="Provide the maximum level of refinement",
                    nargs='?', default=10)
parser.add_argument("-timeList", type=int, help="Provide the list of step, please",
                    nargs='+', default=None)

args = parser.parse_args()
print("args:", args)
filename = args.pvdfile
maxlevel = args.maxlevel
timeList = args.timeList
path = os.path.abspath(os.getcwd())

prefix = os.path.splitext(filename)[0]
print("Current PATH=" + path, ' prefix', prefix)

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1416, 1158]
renderView1.InteractionMode = '2D'
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [0, 0, 10000.0]
renderView1.CameraFocalPoint = [0.0, 0.0, 0.0]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 0.3838091046374134
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.OSPRayMaterialLibrary = materialLibrary1

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------
# get animation scene
animationScene1 = GetAnimationScene()
# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.AssignView(0, renderView1)
layout1.SetSize(1416, 1158)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'PVD Reader'
input = PVDReader(registrationName='myinput', FileName=os.path.join(path, filename))
input.CellArrays = ['f', 'fs', 'p', 'u.x']
timeKeeper1 = GetTimeKeeper()
timesteps = timeKeeper1.TimestepValues # 0, 0.1, 0.2 ...
print('all timesteps:', timesteps)
timeList = []
timeList[:] = timesteps[:]
# timeList = [0.2, 0.4, 1.2, 1.4, 1.5, 1.6, 1.65, 2.25, 15.6] + [0.05*p for p in range(0, 100)]+ [4+0.1*p for p in range(0, 100)]
timeList = [0.1*p for p in range(0, 200)]

timeList = sorted(timeList)
timeList = list(set(timeList))
# timeList = range(0, timesteps[-1], 0.1)
print('timeList:', timeList)

for time in timesteps:
    flag_miss = True
    for t in timeList:
        if abs(t - time) < 1e-3:
            flag_miss = False
            break
    if flag_miss:
        continue
    print("time:" + str(time))
    # Properties modified on animationScene1
    animationScene1.AnimationTime = time
    # Properties modified on timeKeeper1
    timeKeeper1.Time = time

    # create a new 'Cell Data to Point Data'
    cellDatatoPointData = CellDatatoPointData(registrationName='CellDatatoPointData', Input=input)
    cellDatatoPointData.CellDataArraytoprocess = ['f', 'fs', 'p', 'u.x']

    # create a new 'Resample To Image'
    resampleToImage = ResampleToImage(registrationName='ResampleToImage', Input=cellDatatoPointData)
    resampleToImage.UseInputBounds = 0
    N = int(2**maxlevel)
    resampleToImage.SamplingDimensions = [N, N, 1]
    resampleToImage.SamplingBounds = [-0.48, 0.48, -0.5, 0.5, 0.0, 0.0]

    # create a new 'Contour'
    contour3 = Contour(registrationName='Contour3', Input=resampleToImage)
    contour3.ContourBy = ['POINTS', 'f']
    contour3.Isosurfaces = [0.5]
    contour3.PointMergeMethod = 'Uniform Binning'

    # create a new 'Threshold'
    threshold4 = Threshold(registrationName='Threshold4', Input=resampleToImage)
    threshold4.Scalars = ['POINTS', 'fs']
    threshold4.ThresholdRange = [0.45, 1.0]

    # create a new 'Threshold'
    threshold3 = Threshold(registrationName='Threshold3', Input=resampleToImage)
    threshold3.Scalars = ['POINTS', 'fs']
    threshold3.ThresholdRange = [0.0, 0.5]

    # create a new 'Contour'
    contour4 = Contour(registrationName='Contour4', Input=resampleToImage)
    contour4.ContourBy = ['POINTS', 'fs']
    contour4.Isosurfaces = [0.5]
    contour4.PointMergeMethod = 'Uniform Binning'

    # get color transfer function/color map for 'f'
    fLUT = GetColorTransferFunction('f')
    fLUT.RGBPoints = [0.0, 1.0, 1.0, 1.0, 0.5, 1.0, 1.0, 1.0, 1.0, 0.4549019607843137, 0.807843137254902, 1.0]
    fLUT.ColorSpace = 'Step'
    fLUT.ScalarRangeInitialized = 1.0

    # get opacity transfer function/opacity map for 'f'
    fPWF = GetOpacityTransferFunction('f')
    fPWF.Points = [0.0, 0.0, 0.5, 0.0, 0.3800622890658208, 0.38717949390411377, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]
    fPWF.ScalarRangeInitialized = 1

    # show data from threshold3
    threshold3Display = Show(threshold3, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    threshold3Display.Representation = 'Surface'
    threshold3Display.ColorArrayName = ['POINTS', 'f']
    threshold3Display.LookupTable = fLUT
    threshold3Display.Ambient = 1.0
    threshold3Display.Diffuse = 0.0
    threshold3Display.SelectTCoordArray = 'None'
    threshold3Display.SelectNormalArray = 'None'
    threshold3Display.SelectTangentArray = 'None'
    threshold3Display.OSPRayScaleArray = 'f'
    threshold3Display.OSPRayScaleFunction = 'PiecewiseFunction'
    threshold3Display.SelectOrientationVectors = 'None'
    threshold3Display.ScaleFactor = 0.1
    threshold3Display.SelectScaleArray = 'None'
    threshold3Display.GlyphType = 'Arrow'
    threshold3Display.GlyphTableIndexArray = 'None'
    threshold3Display.GaussianRadius = 0.005
    threshold3Display.SetScaleArray = ['POINTS', 'f']
    threshold3Display.ScaleTransferFunction = 'PiecewiseFunction'
    threshold3Display.OpacityArray = ['POINTS', 'f']
    threshold3Display.OpacityTransferFunction = 'PiecewiseFunction'
    threshold3Display.DataAxesGrid = 'GridAxesRepresentation'
    threshold3Display.PolarAxes = 'PolarAxesRepresentation'
    threshold3Display.ScalarOpacityFunction = fPWF
    threshold3Display.ScalarOpacityUnitDistance = 0.021784972811574
    threshold3Display.OpacityArrayName = ['POINTS', 'f']

    # show data from threshold4
    threshold4Display = Show(threshold4, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    threshold4Display.Representation = 'Surface'
    threshold4Display.AmbientColor = [0.6196078431372549, 0.6196078431372549, 0.6196078431372549]
    threshold4Display.ColorArrayName = ['POINTS', '']
    threshold4Display.DiffuseColor = [0.6196078431372549, 0.6196078431372549, 0.6196078431372549]
    threshold4Display.Ambient = 1.0
    threshold4Display.Diffuse = 0.0
    threshold4Display.SelectTCoordArray = 'None'
    threshold4Display.SelectNormalArray = 'None'
    threshold4Display.SelectTangentArray = 'None'
    threshold4Display.OSPRayScaleArray = 'f'
    threshold4Display.OSPRayScaleFunction = 'PiecewiseFunction'
    threshold4Display.SelectOrientationVectors = 'None'
    threshold4Display.ScaleFactor = 0.012328767031431199
    threshold4Display.SelectScaleArray = 'None'
    threshold4Display.GlyphType = 'Arrow'
    threshold4Display.GlyphTableIndexArray = 'None'
    threshold4Display.GaussianRadius = 0.0006164383515715599
    threshold4Display.SetScaleArray = ['POINTS', 'f']
    threshold4Display.ScaleTransferFunction = 'PiecewiseFunction'
    threshold4Display.OpacityArray = ['POINTS', 'f']
    threshold4Display.OpacityTransferFunction = 'PiecewiseFunction'
    threshold4Display.DataAxesGrid = 'GridAxesRepresentation'
    threshold4Display.PolarAxes = 'PolarAxesRepresentation'
    threshold4Display.ScalarOpacityUnitDistance = 0.01173782462695837
    threshold4Display.OpacityArrayName = ['POINTS', 'f']

    # show data from contour3
    contour3Display = Show(contour3, renderView1, 'GeometryRepresentation') # Error occurs
    # trace defaults for the display properties.
    contour3Display.Representation = 'Surface'
    contour3Display.AmbientColor = [0.0, 0.0, 1.0]
    contour3Display.ColorArrayName = ['POINTS', '']
    contour3Display.DiffuseColor = [0.0, 0.0, 1.0]
    contour3Display.LineWidth = 3.0
    contour3Display.SelectTCoordArray = 'None'
    contour3Display.SelectNormalArray = 'None'
    contour3Display.SelectTangentArray = 'None'
    contour3Display.OSPRayScaleArray = 'f'
    contour3Display.OSPRayScaleFunction = 'PiecewiseFunction'
    contour3Display.SelectOrientationVectors = 'None'
    contour3Display.ScaleFactor = 0.025049885362386705
    contour3Display.SelectScaleArray = 'f'
    contour3Display.GlyphType = 'Arrow'
    contour3Display.GlyphTableIndexArray = 'f'
    contour3Display.GaussianRadius = 0.0012524942681193352
    contour3Display.SetScaleArray = ['POINTS', 'f']
    contour3Display.ScaleTransferFunction = 'PiecewiseFunction'
    contour3Display.OpacityArray = ['POINTS', 'f']
    contour3Display.OpacityTransferFunction = 'PiecewiseFunction'
    contour3Display.DataAxesGrid = 'GridAxesRepresentation'
    contour3Display.PolarAxes = 'PolarAxesRepresentation'
    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    contour3Display.ScaleTransferFunction.Points = [0.5, 0.0, 0.5, 0.0, 0.5001220703125, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    contour3Display.OpacityTransferFunction.Points = [0.5, 0.0, 0.5, 0.0, 0.5001220703125, 1.0, 0.5, 0.0]

    # show data from contour4
    contour4Display = Show(contour4, renderView1, 'GeometryRepresentation') #Error occurs

    # trace defaults for the display properties.
    contour4Display.Representation = 'Surface'
    contour4Display.AmbientColor = [0.0, 0.0, 0.0]
    contour4Display.ColorArrayName = ['POINTS', '']
    contour4Display.DiffuseColor = [0.0, 0.0, 0.0]
    contour4Display.LineWidth = 3.0
    contour4Display.SelectTCoordArray = 'None'
    contour4Display.SelectNormalArray = 'None'
    contour4Display.SelectTangentArray = 'None'
    contour4Display.OSPRayScaleArray = 'fs'
    contour4Display.OSPRayScaleFunction = 'PiecewiseFunction'
    contour4Display.SelectOrientationVectors = 'None'
    contour4Display.ScaleFactor = 0.012521313130855562
    contour4Display.SelectScaleArray = 'fs'
    contour4Display.GlyphType = 'Arrow'
    contour4Display.GlyphTableIndexArray = 'fs'
    contour4Display.GaussianRadius = 0.000626065656542778
    contour4Display.SetScaleArray = ['POINTS', 'fs']
    contour4Display.ScaleTransferFunction = 'PiecewiseFunction'
    contour4Display.OpacityArray = ['POINTS', 'fs']
    contour4Display.OpacityTransferFunction = 'PiecewiseFunction'
    contour4Display.DataAxesGrid = 'GridAxesRepresentation'
    contour4Display.PolarAxes = 'PolarAxesRepresentation'

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    contour4Display.ScaleTransferFunction.Points = [0.5, 0.0, 0.5, 0.0, 0.5001220703125, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    contour4Display.OpacityTransferFunction.Points = [0.5, 0.0, 0.5, 0.0, 0.5001220703125, 1.0, 0.5, 0.0]

    # ----------------------------------------------------------------
    # setup color maps and opacity mapes used in the visualization
    # note: the Get..() functions create a new object, if needed
    # ----------------------------------------------------------------

    # ----------------------------------------------------------------
    # restore active source
    SetActiveSource(input)
    # ----------------------------------------------------------------



    fn = path + "/" + prefix +  '_t=' + str(time) +'.png'
    SaveScreenshot( fn, renderView1,
                    TransparentBackground=0,
                    CompressionLevel='2' )
    print('File=' + fn + ' generated succesfully')
