#!/usr/bin/python
# trace generated using paraview version 5.7.0
#
# To ensure correct image size when batch processing, please search
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *
import glob, os, sys
import logging

logging.basicConfig(format='%(message)s')
log = logging.getLogger(__name__)

def eprint(var):
    log.warning(var)

#os.chdir("/haha")
path = os.path.abspath(os.getcwd())
print(path)

# Find files with *.pvtu extension
filenames = []
for file in glob.glob("*.pvtu"):
    filenames.append(path + "/" + file)
    #print(file)
eprint(filenames)


#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'XML Partitioned Unstructured Grid Reader'
rk_0 = XMLPartitionedUnstructuredGridReader(FileName=[filenames])
rk_0.CellArrayStatus = ['fs', 'f', 'omega', 'rhov', 'p', 'l', 'divu', 'my_kappa', 'u.x', 'g.x', 'a.x', 'dbp.x', 'total_rhs.x', 'mapped_data_lower.x', 'mapped_data_upper.x', 'fs_lower.x', 'fs_upper.x']

# get animation scene
animationScene1 = GetAnimationScene()

# get the time-keeper
timeKeeper1 = GetTimeKeeper()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [912, 1318]

# show data in view
rk_0Display = Show(rk_0, renderView1)

# trace defaults for the display properties.
rk_0Display.Representation = 'Surface'
rk_0Display.ColorArrayName = [None, '']
rk_0Display.Opacity = 0.85
rk_0Display.OSPRayScaleArray = 'a.x'
rk_0Display.OSPRayScaleFunction = 'PiecewiseFunction'
rk_0Display.SelectOrientationVectors = 'None'
rk_0Display.ScaleFactor = 0.09375
rk_0Display.SelectScaleArray = 'None'
rk_0Display.GlyphType = 'Arrow'
rk_0Display.GlyphTableIndexArray = 'None'
rk_0Display.GaussianRadius = 0.0046875
rk_0Display.SetScaleArray = ['POINTS', 'a.x']
rk_0Display.ScaleTransferFunction = 'PiecewiseFunction'
rk_0Display.OpacityArray = ['POINTS', 'a.x']
rk_0Display.OpacityTransferFunction = 'PiecewiseFunction'
rk_0Display.DataAxesGrid = 'GridAxesRepresentation'
rk_0Display.PolarAxes = 'PolarAxesRepresentation'
rk_0Display.ScalarOpacityUnitDistance = 0.050871044069041346

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
rk_0Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
rk_0Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
rk_0Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

# reset view to fit data
renderView1.ResetCamera()

# get the material library
materialLibrary1 = GetMaterialLibrary()

# update the view to ensure updated data information
renderView1.Update()

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [-0.03125, -0.03125, 3.6363339075232113]
renderView1.CameraFocalPoint = [-0.03125, -0.03125, 0.0]
renderView1.CameraParallelScale = 0.9580250180960831

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
