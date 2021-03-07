# trace generated using paraview version 5.8.0
#
# To ensure correct image size when batch processing, please search
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'PVD Reader'
tube_bppvd = PVDReader(FileName='/Users/weugene/basilisk/work/tube/tube_bp.pvd')

# get animation scene
animationScene1 = GetAnimationScene()

# get the time-keeper
timeKeeper1 = GetTimeKeeper()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [2120, 1134]

# get layout
layout1 = GetLayout()

# show data in view
tube_bppvdDisplay = Show(tube_bppvd, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
tube_bppvdDisplay.Representation = 'Surface'
tube_bppvdDisplay.ColorArrayName = [None, '']
tube_bppvdDisplay.Opacity = 0.85
tube_bppvdDisplay.OSPRayScaleArray = 'f'
tube_bppvdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
tube_bppvdDisplay.SelectOrientationVectors = 'None'
tube_bppvdDisplay.ScaleFactor = 3.0
tube_bppvdDisplay.SelectScaleArray = 'None'
tube_bppvdDisplay.GlyphType = 'Arrow'
tube_bppvdDisplay.GlyphTableIndexArray = 'None'
tube_bppvdDisplay.GaussianRadius = 0.15
tube_bppvdDisplay.SetScaleArray = ['POINTS', 'f']
tube_bppvdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
tube_bppvdDisplay.OpacityArray = ['POINTS', 'f']
tube_bppvdDisplay.OpacityTransferFunction = 'PiecewiseFunction'
tube_bppvdDisplay.DataAxesGrid = 'GridAxesRepresentation'
tube_bppvdDisplay.PolarAxes = 'PolarAxesRepresentation'
tube_bppvdDisplay.ScalarOpacityUnitDistance = 0.4857920340157689

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
tube_bppvdDisplay.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# reset view to fit data
renderView1.ResetCamera()

# get the material library
materialLibrary1 = GetMaterialLibrary()

# update the view to ensure updated data information
renderView1.Update()

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.CameraPosition = [15.0, 0.0, 100.38195644853695]
renderView1.CameraFocalPoint = [15.0, 0.0, 0.0]
renderView1.CameraParallelScale = 25.98076211353316

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).