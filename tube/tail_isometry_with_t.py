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
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [597, 489]
renderView1.InteractionMode = '2D'
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.CenterOfRotation = [-0.21472561038020577, 0.00010505318641662598, -4.470348358154297e-08]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [-1.7608661667636267, -2.0910945450148724, -0.1012756454524226]
renderView1.CameraFocalPoint = [-0.19566413834533514, -0.014517077240671314, 0.004977968400705676]
renderView1.CameraViewUp = [-0.7889291360994015, 0.6011603861785902, -0.12726746757541868]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 0.46007255351813975
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

# create a new 'PVD Reader'
save_isosurface_tailpvd = PVDReader(FileName='/home/e.sharaborin/wbasilisk/tube/res22_adaptffsu_new/dumps/save_isosurface_tail.pvd')
save_isosurface_tailpvd.CellArrays = ['vtkGhostType']
save_isosurface_tailpvd.PointArrays = ['u.x', 'absOmega', 'f', 'fs', 'l', 'l2', 'omega', 'p', 'vtkValidPointMask', 'vtkGhostType', 'Result']

# create a new 'Annotate Time'
annotateTime1 = AnnotateTime()
annotateTime1.Format = 't=%f'

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from save_isosurface_tailpvd
save_isosurface_tailpvdDisplay = Show(save_isosurface_tailpvd, renderView1, 'GeometryRepresentation')

# get color transfer function/color map for 'ux'
uxLUT = GetColorTransferFunction('ux')
uxLUT.RGBPoints = [0.0, 0.0, 0.0, 0.5625, 0.2777775000000001, 0.0, 0.0, 1.0, 0.9126987500000001, 0.0, 1.0, 1.0, 1.23015875, 0.5, 1.0, 0.5, 1.5476187500000005, 1.0, 1.0, 0.0, 2.18254, 1.0, 0.0, 0.0, 2.5, 0.5, 0.0, 0.0]
uxLUT.ColorSpace = 'RGB'
uxLUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
save_isosurface_tailpvdDisplay.Representation = 'Surface'
save_isosurface_tailpvdDisplay.ColorArrayName = ['POINTS', 'u.x']
save_isosurface_tailpvdDisplay.LookupTable = uxLUT
save_isosurface_tailpvdDisplay.OSPRayScaleArray = 'Result'
save_isosurface_tailpvdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
save_isosurface_tailpvdDisplay.SelectOrientationVectors = 'None'
save_isosurface_tailpvdDisplay.ScaleFactor = 0.09029629230499268
save_isosurface_tailpvdDisplay.SelectScaleArray = 'Result'
save_isosurface_tailpvdDisplay.GlyphType = 'Arrow'
save_isosurface_tailpvdDisplay.GlyphTableIndexArray = 'Result'
save_isosurface_tailpvdDisplay.GaussianRadius = 0.004514814615249634
save_isosurface_tailpvdDisplay.SetScaleArray = ['POINTS', 'Result']
save_isosurface_tailpvdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
save_isosurface_tailpvdDisplay.OpacityArray = ['POINTS', 'Result']
save_isosurface_tailpvdDisplay.OpacityTransferFunction = 'PiecewiseFunction'
save_isosurface_tailpvdDisplay.DataAxesGrid = 'GridAxesRepresentation'
save_isosurface_tailpvdDisplay.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
save_isosurface_tailpvdDisplay.ScaleTransferFunction.Points = [0.005237827916105199, 0.0, 0.5, 0.0, 0.4522067015627967, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
save_isosurface_tailpvdDisplay.OpacityTransferFunction.Points = [0.005237827916105199, 0.0, 0.5, 0.0, 0.4522067015627967, 1.0, 0.5, 0.0]

# show data from annotateTime1
annotateTime1Display = Show(annotateTime1, renderView1, 'TextSourceRepresentation')

# trace defaults for the display properties.
annotateTime1Display.FontFamily = 'Times'
annotateTime1Display.FontSize = 25

# setup the color legend parameters for each legend in this view

# get color legend/bar for uxLUT in view renderView1
uxLUTColorBar = GetScalarBar(uxLUT, renderView1)
uxLUTColorBar.WindowLocation = 'AnyLocation'
uxLUTColorBar.Position = [0.8589743589743588, 0.6319018404907977]
uxLUTColorBar.Title = '|u|'
uxLUTColorBar.ComponentTitle = ''
uxLUTColorBar.RangeLabelFormat = '%-#6.2g'
uxLUTColorBar.ScalarBarLength = 0.32999999999999996

# set color bar visibility
uxLUTColorBar.Visibility = 1

# show color legend
save_isosurface_tailpvdDisplay.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get opacity transfer function/opacity map for 'ux'
uxPWF = GetOpacityTransferFunction('ux')
uxPWF.Points = [0.0, 0.0, 0.5, 0.0, 2.5, 1.0, 0.5, 0.0]
uxPWF.ScalarRangeInitialized = 1

# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(annotateTime1)
# ----------------------------------------------------------------