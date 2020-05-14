# trace generated using paraview version 5.6.2
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# find source
porous_BP_0 = GetActiveSource() #FindSource('rk.pvd')

# create a new 'Calculator'
calculator2 = Calculator(Input=porous_BP_0)
calculator2.Function = ''

# Properties modified on calculator2
calculator2.AttributeType = 'Cell Data'
calculator2.ResultArrayName = 'Solid&L1&L2'
calculator2.Function = 'f+fs'

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1225, 1012]

# show data in view
calculator2Display = Show(calculator2, renderView1)

# get color transfer function/color map for 'SolidL1L2'
solidL1L2LUT = GetColorTransferFunction('SolidL1L2')
solidL1L2LUT.RGBPoints = [0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.14901960784313725, 1.0, 2.0, 0.23529411764705882, 0.00392156862745098, 0.047058823529411764]
solidL1L2LUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'SolidL1L2'
solidL1L2PWF = GetOpacityTransferFunction('SolidL1L2')
solidL1L2PWF.Points = [0.0, 0.0, 0.5, 0.0, 2.0, 1.0, 0.5, 0.0]
solidL1L2PWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
calculator2Display.Representation = 'Surface'
calculator2Display.ColorArrayName = ['CELLS', 'Solid&L1&L2']
calculator2Display.LookupTable = solidL1L2LUT
calculator2Display.OSPRayScaleFunction = 'PiecewiseFunction'
calculator2Display.SelectOrientationVectors = 'None'
calculator2Display.ScaleFactor = 0.09375
calculator2Display.SelectScaleArray = 'Solid&L1&L2'
calculator2Display.GlyphType = 'Arrow'
calculator2Display.GlyphTableIndexArray = 'Solid&L1&L2'
calculator2Display.GaussianRadius = 0.0046875
calculator2Display.SetScaleArray = [None, '']
calculator2Display.ScaleTransferFunction = 'PiecewiseFunction'
calculator2Display.OpacityArray = [None, '']
calculator2Display.OpacityTransferFunction = 'PiecewiseFunction'
calculator2Display.DataAxesGrid = 'GridAxesRepresentation'
calculator2Display.SelectionCellLabelFontFile = ''
calculator2Display.SelectionPointLabelFontFile = ''
calculator2Display.PolarAxes = 'PolarAxesRepresentation'
calculator2Display.ScalarOpacityFunction = solidL1L2PWF
calculator2Display.ScalarOpacityUnitDistance = 0.033594269600637446

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
calculator2Display.DataAxesGrid.XTitleFontFile = ''
calculator2Display.DataAxesGrid.YTitleFontFile = ''
calculator2Display.DataAxesGrid.ZTitleFontFile = ''
calculator2Display.DataAxesGrid.XLabelFontFile = ''
calculator2Display.DataAxesGrid.YLabelFontFile = ''
calculator2Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
calculator2Display.PolarAxes.PolarAxisTitleFontFile = ''
calculator2Display.PolarAxes.PolarAxisLabelFontFile = ''
calculator2Display.PolarAxes.LastRadialAxisTextFontFile = ''
calculator2Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# hide data in view
Hide(porous_BP_0, renderView1)

# show color bar/color legend
calculator2Display.SetScalarBarVisibility(renderView1, True)

# find source
calculator1 = FindSource('Calculator1')

# update the view to ensure updated data information
renderView1.Update()

# hide data in view
# Hide(calculator1, renderView1)

# Hide orientation axes
renderView1.OrientationAxesVisibility = 0

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# uncomment following to set a specific view size
# renderView1.ViewSize = [1152, 1150]

# hide data in view
Hide(porous_BP_0, renderView1)

# find source
calculator1 = FindSource('Calculator1')

# set active source
SetActiveSource(calculator1)

# Properties modified on solidL1L2LUT
solidL1L2LUT.RGBPoints = [0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 2.0, 0.23529411764705882, 0.00392156862745098, 0.047058823529411764]

# Properties modified on solidL1L2LUT
solidL1L2LUT.RGBPoints = [0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 2.0, 0.2196078431372549, 0.0, 0.043137254901960784]

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.0, 0.0, 10000.0]
renderView1.CameraParallelScale = 0.5478616589771803

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).