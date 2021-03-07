# trace generated using paraview version 5.6.2
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# find source
porous_BP_0 = FindSource('porous_BP_0*')

# create a new 'Calculator'
calculator1 = Calculator(Input=porous_BP_0)
calculator1.Function = ''

# Properties modified on calculator1
calculator1.AttributeType = 'Cell Data'
calculator1.ResultArrayName = 'Solid&L1&L2'
calculator1.Function = 'f+2*fs'

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1783, 1012]

# show data in view
calculator1Display = Show(calculator1, renderView1)

# get color transfer function/color map for 'SolidL1L2'
solidL1L2LUT = GetColorTransferFunction('SolidL1L2')
solidL1L2LUT.RGBPoints = [0.0, 0.231373, 0.298039, 0.752941, 1.0, 0.865003, 0.865003, 0.865003, 2.0, 0.705882, 0.0156863, 0.14902]
solidL1L2LUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'SolidL1L2'
solidL1L2PWF = GetOpacityTransferFunction('SolidL1L2')
solidL1L2PWF.Points = [0.0, 0.0, 0.5, 0.0, 2.0, 1.0, 0.5, 0.0]
solidL1L2PWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
calculator1Display.Representation = 'Surface'
calculator1Display.ColorArrayName = ['CELLS', 'Solid&L1&L2']
calculator1Display.LookupTable = solidL1L2LUT
calculator1Display.OSPRayScaleFunction = 'PiecewiseFunction'
calculator1Display.SelectOrientationVectors = 'None'
calculator1Display.ScaleFactor = 0.09375
calculator1Display.SelectScaleArray = 'Solid&L1&L2'
calculator1Display.GlyphType = 'Arrow'
calculator1Display.GlyphTableIndexArray = 'Solid&L1&L2'
calculator1Display.GaussianRadius = 0.0046875
calculator1Display.SetScaleArray = [None, '']
calculator1Display.ScaleTransferFunction = 'PiecewiseFunction'
calculator1Display.OpacityArray = [None, '']
calculator1Display.OpacityTransferFunction = 'PiecewiseFunction'
calculator1Display.DataAxesGrid = 'GridAxesRepresentation'
calculator1Display.SelectionCellLabelFontFile = ''
calculator1Display.SelectionPointLabelFontFile = ''
calculator1Display.PolarAxes = 'PolarAxesRepresentation'
calculator1Display.ScalarOpacityFunction = solidL1L2PWF
calculator1Display.ScalarOpacityUnitDistance = 0.033594269600637446

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
calculator1Display.DataAxesGrid.XTitleFontFile = ''
calculator1Display.DataAxesGrid.YTitleFontFile = ''
calculator1Display.DataAxesGrid.ZTitleFontFile = ''
calculator1Display.DataAxesGrid.XLabelFontFile = ''
calculator1Display.DataAxesGrid.YLabelFontFile = ''
calculator1Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
calculator1Display.PolarAxes.PolarAxisTitleFontFile = ''
calculator1Display.PolarAxes.PolarAxisLabelFontFile = ''
calculator1Display.PolarAxes.LastRadialAxisTextFontFile = ''
calculator1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# hide data in view
Hide(porous_BP_0, renderView1)

# show color bar/color legend
calculator1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# hide color bar/color legend
calculator1Display.SetScalarBarVisibility(renderView1, False)

# Properties modified on solidL1L2LUT
solidL1L2LUT.RGBPoints = [0.0, 1.0, 1.0, 1.0, 1.0, 0.865003, 0.865003, 0.865003, 2.0, 0.705882, 0.0156863, 0.14902]

# Properties modified on solidL1L2LUT
solidL1L2LUT.RGBPoints = [0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.14901960784313725, 1.0, 2.0, 0.705882, 0.0156863, 0.14902]

# Properties modified on solidL1L2LUT
solidL1L2LUT.RGBPoints = [0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.14901960784313725, 1.0, 2.0, 0.23529411764705882, 0.00392156862745098, 0.047058823529411764]

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.0, 0.0, 10000.0]
renderView1.CameraParallelScale = 0.5478616589771803

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).