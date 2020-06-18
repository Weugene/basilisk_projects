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

# create a new 'Threshold'
threshold1 = Threshold(Input=porous_BP_0)
threshold1.Scalars = ['CELLS', 'f']
threshold1.ThresholdRange = [0.0, 1.0]

# Properties modified on threshold1
threshold1.Scalars = ['CELLS', 'fs']
threshold1.ThresholdRange = [0.0, 0.99]

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1783, 1012]

# show data in view
threshold1Display = Show(threshold1, renderView1)

# trace defaults for the display properties.
threshold1Display.Representation = 'Surface'
threshold1Display.ColorArrayName = [None, '']
threshold1Display.OSPRayScaleFunction = 'PiecewiseFunction'
threshold1Display.SelectOrientationVectors = 'None'
threshold1Display.ScaleFactor = 0.09375
threshold1Display.SelectScaleArray = 'None'
threshold1Display.GlyphType = 'Arrow'
threshold1Display.GlyphTableIndexArray = 'None'
threshold1Display.GaussianRadius = 0.0046875
threshold1Display.SetScaleArray = [None, '']
threshold1Display.ScaleTransferFunction = 'PiecewiseFunction'
threshold1Display.OpacityArray = [None, '']
threshold1Display.OpacityTransferFunction = 'PiecewiseFunction'
threshold1Display.DataAxesGrid = 'GridAxesRepresentation'
threshold1Display.SelectionCellLabelFontFile = ''
threshold1Display.SelectionPointLabelFontFile = ''
threshold1Display.PolarAxes = 'PolarAxesRepresentation'
threshold1Display.ScalarOpacityUnitDistance = 0.03465944255666103

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
threshold1Display.DataAxesGrid.XTitleFontFile = ''
threshold1Display.DataAxesGrid.YTitleFontFile = ''
threshold1Display.DataAxesGrid.ZTitleFontFile = ''
threshold1Display.DataAxesGrid.XLabelFontFile = ''
threshold1Display.DataAxesGrid.YLabelFontFile = ''
threshold1Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
threshold1Display.PolarAxes.PolarAxisTitleFontFile = ''
threshold1Display.PolarAxes.PolarAxisLabelFontFile = ''
threshold1Display.PolarAxes.LastRadialAxisTextFontFile = ''
threshold1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# hide data in view
Hide(porous_BP_0, renderView1)

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(threshold1Display, ('CELLS', 'f'))

# rescale color and/or opacity maps used to include current data range
threshold1Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
threshold1Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'f'
fLUT = GetColorTransferFunction('f')
fLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'f'
fPWF = GetOpacityTransferFunction('f')
fPWF.ScalarRangeInitialized = 1

# Properties modified on fLUT
fLUT.RGBPoints = [0.0, 1.0, 1.0, 1.0, 0.5, 0.865003, 0.865003, 0.865003, 1.0, 0.705882, 0.0156863, 0.14902]

# Properties modified on fLUT
fLUT.RGBPoints = [0.0, 1.0, 1.0, 1.0, 0.5, 0.865003, 0.865003, 0.865003, 1.0, 0.0, 0.09803921568627451, 1.0]

# hide color bar/color legend
threshold1Display.SetScalarBarVisibility(renderView1, False)

# Hide orientation axes
renderView1.OrientationAxesVisibility = 0

# get the material library
materialLibrary1 = GetMaterialLibrary()

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.0, 0.0, 10000.0]
renderView1.CameraParallelScale = 0.5478616589771803

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).