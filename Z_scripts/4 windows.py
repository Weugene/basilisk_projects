# trace generated using paraview version 5.6.0
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# find source
rk_ = FindSource('win1')

# set active source
SetActiveSource(rk_)

# find source
rk_0 = FindSource('win2')

# set active source
SetActiveSource(rk_0)

# find source
rk_0_1 = FindSource('win3')

# set active source
SetActiveSource(rk_0_1)

# find source
rk_0_2 = FindSource('win4')

# set active source
SetActiveSource(rk_0_2)

# set active source
SetActiveSource(rk_)

# create a new 'Calculator'
calculator1 = Calculator(Input=None)
calculator1.Function = ''

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [2150, 1366]

# show data in view
calculator1Display = Show(calculator1, renderView1)

# get color transfer function/color map for 'SolidL1L2'
solidL1L2LUT = GetColorTransferFunction('SolidL1L2')
solidL1L2LUT.RGBPoints = [0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 2.0, 0.2196078431372549, 0.0, 0.043137254901960784]
solidL1L2LUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'SolidL1L2'
solidL1L2PWF = GetOpacityTransferFunction('SolidL1L2')
solidL1L2PWF.Points = [0.0, 0.0, 0.5, 0.0, 2.0, 1.0, 0.5, 0.0]
solidL1L2PWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
calculator1Display.Representation = 'Surface'
calculator1Display.ColorArrayName = ['CELLS', 'Solid&L1&L2']
calculator1Display.LookupTable = solidL1L2LUT
calculator1Display.Opacity = 0.85
calculator1Display.OSPRayScaleArray = 'Solid&L1&L2'
calculator1Display.OSPRayScaleFunction = 'PiecewiseFunction'
calculator1Display.SelectOrientationVectors = 'None'
calculator1Display.ScaleFactor = 0.09843750000000001
calculator1Display.SelectScaleArray = 'Solid&L1&L2'
calculator1Display.GlyphType = 'Arrow'
calculator1Display.GlyphTableIndexArray = 'Solid&L1&L2'
calculator1Display.GaussianRadius = 0.004921875
calculator1Display.SetScaleArray = ['POINTS', 'Solid&L1&L2']
calculator1Display.ScaleTransferFunction = 'PiecewiseFunction'
calculator1Display.OpacityArray = ['POINTS', 'Solid&L1&L2']
calculator1Display.OpacityTransferFunction = 'PiecewiseFunction'
calculator1Display.DataAxesGrid = 'GridAxesRepresentation'
calculator1Display.SelectionCellLabelFontFile = ''
calculator1Display.SelectionPointLabelFontFile = ''
calculator1Display.PolarAxes = 'PolarAxesRepresentation'
calculator1Display.ScalarOpacityFunction = solidL1L2PWF
calculator1Display.ScalarOpacityUnitDistance = 0.05599592453661358

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
calculator1Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
calculator1Display.ScaleTransferFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
calculator1Display.OpacityTransferFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

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
Hide(calculator1, renderView1)

# show color bar/color legend
calculator1Display.SetScalarBarVisibility(renderView1, True)

# hide data in view
Hide(calculator1, renderView1)

# set active source
SetActiveSource(rk_0)

# create a new 'Calculator'
calculator2 = Calculator(Input=None)
calculator2.Function = ''

# show data in view
calculator2Display = Show(calculator2, renderView1)

# trace defaults for the display properties.
calculator2Display.Representation = 'Surface'
calculator2Display.ColorArrayName = ['CELLS', 'Solid&L1&L2']
calculator2Display.LookupTable = solidL1L2LUT
calculator2Display.Opacity = 0.85
calculator2Display.OSPRayScaleArray = 'Solid&L1&L2'
calculator2Display.OSPRayScaleFunction = 'PiecewiseFunction'
calculator2Display.SelectOrientationVectors = 'None'
calculator2Display.ScaleFactor = 0.09843750000000001
calculator2Display.SelectScaleArray = 'Solid&L1&L2'
calculator2Display.GlyphType = 'Arrow'
calculator2Display.GlyphTableIndexArray = 'Solid&L1&L2'
calculator2Display.GaussianRadius = 0.004921875
calculator2Display.SetScaleArray = ['POINTS', 'Solid&L1&L2']
calculator2Display.ScaleTransferFunction = 'PiecewiseFunction'
calculator2Display.OpacityArray = ['POINTS', 'Solid&L1&L2']
calculator2Display.OpacityTransferFunction = 'PiecewiseFunction'
calculator2Display.DataAxesGrid = 'GridAxesRepresentation'
calculator2Display.SelectionCellLabelFontFile = ''
calculator2Display.SelectionPointLabelFontFile = ''
calculator2Display.PolarAxes = 'PolarAxesRepresentation'
calculator2Display.ScalarOpacityFunction = solidL1L2PWF
calculator2Display.ScalarOpacityUnitDistance = 0.05599592453661358

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
calculator2Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
calculator2Display.ScaleTransferFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
calculator2Display.OpacityTransferFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

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
Hide(calculator2, renderView1)

# show color bar/color legend
calculator2Display.SetScalarBarVisibility(renderView1, True)

# hide data in view
Hide(calculator1, renderView1)

# set active source
SetActiveSource(rk_0_1)

# create a new 'Calculator'
calculator3 = Calculator(Input=None)
calculator3.Function = ''

# show data in view
calculator3Display = Show(calculator3, renderView1)

# trace defaults for the display properties.
calculator3Display.Representation = 'Surface'
calculator3Display.ColorArrayName = ['CELLS', 'Solid&L1&L2']
calculator3Display.LookupTable = solidL1L2LUT
calculator3Display.Opacity = 0.85
calculator3Display.OSPRayScaleArray = 'Solid&L1&L2'
calculator3Display.OSPRayScaleFunction = 'PiecewiseFunction'
calculator3Display.SelectOrientationVectors = 'None'
calculator3Display.ScaleFactor = 0.09843750000000001
calculator3Display.SelectScaleArray = 'Solid&L1&L2'
calculator3Display.GlyphType = 'Arrow'
calculator3Display.GlyphTableIndexArray = 'Solid&L1&L2'
calculator3Display.GaussianRadius = 0.004921875
calculator3Display.SetScaleArray = ['POINTS', 'Solid&L1&L2']
calculator3Display.ScaleTransferFunction = 'PiecewiseFunction'
calculator3Display.OpacityArray = ['POINTS', 'Solid&L1&L2']
calculator3Display.OpacityTransferFunction = 'PiecewiseFunction'
calculator3Display.DataAxesGrid = 'GridAxesRepresentation'
calculator3Display.SelectionCellLabelFontFile = ''
calculator3Display.SelectionPointLabelFontFile = ''
calculator3Display.PolarAxes = 'PolarAxesRepresentation'
calculator3Display.ScalarOpacityFunction = solidL1L2PWF
calculator3Display.ScalarOpacityUnitDistance = 0.05599592453661358

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
calculator3Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
calculator3Display.ScaleTransferFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
calculator3Display.OpacityTransferFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
calculator3Display.DataAxesGrid.XTitleFontFile = ''
calculator3Display.DataAxesGrid.YTitleFontFile = ''
calculator3Display.DataAxesGrid.ZTitleFontFile = ''
calculator3Display.DataAxesGrid.XLabelFontFile = ''
calculator3Display.DataAxesGrid.YLabelFontFile = ''
calculator3Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
calculator3Display.PolarAxes.PolarAxisTitleFontFile = ''
calculator3Display.PolarAxes.PolarAxisLabelFontFile = ''
calculator3Display.PolarAxes.LastRadialAxisTextFontFile = ''
calculator3Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# hide data in view
Hide(calculator3, renderView1)

# show color bar/color legend
calculator3Display.SetScalarBarVisibility(renderView1, True)

# hide data in view
Hide(calculator1, renderView1)

# set active source
SetActiveSource(rk_0_2)

# create a new 'Calculator'
calculator4 = Calculator(Input=None)
calculator4.Function = ''

# show data in view
calculator4Display = Show(calculator4, renderView1)

# trace defaults for the display properties.
calculator4Display.Representation = 'Surface'
calculator4Display.ColorArrayName = ['CELLS', 'Solid&L1&L2']
calculator4Display.LookupTable = solidL1L2LUT
calculator4Display.Opacity = 0.85
calculator4Display.OSPRayScaleArray = 'Solid&L1&L2'
calculator4Display.OSPRayScaleFunction = 'PiecewiseFunction'
calculator4Display.SelectOrientationVectors = 'None'
calculator4Display.ScaleFactor = 0.09843750000000001
calculator4Display.SelectScaleArray = 'Solid&L1&L2'
calculator4Display.GlyphType = 'Arrow'
calculator4Display.GlyphTableIndexArray = 'Solid&L1&L2'
calculator4Display.GaussianRadius = 0.004921875
calculator4Display.SetScaleArray = ['POINTS', 'Solid&L1&L2']
calculator4Display.ScaleTransferFunction = 'PiecewiseFunction'
calculator4Display.OpacityArray = ['POINTS', 'Solid&L1&L2']
calculator4Display.OpacityTransferFunction = 'PiecewiseFunction'
calculator4Display.DataAxesGrid = 'GridAxesRepresentation'
calculator4Display.SelectionCellLabelFontFile = ''
calculator4Display.SelectionPointLabelFontFile = ''
calculator4Display.PolarAxes = 'PolarAxesRepresentation'
calculator4Display.ScalarOpacityFunction = solidL1L2PWF
calculator4Display.ScalarOpacityUnitDistance = 0.05599592453661358

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
calculator4Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
calculator4Display.ScaleTransferFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
calculator4Display.OpacityTransferFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
calculator4Display.DataAxesGrid.XTitleFontFile = ''
calculator4Display.DataAxesGrid.YTitleFontFile = ''
calculator4Display.DataAxesGrid.ZTitleFontFile = ''
calculator4Display.DataAxesGrid.XLabelFontFile = ''
calculator4Display.DataAxesGrid.YLabelFontFile = ''
calculator4Display.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
calculator4Display.PolarAxes.PolarAxisTitleFontFile = ''
calculator4Display.PolarAxes.PolarAxisLabelFontFile = ''
calculator4Display.PolarAxes.LastRadialAxisTextFontFile = ''
calculator4Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# hide data in view
Hide(calculator4, renderView1)

# show color bar/color legend
calculator4Display.SetScalarBarVisibility(renderView1, True)

# hide data in view
Hide(calculator1, renderView1)

# set active source
SetActiveSource(calculator1)

# set active source
SetActiveSource(calculator2)

# Properties modified on calculator2Display
calculator2Display.Position = [1.0, 0.0, 0.0]

# Properties modified on calculator2Display.DataAxesGrid
calculator2Display.DataAxesGrid.Position = [1.0, 0.0, 0.0]

# Properties modified on calculator2Display.PolarAxes
calculator2Display.PolarAxes.Translation = [1.0, 0.0, 0.0]

# set active source
SetActiveSource(calculator2)

# show data in view
calculator2Display = Show(calculator2, renderView1)

# show color bar/color legend
calculator2Display.SetScalarBarVisibility(renderView1, True)

# hide data in view
Hide(rk_0_2, renderView1)

# hide data in view
Hide(rk_0_1, renderView1)

# hide data in view
Hide(rk_0, renderView1)

# hide data in view
Hide(rk_, renderView1)

# set active source
SetActiveSource(calculator1)

# show data in view
calculator1Display = Show(calculator1, renderView1)

# show color bar/color legend
calculator1Display.SetScalarBarVisibility(renderView1, True)

# hide data in view
Hide(calculator1, renderView1)

# show data in view
calculator1Display = Show(calculator1, renderView1)

# show color bar/color legend
calculator1Display.SetScalarBarVisibility(renderView1, True)

# set active source
SetActiveSource(calculator3)

# Properties modified on calculator3Display
calculator3Display.Position = [0.0, -1.0, 0.0]

# Properties modified on calculator3Display.DataAxesGrid
calculator3Display.DataAxesGrid.Position = [0.0, -1.0, 0.0]

# Properties modified on calculator3Display.PolarAxes
calculator3Display.PolarAxes.Translation = [0.0, -1.0, 0.0]

# set active source
SetActiveSource(calculator3)

# show data in view
calculator3Display = Show(calculator3, renderView1)

# show color bar/color legend
calculator3Display.SetScalarBarVisibility(renderView1, True)

# set active source
SetActiveSource(calculator4)

# Properties modified on calculator4Display
calculator4Display.Position = [1.0, 0.0, 0.0]

# Properties modified on calculator4Display.DataAxesGrid
calculator4Display.DataAxesGrid.Position = [1.0, 0.0, 0.0]

# Properties modified on calculator4Display.PolarAxes
calculator4Display.PolarAxes.Translation = [1.0, 0.0, 0.0]

# Properties modified on calculator4Display
calculator4Display.Position = [1.0, -1.0, 0.0]

# Properties modified on calculator4Display.DataAxesGrid
calculator4Display.DataAxesGrid.Position = [1.0, -1.0, 0.0]

# Properties modified on calculator4Display.PolarAxes
calculator4Display.PolarAxes.Translation = [1.0, -1.0, 0.0]

# set active source
SetActiveSource(calculator4)

# show data in view
calculator4Display = Show(calculator4, renderView1)

# show color bar/color legend
calculator4Display.SetScalarBarVisibility(renderView1, True)

# hide color bar/color legend
calculator4Display.SetScalarBarVisibility(renderView1, False)

# Hide orientation axes
renderView1.OrientationAxesVisibility = 0

# get the material library
materialLibrary1 = GetMaterialLibrary()

# create a new 'Text'
text1 = Text()

# Properties modified on text1
text1.Text = 'Ca=1000 Rb/Rs=1 Rrho=1 Rmu=1 maxlevel=10'

# show data in view
text1Display = Show(text1, renderView1)

# trace defaults for the display properties.
text1Display.FontFile = ''

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on text1Display
text1Display.WindowLocation = 'UpperCenter'

# Properties modified on text1Display
text1Display.FontSize = 17

# Properties modified on text1Display
text1Display.FontSize = 16

# Properties modified on text1Display
text1Display.FontSize = 15

# Properties modified on text1Display
text1Display.FontSize = 14

# Properties modified on text1Display
text1Display.FontSize = 13

# Properties modified on text1Display
text1Display.FontSize = 12

# Properties modified on text1Display
text1Display.FontSize = 11

# Properties modified on text1Display
text1Display.FontSize = 10

# Properties modified on text1Display
text1Display.FontSize = 9

# Properties modified on text1Display
text1Display.FontSize = 8

# Properties modified on text1Display
text1Display.FontSize = 7

# Properties modified on text1Display
text1Display.FontSize = 6

# create a new 'Text'
text2 = Text()

# Properties modified on text2
text2.Text = 'Re=1'

# show data in view
text2Display = Show(text2, renderView1)

# trace defaults for the display properties.
text2Display.FontFile = ''

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(text1)

# set active source
SetActiveSource(text2)

# Properties modified on text2Display
text2Display.FontSize = 6

# Properties modified on text2Display
text2Display.WindowLocation = 'AnyLocation'

# create a new 'Text'
text3 = Text()

# set active source
SetActiveSource(text2)

# set active source
SetActiveSource(text3)

# show data in view
text3Display = Show(text3, renderView1)

# trace defaults for the display properties.
text3Display.FontFile = ''

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(text2)

# set active source
SetActiveSource(text3)

# Properties modified on text3Display
text3Display.FontSize = 6

# Properties modified on text3
text3.Text = 'Re=5'

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on text3Display
text3Display.WindowLocation = 'AnyLocation'
text3Display.Position = [0.01, 0.952281]

# create a new 'Text'
text4 = Text()

# Properties modified on text4
text4.Text = 'Re=100'

# show data in view
text4Display = Show(text4, renderView1)

# trace defaults for the display properties.
text4Display.FontFile = ''
text4Display.FontSize = 6

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on text4Display
text4Display.WindowLocation = 'AnyLocation'

# create a new 'Text'
text5 = Text()
text5.Text = 'Re=100'

# Properties modified on text5
text5.Text = 'Re=1000'

# show data in view
text5Display = Show(text5, renderView1)

# trace defaults for the display properties.
text5Display.FontFile = ''
text5Display.FontSize = 6

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on text5Display
text5Display.WindowLocation = 'AnyLocation'

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.4717849925357097, -0.461364534376025, 10000.000001073184]
renderView1.CameraFocalPoint = [0.4639724925357097, -0.492614534376025, 1.0731821445065543e-06]
renderView1.CameraViewUp = [-2.441406249986595e-12, 0.9999999999951172, -3.1249999999828333e-06]
renderView1.CameraParallelScale = 1.0799618091580863

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).