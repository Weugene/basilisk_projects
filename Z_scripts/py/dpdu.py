# trace generated using paraview version 5.7.0
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# find source
vortex_BP2pvd = FindSource('vortex_BP2.pvd')

# create a new 'Calculator'
calculator3 = Calculator(Input=vortex_BP2pvd)
calculator3.AttributeType = 'Cell Data'
calculator3.Function = ''

# find view
renderView1 = FindViewOrCreate('RenderView1', viewtype='RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [2112, 1366]

# set active view
SetActiveView(renderView1)

# find source
calculator2 = FindSource('Calculator2')

# find source
calculator1 = FindSource('Calculator1')

# find source
calculator1_dp = FindSource('Calculator1_dp')

# find source
plotOverLine1 = FindSource('PlotOverLine1')

# find source
vortex_Popinetpvd = FindSource('vortex_Popinet.pvd')

# find source
threshold2 = FindSource('Threshold2')

# Properties modified on calculator3
calculator3.ResultArrayName = 'Resultu'
calculator3.Function = 'mag(u.x-uexact.x)'

# show data in view
calculator3Display = Show(calculator3, renderView1)

# get color transfer function/color map for 'Resultu'
resultuLUT = GetColorTransferFunction('Resultu')

# get opacity transfer function/opacity map for 'Resultu'
resultuPWF = GetOpacityTransferFunction('Resultu')

# trace defaults for the display properties.
calculator3Display.Representation = 'Surface'
calculator3Display.ColorArrayName = ['CELLS', 'Resultu']
calculator3Display.LookupTable = resultuLUT
calculator3Display.Opacity = 0.85
calculator3Display.OSPRayScaleArray = 'Resultu'
calculator3Display.OSPRayScaleFunction = 'PiecewiseFunction'
calculator3Display.SelectOrientationVectors = 'None'
calculator3Display.ScaleFactor = 0.30000000000000004
calculator3Display.SelectScaleArray = 'Resultu'
calculator3Display.GlyphType = 'Arrow'
calculator3Display.GlyphTableIndexArray = 'Resultu'
calculator3Display.GaussianRadius = 0.015
calculator3Display.SetScaleArray = ['POINTS', 'Resultu']
calculator3Display.ScaleTransferFunction = 'PiecewiseFunction'
calculator3Display.OpacityArray = ['POINTS', 'Resultu']
calculator3Display.OpacityTransferFunction = 'PiecewiseFunction'
calculator3Display.DataAxesGrid = 'GridAxesRepresentation'
calculator3Display.PolarAxes = 'PolarAxesRepresentation'
calculator3Display.ScalarOpacityFunction = resultuPWF
calculator3Display.ScalarOpacityUnitDistance = 0.10535958712446956

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
calculator3Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
calculator3Display.ScaleTransferFunction.Points = [3.984961880906345e-07, 0.0, 0.5, 0.0, 0.006474512918947373, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
calculator3Display.OpacityTransferFunction.Points = [3.984961880906345e-07, 0.0, 0.5, 0.0, 0.006474512918947373, 1.0, 0.5, 0.0]

# show color bar/color legend
calculator3Display.SetScalarBarVisibility(renderView1, True)

# find source
threshold1 = FindSource('Threshold1')

# find source
calculator2_du = FindSource('Calculator2_du')

# find source
vortex_BPpvd = FindSource('vortex_BP.pvd')

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Calculator'
calculator4 = Calculator(Input=calculator3)
calculator4.AttributeType = 'Cell Data'
calculator4.Function = ''

# Properties modified on calculator4
calculator4.ResultArrayName = 'Resultp'
calculator4.Function = 'abs(p-pexact)'

# show data in view
calculator4Display = Show(calculator4, renderView1)

# get color transfer function/color map for 'Resultp'
resultpLUT = GetColorTransferFunction('Resultp')

# get opacity transfer function/opacity map for 'Resultp'
resultpPWF = GetOpacityTransferFunction('Resultp')

# trace defaults for the display properties.
calculator4Display.Representation = 'Surface'
calculator4Display.ColorArrayName = ['CELLS', 'Resultp']
calculator4Display.LookupTable = resultpLUT
calculator4Display.Opacity = 0.85
calculator4Display.OSPRayScaleArray = 'Resultp'
calculator4Display.OSPRayScaleFunction = 'PiecewiseFunction'
calculator4Display.SelectOrientationVectors = 'None'
calculator4Display.ScaleFactor = 0.30000000000000004
calculator4Display.SelectScaleArray = 'Resultp'
calculator4Display.GlyphType = 'Arrow'
calculator4Display.GlyphTableIndexArray = 'Resultp'
calculator4Display.GaussianRadius = 0.015
calculator4Display.SetScaleArray = ['POINTS', 'Resultp']
calculator4Display.ScaleTransferFunction = 'PiecewiseFunction'
calculator4Display.OpacityArray = ['POINTS', 'Resultp']
calculator4Display.OpacityTransferFunction = 'PiecewiseFunction'
calculator4Display.DataAxesGrid = 'GridAxesRepresentation'
calculator4Display.PolarAxes = 'PolarAxesRepresentation'
calculator4Display.ScalarOpacityFunction = resultpPWF
calculator4Display.ScalarOpacityUnitDistance = 0.10535958712446956

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
calculator4Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
calculator4Display.ScaleTransferFunction.Points = [0.00034722218441315955, 0.0, 0.5, 0.0, 35.07826541445242, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
calculator4Display.OpacityTransferFunction.Points = [0.00034722218441315955, 0.0, 0.5, 0.0, 35.07826541445242, 1.0, 0.5, 0.0]

# hide data in view
Hide(calculator3, renderView1)

# show color bar/color legend
calculator4Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Threshold'
threshold3 = Threshold(Input=calculator4)
threshold3.Scalars = ['POINTS', 'Resultp']
threshold3.ThresholdRange = [0.00034722218441315955, 35.07826541445242]

# set active source
SetActiveSource(calculator4)

# destroy threshold3
Delete(threshold3)
del threshold3

# create a new 'Plot Over Line'
plotOverLine2 = PlotOverLine(Input=calculator4,
    Source='High Resolution Line Source')

# init the 'High Resolution Line Source' selected for 'Source'
plotOverLine2.Source.Point1 = [-1.5, -1.5, 0.0]
plotOverLine2.Source.Point2 = [1.5, 1.5, 0.0]

# show data in view
plotOverLine2Display = Show(plotOverLine2, renderView1)

# trace defaults for the display properties.
plotOverLine2Display.Representation = 'Surface'
plotOverLine2Display.ColorArrayName = [None, '']
plotOverLine2Display.OSPRayScaleArray = 'Resultp'
plotOverLine2Display.OSPRayScaleFunction = 'PiecewiseFunction'
plotOverLine2Display.SelectOrientationVectors = 'None'
plotOverLine2Display.ScaleFactor = 0.30000000000000004
plotOverLine2Display.SelectScaleArray = 'None'
plotOverLine2Display.GlyphType = 'Arrow'
plotOverLine2Display.GlyphTableIndexArray = 'None'
plotOverLine2Display.GaussianRadius = 0.015
plotOverLine2Display.SetScaleArray = ['POINTS', 'Resultp']
plotOverLine2Display.ScaleTransferFunction = 'PiecewiseFunction'
plotOverLine2Display.OpacityArray = ['POINTS', 'Resultp']
plotOverLine2Display.OpacityTransferFunction = 'PiecewiseFunction'
plotOverLine2Display.DataAxesGrid = 'GridAxesRepresentation'
plotOverLine2Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
plotOverLine2Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
plotOverLine2Display.ScaleTransferFunction.Points = [0.026302532885149013, 0.0, 0.5, 0.0, 35.07826541445235, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
plotOverLine2Display.OpacityTransferFunction.Points = [0.026302532885149013, 0.0, 0.5, 0.0, 35.07826541445235, 1.0, 0.5, 0.0]

# Create a new 'Line Chart View'
lineChartView2 = CreateView('XYChartView')
# uncomment following to set a specific view size
# lineChartView2.ViewSize = [400, 400]

# show data in view
plotOverLine2Display_1 = Show(plotOverLine2, lineChartView2)

# trace defaults for the display properties.
plotOverLine2Display_1.CompositeDataSetIndex = [0]
plotOverLine2Display_1.UseIndexForXAxis = 0
plotOverLine2Display_1.XArrayName = 'arc_length'
plotOverLine2Display_1.SeriesVisibility = ['fs', 'l', 'omega', 'p', 'pexact', 'Resultp', 'Resultu', 'u.x_Magnitude', 'uexact.x_Magnitude']
plotOverLine2Display_1.SeriesLabel = ['arc_length', 'arc_length', 'fs', 'fs', 'l', 'l', 'omega', 'omega', 'p', 'p', 'pexact', 'pexact', 'Resultp', 'Resultp', 'Resultu', 'Resultu', 'u.x_X', 'u.x_X', 'u.x_Y', 'u.x_Y', 'u.x_Z', 'u.x_Z', 'u.x_Magnitude', 'u.x_Magnitude', 'uexact.x_X', 'uexact.x_X', 'uexact.x_Y', 'uexact.x_Y', 'uexact.x_Z', 'uexact.x_Z', 'uexact.x_Magnitude', 'uexact.x_Magnitude', 'vtkValidPointMask', 'vtkValidPointMask', 'Points_X', 'Points_X', 'Points_Y', 'Points_Y', 'Points_Z', 'Points_Z', 'Points_Magnitude', 'Points_Magnitude']
plotOverLine2Display_1.SeriesColor = ['arc_length', '0', '0', '0', 'fs', '0.89', '0.1', '0.11', 'l', '0.22', '0.49', '0.72', 'omega', '0.3', '0.69', '0.29', 'p', '0.6', '0.31', '0.64', 'pexact', '1', '0.5', '0', 'Resultp', '0.65', '0.34', '0.16', 'Resultu', '0', '0', '0', 'u.x_X', '0.89', '0.1', '0.11', 'u.x_Y', '0.22', '0.49', '0.72', 'u.x_Z', '0.3', '0.69', '0.29', 'u.x_Magnitude', '0.6', '0.31', '0.64', 'uexact.x_X', '1', '0.5', '0', 'uexact.x_Y', '0.65', '0.34', '0.16', 'uexact.x_Z', '0', '0', '0', 'uexact.x_Magnitude', '0.89', '0.1', '0.11', 'vtkValidPointMask', '0.22', '0.49', '0.72', 'Points_X', '0.3', '0.69', '0.29', 'Points_Y', '0.6', '0.31', '0.64', 'Points_Z', '1', '0.5', '0', 'Points_Magnitude', '0.65', '0.34', '0.16']
plotOverLine2Display_1.SeriesPlotCorner = ['arc_length', '0', 'fs', '0', 'l', '0', 'omega', '0', 'p', '0', 'pexact', '0', 'Resultp', '0', 'Resultu', '0', 'u.x_X', '0', 'u.x_Y', '0', 'u.x_Z', '0', 'u.x_Magnitude', '0', 'uexact.x_X', '0', 'uexact.x_Y', '0', 'uexact.x_Z', '0', 'uexact.x_Magnitude', '0', 'vtkValidPointMask', '0', 'Points_X', '0', 'Points_Y', '0', 'Points_Z', '0', 'Points_Magnitude', '0']
plotOverLine2Display_1.SeriesLabelPrefix = ''
plotOverLine2Display_1.SeriesLineStyle = ['arc_length', '1', 'fs', '1', 'l', '1', 'omega', '1', 'p', '1', 'pexact', '1', 'Resultp', '1', 'Resultu', '1', 'u.x_X', '1', 'u.x_Y', '1', 'u.x_Z', '1', 'u.x_Magnitude', '1', 'uexact.x_X', '1', 'uexact.x_Y', '1', 'uexact.x_Z', '1', 'uexact.x_Magnitude', '1', 'vtkValidPointMask', '1', 'Points_X', '1', 'Points_Y', '1', 'Points_Z', '1', 'Points_Magnitude', '1']
plotOverLine2Display_1.SeriesLineThickness = ['arc_length', '2', 'fs', '2', 'l', '2', 'omega', '2', 'p', '2', 'pexact', '2', 'Resultp', '2', 'Resultu', '2', 'u.x_X', '2', 'u.x_Y', '2', 'u.x_Z', '2', 'u.x_Magnitude', '2', 'uexact.x_X', '2', 'uexact.x_Y', '2', 'uexact.x_Z', '2', 'uexact.x_Magnitude', '2', 'vtkValidPointMask', '2', 'Points_X', '2', 'Points_Y', '2', 'Points_Z', '2', 'Points_Magnitude', '2']
plotOverLine2Display_1.SeriesMarkerStyle = ['arc_length', '0', 'fs', '0', 'l', '0', 'omega', '0', 'p', '0', 'pexact', '0', 'Resultp', '0', 'Resultu', '0', 'u.x_X', '0', 'u.x_Y', '0', 'u.x_Z', '0', 'u.x_Magnitude', '0', 'uexact.x_X', '0', 'uexact.x_Y', '0', 'uexact.x_Z', '0', 'uexact.x_Magnitude', '0', 'vtkValidPointMask', '0', 'Points_X', '0', 'Points_Y', '0', 'Points_Z', '0', 'Points_Magnitude', '0']

# get layout
layout1 = GetLayoutByName("Layout #1")

# add view to a layout so it's visible in UI
AssignViewToLayout(view=lineChartView2, layout=layout1, hint=0)

# destroy lineChartView2
Delete(lineChartView2)
del lineChartView2

# close an empty frame
layout1.Collapse(2)

# set active view
SetActiveView(renderView1)

CreateLayout('Layout #4')

# set active view
SetActiveView(None)

# set active source
SetActiveSource(plotOverLine2)

# Create a new 'Line Chart View'
lineChartView2 = CreateView('XYChartView')
# uncomment following to set a specific view size
# lineChartView2.ViewSize = [400, 400]

# show data in view
plotOverLine2Display_1 = Show(plotOverLine2, lineChartView2)

# trace defaults for the display properties.
plotOverLine2Display_1.CompositeDataSetIndex = [0]
plotOverLine2Display_1.UseIndexForXAxis = 0
plotOverLine2Display_1.XArrayName = 'arc_length'
plotOverLine2Display_1.SeriesVisibility = ['fs', 'l', 'omega', 'p', 'pexact', 'Resultp', 'Resultu', 'u.x_Magnitude', 'uexact.x_Magnitude']
plotOverLine2Display_1.SeriesLabel = ['arc_length', 'arc_length', 'fs', 'fs', 'l', 'l', 'omega', 'omega', 'p', 'p', 'pexact', 'pexact', 'Resultp', 'Resultp', 'Resultu', 'Resultu', 'u.x_X', 'u.x_X', 'u.x_Y', 'u.x_Y', 'u.x_Z', 'u.x_Z', 'u.x_Magnitude', 'u.x_Magnitude', 'uexact.x_X', 'uexact.x_X', 'uexact.x_Y', 'uexact.x_Y', 'uexact.x_Z', 'uexact.x_Z', 'uexact.x_Magnitude', 'uexact.x_Magnitude', 'vtkValidPointMask', 'vtkValidPointMask', 'Points_X', 'Points_X', 'Points_Y', 'Points_Y', 'Points_Z', 'Points_Z', 'Points_Magnitude', 'Points_Magnitude']
plotOverLine2Display_1.SeriesColor = ['arc_length', '0', '0', '0', 'fs', '0.89', '0.1', '0.11', 'l', '0.22', '0.49', '0.72', 'omega', '0.3', '0.69', '0.29', 'p', '0.6', '0.31', '0.64', 'pexact', '1', '0.5', '0', 'Resultp', '0.65', '0.34', '0.16', 'Resultu', '0', '0', '0', 'u.x_X', '0.89', '0.1', '0.11', 'u.x_Y', '0.22', '0.49', '0.72', 'u.x_Z', '0.3', '0.69', '0.29', 'u.x_Magnitude', '0.6', '0.31', '0.64', 'uexact.x_X', '1', '0.5', '0', 'uexact.x_Y', '0.65', '0.34', '0.16', 'uexact.x_Z', '0', '0', '0', 'uexact.x_Magnitude', '0.89', '0.1', '0.11', 'vtkValidPointMask', '0.22', '0.49', '0.72', 'Points_X', '0.3', '0.69', '0.29', 'Points_Y', '0.6', '0.31', '0.64', 'Points_Z', '1', '0.5', '0', 'Points_Magnitude', '0.65', '0.34', '0.16']
plotOverLine2Display_1.SeriesPlotCorner = ['arc_length', '0', 'fs', '0', 'l', '0', 'omega', '0', 'p', '0', 'pexact', '0', 'Resultp', '0', 'Resultu', '0', 'u.x_X', '0', 'u.x_Y', '0', 'u.x_Z', '0', 'u.x_Magnitude', '0', 'uexact.x_X', '0', 'uexact.x_Y', '0', 'uexact.x_Z', '0', 'uexact.x_Magnitude', '0', 'vtkValidPointMask', '0', 'Points_X', '0', 'Points_Y', '0', 'Points_Z', '0', 'Points_Magnitude', '0']
plotOverLine2Display_1.SeriesLabelPrefix = ''
plotOverLine2Display_1.SeriesLineStyle = ['arc_length', '1', 'fs', '1', 'l', '1', 'omega', '1', 'p', '1', 'pexact', '1', 'Resultp', '1', 'Resultu', '1', 'u.x_X', '1', 'u.x_Y', '1', 'u.x_Z', '1', 'u.x_Magnitude', '1', 'uexact.x_X', '1', 'uexact.x_Y', '1', 'uexact.x_Z', '1', 'uexact.x_Magnitude', '1', 'vtkValidPointMask', '1', 'Points_X', '1', 'Points_Y', '1', 'Points_Z', '1', 'Points_Magnitude', '1']
plotOverLine2Display_1.SeriesLineThickness = ['arc_length', '2', 'fs', '2', 'l', '2', 'omega', '2', 'p', '2', 'pexact', '2', 'Resultp', '2', 'Resultu', '2', 'u.x_X', '2', 'u.x_Y', '2', 'u.x_Z', '2', 'u.x_Magnitude', '2', 'uexact.x_X', '2', 'uexact.x_Y', '2', 'uexact.x_Z', '2', 'uexact.x_Magnitude', '2', 'vtkValidPointMask', '2', 'Points_X', '2', 'Points_Y', '2', 'Points_Z', '2', 'Points_Magnitude', '2']
plotOverLine2Display_1.SeriesMarkerStyle = ['arc_length', '0', 'fs', '0', 'l', '0', 'omega', '0', 'p', '0', 'pexact', '0', 'Resultp', '0', 'Resultu', '0', 'u.x_X', '0', 'u.x_Y', '0', 'u.x_Z', '0', 'u.x_Magnitude', '0', 'uexact.x_X', '0', 'uexact.x_Y', '0', 'uexact.x_Z', '0', 'uexact.x_Magnitude', '0', 'vtkValidPointMask', '0', 'Points_X', '0', 'Points_Y', '0', 'Points_Z', '0', 'Points_Magnitude', '0']

# get layout
layout4 = GetLayoutByName("Layout #4")

# add view to a layout so it's visible in UI
AssignViewToLayout(view=lineChartView2, layout=layout4, hint=0)

# Properties modified on plotOverLine2Display_1
plotOverLine2Display_1.SeriesVisibility = ['fs', 'omega', 'p', 'pexact', 'Resultp', 'Resultu', 'u.x_Magnitude', 'uexact.x_Magnitude']
plotOverLine2Display_1.SeriesColor = ['arc_length', '0', '0', '0', 'fs', '0.889998', '0.100008', '0.110002', 'l', '0.220005', '0.489998', '0.719997', 'omega', '0.300008', '0.689998', '0.289998', 'p', '0.6', '0.310002', '0.639994', 'pexact', '1', '0.500008', '0', 'Resultp', '0.650004', '0.340002', '0.160006', 'Resultu', '0', '0', '0', 'u.x_X', '0.889998', '0.100008', '0.110002', 'u.x_Y', '0.220005', '0.489998', '0.719997', 'u.x_Z', '0.300008', '0.689998', '0.289998', 'u.x_Magnitude', '0.6', '0.310002', '0.639994', 'uexact.x_X', '1', '0.500008', '0', 'uexact.x_Y', '0.650004', '0.340002', '0.160006', 'uexact.x_Z', '0', '0', '0', 'uexact.x_Magnitude', '0.889998', '0.100008', '0.110002', 'vtkValidPointMask', '0.220005', '0.489998', '0.719997', 'Points_X', '0.300008', '0.689998', '0.289998', 'Points_Y', '0.6', '0.310002', '0.639994', 'Points_Z', '1', '0.500008', '0', 'Points_Magnitude', '0.650004', '0.340002', '0.160006']
plotOverLine2Display_1.SeriesPlotCorner = ['Points_Magnitude', '0', 'Points_X', '0', 'Points_Y', '0', 'Points_Z', '0', 'Resultp', '0', 'Resultu', '0', 'arc_length', '0', 'fs', '0', 'l', '0', 'omega', '0', 'p', '0', 'pexact', '0', 'u.x_Magnitude', '0', 'u.x_X', '0', 'u.x_Y', '0', 'u.x_Z', '0', 'uexact.x_Magnitude', '0', 'uexact.x_X', '0', 'uexact.x_Y', '0', 'uexact.x_Z', '0', 'vtkValidPointMask', '0']
plotOverLine2Display_1.SeriesLineStyle = ['Points_Magnitude', '1', 'Points_X', '1', 'Points_Y', '1', 'Points_Z', '1', 'Resultp', '1', 'Resultu', '1', 'arc_length', '1', 'fs', '1', 'l', '1', 'omega', '1', 'p', '1', 'pexact', '1', 'u.x_Magnitude', '1', 'u.x_X', '1', 'u.x_Y', '1', 'u.x_Z', '1', 'uexact.x_Magnitude', '1', 'uexact.x_X', '1', 'uexact.x_Y', '1', 'uexact.x_Z', '1', 'vtkValidPointMask', '1']
plotOverLine2Display_1.SeriesLineThickness = ['Points_Magnitude', '2', 'Points_X', '2', 'Points_Y', '2', 'Points_Z', '2', 'Resultp', '2', 'Resultu', '2', 'arc_length', '2', 'fs', '2', 'l', '2', 'omega', '2', 'p', '2', 'pexact', '2', 'u.x_Magnitude', '2', 'u.x_X', '2', 'u.x_Y', '2', 'u.x_Z', '2', 'uexact.x_Magnitude', '2', 'uexact.x_X', '2', 'uexact.x_Y', '2', 'uexact.x_Z', '2', 'vtkValidPointMask', '2']
plotOverLine2Display_1.SeriesMarkerStyle = ['Points_Magnitude', '0', 'Points_X', '0', 'Points_Y', '0', 'Points_Z', '0', 'Resultp', '0', 'Resultu', '0', 'arc_length', '0', 'fs', '0', 'l', '0', 'omega', '0', 'p', '0', 'pexact', '0', 'u.x_Magnitude', '0', 'u.x_X', '0', 'u.x_Y', '0', 'u.x_Z', '0', 'uexact.x_Magnitude', '0', 'uexact.x_X', '0', 'uexact.x_Y', '0', 'uexact.x_Z', '0', 'vtkValidPointMask', '0']

# Properties modified on plotOverLine2Display_1
plotOverLine2Display_1.SeriesVisibility = ['fs', 'p', 'pexact', 'Resultp', 'Resultu', 'u.x_Magnitude', 'uexact.x_Magnitude']

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.3846166925942942, 0.2797012832023967, 5.506789909006018]
renderView1.CameraFocalPoint = [0.3846166925942942, 0.2797012832023967, 0.0]
renderView1.CameraParallelScale = 1.425262105829132

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).