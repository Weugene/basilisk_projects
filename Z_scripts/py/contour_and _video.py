# trace generated using paraview version 5.8.0
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# find source
my_source = FindSource('tube_bp_0_0024.pvtu')

# create a new 'Slice'
slice1 = Slice(Input=my_source)
slice1.SliceType = 'Plane'
slice1.HyperTreeGridSlicer = 'Plane'
slice1.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Origin = [15.0, 0.0, 0.0]

# init the 'Plane' selected for 'HyperTreeGridSlicer'
slice1.HyperTreeGridSlicer.Origin = [15.0, 0.0, 0.0]

# Properties modified on slice1.SliceType
slice1.SliceType.Normal = [0.0, 0.0, 1.0]

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
renderView1.ViewSize = [2120, 1134]

# get layout
layout1 = GetLayout()

# show data in view
slice1Display = Show(slice1, renderView1, 'GeometryRepresentation')

# get color transfer function/color map for 'l'
lLUT = GetColorTransferFunction('l')
lLUT.RGBPoints = [5.0, 0.23137254902, 0.298039215686, 0.752941176471, 7.5, 0.865, 0.865, 0.865, 10.0, 0.705882352941, 0.0156862745098, 0.149019607843]
lLUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
slice1Display.Representation = 'Surface'
slice1Display.ColorArrayName = ['CELLS', 'l']
slice1Display.LookupTable = lLUT
slice1Display.OSPRayScaleArray = 'f'
slice1Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice1Display.SelectOrientationVectors = 'None'
slice1Display.ScaleFactor = 3.0
slice1Display.SelectScaleArray = 'None'
slice1Display.GlyphType = 'Arrow'
slice1Display.GlyphTableIndexArray = 'None'
slice1Display.GaussianRadius = 0.15
slice1Display.SetScaleArray = ['POINTS', 'f']
slice1Display.ScaleTransferFunction = 'PiecewiseFunction'
slice1Display.OpacityArray = ['POINTS', 'f']
slice1Display.OpacityTransferFunction = 'PiecewiseFunction'
slice1Display.DataAxesGrid = 'GridAxesRepresentation'
slice1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
slice1Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# show color bar/color legend
slice1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# get opacity transfer function/opacity map for 'l'
lPWF = GetOpacityTransferFunction('l')
lPWF.Points = [5.0, 0.0, 0.5, 0.0, 10.0, 1.0, 0.5, 0.0]
lPWF.ScalarRangeInitialized = 1

# set active source
SetActiveSource(my_source)

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=slice1.SliceType)

# create a new 'Contour'
contour1 = Contour(Input=my_source)
contour1.ContourBy = ['POINTS', 'f']
contour1.Isosurfaces = [0.5]
contour1.PointMergeMethod = 'Uniform Binning'

# show data in view
contour1Display = Show(contour1, renderView1, 'GeometryRepresentation')

# get color transfer function/color map for 'f'
fLUT = GetColorTransferFunction('f')
fLUT.RGBPoints = [0.0, 0.23137254902, 0.298039215686, 0.752941176471, 0.5, 0.865, 0.865, 0.865, 1.0, 0.705882352941, 0.0156862745098, 0.149019607843]
fLUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
contour1Display.Representation = 'Surface'
contour1Display.ColorArrayName = ['POINTS', 'f']
contour1Display.LookupTable = fLUT
contour1Display.OSPRayScaleArray = 'f'
contour1Display.OSPRayScaleFunction = 'PiecewiseFunction'
contour1Display.SelectOrientationVectors = 'None'
contour1Display.ScaleFactor = 0.26233772646968995
contour1Display.SelectScaleArray = 'f'
contour1Display.GlyphType = 'Arrow'
contour1Display.GlyphTableIndexArray = 'f'
contour1Display.GaussianRadius = 0.013116886323484497
contour1Display.SetScaleArray = ['POINTS', 'f']
contour1Display.ScaleTransferFunction = 'PiecewiseFunction'
contour1Display.OpacityArray = ['POINTS', 'f']
contour1Display.OpacityTransferFunction = 'PiecewiseFunction'
contour1Display.DataAxesGrid = 'GridAxesRepresentation'
contour1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
contour1Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
contour1Display.ScaleTransferFunction.Points = [0.5, 0.0, 0.5, 0.0, 0.5001220703125, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
contour1Display.OpacityTransferFunction.Points = [0.5, 0.0, 0.5, 0.0, 0.5001220703125, 1.0, 0.5, 0.0]

# show color bar/color legend
contour1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# get opacity transfer function/opacity map for 'f'
fPWF = GetOpacityTransferFunction('f')
fPWF.ScalarRangeInitialized = 1

# set active source
SetActiveSource(my_source)

# create a new 'Cylinder'
cylinder1 = Cylinder()

# Properties modified on cylinder1
cylinder1.Resolution = 60
cylinder1.Height = 30.0

# show data in view
cylinder1Display = Show(cylinder1, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
cylinder1Display.Representation = 'Surface'
cylinder1Display.ColorArrayName = [None, '']
cylinder1Display.OSPRayScaleArray = 'Normals'
cylinder1Display.OSPRayScaleFunction = 'PiecewiseFunction'
cylinder1Display.SelectOrientationVectors = 'None'
cylinder1Display.ScaleFactor = 3.0
cylinder1Display.SelectScaleArray = 'None'
cylinder1Display.GlyphType = 'Arrow'
cylinder1Display.GlyphTableIndexArray = 'None'
cylinder1Display.GaussianRadius = 0.15
cylinder1Display.SetScaleArray = ['POINTS', 'Normals']
cylinder1Display.ScaleTransferFunction = 'PiecewiseFunction'
cylinder1Display.OpacityArray = ['POINTS', 'Normals']
cylinder1Display.OpacityTransferFunction = 'PiecewiseFunction'
cylinder1Display.DataAxesGrid = 'GridAxesRepresentation'
cylinder1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
cylinder1Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
cylinder1Display.ScaleTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
cylinder1Display.OpacityTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on cylinder1
cylinder1.Capping = 0

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on cylinder1Display
cylinder1Display.Position = [15.0, 0.0, 0.0]

# Properties modified on cylinder1Display.DataAxesGrid
cylinder1Display.DataAxesGrid.Position = [15.0, 0.0, 0.0]

# Properties modified on cylinder1Display.PolarAxes
cylinder1Display.PolarAxes.Translation = [15.0, 0.0, 0.0]

# Properties modified on cylinder1Display
cylinder1Display.Orientation = [0.0, 0.0, 90.0]

# Properties modified on cylinder1Display.PolarAxes
cylinder1Display.PolarAxes.Orientation = [0.0, 0.0, 90.0]

# reset view to fit data bounds
renderView1.ResetCamera(0.0, 30.0, -0.5, 0.5, -0.5, 0.5)

# hide data in view
Hide(my_source, renderView1)

# hide data in view
Hide(slice1, renderView1)

# reset view to fit data bounds
renderView1.ResetCamera(0.0, 30.0, -0.5, 0.5, -0.5, 0.5)

# reset view to fit data
renderView1.ResetCamera()

# reset view to fit data
renderView1.ResetCamera()

# reset view to fit data bounds
renderView1.ResetCamera(0.0, 30.0, -0.5, 0.5, -0.5, 0.5)

# reset view to fit data bounds
renderView1.ResetCamera(0.0, 30.0, -0.5, 0.5, -0.5, 0.5)

# reset view to fit data
renderView1.ResetCamera()

# Properties modified on cylinder1Display
cylinder1Display.Opacity = 0.23

# set active source
SetActiveSource(contour1)

# turn off scalar coloring
ColorBy(contour1Display, None)

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(fLUT, renderView1)

# change solid color
contour1Display.AmbientColor = [0.027450980392156862, 0.07450980392156863, 1.0]
contour1Display.DiffuseColor = [0.027450980392156862, 0.07450980392156863, 1.0]

# Properties modified on contour1Display
contour1Display.Specular = 0.36

# Properties modified on contour1Display
contour1Display.Specular = 0.35

# Properties modified on contour1Display
contour1Display.Specular = 0.75

# Properties modified on contour1Display
contour1Display.Specular = 0.83

# Properties modified on contour1Display
contour1Display.Specular = 0.9

# set active source
SetActiveSource(my_source)

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [8.925681298959862, -0.30823854884808866, 58.01990889709437]
renderView1.CameraFocalPoint = [8.925681298959862, -0.30823854884808866, 0.0]
renderView1.CameraParallelScale = 3.1415802998357254

# save animation
SaveAnimation('/Users/weugene/basilisk/work/tube/bubble_test.avi', renderView1, ImageResolution=[2120, 780],
    FrameRate=10,
    FrameWindow=[0, 9])

# hide data in view
Hide(cylinder1, renderView1)

# hide data in view
Hide(contour1, renderView1)

# set active source
SetActiveSource(slice1)

# toggle 3D widget visibility (only when running from the GUI)
Show3DWidgets(proxy=slice1.SliceType)

# show data in view
slice1Display = Show(slice1, renderView1, 'GeometryRepresentation')

# show color bar/color legend
slice1Display.SetScalarBarVisibility(renderView1, True)

# reset view to fit data
renderView1.ResetCamera()

# create a new 'Contour'
contour2 = Contour(Input=slice1)
contour2.ContourBy = ['POINTS', 'f']
contour2.Isosurfaces = [0.5]
contour2.PointMergeMethod = 'Uniform Binning'

# show data in view
contour2Display = Show(contour2, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
contour2Display.Representation = 'Surface'
contour2Display.ColorArrayName = ['POINTS', 'f']
contour2Display.LookupTable = fLUT
contour2Display.OSPRayScaleArray = 'f'
contour2Display.OSPRayScaleFunction = 'PiecewiseFunction'
contour2Display.SelectOrientationVectors = 'None'
contour2Display.ScaleFactor = 0.2622632263474312
contour2Display.SelectScaleArray = 'f'
contour2Display.GlyphType = 'Arrow'
contour2Display.GlyphTableIndexArray = 'f'
contour2Display.GaussianRadius = 0.01311316131737156
contour2Display.SetScaleArray = ['POINTS', 'f']
contour2Display.ScaleTransferFunction = 'PiecewiseFunction'
contour2Display.OpacityArray = ['POINTS', 'f']
contour2Display.OpacityTransferFunction = 'PiecewiseFunction'
contour2Display.DataAxesGrid = 'GridAxesRepresentation'
contour2Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
contour2Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
contour2Display.ScaleTransferFunction.Points = [0.5, 0.0, 0.5, 0.0, 0.5001220703125, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
contour2Display.OpacityTransferFunction.Points = [0.5, 0.0, 0.5, 0.0, 0.5001220703125, 1.0, 0.5, 0.0]

# hide data in view
Hide(slice1, renderView1)

# show color bar/color legend
contour2Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# reset view to fit data bounds
renderView1.ResetCamera(4.30888399172, 6.93151625519, -0.458760922495, 0.457681992641, 0.0, 0.0)

# hide color bar/color legend
contour2Display.SetScalarBarVisibility(renderView1, False)

CreateLayout('Layout #2')

# set active view
SetActiveView(None)

# Create a new 'SpreadSheet View'
spreadSheetView1 = CreateView('SpreadSheetView')
spreadSheetView1.ColumnToSort = ''
spreadSheetView1.BlockSize = 1024L
# uncomment following to set a specific view size
# spreadSheetView1.ViewSize = [400, 400]

# show data in view
contour2Display_1 = Show(contour2, spreadSheetView1, 'SpreadSheetRepresentation')

# get layout
layout2 = GetLayoutByName("Layout #2")

# assign view to a particular cell in the layout
AssignViewToLayout(view=spreadSheetView1, layout=layout2, hint=0)

# export view
ExportView('/Users/weugene/basilisk/work/tube/contour_data.csv', view=spreadSheetView1)

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [5.620200123456031, -0.0005394649265466622, 5.366956348088599]
renderView1.CameraFocalPoint = [5.620200123456031, -0.0005394649265466622, 0.0]
renderView1.CameraParallelScale = 0.9487538536433321

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).