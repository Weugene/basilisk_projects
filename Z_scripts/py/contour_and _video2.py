# trace generated using paraview version 5.8.0
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# find source
tube_bppvd = FindSource('tube_bp.pvd')

# create a new 'Slice'
slice1 = Slice(Input=tube_bppvd)
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
# renderView1.ViewSize = [2150, 466]

# get layout
layout1 = GetLayout()

# show data in view
slice1Display = Show(slice1, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
slice1Display.Representation = 'Surface'
slice1Display.ColorArrayName = [None, '']
slice1Display.OSPRayScaleArray = 'fs'
slice1Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice1Display.SelectOrientationVectors = 'None'
slice1Display.ScaleFactor = 3.0
slice1Display.SelectScaleArray = 'None'
slice1Display.GlyphType = 'Arrow'
slice1Display.GlyphTableIndexArray = 'None'
slice1Display.GaussianRadius = 0.15
slice1Display.SetScaleArray = ['POINTS', 'fs']
slice1Display.ScaleTransferFunction = 'PiecewiseFunction'
slice1Display.OpacityArray = ['POINTS', 'fs']
slice1Display.OpacityTransferFunction = 'PiecewiseFunction'
slice1Display.DataAxesGrid = 'GridAxesRepresentation'
slice1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
slice1Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# hide data in view
Hide(tube_bppvd, renderView1)

# set active source
SetActiveSource(tube_bppvd)

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=slice1.SliceType)

# create a new 'Contour'
contour1 = Contour(Input=tube_bppvd)
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
contour1Display.ScaleFactor = 0.336715683006052
contour1Display.SelectScaleArray = 'f'
contour1Display.GlyphType = 'Arrow'
contour1Display.GlyphTableIndexArray = 'f'
contour1Display.GaussianRadius = 0.016835784150302596
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

# create a new 'Cylinder'
cylinder1 = Cylinder()

# Properties modified on cylinder1
cylinder1.Resolution = 55
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
cylinder1Display.ScaleTransferFunction.Points = [-0.9983690977096558, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
cylinder1Display.OpacityTransferFunction.Points = [-0.9983690977096558, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# hide data in view
Hide(slice1, renderView1)

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

# Properties modified on cylinder1Display
cylinder1Display.Specular = 0.82

# Properties modified on cylinder1Display
cylinder1Display.Specular = 0.49

# Properties modified on cylinder1Display
cylinder1Display.Ambient = 0.14

# Properties modified on cylinder1Display
cylinder1Display.Opacity = 0.19

# set active source
SetActiveSource(contour1)

# turn off scalar coloring
ColorBy(contour1Display, None)

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(fLUT, renderView1)

# change solid color
contour1Display.AmbientColor = [0.00784313725490196, 0.00784313725490196, 1.0]
contour1Display.DiffuseColor = [0.00784313725490196, 0.00784313725490196, 1.0]

# Properties modified on contour1Display
contour1Display.Specular = 0.8

# Properties modified on contour1Display
contour1Display.Specular = 0.87

# Properties modified on contour1Display
contour1Display.Specular = 0.89

# Properties modified on contour1Display
contour1Display.Specular = 1.0

# get animation scene
animationScene1 = GetAnimationScene()

# get the time-keeper
timeKeeper1 = GetTimeKeeper()

animationScene1.GoToLast()

# set active source
SetActiveSource(slice1)

# toggle 3D widget visibility (only when running from the GUI)
Show3DWidgets(proxy=slice1.SliceType)

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [9.193714171767704, 0.009979232717995515, 58.01980396939221]
renderView1.CameraFocalPoint = [9.193714171767704, 0.009979232717995515, 0.0]
renderView1.CameraParallelScale = 1.9680297154517084

# save animation
SaveAnimation('/Users/weugene/basilisk/work/tube/bubble_test.avi', renderView1, ImageResolution=[2148, 464],
    FrameRate=10,
    FrameWindow=[0, 101])

# hide data in view
Hide(cylinder1, renderView1)

# hide data in view
Hide(contour1, renderView1)

# create a new 'Contour'
contour2 = Contour(Input=slice1)
contour2.ContourBy = ['POINTS', 'fs']
contour2.Isosurfaces = [0.5]
contour2.PointMergeMethod = 'Uniform Binning'

# show data in view
contour2Display = Show(contour2, renderView1, 'GeometryRepresentation')

# get color transfer function/color map for 'fs'
fsLUT = GetColorTransferFunction('fs')
fsLUT.RGBPoints = [0.0, 0.23137254902, 0.298039215686, 0.752941176471, 0.5, 0.865, 0.865, 0.865, 1.0, 0.705882352941, 0.0156862745098, 0.149019607843]
fsLUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
contour2Display.Representation = 'Surface'
contour2Display.ColorArrayName = ['POINTS', 'fs']
contour2Display.LookupTable = fsLUT
contour2Display.OSPRayScaleArray = 'fs'
contour2Display.OSPRayScaleFunction = 'PiecewiseFunction'
contour2Display.SelectOrientationVectors = 'None'
contour2Display.ScaleFactor = 3.0
contour2Display.SelectScaleArray = 'fs'
contour2Display.GlyphType = 'Arrow'
contour2Display.GlyphTableIndexArray = 'fs'
contour2Display.GaussianRadius = 0.15
contour2Display.SetScaleArray = ['POINTS', 'fs']
contour2Display.ScaleTransferFunction = 'PiecewiseFunction'
contour2Display.OpacityArray = ['POINTS', 'fs']
contour2Display.OpacityTransferFunction = 'PiecewiseFunction'
contour2Display.DataAxesGrid = 'GridAxesRepresentation'
contour2Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
contour2Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
contour2Display.ScaleTransferFunction.Points = [0.5, 0.0, 0.5, 0.0, 0.5001220703125, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
contour2Display.OpacityTransferFunction.Points = [0.5, 0.0, 0.5, 0.0, 0.5001220703125, 1.0, 0.5, 0.0]

# reset view to fit data
renderView1.ResetCamera()

#changing interaction mode based on data extents
renderView1.CameraPosition = [15.0, 0.0, 10000.0]
renderView1.CameraFocalPoint = [15.0, 0.0, 0.0]

# get the material library
materialLibrary1 = GetMaterialLibrary()

# hide data in view
Hide(slice1, renderView1)

# show color bar/color legend
contour2Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# get opacity transfer function/opacity map for 'fs'
fsPWF = GetOpacityTransferFunction('fs')
fsPWF.ScalarRangeInitialized = 1

# Properties modified on contour2
contour2.ContourBy = ['POINTS', 'f']

# update the view to ensure updated data information
renderView1.Update()

# reset view to fit data bounds
renderView1.ResetCamera(14.4578813613, 17.838004069, -0.442199892112, 0.445704945132, 0.0, 0.0)

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
renderView1.CameraPosition = [16.147942715187742, 0.0017525265100293286, 6.751428930936986]
renderView1.CameraFocalPoint = [16.147942715187742, 0.0017525265100293286, 0.0]
renderView1.CameraParallelScale = 1.7473983889826432

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).