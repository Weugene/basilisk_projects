# trace generated using paraview version 5.8.0
#
# To ensure correct image size when batch processing, please search
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()
path="/Volumes/GoogleDrive/My Drive/Skoltech/Dumps/tube/10e_coarse/"
# create a new 'PVD Reader'
tube_bp_from_dumppvd = PVDReader(FileName=path+'tube_bp_from_dump.pvd')
tube_bp_from_dumppvd.CellArrays = ['p', 'fs', 'f', 'l', 'residual_of_p', 'u.x']

# get animation scene
animationScene1 = GetAnimationScene()

# get the time-keeper
timeKeeper1 = GetTimeKeeper()
animationScene1.GoToLast()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [2120, 1134]

# get layout
layout1 = GetLayout()

# show data in view
tube_bp_from_dumppvdDisplay = Show(tube_bp_from_dumppvd, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
tube_bp_from_dumppvdDisplay.Representation = 'Surface'
tube_bp_from_dumppvdDisplay.ColorArrayName = [None, '']
tube_bp_from_dumppvdDisplay.Opacity = 0.85
tube_bp_from_dumppvdDisplay.OSPRayScaleArray = 'f'
tube_bp_from_dumppvdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
tube_bp_from_dumppvdDisplay.SelectOrientationVectors = 'None'
tube_bp_from_dumppvdDisplay.ScaleFactor = 3.0
tube_bp_from_dumppvdDisplay.SelectScaleArray = 'None'
tube_bp_from_dumppvdDisplay.GlyphType = 'Arrow'
tube_bp_from_dumppvdDisplay.GlyphTableIndexArray = 'None'
tube_bp_from_dumppvdDisplay.GaussianRadius = 0.15
tube_bp_from_dumppvdDisplay.SetScaleArray = ['POINTS', 'f']
tube_bp_from_dumppvdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
tube_bp_from_dumppvdDisplay.OpacityArray = ['POINTS', 'f']
tube_bp_from_dumppvdDisplay.OpacityTransferFunction = 'PiecewiseFunction'
tube_bp_from_dumppvdDisplay.DataAxesGrid = 'GridAxesRepresentation'
tube_bp_from_dumppvdDisplay.PolarAxes = 'PolarAxesRepresentation'
tube_bp_from_dumppvdDisplay.ScalarOpacityUnitDistance = 0.48577906707444607

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
tube_bp_from_dumppvdDisplay.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# reset view to fit data
renderView1.ResetCamera()

# get the material library
materialLibrary1 = GetMaterialLibrary()

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Cell Data to Point Data'
cellDatatoPointData1 = CellDatatoPointData(Input=tube_bp_from_dumppvd)
cellDatatoPointData1.CellDataArraytoprocess = ['f', 'fs', 'l', 'p', 'residual_of_p', 'u.x']

# show data in view
cellDatatoPointData1Display = Show(cellDatatoPointData1, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
cellDatatoPointData1Display.Representation = 'Surface'
cellDatatoPointData1Display.ColorArrayName = [None, '']
cellDatatoPointData1Display.Opacity = 0.85
cellDatatoPointData1Display.OSPRayScaleArray = 'f'
cellDatatoPointData1Display.OSPRayScaleFunction = 'PiecewiseFunction'
cellDatatoPointData1Display.SelectOrientationVectors = 'None'
cellDatatoPointData1Display.ScaleFactor = 3.0
cellDatatoPointData1Display.SelectScaleArray = 'None'
cellDatatoPointData1Display.GlyphType = 'Arrow'
cellDatatoPointData1Display.GlyphTableIndexArray = 'None'
cellDatatoPointData1Display.GaussianRadius = 0.15
cellDatatoPointData1Display.SetScaleArray = ['POINTS', 'f']
cellDatatoPointData1Display.ScaleTransferFunction = 'PiecewiseFunction'
cellDatatoPointData1Display.OpacityArray = ['POINTS', 'f']
cellDatatoPointData1Display.OpacityTransferFunction = 'PiecewiseFunction'
cellDatatoPointData1Display.DataAxesGrid = 'GridAxesRepresentation'
cellDatatoPointData1Display.PolarAxes = 'PolarAxesRepresentation'
cellDatatoPointData1Display.ScalarOpacityUnitDistance = 0.48577906707444607

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
cellDatatoPointData1Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# hide data in view
Hide(tube_bp_from_dumppvd, renderView1)

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Threshold'
threshold1 = Threshold(Input=cellDatatoPointData1)
threshold1.Scalars = ['POINTS', 'f']
threshold1.ThresholdRange = [0.0, 1.0]

# Properties modified on threshold1
threshold1.Scalars = ['POINTS', 'fs']
threshold1.ThresholdRange = [0.0, 0.95]

# show data in view
threshold1Display = Show(threshold1, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
threshold1Display.Representation = 'Surface'
threshold1Display.ColorArrayName = [None, '']
threshold1Display.Opacity = 0.85
threshold1Display.OSPRayScaleArray = 'f'
threshold1Display.OSPRayScaleFunction = 'PiecewiseFunction'
threshold1Display.SelectOrientationVectors = 'None'
threshold1Display.ScaleFactor = 3.0
threshold1Display.SelectScaleArray = 'None'
threshold1Display.GlyphType = 'Arrow'
threshold1Display.GlyphTableIndexArray = 'None'
threshold1Display.GaussianRadius = 0.15
threshold1Display.SetScaleArray = ['POINTS', 'f']
threshold1Display.ScaleTransferFunction = 'PiecewiseFunction'
threshold1Display.OpacityArray = ['POINTS', 'f']
threshold1Display.OpacityTransferFunction = 'PiecewiseFunction'
threshold1Display.DataAxesGrid = 'GridAxesRepresentation'
threshold1Display.PolarAxes = 'PolarAxesRepresentation'
threshold1Display.ScalarOpacityUnitDistance = 0.30986217978916686

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
threshold1Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# hide data in view
Hide(cellDatatoPointData1, renderView1)

# create a new 'Clip'
clip1 = Clip(Input=threshold1)
clip1.ClipType = 'Plane'
clip1.HyperTreeGridClipper = 'Plane'
clip1.Scalars = ['POINTS', 'f']
clip1.Value = 0.5

# init the 'Plane' selected for 'ClipType'
clip1.ClipType.Origin = [15.0, 0.0, 0.0]

# init the 'Plane' selected for 'HyperTreeGridClipper'
clip1.HyperTreeGridClipper.Origin = [15.0, 0.0, 0.0]

# show data in view
clip1Display = Show(clip1, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
clip1Display.Representation = 'Surface'
clip1Display.ColorArrayName = [None, '']
clip1Display.Opacity = 0.85
clip1Display.OSPRayScaleArray = 'f'
clip1Display.OSPRayScaleFunction = 'PiecewiseFunction'
clip1Display.SelectOrientationVectors = 'None'
clip1Display.ScaleFactor = 1.5
clip1Display.SelectScaleArray = 'None'
clip1Display.GlyphType = 'Arrow'
clip1Display.GlyphTableIndexArray = 'None'
clip1Display.GaussianRadius = 0.075
clip1Display.SetScaleArray = ['POINTS', 'f']
clip1Display.ScaleTransferFunction = 'PiecewiseFunction'
clip1Display.OpacityArray = ['POINTS', 'f']
clip1Display.OpacityTransferFunction = 'PiecewiseFunction'
clip1Display.DataAxesGrid = 'GridAxesRepresentation'
clip1Display.PolarAxes = 'PolarAxesRepresentation'
clip1Display.ScalarOpacityUnitDistance = 0.19649440830229037

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
clip1Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# hide data in view
Hide(threshold1, renderView1)

# Properties modified on clip1.ClipType
clip1.ClipType.Normal = [-1.0, 0.0, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Clip'
clip2 = Clip(Input=clip1)
clip2.ClipType = 'Plane'
clip2.HyperTreeGridClipper = 'Plane'
clip2.Scalars = ['POINTS', 'f']
clip2.Value = 0.5

# init the 'Plane' selected for 'ClipType'
clip2.ClipType.Origin = [22.5, 0.0, 0.0]

# init the 'Plane' selected for 'HyperTreeGridClipper'
clip2.HyperTreeGridClipper.Origin = [22.5, 0.0, 0.0]

# Properties modified on clip2.ClipType
clip2.ClipType.Origin = [17.0, 0.0, 0.0]
clip2.ClipType.Normal = [-1.0, 0.0, 0.0]

# show data in view
clip2Display = Show(clip2, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
clip2Display.Representation = 'Surface'
clip2Display.ColorArrayName = [None, '']
clip2Display.Opacity = 0.85
clip2Display.OSPRayScaleArray = 'f'
clip2Display.OSPRayScaleFunction = 'PiecewiseFunction'
clip2Display.SelectOrientationVectors = 'None'
clip2Display.ScaleFactor = 1.3
clip2Display.SelectScaleArray = 'None'
clip2Display.GlyphType = 'Arrow'
clip2Display.GlyphTableIndexArray = 'None'
clip2Display.GaussianRadius = 0.065
clip2Display.SetScaleArray = ['POINTS', 'f']
clip2Display.ScaleTransferFunction = 'PiecewiseFunction'
clip2Display.OpacityArray = ['POINTS', 'f']
clip2Display.OpacityTransferFunction = 'PiecewiseFunction'
clip2Display.DataAxesGrid = 'GridAxesRepresentation'
clip2Display.PolarAxes = 'PolarAxesRepresentation'
clip2Display.ScalarOpacityUnitDistance = 0.17615700667391

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
clip2Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# hide data in view
Hide(clip1, renderView1)

# set active source
SetActiveSource(clip1)

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=clip2.ClipType)

# toggle 3D widget visibility (only when running from the GUI)
Show3DWidgets(proxy=clip1.ClipType)

# set active source
SetActiveSource(clip2)

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=clip1.ClipType)

# toggle 3D widget visibility (only when running from the GUI)
Show3DWidgets(proxy=clip2.ClipType)

# Properties modified on clip2.ClipType
clip2.ClipType.Normal = [1.0, 0.0, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Delaunay 3D'
delaunay3D1 = Delaunay3D(Input=clip2)

# show data in view
delaunay3D1Display = Show(delaunay3D1, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
delaunay3D1Display.Representation = 'Surface'
delaunay3D1Display.ColorArrayName = [None, '']
delaunay3D1Display.Opacity = 0.85
delaunay3D1Display.OSPRayScaleArray = 'f'
delaunay3D1Display.OSPRayScaleFunction = 'PiecewiseFunction'
delaunay3D1Display.SelectOrientationVectors = 'None'
delaunay3D1Display.ScaleFactor = 0.2
delaunay3D1Display.SelectScaleArray = 'None'
delaunay3D1Display.GlyphType = 'Arrow'
delaunay3D1Display.GlyphTableIndexArray = 'None'
delaunay3D1Display.GaussianRadius = 0.01
delaunay3D1Display.SetScaleArray = ['POINTS', 'f']
delaunay3D1Display.ScaleTransferFunction = 'PiecewiseFunction'
delaunay3D1Display.OpacityArray = ['POINTS', 'f']
delaunay3D1Display.OpacityTransferFunction = 'PiecewiseFunction'
delaunay3D1Display.DataAxesGrid = 'GridAxesRepresentation'
delaunay3D1Display.PolarAxes = 'PolarAxesRepresentation'
delaunay3D1Display.ScalarOpacityUnitDistance = 0.0351274886866612

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
delaunay3D1Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# hide data in view
Hide(clip2, renderView1)

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Contour'
contour1 = Contour(Input=delaunay3D1)
contour1.ContourBy = ['POINTS', 'f']
contour1.Isosurfaces = [0.5]
contour1.PointMergeMethod = 'Uniform Binning'

# show data in view
contour1Display = Show(contour1, renderView1, 'GeometryRepresentation')

# get color transfer function/color map for 'f'
fLUT = GetColorTransferFunction('f')
fLUT.RGBPoints = [0.5, 0.23137254902, 0.298039215686, 0.752941176471, 0.50006103515625, 0.865, 0.865, 0.865, 0.5001220703125, 0.705882352941, 0.0156862745098, 0.149019607843]
fLUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
contour1Display.Representation = 'Surface'
contour1Display.ColorArrayName = ['POINTS', 'f']
contour1Display.LookupTable = fLUT
contour1Display.OSPRayScaleArray = 'f'
contour1Display.OSPRayScaleFunction = 'PiecewiseFunction'
contour1Display.SelectOrientationVectors = 'None'
contour1Display.ScaleFactor = 0.2
contour1Display.SelectScaleArray = 'f'
contour1Display.GlyphType = 'Arrow'
contour1Display.GlyphTableIndexArray = 'f'
contour1Display.GaussianRadius = 0.01
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
fPWF.Points = [0.5, 0.0, 0.5, 0.0, 0.5001220703125, 1.0, 0.5, 0.0]
fPWF.ScalarRangeInitialized = 1

# hide data in view
Hide(delaunay3D1, renderView1)

# reset view to fit data bounds
renderView1.ResetCamera(15.0, 17.0, -0.427789067491, 0.425649342926, -0.427947942457, 0.427292957237)

# current camera placement for renderView1
renderView1.CameraPosition = [16.0, -0.0010698622826302018, 4.51367616745501]
renderView1.CameraFocalPoint = [16.0, -0.0010698622826302018, -0.000327492610103286]
renderView1.CameraParallelScale = 1.1683101168873362

# save screenshot
SaveScreenshot(path+'mm1.png', renderView1, ImageResolution=[2120, 1134],
    TransparentBackground=0)

# Hide orientation axes
renderView1.OrientationAxesVisibility = 0

# set scalar coloring
ColorBy(contour1Display, ('POINTS', 'u.x', 'Magnitude'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(fLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
contour1Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
contour1Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'ux'
uxLUT = GetColorTransferFunction('ux')
uxLUT.RGBPoints = [2.2554501288418073e-07, 0.23137254902, 0.298039215686, 0.752941176471, 1.0513483563140162, 0.865, 0.865, 0.865, 2.1026964870830196, 0.705882352941, 0.0156862745098, 0.149019607843]
uxLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'ux'
uxPWF = GetOpacityTransferFunction('ux')
uxPWF.Points = [2.2554501288418073e-07, 0.0, 0.5, 0.0, 2.1026964870830196, 1.0, 0.5, 0.0]
uxPWF.ScalarRangeInitialized = 1

# hide color bar/color legend
contour1Display.SetScalarBarVisibility(renderView1, False)

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
uxLUT.ApplyPreset('jet', True)

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
uxLUT.ApplyPreset('jet', True)

# reset view to fit data bounds
renderView1.ResetCamera(15.0, 17.0, -0.427789067491, 0.425649342926, -0.427947942457, 0.427292957237)

# reset view to fit data
renderView1.ResetCamera()

# current camera placement for renderView1
renderView1.CameraPosition = [16.0, -0.0010698622826302018, 2.4798274899700514]
renderView1.CameraFocalPoint = [16.0, -0.0010698622826302018, -2.0341761700950607]
renderView1.CameraParallelScale = 1.1683101168873362

# save screenshot
SaveScreenshot(path+'mm2.png', renderView1, ImageResolution=[2120, 1134],
    TransparentBackground=0)

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.CameraPosition = [16.0, -0.0010698622826302018, 2.4798274899700514]
renderView1.CameraFocalPoint = [16.0, -0.0010698622826302018, -2.0341761700950607]
renderView1.CameraParallelScale = 1.1683101168873362

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).