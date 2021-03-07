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
renderView1.ViewSize = [980, 1132]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.CenterOfRotation = [15.5, 0.0, 0.0]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [12.371051723102436, 0.040124162100615104, 0.0]
renderView1.CameraFocalPoint = [15.5, 0.0, 0.0]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 0.8167804860236447
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.OSPRayMaterialLibrary = materialLibrary1

# Create a new 'Render View'
renderView2 = CreateView('RenderView')
renderView2.ViewSize = [980, 1132]
renderView2.AxesGrid = 'GridAxes3DActor'
renderView2.OrientationAxesVisibility = 0
renderView2.CenterOfRotation = [16.13227939605713, 0.002989262342453003, 0.0002599954605102539]
renderView2.StereoType = 'Crystal Eyes'
renderView2.CameraPosition = [16.13227939605713, 0.002989262342453003, 7.966163198512591]
renderView2.CameraFocalPoint = [16.13227939605713, 0.002989262342453003, 0.0002599954605102539]
renderView2.CameraFocalDisk = 1.0
renderView2.CameraParallelScale = 2.0792479820033907
renderView2.BackEnd = 'OSPRay raycaster'
renderView2.OSPRayMaterialLibrary = materialLibrary1

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.SplitHorizontal(0, 0.500000)
layout1.AssignView(1, renderView1)
layout1.AssignView(2, renderView2)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'Text'
text1 = Text()
text1.Text = 'Text Hahah'

# create a new 'PVD Reader'
dump2pvd_compressedpvd = PVDReader(FileName='/Users/weugene/basilisk/work/tube/dump2pvd_compressed.pvd')
dump2pvd_compressedpvd.CellArrays = ['fs', 'f', 'l', 'l2', 'omega', 'u.x']

# create a new 'Clip'
clip1 = Clip(Input=dump2pvd_compressedpvd)
clip1.ClipType = 'Box'
clip1.HyperTreeGridClipper = 'Plane'
clip1.Scalars = ['POINTS', 'f']
clip1.Value = 0.5

# init the 'Box' selected for 'ClipType'
clip1.ClipType.Position = [12.0, -0.6, -0.6]
clip1.ClipType.Length = [7.0, 1.2, 1.2]

# init the 'Plane' selected for 'HyperTreeGridClipper'
clip1.HyperTreeGridClipper.Origin = [16.132652442407533, 0.0012683107921633574, 0.0013393794089324729]

# create a new 'Cell Data to Point Data'
cellDatatoPointData1 = CellDatatoPointData(Input=clip1)
cellDatatoPointData1.CellDataArraytoprocess = ['f', 'f', 'fs', 'l', 'l2', 'omega', 'u.x']

# create a new 'Resample To Image'
resampleToImage1 = ResampleToImage(Input=cellDatatoPointData1)
resampleToImage1.UseInputBounds = 0
resampleToImage1.SamplingBounds = [13.0, 18.0, -0.5, 0.5, -0.5, 0.5]

# create a new 'Slice'
slice1 = Slice(Input=resampleToImage1)
slice1.SliceType = 'Plane'
slice1.HyperTreeGridSlicer = 'Plane'
slice1.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Origin = [15.5, 0.0, 0.0]

# init the 'Plane' selected for 'HyperTreeGridSlicer'
slice1.HyperTreeGridSlicer.Origin = [15.5, 0.0, 0.0]

# create a new 'Contour'
contour1 = Contour(Input=dump2pvd_compressedpvd)
contour1.ContourBy = ['POINTS', 'f']
contour1.Isosurfaces = [0.5]
contour1.PointMergeMethod = 'Uniform Binning'

# create a new 'Contour'
contour2 = Contour(Input=resampleToImage1)
contour2.ContourBy = ['POINTS', 'f']
contour2.Isosurfaces = [0.5]
contour2.PointMergeMethod = 'Uniform Binning'

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from text1
text1Display = Show(text1, renderView1, 'TextSourceRepresentation')

# trace defaults for the display properties.
# text1Display.WindowLocation = 'AnyLocation'
# text1Display.Position = [0.26860647289218714, 0.8672083886925798]
text1Display.WindowLocation = 'UpperCenter'

# show data from slice1
slice1Display = Show(slice1, renderView1, 'GeometryRepresentation')

# get color transfer function/color map for 'omega'
omegaLUT = GetColorTransferFunction('omega')
omegaLUT.RGBPoints = [-27.390586249938153, 0.23137254902, 0.298039215686, 0.752941176471, -1.1160434694887797, 0.865, 0.865, 0.865, 25.158499310960593, 0.705882352941, 0.0156862745098, 0.149019607843]
omegaLUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
slice1Display.Representation = 'Surface'
slice1Display.ColorArrayName = ['POINTS', 'omega']
slice1Display.LookupTable = omegaLUT
slice1Display.OSPRayScaleArray = 'f'
slice1Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice1Display.SelectOrientationVectors = 'None'
slice1Display.ScaleFactor = 0.1
slice1Display.SelectScaleArray = 'None'
slice1Display.GlyphType = 'Arrow'
slice1Display.GlyphTableIndexArray = 'None'
slice1Display.GaussianRadius = 0.005
slice1Display.SetScaleArray = ['POINTS', 'f']
slice1Display.ScaleTransferFunction = 'PiecewiseFunction'
slice1Display.OpacityArray = ['POINTS', 'f']
slice1Display.OpacityTransferFunction = 'PiecewiseFunction'
slice1Display.DataAxesGrid = 'GridAxesRepresentation'
slice1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
slice1Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
slice1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0000000000000002, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
slice1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0000000000000002, 1.0, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for omegaLUT in view renderView1
omegaLUTColorBar = GetScalarBar(omegaLUT, renderView1)
omegaLUTColorBar.Orientation = 'Horizontal'
omegaLUTColorBar.WindowLocation = 'AnyLocation'
omegaLUTColorBar.Position = [0.29346938775510184, 0.02120141342756185]
omegaLUTColorBar.Title = 'omega'
omegaLUTColorBar.ComponentTitle = ''
omegaLUTColorBar.LabelFormat = '%-#6.2g'
omegaLUTColorBar.RangeLabelFormat = '%6.2g'
omegaLUTColorBar.ScalarBarLength = 0.33000000000000035

# set color bar visibility
omegaLUTColorBar.Visibility = 1

# show color legend
slice1Display.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup the visualization in view 'renderView2'
# ----------------------------------------------------------------

# show data from contour2
contour2Display = Show(contour2, renderView2, 'GeometryRepresentation')

# get color transfer function/color map for 'f'
fLUT = GetColorTransferFunction('f')
fLUT.RGBPoints = [0.0, 0.23137254902, 0.298039215686, 0.752941176471, 0.5000000000000002, 0.865, 0.865, 0.865, 1.0000000000000004, 0.705882352941, 0.0156862745098, 0.149019607843]
fLUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
contour2Display.Representation = 'Surface'
contour2Display.ColorArrayName = ['POINTS', 'f']
contour2Display.LookupTable = fLUT
contour2Display.OSPRayScaleArray = 'f'
contour2Display.OSPRayScaleFunction = 'PiecewiseFunction'
contour2Display.SelectOrientationVectors = 'None'
contour2Display.ScaleFactor = 0.3378744125366211
contour2Display.SelectScaleArray = 'f'
contour2Display.GlyphType = 'Arrow'
contour2Display.GlyphTableIndexArray = 'f'
contour2Display.GaussianRadius = 0.016893720626831057
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

# show data from slice1
slice1Display_1 = Show(slice1, renderView2, 'GeometryRepresentation')

# trace defaults for the display properties.
slice1Display_1.Representation = 'Surface'
slice1Display_1.ColorArrayName = [None, '']
slice1Display_1.OSPRayScaleArray = 'f'
slice1Display_1.OSPRayScaleFunction = 'PiecewiseFunction'
slice1Display_1.SelectOrientationVectors = 'None'
slice1Display_1.ScaleFactor = 0.1
slice1Display_1.SelectScaleArray = 'None'
slice1Display_1.GlyphType = 'Arrow'
slice1Display_1.GlyphTableIndexArray = 'None'
slice1Display_1.GaussianRadius = 0.005
slice1Display_1.SetScaleArray = ['POINTS', 'f']
slice1Display_1.ScaleTransferFunction = 'PiecewiseFunction'
slice1Display_1.OpacityArray = ['POINTS', 'f']
slice1Display_1.OpacityTransferFunction = 'PiecewiseFunction'
slice1Display_1.DataAxesGrid = 'GridAxesRepresentation'
slice1Display_1.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
slice1Display_1.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
slice1Display_1.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0000000000000002, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
slice1Display_1.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0000000000000002, 1.0, 0.5, 0.0]

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get opacity transfer function/opacity map for 'omega'
omegaPWF = GetOpacityTransferFunction('omega')
omegaPWF.Points = [-27.390586249938153, 0.0, 0.5, 0.0, 25.158499310960593, 1.0, 0.5, 0.0]
omegaPWF.ScalarRangeInitialized = 1

# get opacity transfer function/opacity map for 'f'
fPWF = GetOpacityTransferFunction('f')
fPWF.Points = [0.0, 0.0, 0.5, 0.0, 1.0000000000000004, 1.0, 0.5, 0.0]
fPWF.ScalarRangeInitialized = 1

# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(slice1)
# ----------------------------------------------------------------