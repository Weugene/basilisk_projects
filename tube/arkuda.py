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
renderView1.ViewSize = [1980, 1132]
renderView1.InteractionMode = '2D'
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.CenterOfRotation = [10.999999761581421, 0.0, 1.3010426069826053e-17]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [10.999999761581421, 0.0, 15.191774263830615]
renderView1.CameraFocalPoint = [10.999999761581421, 0.0, 1.3010426069826053e-17]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 2.6855546126478305
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.OSPRayMaterialLibrary = materialLibrary1

# Create a new 'Render View'
renderView2 = CreateView('RenderView')
renderView2.ViewSize = [1980, 1132]
renderView2.AxesGrid = 'GridAxes3DActor'
renderView2.CenterOfRotation = [12.0, 0.0, 0.0]
renderView2.StereoType = 'Crystal Eyes'
renderView2.CameraPosition = [13.67836732014027, 2.74373007671147, 0.7636881428680581]
renderView2.CameraFocalPoint = [12.0, 0.0, 0.0]
renderView2.CameraViewUp = [-0.12054918240887771, -0.19706914781358603, 0.9729499707593285]
renderView2.CameraFocalDisk = 1.0
renderView2.CameraParallelScale = 0.7071067811865476
renderView2.BackEnd = 'OSPRay raycaster'
renderView2.OSPRayMaterialLibrary = materialLibrary1

# Create a new 'SpreadSheet View'
spreadSheetView1 = CreateView('SpreadSheetView')
spreadSheetView1.ColumnToSort = ''
spreadSheetView1.BlockSize = 1024L
# uncomment following to set a specific view size
# spreadSheetView1.ViewSize = [400, 400]

# Create a new 'SpreadSheet View'
spreadSheetView2 = CreateView('SpreadSheetView')
spreadSheetView2.ColumnToSort = ''
spreadSheetView2.BlockSize = 1024L
# uncomment following to set a specific view size
# spreadSheetView2.ViewSize = [400, 400]

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.AssignView(0, renderView1)

# create new layout object 'Layout #2'
layout2 = CreateLayout(name='Layout #2')
layout2.AssignView(0, spreadSheetView1)

# create new layout object 'Layout #3'
layout3 = CreateLayout(name='Layout #3')
layout3.AssignView(0, spreadSheetView2)

# create new layout object 'Layout #4'
layout4 = CreateLayout(name='Layout #4')
layout4.AssignView(0, renderView2)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView2)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'Cylinder'
cylinder1 = Cylinder()
cylinder1.Height = 30.0
cylinder1.Capping = 0

# create a new 'PVD Reader'
dump2pvd_compressedpvd = PVDReader(FileName='/home-local/e.sharaborin/res_wbasilisk/tube/res22_adaptffsu/dump2pvd_compressed.pvd')
dump2pvd_compressedpvd.CellArrays = ['p', 'fs', 'f', 'l', 'l2', 'omega', 'u.x']

# create a new 'Contour'
contour1 = Contour(Input=dump2pvd_compressedpvd)
contour1.ContourBy = ['POINTS', 'f']
contour1.Isosurfaces = [0.5]
contour1.PointMergeMethod = 'Uniform Binning'

# create a new 'Clip'
clip4 = Clip(Input=dump2pvd_compressedpvd)
clip4.ClipType = 'Plane'
clip4.HyperTreeGridClipper = 'Plane'
clip4.Scalars = ['POINTS', 'f']
clip4.Value = 0.5

# init the 'Plane' selected for 'ClipType'
clip4.ClipType.Origin = [15.0, 0.0, 0.0]
clip4.ClipType.Normal = [0.0, 0.0, 1.0]

# init the 'Plane' selected for 'HyperTreeGridClipper'
clip4.HyperTreeGridClipper.Origin = [15.0, 0.0, 0.0]

# create a new 'Clip'
clip1 = Clip(Input=dump2pvd_compressedpvd)
clip1.ClipType = 'Box'
clip1.HyperTreeGridClipper = 'Plane'
clip1.Scalars = ['POINTS', 'f']
clip1.Value = 0.5

# init the 'Box' selected for 'ClipType'
clip1.ClipType.Position = [7.0, -0.52, -0.52]
clip1.ClipType.Length = [8.0, 1.04, 1.04]

# init the 'Plane' selected for 'HyperTreeGridClipper'
clip1.HyperTreeGridClipper.Origin = [12.383253908307324, 0.0003794258290628072, 6.398686387718011e-05]

# create a new 'Cell Data to Point Data'
cellDatatoPointData1 = CellDatatoPointData(Input=clip1)
cellDatatoPointData1.CellDataArraytoprocess = ['f', 'f', 'fs', 'l', 'l2', 'omega', 'p', 'u.x']

# create a new 'Resample To Image'
resampleToImage1 = ResampleToImage(Input=cellDatatoPointData1)
resampleToImage1.UseInputBounds = 0
resampleToImage1.SamplingDimensions = [1000, 273, 273]
resampleToImage1.SamplingBounds = [7.1, 14.9, -0.5, 0.5, -0.5, 0.5]

# create a new 'Contour'
contour3 = Contour(Input=resampleToImage1)
contour3.ContourBy = ['POINTS', 'l2']
contour3.Isosurfaces = [-1.0, -2.0]
contour3.PointMergeMethod = 'Uniform Binning'

# create a new 'Stream Tracer'
streamTracer1 = StreamTracer(Input=resampleToImage1,
    SeedType='High Resolution Line Source')
streamTracer1.Vectors = ['POINTS', 'u.x']
streamTracer1.SurfaceStreamlines = 1
streamTracer1.InitialStepLength = 0.05
streamTracer1.MaximumStreamlineLength = 7.800000000000001
streamTracer1.TerminalSpeed = 0.0001

# init the 'High Resolution Line Source' selected for 'SeedType'
streamTracer1.SeedType.Point1 = [14.5, -0.49, 0.0]
streamTracer1.SeedType.Point2 = [14.5, 0.49, 0.0]
streamTracer1.SeedType.Resolution = 50

# create a new 'Slice'
slice4 = Slice(Input=resampleToImage1)
slice4.SliceType = 'Plane'
slice4.HyperTreeGridSlicer = 'Plane'
slice4.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice4.SliceType.Origin = [11.0, 0.0, 0.0]
slice4.SliceType.Normal = [0.0, 0.0, 1.0]

# init the 'Plane' selected for 'HyperTreeGridSlicer'
slice4.HyperTreeGridSlicer.Origin = [11.0, 0.0, 0.0]

# create a new 'Clip'
clip2 = Clip(Input=cylinder1)
clip2.ClipType = 'Plane'
clip2.HyperTreeGridClipper = 'Plane'
clip2.Scalars = ['POINTS', 'Normals_Magnitude']
clip2.Value = 0.9999999932694132

# init the 'Plane' selected for 'ClipType'
clip2.ClipType.Normal = [0.0, 0.0, 1.0]

# create a new 'Transform'
transform1 = Transform(Input=clip2)
transform1.Transform = 'Transform'

# init the 'Transform' selected for 'Transform'
transform1.Transform.Translate = [15.0, 0.0, 0.0]
transform1.Transform.Rotate = [0.0, 0.0, 90.0]

# create a new 'Slice'
slice2 = Slice(Input=transform1)
slice2.SliceType = 'Plane'
slice2.HyperTreeGridSlicer = 'Plane'
slice2.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice2.SliceType.Origin = [15.0, 0.0, -0.0001]
slice2.SliceType.Normal = [0.0, 0.0, 1.0]

# init the 'Plane' selected for 'HyperTreeGridSlicer'
slice2.HyperTreeGridSlicer.Origin = [15.0, 0.0, -0.21650634706020355]

# create a new 'Slice'
slice3 = Slice(Input=resampleToImage1)
slice3.SliceType = 'Plane'
slice3.HyperTreeGridSlicer = 'Plane'
slice3.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice3.SliceType.Origin = [13.0, 0.0, 0.0]

# init the 'Plane' selected for 'HyperTreeGridSlicer'
slice3.HyperTreeGridSlicer.Origin = [11.0, 0.0, 0.0]

# create a new 'Contour'
contour2 = Contour(Input=resampleToImage1)
contour2.ContourBy = ['POINTS', 'f']
contour2.Isosurfaces = [0.5]
contour2.PointMergeMethod = 'Uniform Binning'

# create a new 'Clip'
clip3 = Clip(Input=contour2)
clip3.ClipType = 'Box'
clip3.HyperTreeGridClipper = 'Plane'
clip3.Scalars = ['POINTS', 'f']
clip3.Value = 0.5

# init the 'Box' selected for 'ClipType'
clip3.ClipType.Position = [10.6, -0.5, -0.5]
clip3.ClipType.Length = [3.5, 1.0, 1.0]

# init the 'Plane' selected for 'HyperTreeGridClipper'
clip3.HyperTreeGridClipper.Origin = [10.999999761581421, 0.0, 0.0]

# create a new 'Slice'
slice1 = Slice(Input=contour2)
slice1.SliceType = 'Plane'
slice1.HyperTreeGridSlicer = 'Plane'
slice1.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Origin = [10.999999761581421, 0.0, 0.0]
slice1.SliceType.Normal = [0.0, 0.0, 1.0]

# init the 'Plane' selected for 'HyperTreeGridSlicer'
slice1.HyperTreeGridSlicer.Origin = [10.999999761581421, 0.0, 0.0]

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from contour2
contour2Display = Show(contour2, renderView1, 'GeometryRepresentation')

# get color transfer function/color map for 'omega'
omegaLUT = GetColorTransferFunction('omega')
omegaLUT.RGBPoints = [-265.5696132211881, 0.23137254902, 0.298039215686, 0.752941176471, -2.4825307889177566, 0.865, 0.865, 0.865, 260.6045516433526, 0.705882352941, 0.0156862745098, 0.149019607843]
omegaLUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
contour2Display.Representation = 'Surface'
contour2Display.AmbientColor = [0.3333333333333333, 0.3333333333333333, 1.0]
contour2Display.ColorArrayName = ['POINTS', 'omega']
contour2Display.DiffuseColor = [0.3333333333333333, 0.3333333333333333, 1.0]
contour2Display.LookupTable = omegaLUT
contour2Display.Specular = 1.0
contour2Display.OSPRayScaleArray = 'f'
contour2Display.OSPRayScaleFunction = 'PiecewiseFunction'
contour2Display.SelectOrientationVectors = 'None'
contour2Display.ScaleFactor = 0.8
contour2Display.SelectScaleArray = 'f'
contour2Display.GlyphType = 'Arrow'
contour2Display.GlyphTableIndexArray = 'f'
contour2Display.GaussianRadius = 0.04
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

# show data from transform1
transform1Display = Show(transform1, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
transform1Display.Representation = 'Surface'
transform1Display.ColorArrayName = [None, '']
transform1Display.Opacity = 0.85
transform1Display.OSPRayScaleArray = 'Normals'
transform1Display.OSPRayScaleFunction = 'PiecewiseFunction'
transform1Display.SelectOrientationVectors = 'None'
transform1Display.ScaleFactor = 3.0
transform1Display.SelectScaleArray = 'None'
transform1Display.GlyphType = 'Arrow'
transform1Display.GlyphTableIndexArray = 'None'
transform1Display.GaussianRadius = 0.15
transform1Display.SetScaleArray = ['POINTS', 'Normals']
transform1Display.ScaleTransferFunction = 'PiecewiseFunction'
transform1Display.OpacityArray = ['POINTS', 'Normals']
transform1Display.OpacityTransferFunction = 'PiecewiseFunction'
transform1Display.DataAxesGrid = 'GridAxesRepresentation'
transform1Display.PolarAxes = 'PolarAxesRepresentation'
transform1Display.ScalarOpacityUnitDistance = 18.91127960708838

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
transform1Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
transform1Display.ScaleTransferFunction.Points = [-2.220446049250313e-16, 0.0, 0.5, 0.0, 2.220446049250313e-16, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
transform1Display.OpacityTransferFunction.Points = [-2.220446049250313e-16, 0.0, 0.5, 0.0, 2.220446049250313e-16, 1.0, 0.5, 0.0]

# show data from slice3
slice3Display = Show(slice3, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
slice3Display.Representation = 'Surface'
slice3Display.ColorArrayName = ['POINTS', 'omega']
slice3Display.LookupTable = omegaLUT
slice3Display.Specular = 1.0
slice3Display.Luminosity = 35.0
slice3Display.OSPRayScaleArray = 'f'
slice3Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice3Display.SelectOrientationVectors = 'None'
slice3Display.ScaleFactor = 0.1
slice3Display.SelectScaleArray = 'None'
slice3Display.GlyphType = 'Arrow'
slice3Display.GlyphTableIndexArray = 'None'
slice3Display.GaussianRadius = 0.005
slice3Display.SetScaleArray = ['POINTS', 'f']
slice3Display.ScaleTransferFunction = 'PiecewiseFunction'
slice3Display.OpacityArray = ['POINTS', 'f']
slice3Display.OpacityTransferFunction = 'PiecewiseFunction'
slice3Display.DataAxesGrid = 'GridAxesRepresentation'
slice3Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
slice3Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
slice3Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0000000000000004, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
slice3Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0000000000000004, 1.0, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for omegaLUT in view renderView1
omegaLUTColorBar = GetScalarBar(omegaLUT, renderView1)
omegaLUTColorBar.WindowLocation = 'AnyLocation'
omegaLUTColorBar.Position = [0.7700336700336701, 0.6145449590256371]
omegaLUTColorBar.Title = 'omega'
omegaLUTColorBar.ComponentTitle = ''
omegaLUTColorBar.LabelFormat = '%-#6.2g'
omegaLUTColorBar.RangeLabelFormat = '%6.2g'
omegaLUTColorBar.ScalarBarLength = 0.3300000000000002

# set color bar visibility
omegaLUTColorBar.Visibility = 1

# show color legend
contour2Display.SetScalarBarVisibility(renderView1, True)

# show color legend
slice3Display.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup the visualization in view 'renderView2'
# ----------------------------------------------------------------

# show data from slice3
slice3Display_1 = Show(slice3, renderView2, 'GeometryRepresentation')

# trace defaults for the display properties.
slice3Display_1.Representation = 'Surface'
slice3Display_1.ColorArrayName = ['POINTS', 'omega']
slice3Display_1.LookupTable = omegaLUT
slice3Display_1.OSPRayScaleArray = 'f'
slice3Display_1.OSPRayScaleFunction = 'PiecewiseFunction'
slice3Display_1.SelectOrientationVectors = 'None'
slice3Display_1.ScaleFactor = 0.1
slice3Display_1.SelectScaleArray = 'None'
slice3Display_1.GlyphType = 'Arrow'
slice3Display_1.GlyphTableIndexArray = 'None'
slice3Display_1.GaussianRadius = 0.005
slice3Display_1.SetScaleArray = ['POINTS', 'f']
slice3Display_1.ScaleTransferFunction = 'PiecewiseFunction'
slice3Display_1.OpacityArray = ['POINTS', 'f']
slice3Display_1.OpacityTransferFunction = 'PiecewiseFunction'
slice3Display_1.DataAxesGrid = 'GridAxesRepresentation'
slice3Display_1.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
slice3Display_1.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
slice3Display_1.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0000000000000004, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
slice3Display_1.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0000000000000004, 1.0, 0.5, 0.0]

# show data from contour2
contour2Display_1 = Show(contour2, renderView2, 'GeometryRepresentation')

# trace defaults for the display properties.
contour2Display_1.Representation = 'Surface'
contour2Display_1.ColorArrayName = ['POINTS', 'omega']
contour2Display_1.LookupTable = omegaLUT
contour2Display_1.OSPRayScaleArray = 'f'
contour2Display_1.OSPRayScaleFunction = 'PiecewiseFunction'
contour2Display_1.SelectOrientationVectors = 'None'
contour2Display_1.ScaleFactor = 0.7799999713897705
contour2Display_1.SelectScaleArray = 'f'
contour2Display_1.GlyphType = 'Arrow'
contour2Display_1.GlyphTableIndexArray = 'f'
contour2Display_1.GaussianRadius = 0.03899999856948853
contour2Display_1.SetScaleArray = ['POINTS', 'f']
contour2Display_1.ScaleTransferFunction = 'PiecewiseFunction'
contour2Display_1.OpacityArray = ['POINTS', 'f']
contour2Display_1.OpacityTransferFunction = 'PiecewiseFunction'
contour2Display_1.DataAxesGrid = 'GridAxesRepresentation'
contour2Display_1.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
contour2Display_1.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
contour2Display_1.ScaleTransferFunction.Points = [0.5, 0.0, 0.5, 0.0, 0.5001220703125, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
contour2Display_1.OpacityTransferFunction.Points = [0.5, 0.0, 0.5, 0.0, 0.5001220703125, 1.0, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for omegaLUT in view renderView2
omegaLUTColorBar_1 = GetScalarBar(omegaLUT, renderView2)
omegaLUTColorBar_1.WindowLocation = 'UpperLeftCorner'
omegaLUTColorBar_1.Position = [0.00202020202020202, 0.6519434628975265]
omegaLUTColorBar_1.Title = 'omega'
omegaLUTColorBar_1.ComponentTitle = ''
omegaLUTColorBar_1.LabelFormat = '%-#6.2g'
omegaLUTColorBar_1.RangeLabelFormat = '%6.2g'
omegaLUTColorBar_1.ScalarBarLength = 0.3300000000000002

# set color bar visibility
omegaLUTColorBar_1.Visibility = 1

# show color legend
slice3Display_1.SetScalarBarVisibility(renderView2, True)

# show color legend
contour2Display_1.SetScalarBarVisibility(renderView2, True)

# ----------------------------------------------------------------
# setup the visualization in view 'spreadSheetView1'
# ----------------------------------------------------------------

# show data from slice1
slice1Display = Show(slice1, spreadSheetView1, 'SpreadSheetRepresentation')

# ----------------------------------------------------------------
# setup the visualization in view 'spreadSheetView2'
# ----------------------------------------------------------------

# show data from streamTracer1
streamTracer1Display = Show(streamTracer1, spreadSheetView2, 'SpreadSheetRepresentation')

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get opacity transfer function/opacity map for 'omega'
omegaPWF = GetOpacityTransferFunction('omega')
omegaPWF.Points = [-265.5696132211881, 0.0, 0.5, 0.0, 260.6045516433526, 1.0, 0.5, 0.0]
omegaPWF.ScalarRangeInitialized = 1

# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(contour2)
# ----------------------------------------------------------------