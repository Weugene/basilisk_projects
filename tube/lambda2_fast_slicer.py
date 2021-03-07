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
renderView1.ViewSize = [1338, 1134]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.CenterOfRotation = [15.0, 0.0, -0.25]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [15.0, 0.0, 47.68039877838649]
renderView1.CameraFocalPoint = [15.0, 0.0, -0.25]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 15.010413052278075
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.OSPRayMaterialLibrary = materialLibrary1

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.AssignView(0, renderView1)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'PVD Reader'
tube_bp_from_dumppvd = PVDReader(FileName='/Users/weugene/basilisk/work/tube/tube_bp_from_dump.pvd')
tube_bp_from_dumppvd.CellArrays = ['p', 'fs', 'f', 'l', 'residual_of_p', 'u.x']

# create a new 'Contour'
contour1 = Contour(Input=tube_bp_from_dumppvd)
contour1.ContourBy = ['POINTS', 'f']
contour1.Isosurfaces = [0.5]
contour1.PointMergeMethod = 'Uniform Binning'

# create a new 'Clip'
clip2 = Clip(Input=tube_bp_from_dumppvd)
clip2.ClipType = 'Box'
clip2.HyperTreeGridClipper = 'Plane'
clip2.Scalars = ['POINTS', 'f']
clip2.Value = 0.5

# init the 'Box' selected for 'ClipType'
clip2.ClipType.Position = [0.3, -0.55, -0.55]
clip2.ClipType.Length = [5.7784823612399006, 1.1, 1.1]

# init the 'Plane' selected for 'HyperTreeGridClipper'
clip2.HyperTreeGridClipper.Origin = [15.0, 0.0, 0.0]

# create a new 'Cell Data to Point Data'
cellDatatoPointData1 = CellDatatoPointData(Input=clip2)

# create a new 'Resample To Image'
resampleToImage1 = ResampleToImage(Input=cellDatatoPointData1)
resampleToImage1.UseInputBounds = 0
resampleToImage1.SamplingDimensions = [800, 100, 100]
resampleToImage1.SamplingBounds = [0.5, 5.8784823612399, -0.5, 0.5, -0.5, 0.5]

# create a new 'Stream Tracer'
streamTracer1 = StreamTracer(Input=resampleToImage1,
    SeedType='High Resolution Line Source')
streamTracer1.Vectors = ['POINTS', 'u.x']
streamTracer1.SurfaceStreamlines = 1
streamTracer1.MaximumStreamlineLength = 8.0

# init the 'High Resolution Line Source' selected for 'SeedType'
streamTracer1.SeedType.Point1 = [5.6784823612399, -0.5, 0.0]
streamTracer1.SeedType.Point2 = [5.6784823612399, 0.5, 0.0]
streamTracer1.SeedType.Resolution = 50

# create a new 'Contour'
contour2 = Contour(Input=resampleToImage1)
contour2.ContourBy = ['POINTS', '']
contour2.Isosurfaces = [0.5]
contour2.PointMergeMethod = 'Uniform Binning'

# create a new 'Contour'
contour3 = Contour(Input=resampleToImage1)
contour3.ContourBy = ['POINTS', '']
contour3.Isosurfaces = [-1.0, -2.0]
contour3.PointMergeMethod = 'Uniform Binning'

# create a new 'Slice'
slice2 = Slice(Input=resampleToImage1)
slice2.SliceType = 'Plane'
slice2.HyperTreeGridSlicer = 'Plane'
slice2.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice2.SliceType.Origin = [3.18924118061995, 0.0, 0.0]
slice2.SliceType.Normal = [0.0, 0.0, 1.0]

# init the 'Plane' selected for 'HyperTreeGridSlicer'
slice2.HyperTreeGridSlicer.Origin = [3.18924118061995, 0.0, 0.0]

# create a new 'Slice'
slice1 = Slice(Input=resampleToImage1)
slice1.SliceType = 'Plane'
slice1.HyperTreeGridSlicer = 'Plane'
slice1.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Origin = [2.0, 0.0, 0.0]

# init the 'Plane' selected for 'HyperTreeGridSlicer'
slice1.HyperTreeGridSlicer.Origin = [3.18924118061995, 0.0, 0.0]

# create a new 'Cylinder'
cylinder1 = Cylinder()
cylinder1.Resolution = 140
cylinder1.Height = 30.0
cylinder1.Capping = 0

# create a new 'Clip'
clip1 = Clip(Input=cylinder1)
clip1.ClipType = 'Plane'
clip1.HyperTreeGridClipper = 'Plane'
clip1.Scalars = ['POINTS', '']
clip1.Value = 1.0

# init the 'Plane' selected for 'ClipType'
clip1.ClipType.Normal = [0.0, 0.0, 1.0]

# create a new 'Transform'
transform1 = Transform(Input=clip1)
transform1.Transform = 'Transform'

# init the 'Transform' selected for 'Transform'
transform1.Transform.Translate = [15.0, 0.0, 0.0]
transform1.Transform.Rotate = [0.0, 0.0, 90.0]

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from transform1
transform1Display = Show(transform1, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
transform1Display.Representation = 'Surface'
transform1Display.AmbientColor = [1.0, 0.7843137254901961, 0.7529411764705882]
transform1Display.ColorArrayName = ['POINTS', '']
transform1Display.DiffuseColor = [1.0, 0.7843137254901961, 0.7529411764705882]
transform1Display.Opacity = 0.85
transform1Display.Specular = 1.0
transform1Display.Luminosity = 35.0
transform1Display.OSPRayUseScaleArray = 1
transform1Display.OSPRayScaleArray = 'Normals'
transform1Display.OSPRayScaleFunction = 'PiecewiseFunction'
transform1Display.OSPRayMaterial = 'copper'
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
transform1Display.ScalarOpacityUnitDistance = 6.650076732513133

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
transform1Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
transform1Display.ScaleTransferFunction.Points = [-2.220446049250313e-16, 0.0, 0.5, 0.0, 2.220446049250313e-16, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
transform1Display.OpacityTransferFunction.Points = [-2.220446049250313e-16, 0.0, 0.5, 0.0, 2.220446049250313e-16, 1.0, 0.5, 0.0]

# show data from slice2
slice2Display = Show(slice2, renderView1, 'GeometryRepresentation')

# get color transfer function/color map for 'ux'
uxLUT = GetColorTransferFunction('ux')
uxLUT.RGBPoints = [0.0, 0.0, 0.0, 0.5625, 1.6589486673251501, 0.0, 0.0, 1.0, 5.450838800773389, 0.0, 1.0, 1.0, 7.346780134859275, 0.5, 1.0, 0.5, 9.24272146894516, 1.0, 1.0, 0.0, 13.0346116023934, 1.0, 0.0, 0.0, 14.930552936479286, 0.5, 0.0, 0.0]
uxLUT.ColorSpace = 'RGB'
uxLUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
slice2Display.Representation = 'Surface'
slice2Display.ColorArrayName = ['POINTS', 'u.x']
slice2Display.LookupTable = uxLUT
slice2Display.OSPRayScaleArray = 'f'
slice2Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice2Display.SelectOrientationVectors = 'None'
slice2Display.ScaleFactor = 0.5378482341766357
slice2Display.SelectScaleArray = 'None'
slice2Display.GlyphType = 'Arrow'
slice2Display.GlyphTableIndexArray = 'None'
slice2Display.GaussianRadius = 0.026892411708831786
slice2Display.SetScaleArray = ['POINTS', 'f']
slice2Display.ScaleTransferFunction = 'PiecewiseFunction'
slice2Display.OpacityArray = ['POINTS', 'f']
slice2Display.OpacityTransferFunction = 'PiecewiseFunction'
slice2Display.DataAxesGrid = 'GridAxesRepresentation'
slice2Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
slice2Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
slice2Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0000000000000002, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
slice2Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0000000000000002, 1.0, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for uxLUT in view renderView1
uxLUTColorBar = GetScalarBar(uxLUT, renderView1)
uxLUTColorBar.Orientation = 'Horizontal'
uxLUTColorBar.WindowLocation = 'AnyLocation'
uxLUTColorBar.Position = [0.20866965620328853, 0.1]
uxLUTColorBar.Title = 'u'
uxLUTColorBar.ComponentTitle = 'Magnitude'
uxLUTColorBar.AutomaticLabelFormat = 0
uxLUTColorBar.LabelFormat = '%-#6.2g'
uxLUTColorBar.UseCustomLabels = 1
uxLUTColorBar.CustomLabels = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5, 12.0, 12.5]
uxLUTColorBar.RangeLabelFormat = '%6.2g'
uxLUTColorBar.ScalarBarLength = 0.7813153961136026

# set color bar visibility
uxLUTColorBar.Visibility = 1

# show color legend
slice2Display.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get opacity transfer function/opacity map for 'ux'
uxPWF = GetOpacityTransferFunction('ux')
uxPWF.Points = [2.079003474908902e-06, 0.0, 0.5, 0.0, 14.930552936479286, 1.0, 0.5, 0.0]
uxPWF.ScalarRangeInitialized = 1

# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(slice2)
# ----------------------------------------------------------------