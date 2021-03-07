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
renderView1.ViewSize = [1554, 584]
renderView1.InteractionMode = '2D'
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.CenterOfRotation = [8.0, 0.0, 0.0]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [3.3231870473370053, 0.10596421054029943, 2.254022822918342]
renderView1.CameraFocalPoint = [9.2336851975808, -0.1460561274791283, -4.340380850018227]
renderView1.CameraViewUp = [0.08105323969224698, 0.9961097802451542, 0.034578577694032595]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 1.1398365538798487
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.OSPRayMaterialLibrary = materialLibrary1

# init the 'GridAxes3DActor' selected for 'AxesGrid'
renderView1.AxesGrid.XAxisUseCustomLabels = 1
renderView1.AxesGrid.XAxisLabels = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0]
renderView1.AxesGrid.YAxisUseCustomLabels = 1
renderView1.AxesGrid.YAxisLabels = [-1.0, 0.0, 0.0]

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

# create a new 'Cylinder'
cylinder1 = Cylinder()
cylinder1.Resolution = 183
cylinder1.Height = 30.0
cylinder1.Capping = 0

# create a new 'PVD Reader'
tube_bp_from_dumppvd = PVDReader(FileName='/Volumes/GoogleDrive/My Drive/Skoltech/Dumps/tube/10e_coarse/tube_bp_from_dump.pvd')
tube_bp_from_dumppvd.CellArrays = ['p', 'fs', 'f', 'l', 'residual_of_p', 'l2', 'omega', 'u.x']

# create a new 'Clip'
clip1 = Clip(Input=tube_bp_from_dumppvd)
clip1.ClipType = 'Box'
clip1.HyperTreeGridClipper = 'Plane'
clip1.Scalars = ['POINTS', 'f']
clip1.Value = 0.5

# init the 'Box' selected for 'ClipType'
clip1.ClipType.Position = [4.0, -0.55, -0.55]
clip1.ClipType.Length = [8.0, 1.1, 1.1]

# init the 'Plane' selected for 'HyperTreeGridClipper'
clip1.HyperTreeGridClipper.Origin = [15.0, 0.0, 0.0]

# create a new 'Cell Data to Point Data'
cellDatatoPointData1 = CellDatatoPointData(Input=clip1)
cellDatatoPointData1.CellDataArraytoprocess = ['f', 'f', 'fs', 'l', 'l2', 'omega', 'p', 'residual_of_p', 'u.x']

# create a new 'Resample To Image'
resampleToImage1 = ResampleToImage(Input=cellDatatoPointData1)
resampleToImage1.UseInputBounds = 0
resampleToImage1.SamplingDimensions = [800, 100, 100]
resampleToImage1.SamplingBounds = [4.3, 11.9, -0.5, 0.5, -0.5, 0.5]

# create a new 'Contour'
contour2 = Contour(Input=resampleToImage1)
contour2.ContourBy = ['POINTS', 'l2']
contour2.Isosurfaces = [-1.0, -2.0]
contour2.PointMergeMethod = 'Uniform Binning'

# create a new 'Stream Tracer'
streamTracer1 = StreamTracer(Input=resampleToImage1,
    SeedType='High Resolution Line Source')
streamTracer1.Vectors = ['POINTS', 'u.x']
streamTracer1.SurfaceStreamlines = 1
streamTracer1.MaximumStreamlineLength = 7.890000000000002

# init the 'High Resolution Line Source' selected for 'SeedType'
streamTracer1.SeedType.Point1 = [9.3, -0.5, 0.0]
streamTracer1.SeedType.Point2 = [9.3, 0.5, 0.0]
streamTracer1.SeedType.Resolution = 50

# create a new 'Contour'
contour1 = Contour(Input=resampleToImage1)
contour1.ContourBy = ['POINTS', 'f']
contour1.Isosurfaces = [0.5]
contour1.PointMergeMethod = 'Uniform Binning'

# create a new 'Contour'
contour3 = Contour(Input=clip1)
contour3.ContourBy = ['POINTS', 'f']
contour3.Isosurfaces = [0.5]
contour3.PointMergeMethod = 'Uniform Binning'

# create a new 'Clip'
clip2 = Clip(Input=cylinder1)
clip2.ClipType = 'Plane'
clip2.HyperTreeGridClipper = 'Plane'
clip2.Scalars = ['POINTS', 'Normals_Magnitude']
clip2.Value = 1.000000001583269

# init the 'Plane' selected for 'ClipType'
clip2.ClipType.Normal = [0.0, 0.0, 1.0]

# init the 'Plane' selected for 'HyperTreeGridClipper'
clip2.HyperTreeGridClipper.Origin = [3.6835670471191406e-05, 0.0, 0.0]

# create a new 'Transform'
transform1 = Transform(Input=clip2)
transform1.Transform = 'Transform'

# init the 'Transform' selected for 'Transform'
transform1.Transform.Translate = [15.0, 0.0, 0.0]
transform1.Transform.Rotate = [0.0, 0.0, 89.99999999999999]
transform1.Transform.Scale = [1.0000000000000027, 1.0000000000000036, 1.0000000000000018]

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from resampleToImage1
resampleToImage1Display = Show(resampleToImage1, renderView1, 'UniformGridRepresentation')

# get color transfer function/color map for 'ux'
uxLUT = GetColorTransferFunction('ux')
uxLUT.RGBPoints = [0.0, 0.0, 0.0, 0.5625, 0.23617269686646322, 0.0, 0.0, 1.0, 0.7759970667680063, 0.0, 1.0, 1.0, 1.0459087203296784, 0.5, 1.0, 0.5, 1.3158203738913508, 1.0, 1.0, 0.0, 1.855644743792894, 1.0, 0.0, 0.0, 2.125556397354566, 0.5, 0.0, 0.0]
uxLUT.ColorSpace = 'RGB'
uxLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'ux'
uxPWF = GetOpacityTransferFunction('ux')
uxPWF.Points = [0.0, 0.0, 0.5, 0.0, 2.125556397354566, 1.0, 0.5, 0.0]
uxPWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
resampleToImage1Display.Representation = 'Volume'
resampleToImage1Display.ColorArrayName = ['POINTS', 'u.x']
resampleToImage1Display.LookupTable = uxLUT
resampleToImage1Display.OSPRayScaleArray = 'f'
resampleToImage1Display.OSPRayScaleFunction = 'PiecewiseFunction'
resampleToImage1Display.SelectOrientationVectors = 'None'
resampleToImage1Display.ScaleFactor = 0.8
resampleToImage1Display.SelectScaleArray = 'None'
resampleToImage1Display.GlyphType = 'Arrow'
resampleToImage1Display.GlyphTableIndexArray = 'None'
resampleToImage1Display.GaussianRadius = 0.04
resampleToImage1Display.SetScaleArray = ['POINTS', 'f']
resampleToImage1Display.ScaleTransferFunction = 'PiecewiseFunction'
resampleToImage1Display.OpacityArray = ['POINTS', 'f']
resampleToImage1Display.OpacityTransferFunction = 'PiecewiseFunction'
resampleToImage1Display.DataAxesGrid = 'GridAxesRepresentation'
resampleToImage1Display.PolarAxes = 'PolarAxesRepresentation'
resampleToImage1Display.ScalarOpacityUnitDistance = 2.324564821638666
resampleToImage1Display.ScalarOpacityFunction = uxPWF
resampleToImage1Display.Shade = 1
resampleToImage1Display.SliceFunction = 'Plane'
resampleToImage1Display.Slice = 49

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
resampleToImage1Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# init the 'Plane' selected for 'SliceFunction'
resampleToImage1Display.SliceFunction.Origin = [8.0, 0.0, 0.0]

# show data from contour1
contour1Display = Show(contour1, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
contour1Display.Representation = 'Surface'
contour1Display.AmbientColor = [0.0392156862745098, 0.00784313725490196, 1.0]
contour1Display.ColorArrayName = ['POINTS', '']
contour1Display.DiffuseColor = [0.0392156862745098, 0.00784313725490196, 1.0]
contour1Display.Opacity = 0.23
contour1Display.Specular = 1.0
contour1Display.SpecularPower = 1.0
contour1Display.Ambient = 0.21
contour1Display.OSPRayScaleArray = 'f'
contour1Display.OSPRayScaleFunction = 'PiecewiseFunction'
contour1Display.SelectOrientationVectors = 'None'
contour1Display.ScaleFactor = 0.8
contour1Display.SelectScaleArray = 'f'
contour1Display.GlyphType = 'Arrow'
contour1Display.GlyphTableIndexArray = 'f'
contour1Display.GaussianRadius = 0.04
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

# show data from transform1
transform1Display = Show(transform1, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
transform1Display.Representation = 'Surface'
transform1Display.AmbientColor = [1.0, 0.7843137254901961, 0.7529411764705882]
transform1Display.ColorArrayName = [None, '']
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

# show data from streamTracer1
streamTracer1Display = Show(streamTracer1, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
streamTracer1Display.Representation = 'Surface'
streamTracer1Display.ColorArrayName = ['POINTS', 'u.x']
streamTracer1Display.LookupTable = uxLUT
streamTracer1Display.OSPRayScaleArray = 'AngularVelocity'
streamTracer1Display.OSPRayScaleFunction = 'PiecewiseFunction'
streamTracer1Display.SelectOrientationVectors = 'Normals'
streamTracer1Display.ScaleFactor = 0.7899999618530273
streamTracer1Display.SelectScaleArray = 'AngularVelocity'
streamTracer1Display.GlyphType = 'Arrow'
streamTracer1Display.GlyphTableIndexArray = 'AngularVelocity'
streamTracer1Display.GaussianRadius = 0.03949999809265137
streamTracer1Display.SetScaleArray = ['POINTS', 'AngularVelocity']
streamTracer1Display.ScaleTransferFunction = 'PiecewiseFunction'
streamTracer1Display.OpacityArray = ['POINTS', 'AngularVelocity']
streamTracer1Display.OpacityTransferFunction = 'PiecewiseFunction'
streamTracer1Display.DataAxesGrid = 'GridAxesRepresentation'
streamTracer1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
streamTracer1Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
streamTracer1Display.ScaleTransferFunction.Points = [-9.447468533227408, 0.0, 0.5, 0.0, 7.903744217703089, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
streamTracer1Display.OpacityTransferFunction.Points = [-9.447468533227408, 0.0, 0.5, 0.0, 7.903744217703089, 1.0, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for uxLUT in view renderView1
uxLUTColorBar = GetScalarBar(uxLUT, renderView1)
uxLUTColorBar.Orientation = 'Horizontal'
uxLUTColorBar.WindowLocation = 'AnyLocation'
uxLUTColorBar.Position = [0.5581150358180519, 0.05196141809040081]
uxLUTColorBar.Title = 'u'
uxLUTColorBar.ComponentTitle = 'Magnitude'
uxLUTColorBar.LabelFormat = '%-#6.2g'
uxLUTColorBar.RangeLabelFormat = '%6.2g'
uxLUTColorBar.ScalarBarLength = 0.3300000000000003

# set color bar visibility
uxLUTColorBar.Visibility = 1

# show color legend
resampleToImage1Display.SetScalarBarVisibility(renderView1, True)

# show color legend
streamTracer1Display.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(contour1)
# ----------------------------------------------------------------