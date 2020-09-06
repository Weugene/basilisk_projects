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
renderView1.ViewSize = [2150, 1306]
renderView1.InteractionMode = '2D'
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.CenterOfRotation = [15.0, 0.0, -0.0016834563011169923]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [8.536845092430669, 2.400563842200906, 5.224711891158429]
renderView1.CameraFocalPoint = [14.505741857352652, -0.13878497124792705, -0.4587552265465531]
renderView1.CameraViewUp = [0.24121031084380692, 0.9548725736708208, -0.17330883991850277]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 1.7889331297787596
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

# create a new 'Cylinder'
cylinder1 = Cylinder()
cylinder1.Resolution = 140
cylinder1.Height = 30.0
cylinder1.Capping = 0

# create a new 'PVD Reader'
dump2pvd_compressedpvd = PVDReader(FileName='/Users/weugene/basilisk/work/tube/dump2pvd_compressed.pvd')
dump2pvd_compressedpvd.CellArrays = ['p', 'fs', 'f', 'l', 'l2', 'omega', 'u.x']

# create a new 'Clip'
clip2 = Clip(Input=cylinder1)
clip2.ClipType = 'Plane'
clip2.HyperTreeGridClipper = 'Plane'
clip2.Scalars = ['POINTS', 'Normals_Magnitude']
clip2.Value = 0.9999999988389006

# init the 'Plane' selected for 'ClipType'
clip2.ClipType.Normal = [0.0, 0.0, 1.0]

# create a new 'Contour'
contour3 = Contour(Input=dump2pvd_compressedpvd)
contour3.ContourBy = ['POINTS', 'f']
contour3.Isosurfaces = [0.5]
contour3.PointMergeMethod = 'Uniform Binning'

# create a new 'Clip'
clip1 = Clip(Input=dump2pvd_compressedpvd)
clip1.ClipType = 'Box'
clip1.HyperTreeGridClipper = 'Plane'
clip1.Scalars = ['POINTS', 'f']
clip1.Value = 0.5

# init the 'Box' selected for 'ClipType'
clip1.ClipType.Position = [11.0, -0.5, -0.5]
clip1.ClipType.Length = [9.0, 1.0, 1.0]

# init the 'Plane' selected for 'HyperTreeGridClipper'
clip1.HyperTreeGridClipper.Origin = [15.0, 0.0, 0.0]

# create a new 'Cell Data to Point Data'
cellDatatoPointData1 = CellDatatoPointData(Input=clip1)
cellDatatoPointData1.CellDataArraytoprocess = ['f', 'fs', 'l', 'l2', 'omega', 'p', 'u.x']

# create a new 'Delaunay 3D'
delaunay3D1 = Delaunay3D(Input=cellDatatoPointData1)
delaunay3D1.AlphaVerts = 1

# create a new 'Contour'
contour1 = Contour(Input=delaunay3D1)
contour1.ContourBy = ['POINTS', 'f']
contour1.Isosurfaces = [0.5]
contour1.PointMergeMethod = 'Uniform Binning'

# create a new 'Contour'
contour2 = Contour(Input=delaunay3D1)
contour2.ContourBy = ['POINTS', 'l2']
contour2.Isosurfaces = [-1.0]
contour2.PointMergeMethod = 'Uniform Binning'

# create a new 'Transform'
transform1 = Transform(Input=clip2)
transform1.Transform = 'Transform'

# init the 'Transform' selected for 'Transform'
transform1.Transform.Translate = [15.0, 0.0, 0.0]
transform1.Transform.Rotate = [0.0, 0.0, 90.0]

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from contour2
contour2Display = Show(contour2, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
contour2Display.Representation = 'Surface'
contour2Display.AmbientColor = [0.5764705882352941, 0.7607843137254902, 1.0]
contour2Display.ColorArrayName = ['POINTS', '']
contour2Display.DiffuseColor = [0.5764705882352941, 0.7607843137254902, 1.0]
contour2Display.Specular = 1.0
contour2Display.Luminosity = 46.0
contour2Display.Ambient = 0.15
contour2Display.OSPRayScaleArray = 'l2'
contour2Display.OSPRayScaleFunction = 'PiecewiseFunction'
contour2Display.SelectOrientationVectors = 'None'
contour2Display.ScaleFactor = 0.6500720447488617
contour2Display.SelectScaleArray = 'l2'
contour2Display.GlyphType = 'Arrow'
contour2Display.GlyphTableIndexArray = 'l2'
contour2Display.GaussianRadius = 0.03250360223744308
contour2Display.SetScaleArray = ['POINTS', 'l2']
contour2Display.ScaleTransferFunction = 'PiecewiseFunction'
contour2Display.OpacityArray = ['POINTS', 'l2']
contour2Display.OpacityTransferFunction = 'PiecewiseFunction'
contour2Display.DataAxesGrid = 'GridAxesRepresentation'
contour2Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
contour2Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
contour2Display.ScaleTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, -0.9998779296875, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
contour2Display.OpacityTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, -0.9998779296875, 1.0, 0.5, 0.0]

# show data from transform1
transform1Display = Show(transform1, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
transform1Display.Representation = 'Surface'
transform1Display.ColorArrayName = [None, '']
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

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
transform1Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
transform1Display.ScaleTransferFunction.Points = [-2.220446049250313e-16, 0.0, 0.5, 0.0, 2.220446049250313e-16, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
transform1Display.OpacityTransferFunction.Points = [-2.220446049250313e-16, 0.0, 0.5, 0.0, 2.220446049250313e-16, 1.0, 0.5, 0.0]

# show data from contour3
contour3Display = Show(contour3, renderView1, 'GeometryRepresentation')

# get color transfer function/color map for 'f'
fLUT = GetColorTransferFunction('f')
fLUT.RGBPoints = [0.0, 0.23137254902, 0.298039215686, 0.752941176471, 0.5, 0.865, 0.865, 0.865, 1.0, 0.705882352941, 0.0156862745098, 0.149019607843]
fLUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
contour3Display.Representation = 'Surface'
contour3Display.ColorArrayName = ['POINTS', 'f']
contour3Display.LookupTable = fLUT
contour3Display.OSPRayScaleArray = 'f'
contour3Display.OSPRayScaleFunction = 'PiecewiseFunction'
contour3Display.SelectOrientationVectors = 'None'
contour3Display.ScaleFactor = 0.3379646740552927
contour3Display.SelectScaleArray = 'f'
contour3Display.GlyphType = 'Arrow'
contour3Display.GlyphTableIndexArray = 'f'
contour3Display.GaussianRadius = 0.016898233702764633
contour3Display.SetScaleArray = ['POINTS', 'f']
contour3Display.ScaleTransferFunction = 'PiecewiseFunction'
contour3Display.OpacityArray = ['POINTS', 'f']
contour3Display.OpacityTransferFunction = 'PiecewiseFunction'
contour3Display.DataAxesGrid = 'GridAxesRepresentation'
contour3Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
contour3Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
contour3Display.ScaleTransferFunction.Points = [0.5, 0.0, 0.5, 0.0, 0.5001220703125, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
contour3Display.OpacityTransferFunction.Points = [0.5, 0.0, 0.5, 0.0, 0.5001220703125, 1.0, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for fLUT in view renderView1
fLUTColorBar = GetScalarBar(fLUT, renderView1)
fLUTColorBar.WindowLocation = 'UpperRightCorner'
fLUTColorBar.Title = 'f'
fLUTColorBar.ComponentTitle = ''
fLUTColorBar.AutomaticLabelFormat = 0
fLUTColorBar.LabelFormat = '%-#6.5g'
fLUTColorBar.RangeLabelFormat = '%6.5g'
fLUTColorBar.ScalarBarLength = 0.3300000000000002

# set color bar visibility
fLUTColorBar.Visibility = 1

# show color legend
contour3Display.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get opacity transfer function/opacity map for 'f'
fPWF = GetOpacityTransferFunction('f')
fPWF.ScalarRangeInitialized = 1

# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(contour2)
# ----------------------------------------------------------------