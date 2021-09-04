# state file generated using paraview version 5.9.0

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1840, 1156]
renderView1.InteractionMode = '2D'
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.CenterOfRotation = [3.9594372510910034, 0.00012353062629699707, 0.0003523975610733032]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [-0.040512771139578935, 1.290158870335943, 5.248223688485809]
renderView1.CameraFocalPoint = [-0.023208475832503923, 1.284848907700422, 5.227816101479906]
renderView1.CameraViewUp = [0.17432562971565893, 0.9788702853248832, -0.10688095869808088]
renderView1.CameraViewAngle = 15.42391304347826
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 1.1256051366537485
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.OSPRayMaterialLibrary = materialLibrary1

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.AssignView(0, renderView1)
layout1.SetSize(1840, 1156)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'Cylinder'
cylinder1 = Cylinder(registrationName='Cylinder1')
cylinder1.Resolution = 28
cylinder1.Height = 30.0
cylinder1.Capping = 0

# create a new 'Transform'
transform1 = Transform(registrationName='Transform1', Input=cylinder1)
transform1.Transform = 'Transform'

# init the 'Transform' selected for 'Transform'
transform1.Transform.Rotate = [0.0, 0.0, 90.0]

# create a new 'Clip'
clip1 = Clip(registrationName='Clip1', Input=transform1)
clip1.ClipType = 'Plane'
clip1.HyperTreeGridClipper = 'Plane'
clip1.Scalars = ['POINTS', '']

# init the 'Plane' selected for 'ClipType'
clip1.ClipType.Normal = [0.0, 0.0, 1.0]

# create a new 'PVD Reader'
save_isosurfacepvd = PVDReader(registrationName='save_isosurface.pvd', FileName='/Users/weugene/basilisk/work/tube/res23/save_isosurface.pvd')
save_isosurfacepvd.CellArrays = []
save_isosurfacepvd.PointArrays = ['Result', 'absOmega', 'f', 'p', 'u.x']

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from save_isosurfacepvd
save_isosurfacepvdDisplay = Show(save_isosurfacepvd, renderView1, 'GeometryRepresentation')

# get color transfer function/color map for 'ux'
uxLUT = GetColorTransferFunction('ux')
uxLUT.RGBPoints = [0.003649621564209849, 0.0, 0.0, 0.5625, 0.24729264204713047, 0.0, 0.0, 1.0, 0.8041920709742089, 0.0, 1.0, 1.0, 1.082641237240404, 0.5, 1.0, 0.5, 1.3610904035065987, 1.0, 1.0, 0.0, 1.9179898324336777, 1.0, 0.0, 0.0, 2.1964389986998722, 0.5, 0.0, 0.0]
uxLUT.ColorSpace = 'RGB'
uxLUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
save_isosurfacepvdDisplay.Representation = 'Surface'
save_isosurfacepvdDisplay.ColorArrayName = ['POINTS', 'u.x']
save_isosurfacepvdDisplay.LookupTable = uxLUT
save_isosurfacepvdDisplay.Specular = 1.0
save_isosurfacepvdDisplay.SelectTCoordArray = 'None'
save_isosurfacepvdDisplay.SelectNormalArray = 'None'
save_isosurfacepvdDisplay.SelectTangentArray = 'None'
save_isosurfacepvdDisplay.OSPRayScaleArray = 'Result'
save_isosurfacepvdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
save_isosurfacepvdDisplay.SelectOrientationVectors = 'None'
save_isosurfacepvdDisplay.ScaleFactor = 0.31534409523010254
save_isosurfacepvdDisplay.SelectScaleArray = 'Result'
save_isosurfacepvdDisplay.GlyphType = 'Arrow'
save_isosurfacepvdDisplay.GlyphTableIndexArray = 'Result'
save_isosurfacepvdDisplay.GaussianRadius = 0.015767204761505126
save_isosurfacepvdDisplay.SetScaleArray = ['POINTS', 'Result']
save_isosurfacepvdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
save_isosurfacepvdDisplay.OpacityArray = ['POINTS', 'Result']
save_isosurfacepvdDisplay.OpacityTransferFunction = 'PiecewiseFunction'
save_isosurfacepvdDisplay.DataAxesGrid = 'GridAxesRepresentation'
save_isosurfacepvdDisplay.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
save_isosurfacepvdDisplay.ScaleTransferFunction.Points = [0.005237827916105199, 0.0, 0.5, 0.0, 0.46907547386229526, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
save_isosurfacepvdDisplay.OpacityTransferFunction.Points = [0.005237827916105199, 0.0, 0.5, 0.0, 0.46907547386229526, 1.0, 0.5, 0.0]

# show data from clip1
clip1Display = Show(clip1, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
clip1Display.Representation = 'Surface'
clip1Display.AmbientColor = [1.0, 0.7254901960784313, 0.6549019607843137]
clip1Display.ColorArrayName = [None, '']
clip1Display.DiffuseColor = [1.0, 0.7254901960784313, 0.6549019607843137]
clip1Display.Opacity = 0.78
clip1Display.Specular = 1.0
clip1Display.Luminosity = 100.0
clip1Display.Ambient = 0.11
clip1Display.SelectTCoordArray = 'TCoords'
clip1Display.SelectNormalArray = 'Normals'
clip1Display.SelectTangentArray = 'None'
clip1Display.OSPRayScaleArray = 'Normals'
clip1Display.OSPRayScaleFunction = 'PiecewiseFunction'
clip1Display.SelectOrientationVectors = 'None'
clip1Display.ScaleFactor = 3.0
clip1Display.SelectScaleArray = 'None'
clip1Display.GlyphType = 'Arrow'
clip1Display.GlyphTableIndexArray = 'None'
clip1Display.GaussianRadius = 0.15
clip1Display.SetScaleArray = ['POINTS', 'Normals']
clip1Display.ScaleTransferFunction = 'PiecewiseFunction'
clip1Display.OpacityArray = ['POINTS', 'Normals']
clip1Display.OpacityTransferFunction = 'PiecewiseFunction'
clip1Display.DataAxesGrid = 'GridAxesRepresentation'
clip1Display.PolarAxes = 'PolarAxesRepresentation'
clip1Display.ScalarOpacityUnitDistance = 12.172848562792895
clip1Display.OpacityArrayName = ['POINTS', 'Normals']

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
clip1Display.ScaleTransferFunction.Points = [-2.220446049250313e-16, 0.0, 0.5, 0.0, 2.220446049250313e-16, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
clip1Display.OpacityTransferFunction.Points = [-2.220446049250313e-16, 0.0, 0.5, 0.0, 2.220446049250313e-16, 1.0, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for uxLUT in view renderView1
uxLUTColorBar = GetScalarBar(uxLUT, renderView1)
uxLUTColorBar.Orientation = 'Horizontal'
uxLUTColorBar.WindowLocation = 'AnyLocation'
uxLUTColorBar.Position = [0.5700679347826086, 0.08650519031141868]
uxLUTColorBar.Title = '|u|'
uxLUTColorBar.ComponentTitle = ''
uxLUTColorBar.ScalarBarLength = 0.33000000000000007
uxLUTColorBar.CustomLabels = np.ceil(np.linspace(0, range_max, 5))
# set color bar visibility
uxLUTColorBar.Visibility = 1

# show color legend
save_isosurfacepvdDisplay.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get opacity transfer function/opacity map for 'ux'
uxPWF = GetOpacityTransferFunction('ux')
uxPWF.Points = [0.003649621564209849, 0.0, 0.5, 0.0, 2.1964389986998722, 1.0, 0.5, 0.0]
uxPWF.ScalarRangeInitialized = 1

# ----------------------------------------------------------------
# restore active source
SetActiveSource(clip1)
# ----------------------------------------------------------------


if __name__ == '__main__':
    # generate extracts
    SaveExtracts(ExtractsOutputDirectory='extracts')