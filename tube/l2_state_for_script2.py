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
import sys
import timeit

start = timeit.default_timer()
# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView2 = CreateView('RenderView')
renderView2.ViewSize = [1832, 928]
renderView2.InteractionMode = '2D'
renderView2.AxesGrid = 'GridAxes3DActor'
renderView2.OrientationAxesVisibility = 0
renderView2.CenterOfRotation = [9.71011229812088, -0.0002944310094490643, -0.0019979823491385884]
renderView2.StereoType = 'Crystal Eyes'
renderView2.CameraPosition = [5.691130662792347, 1.317963730154289, 3.1554448472414376]
renderView2.CameraFocalPoint = [10.591651496465168, -0.5950607325808129, -1.2580984092708607]
renderView2.CameraViewUp = [0.14111261524587124, 0.9558815502523825, -0.2576386844120518]
renderView2.CameraFocalDisk = 1.0
renderView2.CameraParallelScale = 1.297832166417713
renderView2.BackEnd = 'OSPRay raycaster'
renderView2.OSPRayMaterialLibrary = materialLibrary1

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'

# create new layout object 'Layout #1'
layout1_1 = CreateLayout(name='Layout #1')
layout1_1.AssignView(0, renderView2)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView2)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'Cylinder'
cylinder1 = Cylinder()
cylinder1.Resolution = 112
cylinder1.Height = 30.0
cylinder1.Capping = 0
print("cylinder1...")
# create a new 'PVD Reader'
print("next tube_bp_from_dumppvd...")
tube_bp_from_dumppvd = PVDReader(FileName='/Volumes/GoogleDrive/My Drive/Skoltech/Dumps/tube/10e_coarse/tube_bp_from_dump.pvd')
tube_bp_from_dumppvd.CellArrays = ['p', 'fs', 'f', 'l', 'residual_of_p', 'l2', 'omega', 'u.x']
print("tube_bp_from_dumppvd...")


animationScene1 = GetAnimationScene()
# get the time-keeper
timeKeeper1 = GetTimeKeeper()
# Properties modified on animationScene1
animationScene1.AnimationTime = 4.989
# Properties modified on timeKeeper1
timeKeeper1.Time = 4.989

# create a new 'Contour'
contour2 = Contour(Input=tube_bp_from_dumppvd)
contour2.ContourBy = ['POINTS', 'f']
contour2.Isosurfaces = [0.5]
contour2.PointMergeMethod = 'Uniform Binning'
print("contour2... next clip1")
# create a new 'Clip'
clip1 = Clip(Input=tube_bp_from_dumppvd)
clip1.ClipType = 'Box'
clip1.HyperTreeGridClipper = 'Plane'
clip1.Scalars = ['POINTS', 'f']
clip1.Value = 0.5

# init the 'Box' selected for 'ClipType'
clip1.ClipType.Position = [4.0, -0.5, -0.5]
clip1.ClipType.Length = [8.0, 1.0, 1.0]

# init the 'Plane' selected for 'HyperTreeGridClipper'
clip1.HyperTreeGridClipper.Origin = [15.0, 0.0, 0.0]
print("clip1... next clip2")
# create a new 'Clip'
clip2 = Clip(Input=cylinder1)
clip2.ClipType = 'Plane'
clip2.HyperTreeGridClipper = 'Plane'
clip2.Scalars = ['POINTS', '']
clip2.Value = 0.9999999973129153

# init the 'Plane' selected for 'ClipType'
clip2.ClipType.Normal = [0.0, 0.0, 1.0]
print("clip2 next is cellDatatoPointData1...")
# create a new 'Cell Data to Point Data'
cellDatatoPointData1 = CellDatatoPointData(Input=clip1)
cellDatatoPointData1.CellDataArraytoprocess = ['f', 'f', 'fs', 'l', 'l2', 'omega', 'p', 'residual_of_p', 'u.x']
print("cellDatatoPointData1 next is delaunay3D1...")
# create a new 'Delaunay 3D'
delaunay3D1 = Delaunay3D(Input=cellDatatoPointData1)
print("delaunay3D1 next contour3...")
# create a new 'Contour'
contour3 = Contour(Input=delaunay3D1)
contour3.ContourBy = ['POINTS', 'l2']
contour3.Isosurfaces = [-1.0]
contour3.PointMergeMethod = 'Uniform Binning'

print("contour3 next contour1...")
# create a new 'Contour'
contour1 = Contour(Input=delaunay3D1)
contour1.ContourBy = ['POINTS', 'f']
contour1.Isosurfaces = [0.5]
contour1.PointMergeMethod = 'Uniform Binning'
print("contour1... next contour1Display")

# trace defaults for the display properties.
# tube_bp_from_dumppvdDisplay.Representation = 'Surface'
# tube_bp_from_dumppvdDisplay.ColorArrayName = [None, '']
# tube_bp_from_dumppvdDisplay.Opacity = 0.85
# tube_bp_from_dumppvdDisplay.OSPRayScaleArray = 'f'
# tube_bp_from_dumppvdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
# tube_bp_from_dumppvdDisplay.SelectOrientationVectors = 'None'
# tube_bp_from_dumppvdDisplay.ScaleFactor = 3.0
# tube_bp_from_dumppvdDisplay.SelectScaleArray = 'None'
# tube_bp_from_dumppvdDisplay.GlyphType = 'Arrow'
# tube_bp_from_dumppvdDisplay.GlyphTableIndexArray = 'None'
# tube_bp_from_dumppvdDisplay.GaussianRadius = 0.15
# tube_bp_from_dumppvdDisplay.SetScaleArray = ['POINTS', 'f']
# tube_bp_from_dumppvdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
# tube_bp_from_dumppvdDisplay.OpacityArray = ['POINTS', 'f']
# tube_bp_from_dumppvdDisplay.OpacityTransferFunction = 'PiecewiseFunction'
# tube_bp_from_dumppvdDisplay.DataAxesGrid = 'GridAxesRepresentation'
# tube_bp_from_dumppvdDisplay.PolarAxes = 'PolarAxesRepresentation'
# tube_bp_from_dumppvdDisplay.ScalarOpacityUnitDistance = 0.48577906707444607

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
# tube_bp_from_dumppvdDisplay.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# ----------------------------------------------------------------
# setup the visualization in view 'renderView2'
# ----------------------------------------------------------------
# print ("tube_bp_from_dumppvdDisplay... next Show contour1Display")
# show data from contour1
contour1Display = Show(contour1, renderView2, 'GeometryRepresentation')

# get color transfer function/color map for 'ux'
uxLUT = GetColorTransferFunction('ux')
uxLUT.RGBPoints = [0.0, 0.0, 0.0, 0.5625, 0.23363270937628142, 0.0, 0.0, 1.0, 0.7676513821560252, 0.0, 1.0, 1.0, 1.0346601928717756, 0.5, 1.0, 0.5, 1.3016690035875256, 1.0, 1.0, 0.0, 1.8356876763672694, 1.0, 0.0, 0.0, 2.1026964870830196, 0.5, 0.0, 0.0]
uxLUT.ColorSpace = 'RGB'
uxLUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
contour1Display.Representation = 'Surface'
contour1Display.ColorArrayName = ['POINTS', 'u.x']
contour1Display.LookupTable = uxLUT
contour1Display.Specular = 1.0
contour1Display.OSPRayScaleArray = 'f'
contour1Display.OSPRayScaleFunction = 'PiecewiseFunction'
contour1Display.SelectOrientationVectors = 'None'
contour1Display.ScaleFactor = 0.3379618165202066
contour1Display.SelectScaleArray = 'f'
contour1Display.GlyphType = 'Arrow'
contour1Display.GlyphTableIndexArray = 'f'
contour1Display.GaussianRadius = 0.01689809082601033
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
print ("contour1Display... next Show clip2Display")

# show data from clip2
clip2Display = Show(clip2, renderView2, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
clip2Display.Representation = 'Surface'
clip2Display.ColorArrayName = ['POINTS', '']
clip2Display.Opacity = 0.85
clip2Display.Position = [15.0, 0.0, 0.0]
clip2Display.Orientation = [0.0, 0.0, 90.0]
clip2Display.OSPRayScaleArray = 'Normals'
clip2Display.OSPRayScaleFunction = 'PiecewiseFunction'
clip2Display.SelectOrientationVectors = 'None'
clip2Display.ScaleFactor = 3.0
clip2Display.SelectScaleArray = 'None'
clip2Display.GlyphType = 'Arrow'
clip2Display.GlyphTableIndexArray = 'None'
clip2Display.GaussianRadius = 0.15
clip2Display.SetScaleArray = ['POINTS', 'Normals']
clip2Display.ScaleTransferFunction = 'PiecewiseFunction'
clip2Display.OpacityArray = ['POINTS', 'Normals']
clip2Display.OpacityTransferFunction = 'PiecewiseFunction'
clip2Display.DataAxesGrid = 'GridAxesRepresentation'
clip2Display.PolarAxes = 'PolarAxesRepresentation'
clip2Display.ScalarOpacityUnitDistance = 7.800654093765651

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
clip2Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
clip2Display.ScaleTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
clip2Display.OpacityTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
clip2Display.PolarAxes.Translation = [15.0, 0.0, 0.0]
clip2Display.PolarAxes.Orientation = [0.0, 0.0, 90.0]

# show data from contour3
contour3Display = Show(contour3, renderView2, 'GeometryRepresentation')

# trace defaults for the display properties.
contour3Display.Representation = 'Surface'
contour3Display.AmbientColor = [0.5725490196078431, 0.6509803921568628, 1.0]
contour3Display.ColorArrayName = ['POINTS', '']
contour3Display.DiffuseColor = [0.5725490196078431, 0.6509803921568628, 1.0]
contour3Display.Specular = 1.0
contour3Display.OSPRayScaleArray = 'l2'
contour3Display.OSPRayScaleFunction = 'PiecewiseFunction'
contour3Display.SelectOrientationVectors = 'None'
contour3Display.ScaleFactor = 0.6063199703188431
contour3Display.SelectScaleArray = 'l2'
contour3Display.GlyphType = 'Arrow'
contour3Display.GlyphTableIndexArray = 'l2'
contour3Display.GaussianRadius = 0.030315998515942157
contour3Display.SetScaleArray = ['POINTS', 'l2']
contour3Display.ScaleTransferFunction = 'PiecewiseFunction'
contour3Display.OpacityArray = ['POINTS', 'l2']
contour3Display.OpacityTransferFunction = 'PiecewiseFunction'
contour3Display.DataAxesGrid = 'GridAxesRepresentation'
contour3Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
contour3Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
contour3Display.ScaleTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, -0.9998779296875, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
contour3Display.OpacityTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, -0.9998779296875, 1.0, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for uxLUT in view renderView2
uxLUTColorBar = GetScalarBar(uxLUT, renderView2)
uxLUTColorBar.Orientation = 'Horizontal'
uxLUTColorBar.WindowLocation = 'AnyLocation'
uxLUTColorBar.Position = [0.48467248908296934, 0.13295880149812733]
uxLUTColorBar.Title = 'u'
uxLUTColorBar.ComponentTitle = 'Magnitude'
uxLUTColorBar.AutomaticLabelFormat = 0
uxLUTColorBar.LabelFormat = '%-#6.2g'
uxLUTColorBar.RangeLabelFormat = '%6.2g'
uxLUTColorBar.ScalarBarLength = 0.4293449781659392

# set color bar visibility
uxLUTColorBar.Visibility = 1

# show color legend
contour1Display.SetScalarBarVisibility(renderView2, True)

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get opacity transfer function/opacity map for 'ux'
uxPWF = GetOpacityTransferFunction('ux')
uxPWF.Points = [0.0, 0.0, 0.5, 0.0, 2.1026964870830196, 1.0, 0.5, 0.0]
uxPWF.ScalarRangeInitialized = 1

# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(contour1)
# ----------------------------------------------------------------
Show()
Render()
SaveScreenshot( "MY_TEST1.png", renderView2,
#      ImageResolution=[2316, 2204],
    TransparentBackground=0,
#     magnification=1,
    # PNG options
    CompressionLevel='2' )
print('File=' + "MY_TEST1.png" + ' generated succesfully')

stop = timeit.default_timer()
total_time = stop - start

# output running time in a nice format.
mins, secs = divmod(total_time, 60)
hours, mins = divmod(mins, 60)

sys.stdout.write("Total running time: %d:%d:%d.\n" % (hours, mins, secs))