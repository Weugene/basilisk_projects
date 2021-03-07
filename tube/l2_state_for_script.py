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
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1832, 1068]
renderView1.InteractionMode = '2D'
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.CenterOfRotation = [16.132638154732103, 0, 0]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [13.184119175984472, 3.133389244584978, 3.8710475868614265]
renderView1.CameraFocalPoint = [15.569861393706825, -0.1333933342774375, -0.21935954198731383]
renderView1.CameraViewUp = [0.22487951548838764, 0.821105210650371, -0.5246097945678463]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 1.5893230593172707
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

tube_bp_from_dumppvd = GetActiveSource()
# create a new 'Cylinder'
cylinder1 = Cylinder()
cylinder1.Resolution = 112
cylinder1.Height = 30.0
cylinder1.Capping = 0

# create a new 'PVD Reader'
#tube_bp_from_dumppvd = PVDReader(FileName='/Volumes/GoogleDrive/My Drive/Skoltech/Dumps/tube/10e_coarse/tube_bp_from_dump.pvd')
#tube_bp_from_dumppvd.CellArrays = add_attribute()
#tube_bp_from_dumppvd.CellArrays = ['p', 'fs', 'f', 'l', 'residual_of_p', 'l2', 'omega', 'u.x']

# create a new 'Clip'
clip1 = Clip(Input=tube_bp_from_dumppvd)
clip1.ClipType = 'Box'
clip1.HyperTreeGridClipper = 'Plane'
clip1.Scalars = ['POINTS', 'f']
clip1.Value = 0.5

# init the 'Box' selected for 'ClipType'
# clip1.ClipType.Position = [12.0, -0.5, -0.5]
# clip1.ClipType.Length = [8.0, 1.0, 1.0]
clip1.ClipType.Position = [4.0, -0.5, -0.5]
clip1.ClipType.Length = [9.0, 1.0, 1.0]

# init the 'Plane' selected for 'HyperTreeGridClipper'
# clip1.HyperTreeGridClipper.Origin = [15.0, 0.0, 0.0]
clip1.HyperTreeGridClipper.Origin = [9.0, 0.0, 0.0]

# create a new 'Cell Data to Point Data'
cellDatatoPointData1 = CellDatatoPointData(Input=clip1)
cellDatatoPointData1.CellDataArraytoprocess = ['f', 'fs', 'l', 'l2', 'omega', 'p', 'residual_of_p', 'u.x']

# create a new 'Delaunay 3D'
delaunay3D1 = Delaunay3D(Input=cellDatatoPointData1)

# create a new 'Contour'
contour2 = Contour(Input=delaunay3D1)
contour2.ContourBy = ['POINTS', 'f']
contour2.Isosurfaces = [0.5]
contour2.PointMergeMethod = 'Uniform Binning'

# create a new 'Contour'
contour1 = Contour(Input=tube_bp_from_dumppvd)
contour1.ContourBy = ['POINTS', 'f']
contour1.Isosurfaces = [0.5]
contour1.PointMergeMethod = 'Uniform Binning'
# contour1Display = Show(contour1, renderView1, 'GeometryRepresentation')
# # hide data in view
# Hide(contour1, renderView1)

# create a new 'Clip'
clip2 = Clip(Input=cylinder1)
clip2.ClipType = 'Plane'
clip2.HyperTreeGridClipper = 'Plane'
clip2.Scalars = ['POINTS', 'Normals_Magnitude']
clip2.Value = 1

# init the 'Plane' selected for 'ClipType'
clip2.ClipType.Normal = [0.0, 0.0, 1.0]

# create a new 'Contour'
contour3 = Contour(Input=delaunay3D1)
contour3.ContourBy = ['POINTS', 'l2']
contour3.Isosurfaces = [-1.0]
contour3.PointMergeMethod = 'Uniform Binning'

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from tube_bp_from_dumppvd
tube_bp_from_dumppvdDisplay = Show(tube_bp_from_dumppvd, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
tube_bp_from_dumppvdDisplay.Representation = 'Outline'
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

# show data from contour1
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

# show data from clip1
clip1Display = Show(clip1, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
clip1Display.Representation = 'Surface'
clip1Display.ColorArrayName = [None, '']
clip1Display.Opacity = 0.85
clip1Display.OSPRayScaleArray = 'f'
clip1Display.OSPRayScaleFunction = 'PiecewiseFunction'
clip1Display.SelectOrientationVectors = 'None'
clip1Display.ScaleFactor = 3.0
clip1Display.SelectScaleArray = 'None'
clip1Display.GlyphType = 'Arrow'
clip1Display.GlyphTableIndexArray = 'None'
clip1Display.GaussianRadius = 0.15
clip1Display.SetScaleArray = ['POINTS', 'f']
clip1Display.ScaleTransferFunction = 'PiecewiseFunction'
clip1Display.OpacityArray = ['POINTS', 'f']
clip1Display.OpacityTransferFunction = 'PiecewiseFunction'
clip1Display.DataAxesGrid = 'GridAxesRepresentation'
clip1Display.PolarAxes = 'PolarAxesRepresentation'
clip1Display.ScalarOpacityUnitDistance = 0.29672181428754485

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
clip1Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# show data from cellDatatoPointData1
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
cellDatatoPointData1Display.ScalarOpacityUnitDistance = 0.29672181428754485

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
cellDatatoPointData1Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# show data from delaunay3D1
delaunay3D1Display = Show(delaunay3D1, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
delaunay3D1Display.Representation = 'Surface'
delaunay3D1Display.ColorArrayName = [None, '']
delaunay3D1Display.Opacity = 0.85
delaunay3D1Display.OSPRayScaleArray = 'f'
delaunay3D1Display.OSPRayScaleFunction = 'PiecewiseFunction'
delaunay3D1Display.SelectOrientationVectors = 'None'
delaunay3D1Display.ScaleFactor = 0.8
delaunay3D1Display.SelectScaleArray = 'None'
delaunay3D1Display.GlyphType = 'Arrow'
delaunay3D1Display.GlyphTableIndexArray = 'None'
delaunay3D1Display.GaussianRadius = 0.04
delaunay3D1Display.SetScaleArray = ['POINTS', 'f']
delaunay3D1Display.ScaleTransferFunction = 'PiecewiseFunction'
delaunay3D1Display.OpacityArray = ['POINTS', 'f']
delaunay3D1Display.OpacityTransferFunction = 'PiecewiseFunction'
delaunay3D1Display.DataAxesGrid = 'GridAxesRepresentation'
delaunay3D1Display.PolarAxes = 'PolarAxesRepresentation'
delaunay3D1Display.ScalarOpacityUnitDistance = 0.0663092149088729

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
delaunay3D1Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# show data from contour2
contour2Display = Show(contour2, renderView1, 'GeometryRepresentation')

# get color transfer function/color map for 'ux'
uxLUT = GetColorTransferFunction('ux')
uxLUT.RGBPoints = [0.0, 0.0, 0.0, 0.5625, 0.23363270937628142, 0.0, 0.0, 1.0, 0.7676513821560252, 0.0, 1.0, 1.0, 1.0346601928717756, 0.5, 1.0, 0.5, 1.3016690035875256, 1.0, 1.0, 0.0, 1.8356876763672694, 1.0, 0.0, 0.0, 2.1026964870830196, 0.5, 0.0, 0.0]
uxLUT.ColorSpace = 'RGB'
uxLUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
contour2Display.Representation = 'Surface'
contour2Display.ColorArrayName = ['POINTS', 'u.x']
contour2Display.LookupTable = uxLUT
contour2Display.Specular = 1.0
contour2Display.OSPRayScaleArray = 'f'
contour2Display.OSPRayScaleFunction = 'PiecewiseFunction'
contour2Display.SelectOrientationVectors = 'None'
contour2Display.ScaleFactor = 0.3379618165202066
contour2Display.SelectScaleArray = 'f'
contour2Display.GlyphType = 'Arrow'
contour2Display.GlyphTableIndexArray = 'f'
contour2Display.GaussianRadius = 0.01689809082601033
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

# show data from cylinder1
cylinder1Display = Show(cylinder1, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
cylinder1Display.Representation = 'Surface'
cylinder1Display.ColorArrayName = [None, '']
cylinder1Display.Opacity = 0.32
cylinder1Display.Position = [15.0, 0.0, 0.0]
cylinder1Display.Orientation = [0.0, 0.0, 90.0]
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

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
cylinder1Display.PolarAxes.Translation = [15.0, 0.0, 0.0]
cylinder1Display.PolarAxes.Orientation = [0.0, 0.0, 90.0]

# show data from clip2
clip2Display = Show(clip2, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
clip2Display.Representation = 'Surface'
clip2Display.ColorArrayName = [None, '']
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
contour3Display = Show(contour3, renderView1, 'GeometryRepresentation')

# get color transfer function/color map for 'l2'
l2LUT = GetColorTransferFunction('l2')
l2LUT.RGBPoints = [-218.18022266231446, 0.23137254902, 0.298039215686, 0.752941176471, -77.53745296274602, 0.865, 0.865, 0.865, 63.10531673682243, 0.705882352941, 0.0156862745098, 0.149019607843]
l2LUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
contour3Display.Representation = 'Surface'
contour3Display.ColorArrayName = ['POINTS', 'l2']
contour3Display.LookupTable = l2LUT
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

# get color legend/bar for fLUT in view renderView1
fLUTColorBar = GetScalarBar(fLUT, renderView1)
fLUTColorBar.WindowLocation = 'AnyLocation'
fLUTColorBar.Title = 'f'
fLUTColorBar.ComponentTitle = ''
fLUTColorBar.AutomaticLabelFormat = 0
fLUTColorBar.LabelFormat = '%-#6.2g'
fLUTColorBar.RangeLabelFormat = '%6.2g'
fLUTColorBar.ScalarBarLength = 0.3300000000000002

# set color bar visibility
fLUTColorBar.Visibility = 0

# get color legend/bar for uxLUT in view renderView1
uxLUTColorBar = GetScalarBar(uxLUT, renderView1)
uxLUTColorBar.WindowLocation = 'AnyLocation'
uxLUTColorBar.Title = 'u.x'
uxLUTColorBar.ComponentTitle = 'Magnitude'
uxLUTColorBar.AutomaticLabelFormat = 0
uxLUTColorBar.LabelFormat = '%-#6.2g'
uxLUTColorBar.RangeLabelFormat = '%6.2g'
uxLUTColorBar.ScalarBarLength = 0.3300000000000002

# set color bar visibility
uxLUTColorBar.Visibility = 1

# get color legend/bar for l2LUT in view renderView1
l2LUTColorBar = GetScalarBar(l2LUT, renderView1)
l2LUTColorBar.WindowLocation = 'UpperRightCorner'
l2LUTColorBar.Title = 'l2'
l2LUTColorBar.ComponentTitle = ''
l2LUTColorBar.AutomaticLabelFormat = 0
l2LUTColorBar.LabelFormat = '%-#6.2g'
l2LUTColorBar.RangeLabelFormat = '%6.2g'
l2LUTColorBar.ScalarBarLength = 0.3300000000000002

# set color bar visibility
l2LUTColorBar.Visibility = 0

# hide data in view
Hide(contour1, renderView1)

# hide data in view
Hide(clip1, renderView1)

# hide data in view
Hide(cellDatatoPointData1, renderView1)

# hide data in view
Hide(delaunay3D1, renderView1)

# show color legend
contour2Display.SetScalarBarVisibility(renderView1, True)

# hide data in view
Hide(cylinder1, renderView1)

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get opacity transfer function/opacity map for 'l2'
l2PWF = GetOpacityTransferFunction('l2')
l2PWF.Points = [-218.18022266231446, 0.0, 0.5, 0.0, 63.10531673682243, 1.0, 0.5, 0.0]
l2PWF.ScalarRangeInitialized = 1

# get opacity transfer function/opacity map for 'f'
fPWF = GetOpacityTransferFunction('f')
fPWF.Points = [0.5, 0.0, 0.5, 0.0, 0.5001220703125, 1.0, 0.5, 0.0]
fPWF.ScalarRangeInitialized = 1

# get opacity transfer function/opacity map for 'ux'
uxPWF = GetOpacityTransferFunction('ux')
uxPWF.Points = [0.0, 0.0, 0.5, 0.0, 2.1026964870830196, 1.0, 0.5, 0.0]
uxPWF.ScalarRangeInitialized = 1

# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(contour3)
# ----------------------------------------------------------------
bounds = contour1.GetDataInformation().GetBounds()
print("bounds contour1=", bounds)

bounds = contour1.GetDataInformation().GetBounds()
print("bounds contour2=", bounds)
SaveScreenshot( "MY_TEST.png", renderView1,
#      ImageResolution=[2316, 2204],
    TransparentBackground=0,
#     magnification=1,
    # PNG options
    CompressionLevel='2' )
print('File=' + "MY_TEST.png" + ' generated succesfully')

