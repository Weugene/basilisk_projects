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
renderView1.ViewSize = [1180, 1124]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [10.57037178924736, 2.868010577215367, 6.785167138329059]
renderView1.CameraFocalPoint = [15.0, 0.0, -0.03538362681865692]
renderView1.CameraViewUp = [0.11977480112232697, 0.9405875943633277, -0.3177246832493743]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 15.015509567667033
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.OSPRayMaterialLibrary = materialLibrary1

# init the 'GridAxes3DActor' selected for 'AxesGrid'
renderView1.AxesGrid.XTitleColor = [0.8199282825970855, 0.9197375448233768, 0.9009994659342336]
renderView1.AxesGrid.XTitleFontSize = 18
renderView1.AxesGrid.XLabelColor = [0.8969863431754025, 0.9574883649958038, 0.8886396581979095]
renderView1.AxesGrid.XLabelFontSize = 18
renderView1.AxesGrid.YLabelFontSize = 14
renderView1.AxesGrid.ZLabelFontSize = 14
renderView1.AxesGrid.CustomBounds = [0.0, 30.0, 0.0, 1.0, 0.0, 1.0]

# Create a new 'SpreadSheet View'
spreadSheetView1 = CreateView('SpreadSheetView')
spreadSheetView1.ColumnToSort = ''
spreadSheetView1.BlockSize = 1024L
# uncomment following to set a specific view size
# spreadSheetView1.ViewSize = [400, 400]

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

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'Cylinder'
cylinder1 = Cylinder()
cylinder1.Resolution = 82
cylinder1.Height = 30.0
cylinder1.Capping = 0

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
clip1.ClipType.Position = [14.0, -0.6, -0.6]
clip1.ClipType.Length = [5.0, 1.2, 1.2]

# init the 'Plane' selected for 'HyperTreeGridClipper'
clip1.HyperTreeGridClipper.Origin = [15.0, 0.0, 0.0]

# create a new 'Cell Data to Point Data'
cellDatatoPointData1 = CellDatatoPointData(Input=clip1)
cellDatatoPointData1.CellDataArraytoprocess = ['f', 'fs', 'l', 'l2', 'omega', 'u.x']

# create a new 'Resample To Image'
resampleToImage1 = ResampleToImage(Input=cellDatatoPointData1)
resampleToImage1.UseInputBounds = 0
resampleToImage1.SamplingBounds = [14.2, 18.8, -0.5, 0.5, -0.5, 0.5]

# create a new 'Iso Volume'
isoVolume1 = IsoVolume(Input=resampleToImage1)
isoVolume1.InputScalars = ['POINTS', 'f']
isoVolume1.ThresholdRange = [0.0, 0.5]

# create a new 'Clip'
clip2 = Clip(Input=cylinder1)
clip2.ClipType = 'Plane'
clip2.HyperTreeGridClipper = 'Plane'
clip2.Scalars = ['POINTS', 'Normals_Magnitude']
clip2.Value = 1.0000000004315743

# init the 'Plane' selected for 'ClipType'
clip2.ClipType.Normal = [0.0, 0.0, 1.0]

# create a new 'Contour'
contour2 = Contour(Input=resampleToImage1)
contour2.ContourBy = ['POINTS', 'f']
contour2.Isosurfaces = [0.5000000000000012]
contour2.PointMergeMethod = 'Uniform Binning'

# create a new 'Contour'
contour1 = Contour(Input=dump2pvd_compressedpvd)
contour1.ContourBy = ['POINTS', 'f']
contour1.Isosurfaces = [0.5]
contour1.PointMergeMethod = 'Uniform Binning'

# create a new 'Contour'
contour5 = Contour(Input=isoVolume1)
contour5.ContourBy = ['POINTS', 'l2']
contour5.Isosurfaces = [-1.0, -0.5]
contour5.PointMergeMethod = 'Uniform Binning'

# create a new 'Transform'
transform1 = Transform(Input=clip2)
transform1.Transform = 'Transform'

# init the 'Transform' selected for 'Transform'
transform1.Transform.Translate = [15.0, 0.0, 0.0]
transform1.Transform.Rotate = [0.0, 0.0, 90.0]

# create a new 'Extract Selection'
extractSelection1 = ExtractSelection(Input=contour2)

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

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
transform1Display.ScalarOpacityUnitDistance = 8.706185806543795

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
transform1Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
transform1Display.ScaleTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
transform1Display.OpacityTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# show data from extractSelection1
extractSelection1Display = Show(extractSelection1, renderView1, 'UnstructuredGridRepresentation')

# get color transfer function/color map for 'f'
fLUT = GetColorTransferFunction('f')
fLUT.RGBPoints = [0.0, 0.23137254902, 0.298039215686, 0.752941176471, 0.5000000000000002, 0.865, 0.865, 0.865, 1.0000000000000004, 0.705882352941, 0.0156862745098, 0.149019607843]
fLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'f'
fPWF = GetOpacityTransferFunction('f')
fPWF.Points = [0.0, 0.0, 0.5, 0.0, 1.0000000000000004, 1.0, 0.5, 0.0]
fPWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
extractSelection1Display.Representation = 'Surface'
extractSelection1Display.ColorArrayName = ['POINTS', 'f']
extractSelection1Display.LookupTable = fLUT
extractSelection1Display.Opacity = 0.85
extractSelection1Display.OSPRayScaleArray = 'f'
extractSelection1Display.OSPRayScaleFunction = 'PiecewiseFunction'
extractSelection1Display.SelectOrientationVectors = 'None'
extractSelection1Display.ScaleFactor = 0.1
extractSelection1Display.SelectScaleArray = 'f'
extractSelection1Display.GlyphType = 'Arrow'
extractSelection1Display.GlyphTableIndexArray = 'f'
extractSelection1Display.GaussianRadius = 0.005
extractSelection1Display.SetScaleArray = ['POINTS', 'f']
extractSelection1Display.ScaleTransferFunction = 'PiecewiseFunction'
extractSelection1Display.OpacityArray = ['POINTS', 'f']
extractSelection1Display.OpacityTransferFunction = 'PiecewiseFunction'
extractSelection1Display.DataAxesGrid = 'GridAxesRepresentation'
extractSelection1Display.PolarAxes = 'PolarAxesRepresentation'
extractSelection1Display.ScalarOpacityFunction = fPWF
extractSelection1Display.ScalarOpacityUnitDistance = 0.0

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
extractSelection1Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
extractSelection1Display.ScaleTransferFunction.Points = [0.5, 0.0, 0.5, 0.0, 0.5001220703125, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
extractSelection1Display.OpacityTransferFunction.Points = [0.5, 0.0, 0.5, 0.0, 0.5001220703125, 1.0, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for fLUT in view renderView1
fLUTColorBar = GetScalarBar(fLUT, renderView1)
fLUTColorBar.WindowLocation = 'UpperRightCorner'
fLUTColorBar.Title = 'f'
fLUTColorBar.ComponentTitle = ''
fLUTColorBar.LabelFormat = '%-#6.2g'
fLUTColorBar.RangeLabelFormat = '%6.2g'
fLUTColorBar.ScalarBarLength = 0.3300000000000002

# set color bar visibility
fLUTColorBar.Visibility = 1

# show color legend
extractSelection1Display.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup the visualization in view 'spreadSheetView1'
# ----------------------------------------------------------------

# show data from contour5
contour5Display = Show(contour5, spreadSheetView1, 'SpreadSheetRepresentation')

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(extractSelection1)
# ----------------------------------------------------------------