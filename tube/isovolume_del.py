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
dump2pvd_compressedpvd = PVDReader(FileName='/Users/weugene/basilisk/work/tube/dump2pvd_compressed.pvd')
#dump2pvd_compressedpvd.CellArrays = ['fs', 'f', 'l', 'l2', 'omega', 'u.x']

# create a new 'Contour'
contour1 = Contour(Input=dump2pvd_compressedpvd)
contour1.ContourBy = ['POINTS', 'f']
contour1.Isosurfaces = [0.5]
contour1.PointMergeMethod = 'Uniform Binning'

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

# create a new 'Contour'
contour2 = Contour(Input=resampleToImage1)
contour2.ContourBy = ['POINTS', 'f']
contour2.Isosurfaces = [0.5000000000000012]
contour2.PointMergeMethod = 'Uniform Binning'

# create a new 'Contour'
contour5 = Contour(Input=isoVolume1)
contour5.ContourBy = ['POINTS', 'l2']
contour5.Isosurfaces = [-1.0, -0.5]
contour5.PointMergeMethod = 'Uniform Binning'

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from contour1
contour1Display = Show(contour1, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
contour1Display.Representation = 'Wireframe'
contour1Display.AmbientColor = [0.0, 0.3333333333333333, 1.0]
contour1Display.ColorArrayName = ['POINTS', '']
contour1Display.DiffuseColor = [0.0, 0.3333333333333333, 1.0]
contour1Display.Opacity = 0.11
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

# show data from contour5
contour5Display = Show(contour5, renderView1, 'GeometryRepresentation')

# get color transfer function/color map for 'l2'
l2LUT = GetColorTransferFunction('l2')
l2LUT.RGBPoints = [-1.0, 0.054901960784313725, 0.9411764705882353, 0.12941176470588237, -0.75, 0.865, 0.865, 0.865, -0.5, 1.0, 1.0, 0.0]
l2LUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
contour5Display.Representation = 'Surface'
contour5Display.ColorArrayName = ['POINTS', 'l2']
contour5Display.LookupTable = l2LUT
contour5Display.Opacity = 0.56
contour5Display.OSPRayScaleArray = 'l2'
contour5Display.OSPRayScaleFunction = 'PiecewiseFunction'
contour5Display.SelectOrientationVectors = 'None'
contour5Display.ScaleFactor = 0.3394383430480957
contour5Display.SelectScaleArray = 'l2'
contour5Display.GlyphType = 'Arrow'
contour5Display.GlyphTableIndexArray = 'l2'
contour5Display.GaussianRadius = 0.016971917152404786
contour5Display.SetScaleArray = ['POINTS', 'l2']
contour5Display.ScaleTransferFunction = 'PiecewiseFunction'
contour5Display.OpacityArray = ['POINTS', 'l2']
contour5Display.OpacityTransferFunction = 'PiecewiseFunction'
contour5Display.DataAxesGrid = 'GridAxesRepresentation'
contour5Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
contour5Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
contour5Display.ScaleTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, -0.9998779296875, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
contour5Display.OpacityTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, -0.9998779296875, 1.0, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for l2LUT in view renderView1
l2LUTColorBar = GetScalarBar(l2LUT, renderView1)
l2LUTColorBar.WindowLocation = 'UpperLeftCorner'
l2LUTColorBar.Title = 'l2'
l2LUTColorBar.ComponentTitle = ''
l2LUTColorBar.LabelFormat = '%-#6.2g'
l2LUTColorBar.RangeLabelFormat = '%6.2g'
l2LUTColorBar.ScalarBarLength = 0.3300000000000002

# set color bar visibility
l2LUTColorBar.Visibility = 1

# show color legend
contour5Display.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup the visualization in view 'spreadSheetView1'
# ----------------------------------------------------------------

# show data from isoVolume1
isoVolume1Display = Show(isoVolume1, spreadSheetView1, 'SpreadSheetRepresentation')

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get opacity transfer function/opacity map for 'l2'
l2PWF = GetOpacityTransferFunction('l2')
l2PWF.Points = [-1.0, 0.0, 0.5, 0.0, -0.5, 1.0, 0.5, 0.0]
l2PWF.ScalarRangeInitialized = 1

# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(isoVolume1)
# ----------------------------------------------------------------
