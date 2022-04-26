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
renderView1.ViewSize = [1644, 1298]
renderView1.InteractionMode = '2D'
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [-0.8185821117462657, 0.08929986673595625, 27.320508075688775]
renderView1.CameraFocalPoint = [-0.8185821117462657, 0.08929986673595625, 0.0]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 4.829634459302967
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.OSPRayMaterialLibrary = materialLibrary1

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.AssignView(0, renderView1)
layout1.SetSize(1644, 1298)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'PVD Reader'
output_respvd = PVDReader(registrationName='output_res.pvd', FileName='/Users/weugene/basilisk/work/Heat_transfer_Polymerization/output_res.pvd')

# create a new 'HyperTreeGrid To Dual Grid'
hyperTreeGridToDualGrid1 = HyperTreeGridToDualGrid(registrationName='HyperTreeGridToDualGrid1', Input=output_respvd)

# create a new 'Contour'
contour2 = Contour(registrationName='Contour2', Input=hyperTreeGridToDualGrid1)
contour2.ContourBy = ['POINTS', 'f']
contour2.Isosurfaces = [0.5]
contour2.PointMergeMethod = 'Uniform Binning'

# create a new 'Contour'
contour1 = Contour(registrationName='Contour1', Input=hyperTreeGridToDualGrid1)
contour1.ContourBy = ['POINTS', 'fs']
contour1.Isosurfaces = [0.5]
contour1.PointMergeMethod = 'Uniform Binning'

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from output_respvd
output_respvdDisplay = Show(output_respvd, renderView1, 'GeometryRepresentation')

# get color transfer function/color map for 'u'
uLUT = GetColorTransferFunction('u')
uLUT.RGBPoints = [0.0, 0.231373, 0.298039, 0.752941, 1.5, 0.865003, 0.865003, 0.865003, 3.0, 0.705882, 0.0156863, 0.14902]
uLUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
output_respvdDisplay.Representation = 'Surface'
output_respvdDisplay.ColorArrayName = ['CELLS', 'u']
output_respvdDisplay.LookupTable = uLUT
output_respvdDisplay.SelectTCoordArray = 'None'
output_respvdDisplay.SelectNormalArray = 'None'
output_respvdDisplay.SelectTangentArray = 'None'
output_respvdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
output_respvdDisplay.SelectOrientationVectors = 'None'
output_respvdDisplay.SelectScaleArray = 'fs_face'
output_respvdDisplay.GlyphType = 'Arrow'
output_respvdDisplay.GlyphTableIndexArray = 'fs_face'
output_respvdDisplay.GaussianRadius = 0.05
output_respvdDisplay.SetScaleArray = [None, '']
output_respvdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
output_respvdDisplay.OpacityArray = [None, '']
output_respvdDisplay.OpacityTransferFunction = 'PiecewiseFunction'
output_respvdDisplay.DataAxesGrid = 'GridAxesRepresentation'
output_respvdDisplay.PolarAxes = 'PolarAxesRepresentation'

# show data from hyperTreeGridToDualGrid1
hyperTreeGridToDualGrid1Display = Show(hyperTreeGridToDualGrid1, renderView1, 'UnstructuredGridRepresentation')

# get opacity transfer function/opacity map for 'u'
uPWF = GetOpacityTransferFunction('u')
uPWF.Points = [0.0, 0.0, 0.5, 0.0, 3.0, 1.0, 0.5, 0.0]
uPWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
hyperTreeGridToDualGrid1Display.Representation = 'Surface'
hyperTreeGridToDualGrid1Display.ColorArrayName = ['POINTS', 'u']
hyperTreeGridToDualGrid1Display.LookupTable = uLUT
hyperTreeGridToDualGrid1Display.SelectTCoordArray = 'None'
hyperTreeGridToDualGrid1Display.SelectNormalArray = 'None'
hyperTreeGridToDualGrid1Display.SelectTangentArray = 'None'
hyperTreeGridToDualGrid1Display.OSPRayScaleArray = 'fs_face'
hyperTreeGridToDualGrid1Display.OSPRayScaleFunction = 'PiecewiseFunction'
hyperTreeGridToDualGrid1Display.SelectOrientationVectors = 'None'
hyperTreeGridToDualGrid1Display.SelectScaleArray = 'fs_face'
hyperTreeGridToDualGrid1Display.GlyphType = 'Arrow'
hyperTreeGridToDualGrid1Display.GlyphTableIndexArray = 'fs_face'
hyperTreeGridToDualGrid1Display.GaussianRadius = 0.05
hyperTreeGridToDualGrid1Display.SetScaleArray = ['POINTS', 'fs_face']
hyperTreeGridToDualGrid1Display.ScaleTransferFunction = 'PiecewiseFunction'
hyperTreeGridToDualGrid1Display.OpacityArray = ['POINTS', 'fs_face']
hyperTreeGridToDualGrid1Display.OpacityTransferFunction = 'PiecewiseFunction'
hyperTreeGridToDualGrid1Display.DataAxesGrid = 'GridAxesRepresentation'
hyperTreeGridToDualGrid1Display.PolarAxes = 'PolarAxesRepresentation'
hyperTreeGridToDualGrid1Display.ScalarOpacityFunction = uPWF
hyperTreeGridToDualGrid1Display.ScalarOpacityUnitDistance = 0.6912841351258522
hyperTreeGridToDualGrid1Display.OpacityArrayName = ['POINTS', 'fs_face']

# show data from contour2
contour2Display = Show(contour2, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
contour2Display.Representation = 'Surface'
contour2Display.AmbientColor = [0.0, 0.0, 1.0]
contour2Display.ColorArrayName = ['POINTS', '']
contour2Display.DiffuseColor = [0.0, 0.0, 1.0]
contour2Display.LineWidth = 5.0
contour2Display.SelectTCoordArray = 'None'
contour2Display.SelectNormalArray = 'None'
contour2Display.SelectTangentArray = 'None'
contour2Display.OSPRayScaleArray = 'f'
contour2Display.OSPRayScaleFunction = 'PiecewiseFunction'
contour2Display.SelectOrientationVectors = 'None'
contour2Display.ScaleFactor = 0.18000437021255494
contour2Display.SelectScaleArray = 'f'
contour2Display.GlyphType = 'Arrow'
contour2Display.GlyphTableIndexArray = 'f'
contour2Display.GaussianRadius = 0.009000218510627747
contour2Display.SetScaleArray = ['POINTS', 'f']
contour2Display.ScaleTransferFunction = 'PiecewiseFunction'
contour2Display.OpacityArray = ['POINTS', 'f']
contour2Display.OpacityTransferFunction = 'PiecewiseFunction'
contour2Display.DataAxesGrid = 'GridAxesRepresentation'
contour2Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
contour2Display.ScaleTransferFunction.Points = [0.5, 0.0, 0.5, 0.0, 0.5001220703125, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
contour2Display.OpacityTransferFunction.Points = [0.5, 0.0, 0.5, 0.0, 0.5001220703125, 1.0, 0.5, 0.0]

# show data from contour1
contour1Display = Show(contour1, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
contour1Display.Representation = 'Surface'
contour1Display.AmbientColor = [0.3333333333333333, 0.0, 0.0]
contour1Display.ColorArrayName = ['POINTS', '']
contour1Display.DiffuseColor = [0.3333333333333333, 0.0, 0.0]
contour1Display.LineWidth = 5.0
contour1Display.SelectTCoordArray = 'None'
contour1Display.SelectNormalArray = 'None'
contour1Display.SelectTangentArray = 'None'
contour1Display.OSPRayScaleArray = 'fs'
contour1Display.OSPRayScaleFunction = 'PiecewiseFunction'
contour1Display.SelectOrientationVectors = 'None'
contour1Display.ScaleFactor = 0.6996097326278687
contour1Display.SelectScaleArray = 'fs'
contour1Display.GlyphType = 'Arrow'
contour1Display.GlyphTableIndexArray = 'fs'
contour1Display.GaussianRadius = 0.034980486631393436
contour1Display.SetScaleArray = ['POINTS', 'fs']
contour1Display.ScaleTransferFunction = 'PiecewiseFunction'
contour1Display.OpacityArray = ['POINTS', 'fs']
contour1Display.OpacityTransferFunction = 'PiecewiseFunction'
contour1Display.DataAxesGrid = 'GridAxesRepresentation'
contour1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
contour1Display.ScaleTransferFunction.Points = [0.5, 0.0, 0.5, 0.0, 0.5001220703125, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
contour1Display.OpacityTransferFunction.Points = [0.5, 0.0, 0.5, 0.0, 0.5001220703125, 1.0, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for uLUT in view renderView1
uLUTColorBar = GetScalarBar(uLUT, renderView1)
uLUTColorBar.WindowLocation = 'UpperLeftCorner'
uLUTColorBar.Position = [0.0069051881724923395, 0.3364396537793629]
uLUTColorBar.Title = '$|\\mathbf{u}|$'
uLUTColorBar.ComponentTitle = ''
uLUTColorBar.TitleFontFamily = 'Times'
uLUTColorBar.LabelFontFamily = 'Times'
uLUTColorBar.ScalarBarLength = 0.32999999999999996

# set color bar visibility
uLUTColorBar.Visibility = 1

# show color legend
output_respvdDisplay.SetScalarBarVisibility(renderView1, True)

# show color legend
hyperTreeGridToDualGrid1Display.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# restore active source
SetActiveSource(output_respvd)
# ----------------------------------------------------------------


if __name__ == '__main__':
    # generate extracts
    SaveExtracts(ExtractsOutputDirectory='extracts')