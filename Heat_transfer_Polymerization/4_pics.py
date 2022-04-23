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
renderView1.ViewSize = [1874, 1298]
renderView1.InteractionMode = '2D'
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.CenterOfRotation = [-0.4666639715433121, -2.2644856572151184, 0.0]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [5.2887722743389265, -5.102964182358942, 4.922480693762108]
renderView1.CameraFocalPoint = [5.2887722743389265, -5.102964182358942, 0.0]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 10.253861216098938
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.OSPRayMaterialLibrary = materialLibrary1

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.AssignView(0, renderView1)
layout1.SetSize(1874, 1298)

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
hyperTreeGridToDualGrid2 = HyperTreeGridToDualGrid(registrationName='HyperTreeGridToDualGrid2', Input=output_respvd)

# create a new 'HyperTreeGrid To Dual Grid'
hyperTreeGridToDualGrid3 = HyperTreeGridToDualGrid(registrationName='HyperTreeGridToDualGrid3', Input=output_respvd)

# create a new 'Contour'
contour5 = Contour(registrationName='Contour5', Input=hyperTreeGridToDualGrid3)
contour5.ContourBy = ['POINTS', 'f']
contour5.Isosurfaces = [0.5]
contour5.PointMergeMethod = 'Uniform Binning'

# create a new 'Contour'
contour4 = Contour(registrationName='Contour4', Input=hyperTreeGridToDualGrid2)
contour4.ContourBy = ['POINTS', 'fs']
contour4.Isosurfaces = [0.5]
contour4.PointMergeMethod = 'Uniform Binning'

# create a new 'HyperTreeGrid To Dual Grid'
hyperTreeGridToDualGrid4 = HyperTreeGridToDualGrid(registrationName='HyperTreeGridToDualGrid4', Input=output_respvd)

# create a new 'Contour'
contour7 = Contour(registrationName='Contour7', Input=hyperTreeGridToDualGrid4)
contour7.ContourBy = ['POINTS', 'f']
contour7.Isosurfaces = [0.5]
contour7.PointMergeMethod = 'Uniform Binning'

# create a new 'Contour'
contour6 = Contour(registrationName='Contour6', Input=hyperTreeGridToDualGrid3)
contour6.ContourBy = ['POINTS', 'fs']
contour6.Isosurfaces = [0.5]
contour6.PointMergeMethod = 'Uniform Binning'

# create a new 'Contour'
contour8 = Contour(registrationName='Contour8', Input=hyperTreeGridToDualGrid4)
contour8.ContourBy = ['POINTS', 'fs']
contour8.Isosurfaces = [0.5]
contour8.PointMergeMethod = 'Uniform Binning'

# create a new 'Contour'
contour3 = Contour(registrationName='Contour3', Input=hyperTreeGridToDualGrid2)
contour3.ContourBy = ['POINTS', 'f']
contour3.Isosurfaces = [0.5]
contour3.PointMergeMethod = 'Uniform Binning'

# create a new 'HyperTreeGrid To Dual Grid'
hyperTreeGridToDualGrid1 = HyperTreeGridToDualGrid(registrationName='HyperTreeGridToDualGrid1', Input=output_respvd)

# create a new 'Contour'
contour1 = Contour(registrationName='Contour1', Input=hyperTreeGridToDualGrid1)
contour1.ContourBy = ['POINTS', 'fs']
contour1.Isosurfaces = [0.5]
contour1.PointMergeMethod = 'Uniform Binning'

# create a new 'Contour'
contour2 = Contour(registrationName='Contour2', Input=hyperTreeGridToDualGrid1)
contour2.ContourBy = ['POINTS', 'f']
contour2.Isosurfaces = [0.5]
contour2.PointMergeMethod = 'Uniform Binning'

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from hyperTreeGridToDualGrid1
hyperTreeGridToDualGrid1Display = Show(hyperTreeGridToDualGrid1, renderView1, 'UnstructuredGridRepresentation')

# get color transfer function/color map for 'u'
uLUT = GetColorTransferFunction('u')
uLUT.RGBPoints = [0.0, 0.231373, 0.298039, 0.752941, 1.5, 0.865003, 0.865003, 0.865003, 3.0, 0.705882, 0.0156863, 0.14902]
uLUT.ScalarRangeInitialized = 1.0

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

# show data from hyperTreeGridToDualGrid2
hyperTreeGridToDualGrid2Display = Show(hyperTreeGridToDualGrid2, renderView1, 'UnstructuredGridRepresentation')

# get color transfer function/color map for 'alpha_doc'
alpha_docLUT = GetColorTransferFunction('alpha_doc')
alpha_docLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'alpha_doc'
alpha_docPWF = GetOpacityTransferFunction('alpha_doc')
alpha_docPWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
hyperTreeGridToDualGrid2Display.Representation = 'Surface'
hyperTreeGridToDualGrid2Display.ColorArrayName = ['POINTS', 'alpha_doc']
hyperTreeGridToDualGrid2Display.LookupTable = alpha_docLUT
hyperTreeGridToDualGrid2Display.SelectTCoordArray = 'None'
hyperTreeGridToDualGrid2Display.SelectNormalArray = 'None'
hyperTreeGridToDualGrid2Display.SelectTangentArray = 'None'
hyperTreeGridToDualGrid2Display.Position = [0.0, -10.2, 0.0]
hyperTreeGridToDualGrid2Display.OSPRayScaleArray = 'fs_face'
hyperTreeGridToDualGrid2Display.OSPRayScaleFunction = 'PiecewiseFunction'
hyperTreeGridToDualGrid2Display.SelectOrientationVectors = 'None'
hyperTreeGridToDualGrid2Display.SelectScaleArray = 'fs_face'
hyperTreeGridToDualGrid2Display.GlyphType = 'Arrow'
hyperTreeGridToDualGrid2Display.GlyphTableIndexArray = 'fs_face'
hyperTreeGridToDualGrid2Display.GaussianRadius = 0.05
hyperTreeGridToDualGrid2Display.SetScaleArray = ['POINTS', 'fs_face']
hyperTreeGridToDualGrid2Display.ScaleTransferFunction = 'PiecewiseFunction'
hyperTreeGridToDualGrid2Display.OpacityArray = ['POINTS', 'fs_face']
hyperTreeGridToDualGrid2Display.OpacityTransferFunction = 'PiecewiseFunction'
hyperTreeGridToDualGrid2Display.DataAxesGrid = 'GridAxesRepresentation'
hyperTreeGridToDualGrid2Display.PolarAxes = 'PolarAxesRepresentation'
hyperTreeGridToDualGrid2Display.ScalarOpacityFunction = alpha_docPWF
hyperTreeGridToDualGrid2Display.ScalarOpacityUnitDistance = 0.6912841351258522
hyperTreeGridToDualGrid2Display.OpacityArrayName = ['POINTS', 'fs_face']

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
hyperTreeGridToDualGrid2Display.PolarAxes.Translation = [0.0, -10.2, 0.0]

# show data from hyperTreeGridToDualGrid3
hyperTreeGridToDualGrid3Display = Show(hyperTreeGridToDualGrid3, renderView1, 'UnstructuredGridRepresentation')

# get color transfer function/color map for 'mu_cell'
mu_cellLUT = GetColorTransferFunction('mu_cell')
mu_cellLUT.RGBPoints = [0.0, 0.231373, 0.298039, 0.752941, 9.570980072021484, 0.865003, 0.865003, 0.865003, 19.14196014404297, 0.705882, 0.0156863, 0.14902]
mu_cellLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'mu_cell'
mu_cellPWF = GetOpacityTransferFunction('mu_cell')
mu_cellPWF.Points = [0.0, 0.0, 0.5, 0.0, 19.14196014404297, 1.0, 0.5, 0.0]
mu_cellPWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
hyperTreeGridToDualGrid3Display.Representation = 'Surface'
hyperTreeGridToDualGrid3Display.ColorArrayName = ['POINTS', 'mu_cell']
hyperTreeGridToDualGrid3Display.LookupTable = mu_cellLUT
hyperTreeGridToDualGrid3Display.SelectTCoordArray = 'None'
hyperTreeGridToDualGrid3Display.SelectNormalArray = 'None'
hyperTreeGridToDualGrid3Display.SelectTangentArray = 'None'
hyperTreeGridToDualGrid3Display.Position = [10.2, -10.2, 0.0]
hyperTreeGridToDualGrid3Display.OSPRayScaleArray = 'fs_face'
hyperTreeGridToDualGrid3Display.OSPRayScaleFunction = 'PiecewiseFunction'
hyperTreeGridToDualGrid3Display.SelectOrientationVectors = 'None'
hyperTreeGridToDualGrid3Display.SelectScaleArray = 'fs_face'
hyperTreeGridToDualGrid3Display.GlyphType = 'Arrow'
hyperTreeGridToDualGrid3Display.GlyphTableIndexArray = 'fs_face'
hyperTreeGridToDualGrid3Display.GaussianRadius = 0.05
hyperTreeGridToDualGrid3Display.SetScaleArray = ['POINTS', 'fs_face']
hyperTreeGridToDualGrid3Display.ScaleTransferFunction = 'PiecewiseFunction'
hyperTreeGridToDualGrid3Display.OpacityArray = ['POINTS', 'fs_face']
hyperTreeGridToDualGrid3Display.OpacityTransferFunction = 'PiecewiseFunction'
hyperTreeGridToDualGrid3Display.DataAxesGrid = 'GridAxesRepresentation'
hyperTreeGridToDualGrid3Display.PolarAxes = 'PolarAxesRepresentation'
hyperTreeGridToDualGrid3Display.ScalarOpacityFunction = mu_cellPWF
hyperTreeGridToDualGrid3Display.ScalarOpacityUnitDistance = 0.6912841351258522
hyperTreeGridToDualGrid3Display.OpacityArrayName = ['POINTS', 'fs_face']

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
hyperTreeGridToDualGrid3Display.PolarAxes.Translation = [10.2, -10.2, 0.0]

# show data from contour3
contour3Display = Show(contour3, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
contour3Display.Representation = 'Surface'
contour3Display.AmbientColor = [0.0, 0.0, 1.0]
contour3Display.ColorArrayName = ['POINTS', '']
contour3Display.DiffuseColor = [0.0, 0.0, 1.0]
contour3Display.LineWidth = 4.0
contour3Display.SelectTCoordArray = 'None'
contour3Display.SelectNormalArray = 'None'
contour3Display.SelectTangentArray = 'None'
contour3Display.Position = [0.0, -10.2, 0.0]
contour3Display.OSPRayScaleArray = 'f'
contour3Display.OSPRayScaleFunction = 'PiecewiseFunction'
contour3Display.SelectOrientationVectors = 'None'
contour3Display.ScaleFactor = 0.18000437021255494
contour3Display.SelectScaleArray = 'f'
contour3Display.GlyphType = 'Arrow'
contour3Display.GlyphTableIndexArray = 'f'
contour3Display.GaussianRadius = 0.009000218510627747
contour3Display.SetScaleArray = ['POINTS', 'f']
contour3Display.ScaleTransferFunction = 'PiecewiseFunction'
contour3Display.OpacityArray = ['POINTS', 'f']
contour3Display.OpacityTransferFunction = 'PiecewiseFunction'
contour3Display.DataAxesGrid = 'GridAxesRepresentation'
contour3Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
contour3Display.ScaleTransferFunction.Points = [0.5, 0.0, 0.5, 0.0, 0.5001220703125, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
contour3Display.OpacityTransferFunction.Points = [0.5, 0.0, 0.5, 0.0, 0.5001220703125, 1.0, 0.5, 0.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
contour3Display.PolarAxes.Translation = [0.0, -10.2, 0.0]

# show data from contour4
contour4Display = Show(contour4, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
contour4Display.Representation = 'Surface'
contour4Display.AmbientColor = [0.3333333333333333, 0.0, 0.0]
contour4Display.ColorArrayName = ['POINTS', '']
contour4Display.DiffuseColor = [0.3333333333333333, 0.0, 0.0]
contour4Display.LineWidth = 4.0
contour4Display.SelectTCoordArray = 'None'
contour4Display.SelectNormalArray = 'None'
contour4Display.SelectTangentArray = 'None'
contour4Display.Position = [0.0, -10.2, 0.0]
contour4Display.OSPRayScaleArray = 'fs'
contour4Display.OSPRayScaleFunction = 'PiecewiseFunction'
contour4Display.SelectOrientationVectors = 'None'
contour4Display.ScaleFactor = 0.6996097326278687
contour4Display.SelectScaleArray = 'fs'
contour4Display.GlyphType = 'Arrow'
contour4Display.GlyphTableIndexArray = 'fs'
contour4Display.GaussianRadius = 0.034980486631393436
contour4Display.SetScaleArray = ['POINTS', 'fs']
contour4Display.ScaleTransferFunction = 'PiecewiseFunction'
contour4Display.OpacityArray = ['POINTS', 'fs']
contour4Display.OpacityTransferFunction = 'PiecewiseFunction'
contour4Display.DataAxesGrid = 'GridAxesRepresentation'
contour4Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
contour4Display.ScaleTransferFunction.Points = [0.5, 0.0, 0.5, 0.0, 0.5001220703125, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
contour4Display.OpacityTransferFunction.Points = [0.5, 0.0, 0.5, 0.0, 0.5001220703125, 1.0, 0.5, 0.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
contour4Display.PolarAxes.Translation = [0.0, -10.2, 0.0]

# show data from contour5
contour5Display = Show(contour5, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
contour5Display.Representation = 'Surface'
contour5Display.AmbientColor = [0.0, 0.0, 1.0]
contour5Display.ColorArrayName = ['POINTS', '']
contour5Display.DiffuseColor = [0.0, 0.0, 1.0]
contour5Display.LineWidth = 4.0
contour5Display.SelectTCoordArray = 'None'
contour5Display.SelectNormalArray = 'None'
contour5Display.SelectTangentArray = 'None'
contour5Display.Position = [10.2, -10.2, 0.0]
contour5Display.OSPRayScaleArray = 'f'
contour5Display.OSPRayScaleFunction = 'PiecewiseFunction'
contour5Display.SelectOrientationVectors = 'None'
contour5Display.ScaleFactor = 0.18000437021255494
contour5Display.SelectScaleArray = 'f'
contour5Display.GlyphType = 'Arrow'
contour5Display.GlyphTableIndexArray = 'f'
contour5Display.GaussianRadius = 0.009000218510627747
contour5Display.SetScaleArray = ['POINTS', 'f']
contour5Display.ScaleTransferFunction = 'PiecewiseFunction'
contour5Display.OpacityArray = ['POINTS', 'f']
contour5Display.OpacityTransferFunction = 'PiecewiseFunction'
contour5Display.DataAxesGrid = 'GridAxesRepresentation'
contour5Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
contour5Display.ScaleTransferFunction.Points = [0.5, 0.0, 0.5, 0.0, 0.5001220703125, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
contour5Display.OpacityTransferFunction.Points = [0.5, 0.0, 0.5, 0.0, 0.5001220703125, 1.0, 0.5, 0.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
contour5Display.PolarAxes.Translation = [10.2, -10.2, 0.0]

# show data from contour6
contour6Display = Show(contour6, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
contour6Display.Representation = 'Surface'
contour6Display.AmbientColor = [0.3333333333333333, 0.0, 0.0]
contour6Display.ColorArrayName = ['POINTS', '']
contour6Display.DiffuseColor = [0.3333333333333333, 0.0, 0.0]
contour6Display.LineWidth = 4.0
contour6Display.SelectTCoordArray = 'None'
contour6Display.SelectNormalArray = 'None'
contour6Display.SelectTangentArray = 'None'
contour6Display.Position = [10.2, -10.2, 0.0]
contour6Display.OSPRayScaleArray = 'fs'
contour6Display.OSPRayScaleFunction = 'PiecewiseFunction'
contour6Display.SelectOrientationVectors = 'None'
contour6Display.ScaleFactor = 0.6996097326278687
contour6Display.SelectScaleArray = 'fs'
contour6Display.GlyphType = 'Arrow'
contour6Display.GlyphTableIndexArray = 'fs'
contour6Display.GaussianRadius = 0.034980486631393436
contour6Display.SetScaleArray = ['POINTS', 'fs']
contour6Display.ScaleTransferFunction = 'PiecewiseFunction'
contour6Display.OpacityArray = ['POINTS', 'fs']
contour6Display.OpacityTransferFunction = 'PiecewiseFunction'
contour6Display.DataAxesGrid = 'GridAxesRepresentation'
contour6Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
contour6Display.ScaleTransferFunction.Points = [0.5, 0.0, 0.5, 0.0, 0.5001220703125, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
contour6Display.OpacityTransferFunction.Points = [0.5, 0.0, 0.5, 0.0, 0.5001220703125, 1.0, 0.5, 0.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
contour6Display.PolarAxes.Translation = [10.2, -10.2, 0.0]

# show data from hyperTreeGridToDualGrid4
hyperTreeGridToDualGrid4Display = Show(hyperTreeGridToDualGrid4, renderView1, 'UnstructuredGridRepresentation')

# get color transfer function/color map for 'T'
tLUT = GetColorTransferFunction('T')
tLUT.RGBPoints = [1.0, 0.231373, 0.298039, 0.752941, 2.0, 0.865003, 0.865003, 0.865003, 3.0, 0.705882, 0.0156863, 0.14902]
tLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'T'
tPWF = GetOpacityTransferFunction('T')
tPWF.Points = [1.0, 0.0, 0.5, 0.0, 3.0, 1.0, 0.5, 0.0]
tPWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
hyperTreeGridToDualGrid4Display.Representation = 'Surface'
hyperTreeGridToDualGrid4Display.ColorArrayName = ['POINTS', 'T']
hyperTreeGridToDualGrid4Display.LookupTable = tLUT
hyperTreeGridToDualGrid4Display.SelectTCoordArray = 'None'
hyperTreeGridToDualGrid4Display.SelectNormalArray = 'None'
hyperTreeGridToDualGrid4Display.SelectTangentArray = 'None'
hyperTreeGridToDualGrid4Display.Position = [10.2, 0.0, 0.0]
hyperTreeGridToDualGrid4Display.OSPRayScaleArray = 'fs_face'
hyperTreeGridToDualGrid4Display.OSPRayScaleFunction = 'PiecewiseFunction'
hyperTreeGridToDualGrid4Display.SelectOrientationVectors = 'None'
hyperTreeGridToDualGrid4Display.SelectScaleArray = 'fs_face'
hyperTreeGridToDualGrid4Display.GlyphType = 'Arrow'
hyperTreeGridToDualGrid4Display.GlyphTableIndexArray = 'fs_face'
hyperTreeGridToDualGrid4Display.GaussianRadius = 0.05
hyperTreeGridToDualGrid4Display.SetScaleArray = ['POINTS', 'fs_face']
hyperTreeGridToDualGrid4Display.ScaleTransferFunction = 'PiecewiseFunction'
hyperTreeGridToDualGrid4Display.OpacityArray = ['POINTS', 'fs_face']
hyperTreeGridToDualGrid4Display.OpacityTransferFunction = 'PiecewiseFunction'
hyperTreeGridToDualGrid4Display.DataAxesGrid = 'GridAxesRepresentation'
hyperTreeGridToDualGrid4Display.PolarAxes = 'PolarAxesRepresentation'
hyperTreeGridToDualGrid4Display.ScalarOpacityFunction = tPWF
hyperTreeGridToDualGrid4Display.ScalarOpacityUnitDistance = 0.6912841351258522
hyperTreeGridToDualGrid4Display.OpacityArrayName = ['POINTS', 'fs_face']

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
hyperTreeGridToDualGrid4Display.PolarAxes.Translation = [10.2, 0.0, 0.0]

# show data from contour7
contour7Display = Show(contour7, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
contour7Display.Representation = 'Surface'
contour7Display.AmbientColor = [0.0, 0.0, 1.0]
contour7Display.ColorArrayName = ['POINTS', '']
contour7Display.DiffuseColor = [0.0, 0.0, 1.0]
contour7Display.LineWidth = 4.0
contour7Display.SelectTCoordArray = 'None'
contour7Display.SelectNormalArray = 'None'
contour7Display.SelectTangentArray = 'None'
contour7Display.Position = [10.2, 0.0, 0.0]
contour7Display.OSPRayScaleArray = 'f'
contour7Display.OSPRayScaleFunction = 'PiecewiseFunction'
contour7Display.SelectOrientationVectors = 'None'
contour7Display.ScaleFactor = 0.18000437021255494
contour7Display.SelectScaleArray = 'f'
contour7Display.GlyphType = 'Arrow'
contour7Display.GlyphTableIndexArray = 'f'
contour7Display.GaussianRadius = 0.009000218510627747
contour7Display.SetScaleArray = ['POINTS', 'f']
contour7Display.ScaleTransferFunction = 'PiecewiseFunction'
contour7Display.OpacityArray = ['POINTS', 'f']
contour7Display.OpacityTransferFunction = 'PiecewiseFunction'
contour7Display.DataAxesGrid = 'GridAxesRepresentation'
contour7Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
contour7Display.ScaleTransferFunction.Points = [0.5, 0.0, 0.5, 0.0, 0.5001220703125, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
contour7Display.OpacityTransferFunction.Points = [0.5, 0.0, 0.5, 0.0, 0.5001220703125, 1.0, 0.5, 0.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
contour7Display.PolarAxes.Translation = [10.2, 0.0, 0.0]

# show data from contour8
contour8Display = Show(contour8, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
contour8Display.Representation = 'Surface'
contour8Display.AmbientColor = [0.3333333333333333, 0.0, 0.0]
contour8Display.ColorArrayName = ['POINTS', '']
contour8Display.DiffuseColor = [0.3333333333333333, 0.0, 0.0]
contour8Display.LineWidth = 4.0
contour8Display.SelectTCoordArray = 'None'
contour8Display.SelectNormalArray = 'None'
contour8Display.SelectTangentArray = 'None'
contour8Display.Position = [10.2, 0.0, 0.0]
contour8Display.OSPRayScaleArray = 'fs'
contour8Display.OSPRayScaleFunction = 'PiecewiseFunction'
contour8Display.SelectOrientationVectors = 'None'
contour8Display.ScaleFactor = 0.6996097326278687
contour8Display.SelectScaleArray = 'fs'
contour8Display.GlyphType = 'Arrow'
contour8Display.GlyphTableIndexArray = 'fs'
contour8Display.GaussianRadius = 0.034980486631393436
contour8Display.SetScaleArray = ['POINTS', 'fs']
contour8Display.ScaleTransferFunction = 'PiecewiseFunction'
contour8Display.OpacityArray = ['POINTS', 'fs']
contour8Display.OpacityTransferFunction = 'PiecewiseFunction'
contour8Display.DataAxesGrid = 'GridAxesRepresentation'
contour8Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
contour8Display.ScaleTransferFunction.Points = [0.5, 0.0, 0.5, 0.0, 0.5001220703125, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
contour8Display.OpacityTransferFunction.Points = [0.5, 0.0, 0.5, 0.0, 0.5001220703125, 1.0, 0.5, 0.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
contour8Display.PolarAxes.Translation = [10.2, 0.0, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for uLUT in view renderView1
uLUTColorBar = GetScalarBar(uLUT, renderView1)
uLUTColorBar.WindowLocation = 'AnyLocation'
uLUTColorBar.Position = [0.017613917013485223, 0.5845136137177296]
uLUTColorBar.Title = '$|\\mathbf{u}|$'
uLUTColorBar.ComponentTitle = ''
uLUTColorBar.TitleFontFamily = 'Times'
uLUTColorBar.LabelFontFamily = 'Times'
uLUTColorBar.ScalarBarLength = 0.3300000000000002

# set color bar visibility
uLUTColorBar.Visibility = 1

# get color legend/bar for tLUT in view renderView1
tLUTColorBar = GetScalarBar(tLUT, renderView1)
tLUTColorBar.WindowLocation = 'AnyLocation'
tLUTColorBar.Position = [0.8419285204897381, 0.5954644459263543]
tLUTColorBar.Title = '$T$'
tLUTColorBar.ComponentTitle = ''
tLUTColorBar.TitleFontFamily = 'Times'
tLUTColorBar.LabelFontFamily = 'Times'
tLUTColorBar.ScalarBarLength = 0.32999999999999996

# set color bar visibility
tLUTColorBar.Visibility = 1

# get color legend/bar for alpha_docLUT in view renderView1
alpha_docLUTColorBar = GetScalarBar(alpha_docLUT, renderView1)
alpha_docLUTColorBar.WindowLocation = 'AnyLocation'
alpha_docLUTColorBar.Position = [0.026857406973930825, 0.08434100411267836]
alpha_docLUTColorBar.Title = '$\\alpha$'
alpha_docLUTColorBar.ComponentTitle = ''
alpha_docLUTColorBar.TitleFontFamily = 'Times'
alpha_docLUTColorBar.LabelFontFamily = 'Times'
alpha_docLUTColorBar.ScalarBarLength = 0.330000000000001

# set color bar visibility
alpha_docLUTColorBar.Visibility = 1

# get color legend/bar for mu_cellLUT in view renderView1
mu_cellLUTColorBar = GetScalarBar(mu_cellLUT, renderView1)
mu_cellLUTColorBar.WindowLocation = 'AnyLocation'
mu_cellLUTColorBar.Position = [0.8418046626106415, 0.09206414595343822]
mu_cellLUTColorBar.Title = '$\\mu(T,\\alpha)$'
mu_cellLUTColorBar.ComponentTitle = ''
mu_cellLUTColorBar.TitleFontFamily = 'Times'
mu_cellLUTColorBar.LabelFontFamily = 'Times'

# set color bar visibility
mu_cellLUTColorBar.Visibility = 1

# show color legend
hyperTreeGridToDualGrid1Display.SetScalarBarVisibility(renderView1, True)

# show color legend
hyperTreeGridToDualGrid2Display.SetScalarBarVisibility(renderView1, True)

# show color legend
hyperTreeGridToDualGrid3Display.SetScalarBarVisibility(renderView1, True)

# show color legend
hyperTreeGridToDualGrid4Display.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# restore active source
SetActiveSource(hyperTreeGridToDualGrid1)
# ----------------------------------------------------------------


if __name__ == '__main__':
    # generate extracts
    SaveExtracts(ExtractsOutputDirectory='extracts')