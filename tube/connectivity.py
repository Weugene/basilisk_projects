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
renderView1.ViewSize = [1096, 1354]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.CenterOfRotation = [16.132638154732103, 0.0012683107921633574, 0.0012655406649076106]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [15.758235520512642, -6.12518755156023, 9.910041743991151]
renderView1.CameraFocalPoint = [16.132638154732103, 0.001268310792163371, 0.0012655406649075302]
renderView1.CameraViewUp = [-0.01275768042689818, 0.8507013271322498, 0.5254945229072842]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 2.0778206447408616
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

# Create a new 'SpreadSheet View'
spreadSheetView2 = CreateView('SpreadSheetView')
spreadSheetView2.ColumnToSort = ''
spreadSheetView2.BlockSize = 1024L
# uncomment following to set a specific view size
# spreadSheetView2.ViewSize = [400, 400]

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.SplitHorizontal(0, 0.521576)
layout1.AssignView(1, renderView1)
layout1.AssignView(2, spreadSheetView2)

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

# create a new 'PVD Reader'
dump2pvd_compressedpvd = PVDReader(FileName='/Users/weugene/basilisk/work/tube/dump2pvd_compressed.pvd')
dump2pvd_compressedpvd.CellArrays = ['fs', 'f', 'l', 'l2', 'omega', 'u.x']

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

# create a new 'Connectivity'
connectivity1 = Connectivity(Input=contour2)
connectivity1.ExtractionMode = 'Extract Largest Region'
connectivity1.ColorRegions = 0
connectivity1.RegionIdAssignmentMode = 'Cell Count Ascending'

# create a new 'Cylinder'
cylinder1 = Cylinder()
cylinder1.Resolution = 82
cylinder1.Height = 30.0
cylinder1.Capping = 0

# create a new 'Extract Selection'
extractSelection1 = ExtractSelection(Input=contour2)

# create a new 'Clip'
clip2 = Clip(Input=cylinder1)
clip2.ClipType = 'Plane'
clip2.HyperTreeGridClipper = 'Plane'
clip2.Scalars = ['POINTS', 'Normals_Magnitude']
clip2.Value = 1.0000000004315743

# init the 'Plane' selected for 'ClipType'
clip2.ClipType.Normal = [0.0, 0.0, 1.0]

# create a new 'Calculator'
calculator1 = Calculator(Input=contour2)
calculator1.ResultArrayName = 'ResultR'
calculator1.Function = 'sqrt(coordsY^2+coordsZ^2)'

# create a new 'Contour'
contour5 = Contour(Input=isoVolume1)
contour5.ContourBy = ['POINTS', 'l2']
contour5.Isosurfaces = [-1.0, -0.5]
contour5.PointMergeMethod = 'Uniform Binning'

# create a new 'Iso Volume'
isoVolume2 = IsoVolume(Input=clip1)
isoVolume2.InputScalars = ['POINTS', 'f']
isoVolume2.ThresholdRange = [0.0, 0.5]

# create a new 'Cell Data to Point Data'
cellDatatoPointData2 = CellDatatoPointData(Input=isoVolume2)
cellDatatoPointData2.CellDataArraytoprocess = ['f', 'f', 'fs', 'l', 'l2', 'omega', 'u.x']

# create a new 'Calculator'
calculator2 = Calculator(Input=cellDatatoPointData2)
calculator2.Function = 'coordsX*(1-f)'

# create a new 'Integrate Variables'
integrateVariables2 = IntegrateVariables(Input=cellDatatoPointData2)

# create a new 'Programmable Filter'
programmableFilter1 = ProgrammableFilter(Input=contour2)
programmableFilter1.Script = """inp = self.GetInput()
numPoints = inp.GetNumberOfPoints()
print numPoints
data=[]
for x in range(numPoints):
    data.append(inp.GetPointData().GetArray('u.x').GetValue(x))
print data """
programmableFilter1.RequestInformationScript = ''
programmableFilter1.RequestUpdateExtentScript = ''
programmableFilter1.PythonPath = ''

# create a new 'Clip'
clip3 = Clip(Input=calculator1)
clip3.ClipType = 'Box'
clip3.HyperTreeGridClipper = 'Plane'
clip3.Scalars = ['POINTS', 'ResultX']
clip3.Value = 0.2257609853976255

# init the 'Box' selected for 'ClipType'
clip3.ClipType.Position = [14.8, -0.4396306872367859, -0.43842387199401855]
clip3.ClipType.Length = [1.5, 0.8837202489376068, 0.8770357668399811]

# init the 'Plane' selected for 'HyperTreeGridClipper'
clip3.HyperTreeGridClipper.Origin = [16.132742404937744, 0.002229437232017517, 9.401142597198486e-05]

# create a new 'Integrate Variables'
integrateVariables1 = IntegrateVariables(Input=clip3)

# create a new 'Transform'
transform1 = Transform(Input=clip2)
transform1.Transform = 'Transform'

# init the 'Transform' selected for 'Transform'
transform1.Transform.Translate = [15.0, 0.0, 0.0]
transform1.Transform.Rotate = [0.0, 0.0, 90.0]

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from connectivity1
connectivity1Display = Show(connectivity1, renderView1, 'UnstructuredGridRepresentation')

# get color transfer function/color map for 'vtkGhostType'
vtkGhostTypeLUT = GetColorTransferFunction('vtkGhostType')
vtkGhostTypeLUT.RGBPoints = [0.0, 0.23137254902, 0.298039215686, 0.752941176471, 5.878906683738906e-39, 0.865, 0.865, 0.865, 1.1757813367477812e-38, 0.705882352941, 0.0156862745098, 0.149019607843]
vtkGhostTypeLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'vtkGhostType'
vtkGhostTypePWF = GetOpacityTransferFunction('vtkGhostType')
vtkGhostTypePWF.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]
vtkGhostTypePWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
connectivity1Display.Representation = 'Surface'
connectivity1Display.ColorArrayName = ['CELLS', 'vtkGhostType']
connectivity1Display.LookupTable = vtkGhostTypeLUT
connectivity1Display.Opacity = 0.73
connectivity1Display.OSPRayScaleArray = 'RegionId'
connectivity1Display.OSPRayScaleFunction = 'PiecewiseFunction'
connectivity1Display.SelectOrientationVectors = 'None'
connectivity1Display.ScaleFactor = 0.3379618165202066
connectivity1Display.SelectScaleArray = 'RegionId'
connectivity1Display.GlyphType = 'Arrow'
connectivity1Display.GlyphTableIndexArray = 'RegionId'
connectivity1Display.GaussianRadius = 0.01689809082601033
connectivity1Display.SetScaleArray = ['POINTS', 'RegionId']
connectivity1Display.ScaleTransferFunction = 'PiecewiseFunction'
connectivity1Display.OpacityArray = ['POINTS', 'RegionId']
connectivity1Display.OpacityTransferFunction = 'PiecewiseFunction'
connectivity1Display.DataAxesGrid = 'GridAxesRepresentation'
connectivity1Display.PolarAxes = 'PolarAxesRepresentation'
connectivity1Display.ScalarOpacityFunction = vtkGhostTypePWF
connectivity1Display.ScalarOpacityUnitDistance = 0.08645112135295362

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
connectivity1Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
connectivity1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 9.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
connectivity1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 9.0, 1.0, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for vtkGhostTypeLUT in view renderView1
vtkGhostTypeLUTColorBar = GetScalarBar(vtkGhostTypeLUT, renderView1)
vtkGhostTypeLUTColorBar.WindowLocation = 'UpperLeftCorner'
vtkGhostTypeLUTColorBar.Title = 'vtkGhostType'
vtkGhostTypeLUTColorBar.ComponentTitle = ''
vtkGhostTypeLUTColorBar.LabelFormat = '%-#6.2g'
vtkGhostTypeLUTColorBar.RangeLabelFormat = '%6.2g'
vtkGhostTypeLUTColorBar.ScalarBarLength = 0.3300000000000002

# set color bar visibility
vtkGhostTypeLUTColorBar.Visibility = 1

# show color legend
connectivity1Display.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup the visualization in view 'spreadSheetView1'
# ----------------------------------------------------------------

# show data from contour5
contour5Display = Show(contour5, spreadSheetView1, 'SpreadSheetRepresentation')

# ----------------------------------------------------------------
# setup the visualization in view 'spreadSheetView2'
# ----------------------------------------------------------------

# show data from integrateVariables1
integrateVariables1Display = Show(integrateVariables1, spreadSheetView2, 'SpreadSheetRepresentation')

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(connectivity1)
# ----------------------------------------------------------------