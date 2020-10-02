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
renderView1.ViewSize = [2048, 2048]
renderView1.InteractionMode = '2D'
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.CenterOfRotation = [15.000000155913408, -3.585828789454126e-08, -0.1626075438383584]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [15.000000155913408, -3.585828789454126e-08, 58.6641342543975]
renderView1.CameraFocalPoint = [15.000000155913408, -3.585828789454126e-08, -0.1626075438383584]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 37.55027056222264
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.OSPRayMaterialLibrary = materialLibrary1

# init the 'GridAxes3DActor' selected for 'AxesGrid'
renderView1.AxesGrid.Visibility = 1
renderView1.AxesGrid.XTitle = 'X'
renderView1.AxesGrid.YTitle = 'Y'
renderView1.AxesGrid.ZTitle = 'Z'
renderView1.AxesGrid.XTitleColor = [0.9, 0.9, 0.9]
renderView1.AxesGrid.XTitleFontSize = 20
renderView1.AxesGrid.YTitleColor = [0.9, 0.9, 0.9]
renderView1.AxesGrid.YTitleFontSize = 20
renderView1.AxesGrid.ZTitleColor = [0.9, 0.9, 0.9]
renderView1.AxesGrid.ZTitleFontSize = 20
renderView1.AxesGrid.XLabelColor = [0.9, 0.9, 0.9]
renderView1.AxesGrid.XLabelFontSize = 18
renderView1.AxesGrid.YLabelColor = [0.9, 0.9, 0.9]
renderView1.AxesGrid.YLabelFontSize = 18
renderView1.AxesGrid.ZLabelColor = [0.9, 0.9, 0.9]
renderView1.AxesGrid.ZLabelFontSize = 18
renderView1.AxesGrid.XAxisUseCustomLabels = 1
renderView1.AxesGrid.XAxisLabels = [0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0, 5.25, 5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 7.0, 7.25, 7.5, 7.75, 8.0, 8.25, 8.5, 8.75, 9.0, 9.25, 9.5, 9.75, 10.0, 10.25, 10.5, 10.75, 11.0, 11.25, 11.5, 11.75, 12.0, 12.25, 12.5, 12.75, 13.0, 13.25, 13.5, 13.75, 14.0, 14.25, 14.5, 14.75, 15.0, 15.25, 15.5, 15.75, 16.0, 16.25, 16.5, 16.75, 17.0, 17.25, 17.5, 17.75, 18.0, 18.25, 18.5, 18.75, 19.0, 19.25, 19.5, 19.75, 20.0, 20.25, 20.5, 20.75, 21.0, 21.25, 21.5, 21.75, 22.0, 22.25, 22.5, 22.75, 23.0, 23.25, 23.5, 23.75, 24.0, 24.25, 24.5, 24.75, 25.0, 25.25, 25.5, 25.75, 26.0, 26.25, 26.5, 26.75, 27.0, 27.25, 27.5, 27.75, 28.0, 28.25, 28.5, 28.75, 29.0, 29.25, 29.5, 29.75, 30.0]
renderView1.AxesGrid.YAxisUseCustomLabels = 1
renderView1.AxesGrid.YAxisLabels = [-0.5, -0.25, 0.25, 0.5]
renderView1.AxesGrid.ZAxisUseCustomLabels = 1
renderView1.AxesGrid.ZAxisLabels = [-0.5, -0.25, 0.25, 0.5]

# Create a new 'Render View'
renderView2 = CreateView('RenderView')
renderView2.ViewSize = [2046, 2048]
renderView2.InteractionMode = '2D'
renderView2.AxesGrid = 'GridAxes3DActor'
renderView2.OrientationAxesVisibility = 0
renderView2.CenterOfRotation = [15.0, 0.0, -0.25]
renderView2.StereoType = 'Crystal Eyes'
renderView2.CameraPosition = [15.0, 0.0, 57.7986784864182]
renderView2.CameraFocalPoint = [15.0, 0.0, -0.25]
renderView2.CameraFocalDisk = 1.0
renderView2.CameraParallelScale = 15.025085987813048
renderView2.BackEnd = 'OSPRay raycaster'
renderView2.OSPRayMaterialLibrary = materialLibrary1

# init the 'GridAxes3DActor' selected for 'AxesGrid'
renderView1.AxesGrid.Visibility = 1
renderView1.AxesGrid.XTitle = 'X'
renderView1.AxesGrid.YTitle = 'Y'
renderView1.AxesGrid.ZTitle = 'Z'
renderView1.AxesGrid.XTitleColor = [0.9, 0.9, 0.9]
renderView1.AxesGrid.XTitleFontSize = 20
renderView1.AxesGrid.YTitleColor = [0.9, 0.9, 0.9]
renderView1.AxesGrid.YTitleFontSize = 20
renderView1.AxesGrid.ZTitleColor = [0.9, 0.9, 0.9]
renderView1.AxesGrid.ZTitleFontSize = 20
renderView1.AxesGrid.XLabelColor = [0.9, 0.9, 0.9]
renderView1.AxesGrid.XLabelFontSize = 18
renderView1.AxesGrid.YLabelColor = [0.9, 0.9, 0.9]
renderView1.AxesGrid.YLabelFontSize = 18
renderView1.AxesGrid.ZLabelColor = [0.9, 0.9, 0.9]
renderView1.AxesGrid.ZLabelFontSize = 18
renderView1.AxesGrid.XAxisUseCustomLabels = 1
renderView1.AxesGrid.XAxisLabels = [0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0, 5.25, 5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 7.0, 7.25, 7.5, 7.75, 8.0, 8.25, 8.5, 8.75, 9.0, 9.25, 9.5, 9.75, 10.0, 10.25, 10.5, 10.75, 11.0, 11.25, 11.5, 11.75, 12.0, 12.25, 12.5, 12.75, 13.0, 13.25, 13.5, 13.75, 14.0, 14.25, 14.5, 14.75, 15.0, 15.25, 15.5, 15.75, 16.0, 16.25, 16.5, 16.75, 17.0, 17.25, 17.5, 17.75, 18.0, 18.25, 18.5, 18.75, 19.0, 19.25, 19.5, 19.75, 20.0, 20.25, 20.5, 20.75, 21.0, 21.25, 21.5, 21.75, 22.0, 22.25, 22.5, 22.75, 23.0, 23.25, 23.5, 23.75, 24.0, 24.25, 24.5, 24.75, 25.0, 25.25, 25.5, 25.75, 26.0, 26.25, 26.5, 26.75, 27.0, 27.25, 27.5, 27.75, 28.0, 28.25, 28.5, 28.75, 29.0, 29.25, 29.5, 29.75, 30.0]
renderView1.AxesGrid.YAxisUseCustomLabels = 1
renderView1.AxesGrid.YAxisLabels = [-0.5, -0.25, 0.25, 0.5]
renderView1.AxesGrid.ZAxisUseCustomLabels = 1
renderView1.AxesGrid.ZAxisLabels = [-0.5, -0.25, 0.25, 0.5]

# Create a new 'Render View'
renderView3 = CreateView('RenderView')
renderView3.ViewSize = [2120, 1124]
renderView3.InteractionMode = '2D'
renderView3.AxesGrid = 'GridAxes3DActor'
renderView3.CenterOfRotation = [14.834758758544922, -0.0007286816835403442, -0.0006128847599029541]
renderView3.StereoType = 'Crystal Eyes'
renderView3.CameraPosition = [17.075357284218406, -0.0007286816835403442, -0.0006128847599029541]
renderView3.CameraFocalPoint = [14.834758758544922, -0.0007286816835403442, -0.0006128847599029541]
renderView3.CameraViewUp = [0.0, 0.0, 1.0]
renderView3.CameraFocalDisk = 1.0
renderView3.CameraParallelScale = 0.5799095708729264
renderView3.BackEnd = 'OSPRay raycaster'
renderView3.OSPRayMaterialLibrary = materialLibrary1

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
layout1.SplitHorizontal(0, 0.500000)
layout1.AssignView(1, renderView1)
layout1.AssignView(2, renderView2)

# create new layout object 'Layout #2'
layout2 = CreateLayout(name='Layout #2')
layout2.AssignView(0, spreadSheetView1)

# create new layout object 'Layout #3'
layout3 = CreateLayout(name='Layout #3')
layout3.AssignView(0, spreadSheetView2)

# create new layout object 'Layout #4'
layout4 = CreateLayout(name='Layout #4')
layout4.AssignView(0, renderView3)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView3)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'Cylinder'
cylinder1 = Cylinder()
cylinder1.Resolution = 100
cylinder1.Height = 30.0
cylinder1.Capping = 0

# create a new 'PVD Reader'
dump2pvd_compressedpvd = PVDReader(FileName='/Users/weugene/basilisk/work/tube/dump2pvd_compressed.pvd')
dump2pvd_compressedpvd.CellArrays = ['fs', 'f', 'l', 'l2', 'omega', 'u.x']

# create a new 'Clip'
clip1 = Clip(Input=cylinder1)
clip1.ClipType = 'Plane'
clip1.HyperTreeGridClipper = 'Plane'
clip1.Scalars = ['POINTS', '']
clip1.Value = 1.0

# init the 'Plane' selected for 'ClipType'
clip1.ClipType.Normal = [0.0, 0.0, 1.0]

# create a new 'Transform'
transform1 = Transform(Input=clip1)
transform1.Transform = 'Transform'

# init the 'Transform' selected for 'Transform'
transform1.Transform.Translate = [15.005360160722741, -0.00080364434197433, -5.876124017079043e-05]
transform1.Transform.Rotate = [1.240091495378975, 0.07381443459238507, 90.0]
transform1.Transform.Scale = [0.9999999999999994, 1.0000000000000033, 0.9999999999999953]

# create a new 'Clip'
clip2 = Clip(Input=dump2pvd_compressedpvd)
clip2.ClipType = 'Box'
clip2.HyperTreeGridClipper = 'Plane'
clip2.Scalars = ['POINTS', '']
clip2.Value = 0.5

# init the 'Box' selected for 'ClipType'
clip2.ClipType.Position = [0.0, -0.6, -0.6]
clip2.ClipType.Length = [30.0, 1.2, 1.2]

# init the 'Plane' selected for 'HyperTreeGridClipper'
clip2.HyperTreeGridClipper.Origin = [15.0, 0.0, 0.0]

# create a new 'Clip'
clip3 = Clip(Input=clip2)
clip3.ClipType = 'Box'
clip3.HyperTreeGridClipper = 'Plane'
clip3.Scalars = ['POINTS', 'f']
clip3.Value = 0.5

# init the 'Box' selected for 'ClipType'
clip3.ClipType.Position = [10.909568023681642, -0.6, -0.6]
clip3.ClipType.Length = [7.599999999999999, 1.2, 1.2]

# init the 'Plane' selected for 'HyperTreeGridClipper'
clip3.HyperTreeGridClipper.Origin = [15.0, 0.0, 0.0]

# create a new 'Cell Data to Point Data'
cellDatatoPointData1 = CellDatatoPointData(Input=clip3)

# create a new 'Calculator'
calculator1 = Calculator(Input=cellDatatoPointData1)
calculator1.ResultArrayName = 'absOmega'
calculator1.Function = 'abs(omega)'

# create a new 'Resample To Image'
resampleToImage1 = ResampleToImage(Input=calculator1)
resampleToImage1.UseInputBounds = 0
resampleToImage1.SamplingDimensions = [1010, 136, 136]
resampleToImage1.SamplingBounds = [11.009568023681641, 18.40956802368164, -0.5, 0.5, -0.5, 0.5]

# create a new 'Iso Volume'
isoVolume1 = IsoVolume(Input=resampleToImage1)
isoVolume1.InputScalars = ['POINTS', 'f']
isoVolume1.ThresholdRange = [0.0, 0.5]

# create a new 'Connectivity'
connectivity1 = Connectivity(Input=isoVolume1)
connectivity1.ExtractionMode = 'Extract Largest Region'
connectivity1.ColorRegions = 0

# create a new 'Contour'
contour2 = Contour(Input=connectivity1)
contour2.ContourBy = ['POINTS', 'l2']
contour2.Isosurfaces = [-1.5]
contour2.PointMergeMethod = 'Uniform Binning'

# create a new 'Integrate Variables'
integrateVariables1 = IntegrateVariables(Input=connectivity1)

# create a new 'PassArrays'
passArrays1 = PassArrays(Input=integrateVariables1)
passArrays1.PointDataArrays = ['u.x']
passArrays1.CellDataArrays = ['absOmega', 'f', 'fs', 'l', 'l2', 'omega', 'u.x', 'vtkGhostType', 'vtkValidPointMask', 'Volume', 'vtkGhostType']

# create a new 'Extract Surface'
f = ExtractSurface(Input=connectivity1)

# create a new 'Slice'
slice6 = Slice(Input=extractSurface1)
slice6.SliceType = 'Plane'
slice6.HyperTreeGridSlicer = 'Plane'
slice6.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice6.SliceType.Origin = [14.83475923538208, -0.00173245370388031, -0.00017940998077392578]

# init the 'Plane' selected for 'HyperTreeGridSlicer'
slice6.HyperTreeGridSlicer.Origin = [14.83475923538208, -0.00173245370388031, -0.00017940998077392578]

# create a new 'PassArrays'
passArrays2 = PassArrays(Input=connectivity1)

# create a new 'Extract Selection'
extractSelection1 = ExtractSelection(Input=connectivity1)

# create a new 'Slice'
slice2 = Slice(Input=connectivity1)
slice2.SliceType = 'Plane'
slice2.HyperTreeGridSlicer = 'Plane'
slice2.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice2.SliceType.Origin = [14.83475923538208, 0.0, 0.0]
slice2.SliceType.Normal = [0.0, 0.0, 1.0]

# init the 'Plane' selected for 'HyperTreeGridSlicer'
slice2.HyperTreeGridSlicer.Origin = [14.83475923538208, -0.00173245370388031, -0.00017940998077392578]

# create a new 'Contour'
contour1 = Contour(Input=resampleToImage1)
contour1.ContourBy = ['POINTS', 'l2']
contour1.Isosurfaces = [-2.0, -1.0]
contour1.PointMergeMethod = 'Uniform Binning'

# create a new 'Calculator'
calculator2 = Calculator(Input=resampleToImage1)
calculator2.ResultArrayName = 'deltaU'
calculator2.Function = 'u.x - iHat*1.34014457095'

# create a new 'Stream Tracer'
streamTracer1 = StreamTracer(Input=calculator2,
    SeedType='High Resolution Line Source')
streamTracer1.Vectors = ['POINTS', 'deltaU']
streamTracer1.SurfaceStreamlines = 1
streamTracer1.InitialStepLength = 0.05
streamTracer1.MaximumSteps = 20000
streamTracer1.MaximumStreamlineLength = 7.386458396911621
streamTracer1.TerminalSpeed = 0.0001

# init the 'High Resolution Line Source' selected for 'SeedType'
streamTracer1.SeedType.Point1 = [17.52798843383789, -0.5, 0.0]
streamTracer1.SeedType.Point2 = [17.52798843383789, 0.5, 0.0]
streamTracer1.SeedType.Resolution = 50

# create a new 'Stream Tracer'
streamTracer2 = StreamTracer(Input=calculator2,
    SeedType='High Resolution Line Source')
streamTracer2.Vectors = ['POINTS', 'deltaU']
streamTracer2.SurfaceStreamlines = 1
streamTracer2.InitialStepLength = 0.01
streamTracer2.MaximumSteps = 20000
streamTracer2.MaximumStreamlineLength = 7.386458396911621
streamTracer2.TerminalSpeed = 0.0001

# init the 'High Resolution Line Source' selected for 'SeedType'
streamTracer2.SeedType.Point1 = [12.14153003692627, -0.5, 0.0]
streamTracer2.SeedType.Point2 = [12.14153003692627, 0.5, 0.0]
streamTracer2.SeedType.Resolution = 50

# create a new 'Slice'
slice3 = Slice(Input=resampleToImage1)
slice3.SliceType = 'Plane'
slice3.HyperTreeGridSlicer = 'Plane'
slice3.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice3.SliceType.Origin = [14.83475923538208, 0.0, 0.0]
slice3.SliceType.Normal = [0.0, 0.0, 1.0]

# init the 'Plane' selected for 'HyperTreeGridSlicer'
slice3.HyperTreeGridSlicer.Origin = [14.70956802368164, 0.0, 0.0]

# create a new 'Stream Tracer'
streamTracer3 = StreamTracer(Input=calculator2,
    SeedType='High Resolution Line Source')
streamTracer3.Vectors = ['POINTS', 'deltaU']
streamTracer3.SurfaceStreamlines = 1
streamTracer3.InitialStepLength = 0.01
streamTracer3.MaximumSteps = 20000
streamTracer3.MaximumStreamlineLength = 7.386458396911621
streamTracer3.TerminalSpeed = 0.0001

# init the 'High Resolution Line Source' selected for 'SeedType'
streamTracer3.SeedType.Point1 = [14.83475923538208, -0.4, 0.0]
streamTracer3.SeedType.Point2 = [14.83475923538208, 0.4, 0.0]
streamTracer3.SeedType.Resolution = 50

# create a new 'Slice'
slice1 = Slice(Input=transform1)
slice1.SliceType = 'Plane'
slice1.HyperTreeGridSlicer = 'Plane'
slice1.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Origin = [15.0, 0.0, -1e-05]
slice1.SliceType.Normal = [0.0, 0.0, 1.0]

# init the 'Plane' selected for 'HyperTreeGridSlicer'
slice1.HyperTreeGridSlicer.Origin = [15.0, 0.0, -0.21650634706020355]

# create a new 'Slice'
slice5 = Slice(Input=connectivity1)
slice5.SliceType = 'Plane'
slice5.HyperTreeGridSlicer = 'Plane'
slice5.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice5.SliceType.Origin = [16.42798843383789, 0.0, 0.0]

# init the 'Plane' selected for 'HyperTreeGridSlicer'
slice5.HyperTreeGridSlicer.Origin = [14.83475923538208, -0.00173245370388031, -0.00017940998077392578]

# create a new 'Slice'
slice4 = Slice(Input=resampleToImage1)
slice4.SliceType = 'Plane'
slice4.HyperTreeGridSlicer = 'Plane'
slice4.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice4.SliceType.Origin = [16.42798843383789, 0.0, 0.0]

# init the 'Plane' selected for 'HyperTreeGridSlicer'
slice4.HyperTreeGridSlicer.Origin = [14.70956802368164, 0.0, 0.0]

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from transform1
transform1Display = Show(transform1, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
transform1Display.Representation = 'Surface'
transform1Display.AmbientColor = [1.0, 0.7843137254901961, 0.7529411764705882]
transform1Display.ColorArrayName = ['POINTS', '']
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

# show data from slice4
slice4Display = Show(slice4, renderView1, 'GeometryRepresentation')

# get separate color transfer function/color map for 'absOmega'
separate_slice4Display_absOmegaLUT = GetColorTransferFunction('absOmega', slice4Display, separate=True)
separate_slice4Display_absOmegaLUT.RGBPoints = [0.0, 0.0, 0.0, 0.5625, 1.198753696686236, 0.0, 0.0, 1.0, 3.938767540651804, 0.0, 1.0, 1.0, 5.3087717654360755, 0.5, 1.0, 0.5, 6.678775990220343, 1.0, 1.0, 0.0, 9.41878983418591, 1.0, 0.0, 0.0, 10.78879405897018, 0.5, 0.0, 0.0]
separate_slice4Display_absOmegaLUT.ColorSpace = 'RGB'
separate_slice4Display_absOmegaLUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
slice4Display.Representation = 'Surface'
slice4Display.ColorArrayName = ['POINTS', 'absOmega']
slice4Display.LookupTable = separate_slice4Display_absOmegaLUT
slice4Display.OSPRayScaleArray = 'f'
slice4Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice4Display.SelectOrientationVectors = 'None'
slice4Display.ScaleFactor = 0.3379646740552927
slice4Display.SelectScaleArray = 'f'
slice4Display.GlyphType = 'Arrow'
slice4Display.GlyphTableIndexArray = 'f'
slice4Display.GaussianRadius = 0.016898233702764633
slice4Display.SetScaleArray = ['POINTS', 'f']
slice4Display.ScaleTransferFunction = 'PiecewiseFunction'
slice4Display.OpacityArray = ['POINTS', 'f']
slice4Display.OpacityTransferFunction = 'PiecewiseFunction'
slice4Display.DataAxesGrid = 'GridAxesRepresentation'
slice4Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
slice4Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
slice4Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0000000000000002, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
slice4Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0000000000000002, 1.0, 0.5, 0.0]

# show data from slice5
slice5Display = Show(slice5, renderView1, 'GeometryRepresentation')

# get color transfer function/color map for 'absOmega'
absOmegaLUT = GetColorTransferFunction('absOmega')
absOmegaLUT.RGBPoints = [0.0, 0.0, 0.0, 0.5625, 3.324178890387746, 0.0, 0.0, 1.0, 10.922317027236845, 0.0, 1.0, 1.0, 14.721378616251418, 0.5, 1.0, 0.5, 18.520440205265977, 1.0, 1.0, 0.0, 26.11857834211507, 1.0, 0.0, 0.0, 29.91763993112964, 0.5, 0.0, 0.0]
absOmegaLUT.ColorSpace = 'RGB'
absOmegaLUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
slice5Display.Representation = 'Surface'
slice5Display.AmbientColor = [0.0, 0.0, 0.0]
slice5Display.ColorArrayName = ['POINTS', '']
slice5Display.DiffuseColor = [0.0, 0.0, 0.0]
slice5Display.LookupTable = absOmegaLUT
slice5Display.LineWidth = 4.0
slice5Display.OSPRayScaleArray = 'f'
slice5Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice5Display.SelectOrientationVectors = 'None'
slice5Display.ScaleFactor = 0.3379646740552927
slice5Display.SelectScaleArray = 'f'
slice5Display.GlyphType = 'Arrow'
slice5Display.GlyphTableIndexArray = 'f'
slice5Display.GaussianRadius = 0.003989625424146652
slice5Display.SetScaleArray = ['POINTS', 'f']
slice5Display.ScaleTransferFunction = 'PiecewiseFunction'
slice5Display.OpacityArray = ['POINTS', 'f']
slice5Display.OpacityTransferFunction = 'PiecewiseFunction'
slice5Display.DataAxesGrid = 'GridAxesRepresentation'
slice5Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
slice5Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
slice5Display.ScaleTransferFunction.Points = [0.06722302855573097, 0.0, 0.5, 0.0, 5.7994661811456565, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
slice5Display.OpacityTransferFunction.Points = [0.06722302855573097, 0.0, 0.5, 0.0, 5.7994661811456565, 1.0, 0.5, 0.0]

# ----------------------------------------------------------------
# setup the visualization in view 'renderView2'
# ----------------------------------------------------------------

# show data from transform1
transform1Display_1 = Show(transform1, renderView2, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
transform1Display_1.Representation = 'Surface'
transform1Display_1.AmbientColor = [1.0, 0.7843137254901961, 0.7529411764705882]
transform1Display_1.ColorArrayName = ['POINTS', '']
transform1Display_1.DiffuseColor = [1.0, 0.7843137254901961, 0.7529411764705882]
transform1Display_1.Opacity = 0.85
transform1Display_1.Specular = 1.0
transform1Display_1.Luminosity = 35.0
transform1Display_1.OSPRayUseScaleArray = 1
transform1Display_1.OSPRayScaleArray = 'Normals'
transform1Display_1.OSPRayScaleFunction = 'PiecewiseFunction'
transform1Display_1.OSPRayMaterial = 'copper'
transform1Display_1.SelectOrientationVectors = 'None'
transform1Display_1.ScaleFactor = 3.0
transform1Display_1.SelectScaleArray = 'None'
transform1Display_1.GlyphType = 'Arrow'
transform1Display_1.GlyphTableIndexArray = 'None'
transform1Display_1.GaussianRadius = 0.15
transform1Display_1.SetScaleArray = ['POINTS', 'Normals']
transform1Display_1.ScaleTransferFunction = 'PiecewiseFunction'
transform1Display_1.OpacityArray = ['POINTS', 'Normals']
transform1Display_1.OpacityTransferFunction = 'PiecewiseFunction'
transform1Display_1.DataAxesGrid = 'GridAxesRepresentation'
transform1Display_1.PolarAxes = 'PolarAxesRepresentation'
transform1Display_1.ScalarOpacityUnitDistance = 6.650076732513133

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
transform1Display_1.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
transform1Display_1.ScaleTransferFunction.Points = [-2.220446049250313e-16, 0.0, 0.5, 0.0, 2.220446049250313e-16, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
transform1Display_1.OpacityTransferFunction.Points = [-2.220446049250313e-16, 0.0, 0.5, 0.0, 2.220446049250313e-16, 1.0, 0.5, 0.0]

# show data from dump2pvd_compressedpvd
dump2pvd_compressedpvdDisplay = Show(dump2pvd_compressedpvd, renderView2, 'UnstructuredGridRepresentation')

# get color transfer function/color map for 'f'
fLUT = GetColorTransferFunction('f')
fLUT.RGBPoints = [0.0, 0.23137254902, 0.298039215686, 0.752941176471, 0.5, 0.865, 0.865, 0.865, 1.0, 0.705882352941, 0.0156862745098, 0.149019607843]
fLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'f'
fPWF = GetOpacityTransferFunction('f')
fPWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
dump2pvd_compressedpvdDisplay.Representation = 'Surface'
dump2pvd_compressedpvdDisplay.ColorArrayName = ['CELLS', 'f']
dump2pvd_compressedpvdDisplay.LookupTable = fLUT
dump2pvd_compressedpvdDisplay.Opacity = 0.85
dump2pvd_compressedpvdDisplay.OSPRayScaleArray = 'f'
dump2pvd_compressedpvdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
dump2pvd_compressedpvdDisplay.SelectOrientationVectors = 'None'
dump2pvd_compressedpvdDisplay.ScaleFactor = 3.0
dump2pvd_compressedpvdDisplay.SelectScaleArray = 'None'
dump2pvd_compressedpvdDisplay.GlyphType = 'Arrow'
dump2pvd_compressedpvdDisplay.GlyphTableIndexArray = 'None'
dump2pvd_compressedpvdDisplay.GaussianRadius = 0.15
dump2pvd_compressedpvdDisplay.SetScaleArray = ['POINTS', 'f']
dump2pvd_compressedpvdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
dump2pvd_compressedpvdDisplay.OpacityArray = ['POINTS', 'f']
dump2pvd_compressedpvdDisplay.OpacityTransferFunction = 'PiecewiseFunction'
dump2pvd_compressedpvdDisplay.DataAxesGrid = 'GridAxesRepresentation'
dump2pvd_compressedpvdDisplay.PolarAxes = 'PolarAxesRepresentation'
dump2pvd_compressedpvdDisplay.ScalarOpacityFunction = fPWF
dump2pvd_compressedpvdDisplay.ScalarOpacityUnitDistance = 0.49414351405732404

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
dump2pvd_compressedpvdDisplay.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for fLUT in view renderView2
fLUTColorBar = GetScalarBar(fLUT, renderView2)
fLUTColorBar.WindowLocation = 'UpperLeftCorner'
fLUTColorBar.Title = 'f'
fLUTColorBar.ComponentTitle = ''
fLUTColorBar.LabelFormat = '%-#6.2g'
fLUTColorBar.RangeLabelFormat = '%6.2g'
fLUTColorBar.ScalarBarLength = 0.3300000000000002

# set color bar visibility
fLUTColorBar.Visibility = 1

# show color legend
dump2pvd_compressedpvdDisplay.SetScalarBarVisibility(renderView2, True)

# ----------------------------------------------------------------
# setup the visualization in view 'renderView3'
# ----------------------------------------------------------------

# show data from slice6
slice6Display = Show(slice6, renderView3, 'GeometryRepresentation')

# trace defaults for the display properties.
slice6Display.Representation = 'Surface'
slice6Display.ColorArrayName = ['POINTS', 'absOmega']
slice6Display.LookupTable = absOmegaLUT
slice6Display.OSPRayScaleArray = 'absOmega'
slice6Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice6Display.SelectOrientationVectors = 'None'
slice6Display.ScaleFactor = 0.08209827542304993
slice6Display.SelectScaleArray = 'absOmega'
slice6Display.GlyphType = 'Arrow'
slice6Display.GlyphTableIndexArray = 'absOmega'
slice6Display.GaussianRadius = 0.004104913771152496
slice6Display.SetScaleArray = ['POINTS', 'absOmega']
slice6Display.ScaleTransferFunction = 'PiecewiseFunction'
slice6Display.OpacityArray = ['POINTS', 'absOmega']
slice6Display.OpacityTransferFunction = 'PiecewiseFunction'
slice6Display.DataAxesGrid = 'GridAxesRepresentation'
slice6Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
slice6Display.OSPRayScaleFunction.Points = [0.001414213562373095, 0.0, 0.5, 0.0, 1.4142135623730951, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
slice6Display.ScaleTransferFunction.Points = [0.5690040558055741, 0.0, 0.5, 0.0, 9.136316321136388, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
slice6Display.OpacityTransferFunction.Points = [0.5690040558055741, 0.0, 0.5, 0.0, 9.136316321136388, 1.0, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for absOmegaLUT in view renderView3
absOmegaLUTColorBar = GetScalarBar(absOmegaLUT, renderView3)
absOmegaLUTColorBar.WindowLocation = 'AnyLocation'
absOmegaLUTColorBar.Position = [0.10660377358490566, 0.6432384341637011]
absOmegaLUTColorBar.Title = 'absOmega'
absOmegaLUTColorBar.ComponentTitle = ''
absOmegaLUTColorBar.LabelFormat = '%-#6.2g'
absOmegaLUTColorBar.RangeLabelFormat = '%6.2g'
absOmegaLUTColorBar.ScalarBarLength = 0.3300000000000002

# set color bar visibility
absOmegaLUTColorBar.Visibility = 1

# show color legend
slice6Display.SetScalarBarVisibility(renderView3, True)

# ----------------------------------------------------------------
# setup the visualization in view 'spreadSheetView1'
# ----------------------------------------------------------------

# show data from slice2
slice2Display = Show(slice2, spreadSheetView1, 'SpreadSheetRepresentation')

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get separate opacity transfer function/opacity map for 'absOmega'
separate_slice4Display_absOmegaPWF = GetOpacityTransferFunction('absOmega', slice4Display, separate=True)
separate_slice4Display_absOmegaPWF.Points = [2.2554501288418073e-07, 0.0, 0.5, 0.0, 9.294432572752617, 1.0, 0.5, 0.0]
separate_slice4Display_absOmegaPWF.ScalarRangeInitialized = 1

# get opacity transfer function/opacity map for 'absOmega'
absOmegaPWF = GetOpacityTransferFunction('absOmega')
absOmegaPWF.Points = [0.0, 0.0, 0.5, 0.0, 23.324664378194903, 1.0, 0.5, 0.0]
absOmegaPWF.ScalarRangeInitialized = 1

# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(extractSurface1)
# ----------------------------------------------------------------