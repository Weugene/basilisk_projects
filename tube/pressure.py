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
renderView1.ViewSize = [1049, 222]
renderView1.InteractionMode = '2D'
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.CenterOfRotation = [15.0, 0.0, 0.0]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [15.0, 0.0, 15.269983825650256]
renderView1.CameraFocalPoint = [15.0, 0.0, 0.0]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 3.266250109492342
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

# create a new 'PVD Reader'
#path = 'G:\\My Drive\\Skoltech\\Dumps\\tube\\10e_coarse\\'
path = '/home/e.sharaborin/basilisk/work/tube/' 
dump2pvd_compressedpvd = PVDReader(FileName=path + 'dump2pvd_compressed.pvd')
dump2pvd_compressedpvd.CellArrays = ['fs', 'f', 'l', 'l2', 'omega', 'u.x']

# create a new 'Slice'
slice1 = Slice(Input=dump2pvd_compressedpvd)
slice1.SliceType = 'Plane'
slice1.HyperTreeGridSlicer = 'Plane'
slice1.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Origin = [15.0, 0.0, 0.0]
slice1.SliceType.Normal = [0.0, 0.0, 1.0]

# init the 'Plane' selected for 'HyperTreeGridSlicer'
slice1.HyperTreeGridSlicer.Origin = [15.0, 0.0, 0.0]

# create a new 'Clip'
clip1 = Clip(Input=slice1)
clip1.ClipType = 'Box'
clip1.HyperTreeGridClipper = 'Plane'
clip1.Scalars = ['CELLS', 'f']
clip1.Value = 0.5

# init the 'Box' selected for 'ClipType'
clip1.ClipType.Position = [0.0, -0.5, -0.5]
clip1.ClipType.Length = [30.0, 1.0, 1.0]

# init the 'Plane' selected for 'HyperTreeGridClipper'
clip1.HyperTreeGridClipper.Origin = [15.0, 0.0, 0.0]

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from dump2pvd_compressedpvd
dump2pvd_compressedpvdDisplay = Show(dump2pvd_compressedpvd, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
dump2pvd_compressedpvdDisplay.Representation = 'Surface'
dump2pvd_compressedpvdDisplay.ColorArrayName = [None, '']
dump2pvd_compressedpvdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
dump2pvd_compressedpvdDisplay.SelectOrientationVectors = 'None'
dump2pvd_compressedpvdDisplay.ScaleFactor = 3.0
dump2pvd_compressedpvdDisplay.SelectScaleArray = 'None'
dump2pvd_compressedpvdDisplay.GlyphType = 'Arrow'
dump2pvd_compressedpvdDisplay.GlyphTableIndexArray = 'None'
dump2pvd_compressedpvdDisplay.GaussianRadius = 0.15
dump2pvd_compressedpvdDisplay.SetScaleArray = [None, '']
dump2pvd_compressedpvdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
dump2pvd_compressedpvdDisplay.OpacityArray = [None, '']
dump2pvd_compressedpvdDisplay.OpacityTransferFunction = 'PiecewiseFunction'
dump2pvd_compressedpvdDisplay.DataAxesGrid = 'GridAxesRepresentation'
dump2pvd_compressedpvdDisplay.PolarAxes = 'PolarAxesRepresentation'
dump2pvd_compressedpvdDisplay.ScalarOpacityUnitDistance = 0.4920658244349408

# show data from slice1
slice1Display = Show(slice1, renderView1, 'GeometryRepresentation')

# get color transfer function/color map for 'f'
fLUT = GetColorTransferFunction('f')
fLUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
slice1Display.Representation = 'Surface'
slice1Display.ColorArrayName = ['CELLS', 'f']
slice1Display.LookupTable = fLUT
slice1Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice1Display.OrientationMode = 'Rotation'
slice1Display.SelectOrientationVectors = 'None'
slice1Display.ScaleFactor = 3.0
slice1Display.SelectScaleArray = 'None'
slice1Display.GlyphType = 'Arrow'
slice1Display.GlyphTableIndexArray = 'None'
slice1Display.GaussianRadius = 0.15
slice1Display.ShaderPreset = 'Plain circle'
slice1Display.SetScaleArray = [None, '']
slice1Display.ScaleTransferFunction = 'PiecewiseFunction'
slice1Display.OpacityArray = [None, '']
slice1Display.OpacityTransferFunction = 'PiecewiseFunction'
slice1Display.DataAxesGrid = 'GridAxesRepresentation'
slice1Display.PolarAxes = 'PolarAxesRepresentation'

# show data from clip1
clip1Display = Show(clip1, renderView1, 'UnstructuredGridRepresentation')

# get opacity transfer function/opacity map for 'f'
fPWF = GetOpacityTransferFunction('f')
fPWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
clip1Display.Representation = 'Surface'
clip1Display.ColorArrayName = ['CELLS', 'f']
clip1Display.LookupTable = fLUT
clip1Display.OSPRayScaleFunction = 'PiecewiseFunction'
clip1Display.SelectOrientationVectors = 'None'
clip1Display.ScaleFactor = 3.0
clip1Display.SelectScaleArray = 'None'
clip1Display.GlyphType = 'Arrow'
clip1Display.GlyphTableIndexArray = 'None'
clip1Display.GaussianRadius = 0.15
clip1Display.SetScaleArray = [None, '']
clip1Display.ScaleTransferFunction = 'PiecewiseFunction'
clip1Display.OpacityArray = [None, '']
clip1Display.OpacityTransferFunction = 'PiecewiseFunction'
clip1Display.DataAxesGrid = 'GridAxesRepresentation'
clip1Display.PolarAxes = 'PolarAxesRepresentation'
clip1Display.ScalarOpacityFunction = fPWF
clip1Display.ScalarOpacityUnitDistance = 1.0136285287300626

# setup the color legend parameters for each legend in this view

# get color legend/bar for fLUT in view renderView1
fLUTColorBar = GetScalarBar(fLUT, renderView1)
fLUTColorBar.Orientation = 'Horizontal'
fLUTColorBar.WindowLocation = 'AnyLocation'
fLUTColorBar.Position = [0.3340807465443276, 0.08746606334841664]
fLUTColorBar.Title = 'f'
fLUTColorBar.ComponentTitle = ''
fLUTColorBar.TitleFontSize = 10
fLUTColorBar.LabelFontSize = 9
fLUTColorBar.ScalarBarLength = 0.32999999999999974

# set color bar visibility
fLUTColorBar.Visibility = 1

# hide data in view
Hide(dump2pvd_compressedpvd, renderView1)

# show color legend
slice1Display.SetScalarBarVisibility(renderView1, True)

# hide data in view
Hide(slice1, renderView1)

# show color legend
clip1Display.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(clip1)
# ----------------------------------------------------------------

fn = path  +  'pressure.png'
SaveScreenshot( fn, renderView1,
	TransparentBackground=0,
	CompressionLevel='2' )
print('File=' + fn + ' generated succesfully')
