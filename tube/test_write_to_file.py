# trace generated using paraview version 5.8.0
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'XML Partitioned Unstructured Grid Reader'
dump2pvd_compressed_0_000 = XMLPartitionedUnstructuredGridReader(FileName=['/home/e.sharaborin/basilisk/work/tube/res/dump2pvd_compressed_0_0000.pvtu', '/home/e.sharaborin/basilisk/work/tube/res/dump2pvd_compressed_0_0001.pvtu', '/home/e.sharaborin/basilisk/work/tube/res/dump2pvd_compressed_0_0002.pvtu', '/home/e.sharaborin/basilisk/work/tube/res/dump2pvd_compressed_0_0003.pvtu', '/home/e.sharaborin/basilisk/work/tube/res/dump2pvd_compressed_0_0004.pvtu', '/home/e.sharaborin/basilisk/work/tube/res/dump2pvd_compressed_0_0005.pvtu', '/home/e.sharaborin/basilisk/work/tube/res/dump2pvd_compressed_0_0006.pvtu', '/home/e.sharaborin/basilisk/work/tube/res/dump2pvd_compressed_0_0007.pvtu', '/home/e.sharaborin/basilisk/work/tube/res/dump2pvd_compressed_0_0008.pvtu', '/home/e.sharaborin/basilisk/work/tube/res/dump2pvd_compressed_0_0009.pvtu'])
dump2pvd_compressed_0_000.CellArrayStatus = ['fs', 'f', 'l', 'l2', 'omega', 'u.x']

# get animation scene
animationScene1 = GetAnimationScene()

# get the time-keeper
timeKeeper1 = GetTimeKeeper()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [914, 491]

# get layout
layout1 = GetLayout()

# show data in view
dump2pvd_compressed_0_000Display = Show(dump2pvd_compressed_0_000, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
dump2pvd_compressed_0_000Display.Representation = 'Surface'
dump2pvd_compressed_0_000Display.ColorArrayName = [None, '']
dump2pvd_compressed_0_000Display.OSPRayScaleFunction = 'PiecewiseFunction'
dump2pvd_compressed_0_000Display.SelectOrientationVectors = 'None'
dump2pvd_compressed_0_000Display.ScaleFactor = 3.0
dump2pvd_compressed_0_000Display.SelectScaleArray = 'None'
dump2pvd_compressed_0_000Display.GlyphType = 'Arrow'
dump2pvd_compressed_0_000Display.GlyphTableIndexArray = 'None'
dump2pvd_compressed_0_000Display.GaussianRadius = 0.15
dump2pvd_compressed_0_000Display.SetScaleArray = [None, '']
dump2pvd_compressed_0_000Display.ScaleTransferFunction = 'PiecewiseFunction'
dump2pvd_compressed_0_000Display.OpacityArray = [None, '']
dump2pvd_compressed_0_000Display.OpacityTransferFunction = 'PiecewiseFunction'
dump2pvd_compressed_0_000Display.DataAxesGrid = 'GridAxesRepresentation'
dump2pvd_compressed_0_000Display.PolarAxes = 'PolarAxesRepresentation'
dump2pvd_compressed_0_000Display.ScalarOpacityUnitDistance = 0.49206582443494074

# reset view to fit data
renderView1.ResetCamera()

# get the material library
materialLibrary1 = GetMaterialLibrary()

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Clip'
clip1 = Clip(Input=dump2pvd_compressed_0_000)
clip1.ClipType = 'Plane'
clip1.HyperTreeGridClipper = 'Plane'
clip1.Scalars = ['CELLS', 'f']
clip1.Value = 0.5

# init the 'Plane' selected for 'ClipType'
clip1.ClipType.Origin = [15.0, 0.0, 0.0]

# init the 'Plane' selected for 'HyperTreeGridClipper'
clip1.HyperTreeGridClipper.Origin = [15.0, 0.0, 0.0]

# toggle 3D widget visibility (only when running from the GUI)
Show3DWidgets(proxy=clip1.ClipType)

# Properties modified on clip1
clip1.ClipType = 'Box'

# Properties modified on clip1.ClipType
clip1.ClipType.Position = [0.0, -0.6, -0.6]
clip1.ClipType.Length = [30.0, 1.2, 1.2]

# show data in view
clip1Display = Show(clip1, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
clip1Display.Representation = 'Surface'
clip1Display.ColorArrayName = [None, '']
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
clip1Display.ScalarOpacityUnitDistance = 0.2853051783454308

# hide data in view
Hide(dump2pvd_compressed_0_000, renderView1)

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Cell Data to Point Data'
cellDatatoPointData1 = CellDatatoPointData(Input=clip1)
cellDatatoPointData1.CellDataArraytoprocess = ['f', 'fs', 'l', 'l2', 'omega', 'u.x']

# show data in view
cellDatatoPointData1Display = Show(cellDatatoPointData1, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
cellDatatoPointData1Display.Representation = 'Surface'
cellDatatoPointData1Display.ColorArrayName = [None, '']
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
cellDatatoPointData1Display.ScalarOpacityUnitDistance = 0.2853051783454308

# hide data in view
Hide(clip1, renderView1)

# update the view to ensure updated data information
renderView1.Update()

# save data
SaveData('/home/e.sharaborin/basilisk/work/tube/res/test.pvd', proxy=cellDatatoPointData1, PointDataArrays=['f', 'fs', 'l', 'l2', 'omega', 'u.x'], WriteTimeSteps=1)

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.CameraPosition = [15.0, 0.0, 100.38195644853695]
renderView1.CameraFocalPoint = [15.0, 0.0, 0.0]
renderView1.CameraParallelScale = 25.98076211353316

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
