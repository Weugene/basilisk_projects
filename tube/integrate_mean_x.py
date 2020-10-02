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
renderView1.ViewSize = [1128, 1124]
renderView1.InteractionMode = 'Selection'
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.CenterOfRotation = [4.618346055841663, -0.00012084120092947792, -0.00014717160817770414]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [11.24488764252369, 59.2309828436138, 49.32020070713154]
renderView1.CameraFocalPoint = [4.618346055841664, -0.00012084120092961803, -0.00014717160817782083]
renderView1.CameraViewUp = [0.1769483238549785, 0.6180233083864424, -0.7659872590167477]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 1.67999793815955
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.OSPRayMaterialLibrary = materialLibrary1

# Create a new 'SpreadSheet View'
spreadSheetView1 = CreateView('SpreadSheetView')
spreadSheetView1.ColumnToSort = ''
spreadSheetView1.BlockSize = 1024L
spreadSheetView1.FieldAssociation = 'Cell Data'
# uncomment following to set a specific view size
# spreadSheetView1.ViewSize = [400, 400]

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.SplitHorizontal(0, 0.541154)
layout1.AssignView(1, renderView1)
layout1.AssignView(2, spreadSheetView1)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'PVD Reader'
pVDReader1 = PVDReader(FileName='/Users/weugene/basilisk/work/tube/dump2pvd_compressed.pvd')
pVDReader1.CellArrays = ['fs', 'f', 'l', 'l2', 'omega', 'u.x']

# create a new 'Cell Data to Point Data'
cellDatatoPointData1 = CellDatatoPointData(Input=pVDReader1)
cellDatatoPointData1.CellDataArraytoprocess = ['f', 'fs', 'l', 'l2', 'omega', 'u.x']

# create a new 'Calculator'
calculator1 = Calculator(Input=cellDatatoPointData1)
calculator1.ResultArrayName = 'meanX'
calculator1.Function = 'coordsX*(1-f)*(1-fs)'

# create a new 'Integrate Variables'
integrateVariables1 = IntegrateVariables(Input=calculator1)

# ----------------------------------------------------------------
# setup the visualization in view 'spreadSheetView1'
# ----------------------------------------------------------------

# show data from integrateVariables1
integrateVariables1Display = Show(integrateVariables1, spreadSheetView1, 'SpreadSheetRepresentation')

# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(calculator1)
# ----------------------------------------------------------------