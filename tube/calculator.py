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
renderView1.ViewSize = [980, 1354]
renderView1.InteractionMode = 'Selection'
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.CenterOfRotation = [2.1022379343433126, -0.00020830861001330803, -0.004714558376574363]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [2.1022379343433126, -0.00020830861001330803, 3.5940055904215447]
renderView1.CameraFocalPoint = [2.1022379343433126, -0.00020830861001330803, -0.004714558376574363]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 0.9466363053263476
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.OSPRayMaterialLibrary = materialLibrary1

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
layout1.SplitHorizontal(0, 0.500000)
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
tube_bp_from_dumppvd = PVDReader(FileName='/Users/weugene/basilisk/work/tube/tube_bp_from_dump.pvd')
tube_bp_from_dumppvd.CellArrays = ['p', 'fs', 'f', 'l', 'residual_of_p', 'u.x']

# create a new 'Cell Data to Point Data'
cellDatatoPointData1 = CellDatatoPointData(Input=tube_bp_from_dumppvd)
cellDatatoPointData1.CellDataArraytoprocess = ['f', 'fs', 'l', 'p', 'residual_of_p', 'u.x']

# create a new 'Calculator'
calculator1 = Calculator(Input=cellDatatoPointData1)
calculator1.Function = 'coordsX*(1-f)/0.26'

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