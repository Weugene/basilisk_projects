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
from paraview.servermanager import Fetch

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')


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
save_isosurfacepvd = PVDReader(FileName='/home/e.sharaborin/basilisk/work/tube/save_isosurface.pvd')
#save_isosurfacepvd.CellArrays = ['vtkGhostType']
#save_isosurfacepvd.PointArrays = ['absOmega', 'f', 'fs', 'l', 'l2', 'omega', 'p', 'u.x', 'vtkValidPointMask', 'vtkGhostType']

# create a new 'Calculator'
calculator1 = Calculator(Input=save_isosurfacepvd)
calculator1.Function = 'sqrt(coordsY^2+coordsZ^2)'

# create a new 'Clip'
clip1 = Clip(Input=calculator1)
clip1.ClipType = 'Box'
clip1.HyperTreeGridClipper = 'Plane'
clip1.Scalars = ['POINTS', 'absOmega']
clip1.Value = 13.388705481122626

# init the 'Box' selected for 'ClipType'
clip1.ClipType.Position = [13.9, -0.6, -0.6]
clip1.ClipType.Length = [1.5, 1.1999999999999984, 1.1999999999999993]

# init the 'Plane' selected for 'HyperTreeGridClipper'
clip1.HyperTreeGridClipper.Origin = [4.61862313747406, -8.401274681091309e-05, -0.00011083483695983887]

# create a new 'Clip'
clip2 = Clip(Input=clip1)
clip2.ClipType = 'Cylinder'
clip2.HyperTreeGridClipper = 'Plane'
clip2.Scalars = ['POINTS', 'absOmega']
clip2.Value = 13.388705481122626
clip2.Invert = 0

# init the 'Cylinder' selected for 'ClipType'
clip2.ClipType.Center = [14.649999618530273, 0.0003166794776916504, -0.00022161006927490234]
clip2.ClipType.Axis = [1.0, 0.0, 0.0]
clip2.ClipType.Radius = 0.4

# init the 'Plane' selected for 'HyperTreeGridClipper'
clip2.HyperTreeGridClipper.Origin = [4.61862313747406, -8.401274681091309e-05, -0.00011083483695983887]

# create a new 'Integrate Variables'
integrateVariables1 = IntegrateVariables(Input=clip2)

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from integrateVariables1
integrateVariables1Display = Show(integrateVariables1, renderView1, 'UnstructuredGridRepresentation')

# get animation scene
animationScene1 = GetAnimationScene()

# get the time-keeper
timeKeeper1 = GetTimeKeeper()
timesteps = timeKeeper1.TimestepValues # 0, 0.1, 0.2 ...

i=3
animationScene1.AnimationTime = timesteps[i]
# Properties modified on timeKeeper1
timeKeeper1.Time = timesteps[i]

# Create a new 'SpreadSheet View'
spreadSheetView1 = CreateView('SpreadSheetView')       

# create a new 'Pass Arrays'
passArrays1 = PassArrays(Input=integrateVariables1)
passArrays1.PointDataArrays = [ 'Result']
passArrays1.CellDataArrays = ['Area']
  

# update the view to ensure updated data information
spreadSheetView1.Update()

ss_data = paraview.servermanager.Fetch(passArrays1)    
area = ss_data.GetCellData().GetArray('Area').GetValue(0)
rad_dv = ss_data.GetPointData().GetArray('Result').GetValue(0)
print('N=', ss_data.GetNumberOfPoints(), 'Area=', area, 'rad_dv=', rad_dv, 'rad_mean=', rad_dv/area)


# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(integrateVariables1)
# ----------------------------------------------------------------
