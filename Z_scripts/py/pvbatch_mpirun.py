from paraview.simple import *

renderView1 = CreateView('RenderView')
renderView1.ViewSize = [5000, 4000]
sp=Sphere()
Show()
Render()
SaveScreenshot( "resolution_test.png", renderView1,
#      ImageResolution=[2316, 2204],
    TransparentBackground=0,
#     magnification=1,
    # PNG options
    CompressionLevel='2' )
writer = XMLPPolyDataWriter(FileName='sphere.pvtp')
writer.UpdatePipeline()
