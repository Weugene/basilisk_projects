import paraview.simple as pv

path = '/Applications/ParaView-5.6.0.for\ Arkuda.app/Contents/examples/'
fn =  'can.ex2'
print(fn)
reader = pv.ExodusIIReader(FileName=fn)
pv.Show(reader); pv.Render()
pv.WriteAnimation('movie0.avi', Compression=1, Quality=2)

pv.SaveAnimation('movie.avi', ImageQuality=50) # <== works (Writer->Quality=1, Writer->Compression=true)
pv.SaveAnimation('movie2.avi', ImageQuality=75) # <== fails (Writer->Quality=2, Writer->Compression=false)

