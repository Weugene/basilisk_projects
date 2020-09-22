# This script demonstrates how to track a moving object.
#
# Each time it is executed, it marks the active view
# and object. It also adds a track to the animation
# that moves the camera with the center of the bounds
# of that object.
#
# To use, setup the view like you want the camera to
# look, run this script and then run the animation.
#
# To disable tracking, select the server object in the
# pipeline browser and rerun this script.

from paraview.simple import *

ans = None
try:
    ans = GetAnimationScene()
except NameError:
    for i in servermanager.ProxyManager().GetProxiesInGroup("animation").values():
        if i.IsA("vtkSMAnimationSceneProxy"):
            ans = i
            break

if ans:
    # Turn of caching so that we can get the bounds for the
    # data for each time step
    #ans.Caching = False
    scene = ans.GetClientSideObject()
    # Remove the previous observer.
    try:
        scene.RemoveObserver(callbackid)
    except NameError:
        pass
    def callback(caller, *args):
        global rv, source, first_bounds, first_focal, first_pos, view_up
        if source:
            #The following line does not work in Paraview 3.10.1 -- the error printed is "AttributeError: GetViewUpdateTime". The script seems to work fine without it.
            #source.UpdatePipeline(rv.GetViewUpdateTime())
            bds = source.GetDataInformation().GetBounds()
            # Center of bounds
            bds = map(lambda x,y : (x+y)/2, bds[0:6:2], bds[1:6:2])
            # Difference between initial bounds and current bounds
            bds_diff = map(lambda x,y : x-y, bds, first_bounds)
            # Move the focal point as much as the center of the bounds
            # moved
            rv.CameraFocalPoint = map(lambda x,y: x+y, first_focal, bds_diff)
            # Move the camera as much as the center of the bounds
            # moved
            rv.CameraPosition = map(lambda x,y: x+y, first_pos, bds_diff)
            # Keep the original view up
            rv.CameraViewUp = view_up

    rv = GetActiveView()

    if rv:
        # Add the observer to the VTK scene with priority 1 so that it happens before
        # anything else
        callbackid = scene.AddObserver("AnimationCueTickEvent", callback, 1.0)

        source = GetActiveSource()
        bds = source.GetDataInformation().GetBounds()
        # Center of bounds
        first_bounds = map(lambda x,y : (x+y)/2, bds[0:6:2], bds[1:6:2])
        first_focal = rv.CameraFocalPoint[:]
        first_pos = rv.CameraPosition[:]
        view_up = rv.CameraViewUp[:]

def removeCameraObserver():
    global callbackid, scene
    scene.RemoveObserver(callbackid)
