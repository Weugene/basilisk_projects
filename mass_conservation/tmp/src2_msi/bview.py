from Tkinter import *
from PIL import Image, ImageTk
from threading import Thread, Lock
from sys import stdin, stdout, platform
from os import fdopen, devnull
from time import sleep, time
from StringIO import StringIO
from math import pi

keymap = {
    "'-'"   : 'view (fov =  0.02, relative = 1);\n',
    "'+'"   : 'view (fov = -0.02, relative = 1);\n',
    "'='"   : 'view (fov = -0.02, relative = 1);\n',
    "'B'"   : 'box();\n',
    "'c'"   : 'clear();\n',
    "'l'"   : 'view (camera = "left");\n',
    "'r'"   : 'view (camera = "right");\n',
    "'t'"   : 'view (camera = "top");\n',
    "'b'"   : 'view (camera = "bottom");\n',
    "'f'"   : 'view (camera = "front");\n',
    "'z'"   : 'view (camera = "back");\n',
    "'i'"   : 'view (camera = "iso");\n',
    "'1'"   : 'view (samples = 1);\n',
    "'4'"   : 'view (samples = 4);\n',
    "'8'"   : 'view (psi = 1.5707963267949, relative = 1);\n',
    "'9'"   : 'view (psi = -1.5707963267949, relative = 1);\n',
    'Right' : 'view (theta = 0.0523598775598299, phi = 0, relative = 1);\n',
    'Left'  : 'view (theta = -0.0523598775598299, phi = 0, relative = 1);\n',
    'Up'    : 'view (theta = 0, phi = 0.0523598775598299, relative = 1);\n',
    'Down'  : 'view (theta = 0, phi = -0.0523598775598299, relative = 1);\n'
}

class BCanvas (Canvas):

    def draw(self, s):
        try:
            self.pipe.write (s)
            self.changed = True
        except IOError, e:
            if e.errno == 32: # Broken pipe
                self.root.destroy()
            else:
                raise

    def button1 (self, event):
        self.beginx, self.beginy = event.x, event.y
        self.motion = True

    def button_release (self, event):
        self.motion = False
        if self.res > 1:
            self.redrawing = False
            self.root.after (1000, self.redraw)

    def zoom (self, event):
        self.draw ("view (fov = %g, relative = 1);\n" %
                   ((event.y - self.beginy)/float(self.winfo_height())))
        self.beginx, self.beginy = event.x, event.y

    def zoomin (self, event):
        self.draw ("view (fov = -0.01, relative = 1);\n")

    def zoomout (self, event):
        self.draw ("view (fov = +0.01, relative = 1);\n")

    def zoomMac (self, event):
        self.draw("view (fov = %g, relative = 1);\n" % (event.delta/100.))

    def move (self, event):
        self.draw ("view (tx = %g, ty = %g, relative = 1);\n" %
                   ((event.x - self.beginx)/float(self.winfo_width()),
                    (self.beginy - event.y)/float(self.winfo_height())))
        self.beginx, self.beginy = event.x, event.y

    def trackball (self, event):
        width, height = self.winfo_reqwidth(), self.winfo_reqheight()
        self.draw ("view (p1x = %g, p1y = %g, p2x = %g, p2y = %g);\n" %
                   ((2.*self.beginx - width)/width,
                    (height - 2.*self.beginy)/height,
                    (2.*event.x - width)/width,
                    (height - 2.*event.y)/height))
        self.beginx, self.beginy = event.x, event.y

    def quit (self):
        try:
            self.pipe.write ("quit();\n")
        except IOError, e:
            if e.errno != 32: # Broken pipe
                raise
        self.root.destroy()

    def key (self, event):
        c = repr (event.char)
        if c == "'<'":
            if self.refresh > 1:
                self.refresh = self.refresh - 1
        elif c == "'>'":
            self.refresh = self.refresh + 1
        elif c == "'q'":
            self.quit()
        elif c in keymap:
            self.draw (keymap[c])
        elif event.keysym in keymap:
            self.draw (keymap[event.keysym])
#        else:
#            print c, event.keysym, event.keycode

    def resize_canvas (self):
        self.draw ('view (width = %d, height = %d);\n' %
              (self.winfo_width(), self.winfo_height()))
        self.resize = False

    def configure (self, event):
        if not self.resize:
            self.resize = True
            self.root.after (1000, self.resize_canvas)

    def listen_stdin (self):
        condition = True;
        while condition:
            input = sys.stdin.readline().split(' ')
            condition = input[0] == "P6"
            if condition:
                width = int(input[1])
                height = int(input[2])
                colors = int(input[3])
                data = sys.stdin.read(width*height*3)

                buf = StringIO()
                buf.write ("P6 %d %d %d\n" % (width, height, colors))
                buf.write (data)
                buf.seek (0)

                with self.new_image_lock:
                    self.new_image = Image.open (buf)

                sleep (self.refresh/1000.) # let the GUI catch up

    def redraw (self):
        if not self.redrawing:
            if self.motion:
                if self.dt > self.reactivity:
                    self.res *= 2
                    if self.res > 500:
                        self.res = 500
                elif self.dt < self.reactivity/2. and self.res > 1:
                    self.res /= 2
                    if self.res < 1:
                        self.res = 1
                self.pipe.write ("view (res = %d);\n" % self.res)
            else:
                self.pipe.write ("view (res = 1);\n")
            self.redrawing = True
            self.start = time()
            self.draw ("display();\n")
            self.changed = False

    def check_new_image (self):
        with self.new_image_lock:
            if self.new_image is not None:
                img = ImageTk.PhotoImage (self.new_image)
                self.config (width = img.width(), height = img.height())
                self.itemconfig (self.image, image = img)
                self.image_tk = img
                self.redrawing = False
                self.root.deiconify()
                t = time()
                self.dt = t - self.start
                self.start = t
#                print "redrawn in", self.dt, "sec"
                self.new_image = None
            else:
                if self.changed:
                    self.redraw()
        self.root.after (self.refresh, self.check_new_image)

    def __init__(self, root, out):
        Canvas.__init__(self, root, bd = 0, highlightthickness = 0)
        self.root = root
        self.refresh = 10
        self.new_image = None
        self.new_image_lock = Lock()

        if out.isatty():
            self.pipe = open (devnull, "w")
        else:
            self.pipe = fdopen (out.fileno(), 'w', 1)

        self.beginx = 0
        self.beginy = 0
        self.changed = False
        self.resize = False
        self.redrawing = False
        self.motion = False
        self.image = self.create_image (0, 0, anchor = NW)
        self.start = time()
        self.dt = 0
        self.res = 10
        self.reactivity = 0.3
        self.focus_set()
        self.bind ("<Key>", self.key)
        self.bind ("<Button-1>", self.button1)
        self.bind ("<Button-2>", self.button1)
        self.bind ("<Button-3>", self.button1)
        self.bind ("<ButtonRelease-1>", self.button_release)
        self.bind ("<ButtonRelease-2>", self.button_release)
        self.bind ("<ButtonRelease-3>", self.button_release)
        self.bind ("<B1-Motion>", self.trackball)
        if platform == "darwin":
            self.bind ("<B3-Motion>", self.zoom)
            self.bind ("<B2-Motion>", self.move)
            self.bind ("<MouseWheel>", self.zoomMac)
        else:
            self.bind ("<B2-Motion>", self.zoom)
            self.bind ("<B3-Motion>", self.move)
            self.bind ("<Button-4>", self.zoomin)
            self.bind ("<Button-5>", self.zoomout)
        self.pack (fill = BOTH, expand = 1)
        self.bind ("<Configure>", self.configure)

        self.root.after (self.refresh, self.check_new_image)
        self.root.protocol ("WM_DELETE_WINDOW", self.quit)
    
        self.t = Thread (target = self.listen_stdin)
        self.t.daemon = True
        self.t.start()
