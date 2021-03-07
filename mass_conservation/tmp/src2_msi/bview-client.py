"""
# A Graphical User Interface for Basilisk View

The interactive version of [Basilisk View](view.h) relies on a
[client/server
model](https://en.wikipedia.org/wiki/Client-server_model). The
[Basilisk View server](bview-server.c) can run on a distant system
(for example the large parallel machine on which the Basilisk runs are
done), or on the local system.

A Basilisk View client running on the local system then sends
"commands" to the server and receives in return a stream of PPM images
containing the updated views.

The commands are just function calls, sent as a stream of text by the
client, using the same syntax as for the
[*load()*](view.h#load-read-drawing-commands-from-a-file-or-buffer)
function.

The following Python code implements a simple Graphical User Interface
(GUI) client using this model.

The client standard input should be connected (typically through [Unix
pipes](https://en.wikipedia.org/wiki/Pipeline_(Unix)) or [named
pipes](https://en.wikipedia.org/wiki/Named_pipe)) to the standard
output of the server. Conversely, the client standard output should be
connected to the server's standard input. For convenience, this is
typically done using the [bview]() shell script.

The client then creates a
[Tkinter](https://wiki.python.org/moin/TkInter) window and waits both
for user interaction (using the mouse or keyboard) and for images sent
by the server on input. This is done concurrently using two threads.

User interaction is converted to Basilisk View commands which are
written to standard output, and thus sent to the server, which in turn
responds with images which are used to refresh the Tkinter window.

## Summary of controls

* Left-mouse button + drag: rotate camera.
* Right-mouse button + drag: translate camera.
* Center-mouse button + drag: zoom.
* Mouse wheel: zoom.
* '+' or '-' keys: zoom.
* 'B' key: calls [*box()*](draw.h#box-displays-box-boundaries-and-axis-coordinates).
* 'c' key: calls [*clear()*](draw.h#clear-removes-all-objects-previously-drawn).
* 'l', 'r', 't', 'b', 'f', 'z', 'i' keys: changes view to 'left', 'right',
  'top', 'bottom', 'front', 'back', 'iso' respectively.
* Right/left arrow keys: rotates right/left by 3 degrees.
* Up/down arrow keys: rotates up/down by 3 degrees.
* 'q' key: quits.
* '<' or '>' keys: decrease or increase minimum delay between screen refreshes.
* '1' and '4': set the number of samples (1 fast/coarse, 4 slow/fine).
* '8' and '9': rotate view clockwise (resp. anti-clockwise).

## User customisation

If it exists, the script will import the file `$HOME/.bview.py`. This
can be used to customise the interface. For example, if one wants to
define a new keyboard mapping, one could use the following in
`$HOME/.bview.py`.

~~~python
from bview import keymap
keymap["'v'"] = 'draw_vof("f")\n'
~~~

## Dependencies

The current version will work with python2.7 but not
python3.x. Besides Tkinter, the program also uses the ImageTk
extension of the Python Imaging Library (PIL).

### Debian-like systems (Debian, Ubuntu etc.)

The required dependencies can be installed easily using:

~~~bash
sudo apt-get install python-pil.imagetk
~~~

or for older Debian versions ($\leq 7$):

~~~bash
sudo apt-get install python-imaging
~~~

### Mac OSX

To check whether Tkinter is correctly installed on your system, do:

~~~bash
python -m Tkinter
~~~

If a window pops up with a click button, then it works. Otherwise you
will need to install it.

It is recommended to reinstall python as well. See:

* [python 2.7](https://www.python.org/downloads/release/python-2714/)
* [Tkinter](https://www.python.org/download/mac/tcltk/)

You may also need to install the PIL (or Pillow) module. This can
be done using:

~~~bash
sudo easy_install pip
pip install pillow
~~~

# Implementation

"""

from Tkinter import Tk, mainloop
from sys import stdout, stderr
from bview import BCanvas

# Try user customization
from os.path import expanduser
import imp
try:
    imp.load_source('userinit', expanduser("~/.bview.py"))
except IOError:
    pass
        
root = Tk()
root.withdraw()
root.title ('bview')
root.canvas = BCanvas (root, stdout)
stdout = stderr

try:
    mainloop()
except (KeyboardInterrupt, SystemExit):
    None
