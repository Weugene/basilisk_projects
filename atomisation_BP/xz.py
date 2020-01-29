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
