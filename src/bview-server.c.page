/**
# The Basilisk View server

This simple server waits for [Basilisk View
commands](view.h#load-read-drawing-commands-from-a-file-or-buffer) on
standard input and returns PPM images on standard output.

Optional file names of either Basilisk View command files (with the
*.bv* extension) or [Basilisk snapshot
files](output.h#dump-basilisk-snapshots) can be given on the command
line.

It can run in parallel using MPI.

It is typically connected to a client such as [bview-client.py](), as
automated by the [bview]() script. 

## Installation

The server is not installed by default and relies on additional
libraries which need to be [installed](gl/INSTALL) first.

Once either `libfb_osmesa.a` or `libfb_glx.a` have been installed, the
`config` file needs to define the *OPENGLIBS* environment variable
accordingly i.e. using either

~~~bash
OPENGLIBS = -lfb_osmesa -lGLU -lOSMesa
~~~

or

~~~bash
OPENGLIBS = -lfb_glx -lGLU -lGLEW -lGL -lX11
~~~

For Mac OSX use:

~~~bash
OPENGLIBS = -L/opt/local/lib/ -lfb_osmesa -lGLU -lOSMesa
~~~

See [config.gcc]() and [config.osx]() for examples. Once this variable
is defined in *config*, the Basilisk servers can be compiled using:

~~~bash
cd $BASILISK
make bview-servers
~~~

For OSX 10.9 and above, you may see some warnings about deprecated
functions. You can safely ignore them.

## Implementation
*/

#include "view.h"

int main (int argc, char * argv[])
{
  Array * history = array_new();
  
  view (samples = 1);

  for (int i = 1; i < argc; i++) {
    if (strlen (argv[i]) >= 3 &&
	!strcmp (&argv[i][strlen(argv[i]) - 3], ".bv")) {
      if (!load (file = argv[i], history = history))
	exit (1);
    }
    else {
      if (!restore (file = argv[i], list = all)) {
	fprintf (ferr, "bview-server: could not restore from '%s'\n", argv[i]);
	exit (1);
      }
      restriction (all);
      fields_stats();
    }
  }

  if (history->len)
    save (fp = stdout);

  char line[256];
  FILE * interactive = stdout;
  do {
    line[0] = '\0';
    if (pid() == 0)
      fgets (line, 256, stdin);
#if _MPI
    MPI_Bcast (line, 256, MPI_BYTE, 0, MPI_COMM_WORLD);
#endif
    if (!strcmp (line, "interactive (true);\n"))
      interactive = stdout, line[0] = '\0';
    else if (!strcmp (line, "interactive (false);\n"))
      interactive = NULL, line[0] = '\0';
  } while (process_line (line, history, interactive));
  array_free (history);
}
