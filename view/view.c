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
