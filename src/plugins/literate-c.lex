%option reentrant noyywrap noinput extra-type="struct MyScanner *"
%{
  #include <stdlib.h>
  #include <stdarg.h>
  #include <sys/wait.h>
  #include <sys/types.h>
  #include <sys/stat.h>
  #include <unistd.h>
  #include <glob.h>

  typedef struct {
    char * error, * warning;
  } Error;
  
  static char C[]      = "\n~~~literatec";
  static char Python[] = "\n~~~python";
  static char Octave[] = "\n~~~matlab";
  static char Bash[] = "\n~~~bash";

  struct MyScanner {
    char * input;
    char * output, * page, * gnuplot, * gnuplot_output, * type;
    int i, len, ncodes, incode, first, line, indent;
    Error * error;
    int nerror, nplots;
  };

  #define output_c(c) output_c1 (yyextra, c)
  static void output_c1 (struct MyScanner * scan, int c) {
    if (scan->i == scan->len) {
      scan->len += 100;
      scan->output = realloc (scan->output, scan->len + 1);
    }
    scan->output[scan->i++] = c;
    scan->output[scan->i] = '\0';
  }

  #define output_s(s) output_s1 (yyextra, s)
  static void output_s1 (struct MyScanner * scan, char * s) {
    while (*s)
      output_c1 (scan, *s++);
  }

  static void error_start (struct MyScanner * scan) {
     output_s1 (scan,
		"<div class=error>"
		"<div id=msg_logo>"
		"<img src=/img/error.png>"
		"</div>"
		"<div id=msg_label>");
  }

  static void warning_start (struct MyScanner * scan) {
     output_s1 (scan,
		"<div class=message>"
		"<div id=msg_logo>"
		"<img src=/img/warning.png>"
		"</div>"
		"<div id=msg_label>");
  }

  static void error_end (struct MyScanner * scan) {
    output_s1 (scan, "</div></div>");
  }

  static void check_error (struct MyScanner * scan) {
    if (scan->line < scan->nerror) {
      Error * e = &scan->error[scan->line];
      if (e->error || e->warning) {
	if (!scan->first && scan->incode)
	  output_s1 (scan, "\n~~~\n");
	if (e->error) {
	  error_start (scan);
	  output_s1 (scan, e->error);
	  free (e->error); e->error = NULL;
	}
	else {
	  warning_start (scan);
	  output_s1 (scan, e->warning);
	  free (e->warning); e->warning = NULL;
	}
	error_end (scan);
	if (!scan->first && scan->incode)
	  output_s1 (scan, scan->type);
      }
    }
  }

  static int spacenb (char * s) {
    int ns = 0;
    while (*s != '\0' && strchr (" \t\v\n\f", *s)) {
      switch (*s) {
      case ' ':  ns++; break;
      case '\t': ns += 8; break;
      case '\n': case '\v': case '\f': ns = 0; break;
      }
      s++;
    }
    return ns;
  }

  static char * acat (const char * s, ...)
  {
    int len = strlen (s) + 1;
    char * c = malloc (len);
    strcpy (c, s);

    va_list ap;
    va_start (ap, s);
    char * s1;
    while ((s1 = va_arg (ap, char *))) {
      len += strlen (s1);
      c = realloc (c, len);
      strcat (c, s1);
    }
    va_end (ap);

    return c;
  }

  static char * stamped_image (char * name, struct MyScanner * scan) {
    char * path = name[0] == '/' ? strdup (&name[1]) : 
                                   acat (scan->page, "/", name, NULL);
    char * link;
    struct stat buf;
    if (!stat (path, &buf)) {
      // create a symbolic link by appending the modification time
      // to force the browser to reload the updated image
      char * dot = path + strlen(path) - 1;
      while (dot != path && *dot != '.') dot--;
      char link1[160];
      if (*dot == '.') {
	*dot++ = '\0';
	sprintf (link1, "%s_%ld.%s", path, buf.st_mtime, dot);
	*--dot = '.';
      }
      else
	sprintf (link1, "%s_%ld", path, buf.st_mtime);
      char cwd[80], old[160];
      getcwd (cwd, 80);
      sprintf (old, "%s/%s", cwd, path);
      symlink (old, link1);
      link = acat ("/", link1, NULL);
    }
    else
      link = strdup (name);
    free (path);
    return link;
  }

  #define YY_INPUT(buf,result,max_size)			      \
  {							      \
    int c = *(yyextra->input);                                \
    if (c == '\0') result = YY_NULL;                          \
    else { buf[0] = c; result = 1; yyextra->input++; }	      \
    if (c == '\n') yyextra->line++;			      \
  }
%}

ID  [a-zA-Z_0-9]
SP  [ \t]
ES  (\\([\'\"\?\\abfnrtv]|[0-7]{1,3}|x[a-fA-F0-9]+))
WS  [ \t\v\n\f]
SCALAR [a-zA-Z_0-9]+[.xyz]*

%%

{WS}*\/[*][*]{SP}* {
  // start of C documentation comment i.e. "/**"
  if (yyextra->ncodes > 0 && !yyextra->first)
    output_s ("\n~~~\n");
  output_c ('\n');
  yyextra->incode = 0;
  yyextra->indent = spacenb (yytext);
  yyextra->type = C;
}

{SP}*[*]\/{SP}* {
  // end of any C comment block i.e. "*/"
  if (yyextra->incode)
    // end of standard comment
    output_s (yytext);
  else {
    // end of documentation comment
    output_s ("\n");
    yyextra->ncodes++;
    yyextra->incode = 1;
    yyextra->first = 1;
  }
}

{WS}*\"\"\"{SP}* {
  if (yyextra->incode) {
    // start of Python documentation comment i.e. """
    if (yyextra->ncodes > 0 && !yyextra->first)
      output_s ("\n~~~\n");
    output_c ('\n');
    yyextra->incode = 0;
    yyextra->indent = spacenb (yytext);
    yyextra->type = Python;
  }
  else {
    // end of Python documentation comment
    output_s ("\n");
    yyextra->ncodes++;
    yyextra->incode = 1;
    yyextra->first = 1;
  }
}

^{WS}*%\{{WS}*$ {
  if (yyextra->incode) {
    // start of Octave documentation comment i.e. "%{"
    if (yyextra->ncodes > 0 && !yyextra->first)
      output_s ("\n~~~\n");
    output_c ('\n');
    yyextra->incode = 0;
    yyextra->indent = spacenb (yytext);
    yyextra->type = Octave;
  }
  else
    REJECT;
}

^{WS}*%\}{WS}*$ {
  if (!yyextra->incode) {
    // end of Octave documentation comment i.e. "%}"
    output_s ("\n");
    yyextra->ncodes++;
    yyextra->incode = 1;
    yyextra->first = 1;
  }
  else
    REJECT;
}

^#!\/bin\/bash{WS}*$ {
  if (yyextra->line > 2)
    REJECT;
}

^:<<'DOC'{WS}*$ {
  if (yyextra->incode) {
    // start of Bash documentation comment i.e. ":<<'DOC'"
    if (yyextra->ncodes > 0 && !yyextra->first)
      output_s ("\n~~~\n");
    output_c ('\n');
    yyextra->incode = 0;
    yyextra->indent = spacenb (yytext);
    yyextra->type = Bash;
  }
  else
    REJECT;
}

^DOC{WS}*$ {
  if (!yyextra->incode) {
    // end of Bash documentation comment i.e. "DOC"
    output_s ("\n");
    yyextra->ncodes++;
    yyextra->incode = 1;
    yyextra->first = 1;
  }
  else
    REJECT;
}

^{SP}*[\v\n\f] {
  // empty line
  check_error (yyextra);  
  if (!yyextra->incode || !yyextra->first)
    output_s (yytext);
}

^{SP}* {
  // spaces at the beginning of a line
  if (yyextra->incode) {
    if (yyextra->first) {
      output_s (yyextra->type);
      output_c ('\n');
      yyextra->first = 0;
    }
    output_s (yytext);
  }
  else {
    int ns = spacenb (yytext) - yyextra->indent;
    while (ns-- > 0)
      output_c (' ');
  }
}

^{SP}*~~~gnuplot.*$ {
  if (yyextra->incode)
    REJECT;
  yyextra->gnuplot = strdup (strstr (yytext, "gnuplot") + 7);
}

^{SP}*~~~pythonplot.*$ {
  if (yyextra->incode)
    REJECT;
  yyextra->gnuplot = strdup (strstr (yytext, "pythonplot") + 10);
}

set{SP}+output{SP}*['"][^'"]+['"] |
savefig{SP}*[(]{SP}*['"][^'"]+['"] {
  if (yyextra->gnuplot) {
    char * s = strchr (yytext, '\'');
    if (!s)
      s = strchr (yytext, '"');
    s++;
    yyextra->gnuplot_output = strdup (s);
    yyextra->gnuplot_output[strlen(s) - 1] = '\0';
  }
  else
    REJECT;
}

^{SP}*~~~{SP}*$ {
  if (!yyextra->incode && yyextra->gnuplot) {
    if (!yyextra->gnuplot_output) {
      char name[30];
      sprintf (name, "_plot%d.svg", yyextra->nplots++);
      yyextra->gnuplot_output = strdup (name);
    }
    output_s ("\n![");
    output_s (yyextra->gnuplot);
    output_s ("](");
    char * file = &yyextra->page[strlen(yyextra->page) + 1];
    char * name = acat (file, "/", yyextra->gnuplot_output, NULL);
    char * link = stamped_image (name, yyextra);
    output_s (link);
    free (name); free (link);
    output_s(")");
    free (yyextra->gnuplot);
    yyextra->gnuplot = NULL;
    free (yyextra->gnuplot_output);
    yyextra->gnuplot_output = NULL;
  }
  else
    output_s (yytext);
}

^{SP}*~~~bib{SP}*$ {
  // bibtex file
#if !STANDALONE
  // only for static HTML generation (otherwise it is handled by gitit plugin)
  REJECT;
#else // STANDALONE: use bibtex2html
  char command[256], file[80];
  strcpy (file, yyextra->page);
  strcat (file, ".bib2html");
  sprintf (command,
	   "bibtex2html -a -use-keys -nodoc -noheader -q | "
	   "sed -e 's|</table>.*|</table>|' -e '/<\\/table>/q' > %s", file);
  FILE * fp = popen (command, "w");
  if (fp == NULL) {
    perror (command);
    REJECT;
  }
  else {
    int b, res, * ab = &b;
    do {
      YY_INPUT (ab, res, 0);
      if (res) {
	if (b == '~') {
	  int n = 0;
	  while (res && n < 3 && b == '~') {
	    YY_INPUT (ab, res, 0);
	    n++;
	  }
	  if (n == 3)
	    break;
	  int i;
	  for (i = 0; i < n; i++)
	    fputc ('~', fp);
	  if (res)
	    fputc (b, fp);
	}
	else
	  fputc (b, fp);
      }
    } while (res);
    fflush (fp);
    pclose (fp);
    fp = fopen (file, "r");
    if (!fp)
      perror (file);
    else {
      output_s ("<div class=\"bibtex\">\n");
      int c;
      while ((c = fgetc (fp)) != EOF)
	output_c (c);
      fclose (fp);
      output_s ("</div>\n");
    }
  }
  remove (file);
#endif // STANDALONE
}
  
!\[[^\]]*\][(][^)]+\/[^)]+\.(png|gif|jpg|mp4|ogv)[)](\([^)]*\))? {
  if (yyextra->incode)
    REJECT;

  // Wiki links to images (only in page directory)
  char * end = strchr (yytext, ']');
  char * name = strchr (end, '(') + 1;

#if !STANDALONE
  char * file = &yyextra->page[strlen(yyextra->page) + 1];
  if (strncmp (file, name, strlen (file)))
    REJECT;
#endif
  
  char * options = strchr (name, ')');
  *options++ = '\0';

#if !STANDALONE
  char * link = stamped_image (name, yyextra);
#else
  char * link = strdup (name);
#endif

  *end = '\0';
  if (!strcmp(link + strlen(link) - 4, ".mp4") ||
      !strcmp(link + strlen(link) - 4, ".ogv")) {
    char * caption = strchr (yytext, '[') + 1;
    output_s ("<div class=\"figure\"><video ");
    if (options[0] == '(') {
      *strchr (++options, ')') = '\0';
      output_s (options);
      output_c (' ');
    }
    output_s ("controls><source src=\"");
    output_s (link);
    output_s ("\" type = \"video/");
    output_s (!strcmp(link + strlen(link) - 4, ".mp4") ? "mp4" : "ogg");
    output_s ("\">Your browser does not support the video tag.</video>");
    if (caption[0] != '\0') {
      output_s ("<p class=\"caption\">");
      output_s (caption);
      output_s ("</p>");
    }
    output_s ("</div>");
  }
  else {
    output_s (yytext);
    output_s ("](");
    output_s (link);
    output_c (')');
  }
  free (link);
}

\n {
  check_error (yyextra);
  if (!yyextra->gnuplot)
    output_s (yytext);
}

. {
  if (yyextra->incode && yyextra->first) {
    output_s (yyextra->type);
    output_c ('\n');
    yyextra->first = 0;
  }
  if (!yyextra->gnuplot)
    output_s (yytext);
}

({SP}?\"([^\"\\\n]|{ES})*\"{WS}*)+  {
  /* STRING_LITERAL */
  if (!yyextra->gnuplot)
    output_s (yytext);
}

%%

static void revert (char * src, char * bak)
{
  if (src) {
    if (bak) {
      rename (bak, src);
      free (bak);
    }
    else
      remove (src);
    free (src);
  }
}

static char * append (char * s, char * s1)
{
  if (s) {
    if (strstr (s, s1))
      return s;
    char * n = acat (s, "<br>", s1, NULL);
    free (s);
    return n;
  }
  return strdup (s1);
}

static int scan_errors (FILE * fp, char * file, struct MyScanner * scan,
			int output)
{
  char * header = acat (file, ".c:", NULL);
  char * line = NULL;
  size_t n = 0, nl = 0;
  while (getline (&line, &n, fp) > 0) {
    if (!strncmp (line, header, strlen (header))) {
      char * lineno = line; while (*lineno != ':') lineno++;
      char * type = ++lineno; while (*type != ':' && *type != '\0') type++;
      if (*type == ':') {
	*type++ = '\0';
	if (*type >= '0' && *type <= '9') {
	  while (*type != ':' && *type != '\0') type++;
	  *type++ = '\0';
	}
	while (strchr(" \t", *type)) type++;
	char * msg = type; while (*msg != ':' && *msg != '\0') msg++;
	if (*msg == ':') {
	  *msg++ = '\0';
	  if (!strcmp (type, "error") || !strcmp(type, "warning")) {
	    int line = atoi (lineno);
	    if (line > 0) {
	      if (line > scan->nerror) {
		scan->error = realloc (scan->error, line*sizeof (Error));
		int i;
		for (i = scan->nerror; i < line; i++)
		  scan->error[i].error = scan->error[i].warning = NULL;
		scan->nerror = line;
	      }
	      if (!strcmp (type, "error"))
		scan->error[line-1].error = 
		  append (scan->error[line-1].error, msg);
	      else
		scan->error[line-1].warning = 
		  append (scan->error[line-1].warning, msg);
	    }
	  }
	  free (line); line = NULL;
	}
	else if (type[0] != '\n')
	  line[strlen(line)] = ':';
	else {
	  free (line); line = NULL;
	}	  
      }
    }
    if (output && line) {
      if (nl < 20) {
	char * s = line;
	while (*s != '\0') {
	  if (*s == '`')
	    *s = '\'';
	  s++;
	}
	output_s1 (scan, line);
	output_s1 (scan, "<br>");
      }
      else if (nl == 20)
	output_s1 (scan, "...<br>");
      nl++;
    }
    free (line); line = NULL;
  }
  free (header);

  return scan->nerror;
}

static int run (char * text, char * dir, char * file, struct MyScanner * scan)
{
  // copy contents of text into file
  char * src = acat (dir, "/", file, ".c.page", NULL);
  FILE * fp = fopen (src, "r");
  char * bak = NULL;
  if (fp) {
    fclose (fp);
    bak = acat (dir, "/", file, ".c.bak", NULL);
    if (rename (src, bak)) {
      perror (bak);
      free (src);
      free (bak);
      return 1;
    }
  }

  // create directory in case it's not there
  char * mkpath = acat ("mkdir -p ", dir, NULL);
  system (mkpath);
  free (mkpath);

  // open new page
  fp = fopen (src, "w");
  if (!fp) {
    perror (src);
    free (src);
    free (bak);
    return 1;
  }
  fwrite (text, 1, strlen (text), fp);
  fclose (fp);

  // create link to Makefile in case it's not there
  mkpath = acat ("test -f ", dir, 
		 "/Makefile || ln -s `pwd`/sandbox/Makefile ", dir, "/Makefile",
		 NULL);
  system (mkpath);
  free (mkpath);

  // cd dir && make file.c.tags && make file.tst > file.log 2>&1
  char * command = acat ("cd ", dir, 
			 " && make ", file, ".c.tags",
			 " && make ", file, ".tst",
			 " > ", file, ".log 2>&1", NULL);
  int status = system (command);
  free (command);
  command = acat ("cd ", dir, " && cat ", file, ".log", NULL);
  system (command);
  free (command);
  if (status == -1 ||
      (WIFSIGNALED (status) && (WTERMSIG (status) == SIGINT || 
				WTERMSIG (status) == SIGQUIT)))
    status = 1;
  else
    status = WEXITSTATUS (status);

  if (status) {
    // scan for errors
    char * log = acat (dir, "/", file, ".log", NULL);
    fp = fopen (log, "r");
    if (!fp)
      perror (log);
    else {
      error_start (scan);
      scan_errors (fp, file, scan, 1);
      error_end (scan);
      fclose (fp);
    }
    free (log);
    revert (src, bak);
    return 1;
  }
  else {
    // scan for warnings
    char * log = acat (dir, "/", file, ".log", NULL);
    fp = fopen (log, "r");
    if (!fp)
      perror (log);
    else {
      scan_errors (fp, file, scan, 0);
      if (scan->nerror) {
	warning_start (scan);
	output_s1 (scan, "Please check the code for warnings.");
	error_end (scan);
      }
      fclose (fp);
    }
    free (log);    
  }

  revert (src, bak);
  return 0;
}

static void check_status (char * dir, char * file, struct MyScanner * scan)
{
  // check for failure
  char * fail = acat (dir, "/", file, "/fail", NULL);
  FILE * fp = fopen (fail, "r");
  free (fail);
  if (fp)
    fclose (fp);
  else {
    // are we still running?
    char * pids = acat (dir, "/", file, "/pid.*", NULL);
    glob_t globbuf;
    if (!glob(pids, GLOB_ERR, NULL, &globbuf))
      output_s1 (scan,
		 "<div class=message><div id=msg_logo>"
		 "<img src=/img/run.png></div>"
		 "<div id=msg_label>"
		 "This simulation is running."
		 " Please reload/preview the page to update results."
		 "</div></div>");
    globfree (&globbuf);
  }
}

static int belongs (const char * user, const char * group)
{
  FILE * fp = fopen (group, "r");
  if (fp) {
    char * line = NULL;
    size_t n = 0;
    while (getline (&line, &n, fp) > 0) {
      line[strlen(line) - 1] = '\0';
      if (!strcmp (user, line)) {
	free (line);
	fclose (fp);
	return 1;
      }
      free (line); line = NULL;
    }
    fclose (fp);
  }
  return 0;
}

static int is_preview (char * s, const char * page)
{
  // diff between buffer and stored page
  char * file = acat (page, ".page", NULL);
  FILE * fp = fopen (file, "r");
  free (file);
  if (!fp)
    return 1;
  int c = fgetc (fp);
  while (*s != '\0' && c != EOF) {
    if ((unsigned char) *s != c) {
      fclose (fp);
      return 1;
    }
    c = fgetc (fp); s++;
  }
  fclose (fp);
  return (*s != '\0' || c != EOF);
}

static int root_is (const char * path, const char * root)
{
  return !strncmp (path, root, strlen(root)) &&
    (path[strlen(root)] == '\0' || path[strlen(root)] == '/');
}

static int tail_is (const char * path, const char * tail)
{
  int plen = strlen (path);
  int tlen = strlen (tail);
  return plen >= tlen && !strcmp (&path[plen - tlen], tail);
}

static void usage (struct MyScanner * scan)
{
  char s[80];
  strcpy (s, scan->page);
  strcat (s, ".itags");
  FILE * fp = fopen (s, "r");
  if (!fp)
    perror (s);
  else {
    char type[10], title[80], file[80], line[10];
    int nf = 0, ne = 0, nt = 0;
    while (fscanf (fp, "%s %[^\t]\t%s %s", type, title, file, line) == 4) {
      nf++;
      if (root_is (file, "/src/test")) nt++;
      if (root_is (file, "/src/examples")) ne++;
    }
    if (nf > 0) {
      output_s1 (scan, "## Usage\n");
      rewind (fp);
      while (fscanf (fp, "%s %[^\t]\t%s %s", type, title, file, line) == 4)
	if (!root_is (file, "/src/test") &&
	    !root_is (file, "/src/examples")) {
	  output_s1 (scan, "* ["); output_s1 (scan, title); 
	  output_s1 (scan, "](");
	  output_s1 (scan, file); output_s1 (scan, ")\n");
	}
      if (ne > 0) {
	output_s1 (scan, "\n### Examples\n");
	rewind (fp);
	while (fscanf (fp, "%s %[^\t]\t%s %s", type, title, file, line) == 4)
	  if (root_is (file, "/src/examples")) {
	    output_s1 (scan, "* ["); output_s1 (scan, title); 
	    output_s1 (scan, "](");
	    output_s1 (scan, file); output_s1 (scan, ")\n");
	  }
      }
      if (nt > 0) {
	output_s1 (scan, "\n### Tests\n");
	rewind (fp);
	while (fscanf (fp, "%s %[^\t]\t%s %s", type, title, file, line) == 4)
	  if (root_is (file, "/src/test")) {
	    output_s1 (scan, "* ["); output_s1 (scan, title); 
	    output_s1 (scan, "](");
	    output_s1 (scan, file); output_s1 (scan, ")\n");
	  }
      }
    }
    fclose (fp);
  }
}

char * literate (char * s, char * page, char * user)
{
  yyscan_t scanner;
  struct MyScanner sdata;
  sdata.input = s;
  sdata.len = strlen(s);
  sdata.output = malloc (sdata.len + 1);
  sdata.i = 0;
  sdata.output[0] = '\0';
  sdata.incode = 1;
  sdata.ncodes = 0;
  sdata.first = 0;
  sdata.error = NULL;
  sdata.nerror = sdata.nplots = 0;
  sdata.line = 0;
  sdata.page = page;
  sdata.gnuplot = sdata.gnuplot_output = NULL;

  int developer = belongs (user, "../developers");
  int runner = developer || belongs (user, "../runners");

  if (root_is (page, "src") && !developer && is_preview (s, page)) {
    error_start (&sdata);
    output_s1 (&sdata, 
	       "Only developers are allowed to modify files in src/.<br>"
	       "You will not be able to save your changes.");
    error_end (&sdata);
  }

  if (tail_is (page, ".c") &&
      (root_is (page, "src/test") ||
       root_is (page, "src/examples") ||
       root_is (page, "sandbox"))) {
    char * file = &page[strlen(page) - 2];
    *file-- = '\0';
    while (file != page && *file != '/')
      file--;
    if (file != page)
      *file++ = '\0';
    if ((root_is (page, "src/test") && developer) ||
	(root_is (page, "src/examples") && developer) ||
	(root_is (page, "sandbox") && runner)) {
      if (!run (s, page, file, &sdata))
	check_status (page, file, &sdata);
    }
    else if (root_is (page, "sandbox")) {
      if (strcmp (user, "???")) {
        warning_start (&sdata);
        output_s1 (&sdata, 
    		   "You do not have the rights to run code in /sandbox.");
        error_end (&sdata);
      }
      check_status (page, file, &sdata);
    }
  }

  yylex_init_extra (&sdata, &scanner);
  FILE * fp = fopen ("/dev/null", "r");
  yyset_in (fp, scanner);
  yyset_out (stderr, scanner);
  yylex (scanner);
  yylex_destroy (scanner);
  fclose (fp);

  if (sdata.ncodes > 0 && sdata.incode && !sdata.first)
    output_s1 (&sdata, "\n~~~\n");

  if (sdata.ncodes > 0 && tail_is (page, ".h"))
    usage (&sdata);

#if 0
  fputc ('"', stderr);
  fputs (sdata.output, stderr);
  fputs ("\"\n", stderr);
#endif

  int i;
  for (i = 0; i < sdata.nerror; i++) {
    free (sdata.error[i].error);
    free (sdata.error[i].warning);
  }
  free (sdata.error);

  return sdata.output;
}

#if STANDALONE
int main (int argc, char * argv[])
{
  if (argc < 3) {
    fprintf (stderr, "usage: ./literate USER FILE\n");
    return 1;
  }
  char name[80];
  sprintf (name, "%s.page", argv[2]);
  FILE * f = fopen (name, "r");
  if (!f) {
    perror (name);
    return 1;
  }

  char basename[256] = "";
  strcat (basename, argv[2]);
  char * dirname = basename + strlen(basename) + 1;
  if (strcmp (basename + strlen(basename) - 2, ".c"))
    strcpy (dirname, "");
  else {
    strcpy (dirname, basename);
    dirname[strlen(dirname) - 2] = '\0';
  }
  
  fseek (f, 0, SEEK_END);
  long fsize = ftell(f);
  fseek (f, 0, SEEK_SET);
  char * string = malloc (fsize + 1);
  fread (string, fsize, 1, f);
  fclose (f);
  string[fsize] = '\0';
  char * output = literate (string, basename, argv[1]);
  puts (output);
  free (output);
  free (string);
  return 0;
}
#endif
