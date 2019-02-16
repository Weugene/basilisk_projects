%option reentrant noyywrap extra-type="struct MyScanner *"
%{
  #include <stdlib.h>

  typedef struct {
    char id[80], file[512], line[80];
  } Tag;

  struct MyScanner {
    char * input, * page;
    char * output;
    Tag * decl, * call, * incl;
    int ndecl, ncall, nincl;
    int i, len, intypedef, intypedefscope, scope;
    int incode;
  };

  typedef struct { char * s, * sub; } Subst;

  static Subst greek[] = {
    {"alpha","α"}, {"beta","β"}, {"gamma","γ"}, {"delta","δ"},
    {"epsilon","ε"}, {"zeta","ζ"}, {"eta","η"}, {"theta","θ"},
    {"iota","ι"}, {"kappa","κ"}, {"lambda","λ"}, {"mu","μ"},
    {"nu","ν"}, {"xi","ξ"}, {"omicron","ο"}, {"pi","π"},
    {"rho","ρ"}, {"sigma","σ"}, {"tau","τ"}, {"upsilon","υ"},
    {"phi","φ"}, {"chi","χ"}, {"psi","ψ"}, {"omega","ω"},
    {"Alpha","Α"}, {"Beta","Β"}, {"Gamma","Γ"}, {"Delta","Δ"},
    {"Epsilon","Ε"}, {"Zeta","Ζ"}, {"Eta","Η"}, {"Theta","Θ"},
    {"Iota","Ι"}, {"Kappa","Κ"}, {"Lambda","Λ"}, {"Mu","Μ"},
    {"Nu","Ν"}, {"Xi","Ξ"}, {"Omicron","Ο"}, {"Pi","Π"},
    {"Rho","Ρ"}, {"Sigma","Σ"}, {"Tau","Τ"}, {"Upsilon","Υ"},
    {"Phi","Φ"}, {"Chi","Χ"}, {"Psi","Ψ"}, {"Omega","Ω"},
    {NULL}
  };

  #define lookup_decl(x) lookup_tag(x, yyextra->decl, yyextra->ndecl);
  #define lookup_call(x) lookup_tag(x, yyextra->call, yyextra->ncall);
  static Tag * lookup_tag (const char * id, Tag * tags, int ntags) {
    int i;
    for (i = 0; i < ntags; i++)
      if (!strcmp (id, tags[i].id))
	return &tags[i];
    return NULL;
  }

  #define output_c(c) output_c1 (yyextra, c)
  static void output_c1 (struct MyScanner * scan, int c) {
    if (scan->i == scan->len) {
      scan->len += 100;
      scan->output = realloc (scan->output, scan->len + 1);
    }
    scan->output[scan->i++] = c;
    scan->output[scan->i] = '\0';
  }

  #define echo() output_s(yytext)
  #define output_s(s) output_s1 (yyextra, s)
  static void output_s1 (struct MyScanner * scan, char * s) {
    while (*s)
      output_c1 (scan, *s++);
  }

#if STANDALONE
  static char * baseurl = "", * ext = NULL;

  char * url (char * s)
  {
    static char s1[256];
    if (s[0] == '/') {
      // assert (strlen(baseurl) + strlen(s) + 1) < 256;
      strcat (strcpy (s1, baseurl), s);
    }
    else
      strcpy (s1, s);
    if (ext) {
      char page[strlen(s1) + strlen (".page") + 1];
      strcpy (page, s1);
      char * anchor = strchr (page, '#');
      if (anchor)
	*anchor = '\0';
      strcat (page, ".page");
      // fprintf (stderr, "testing '%s'\n", page);
      FILE * fp = fopen (page, "r");
      if (fp) {
	if (anchor) {
	  strcpy (page, s1);
	  *strchr (s1, '#') = '\0';
	  strcat (s1, ".html");
	  strcat (s1, anchor);
	}
	else
	  strcat (s1, ".html");
	fclose (fp);
      }
    }
    return s1;
  }
  
# define ESCAPE  ""
#else  // for gitit module
# define ESCAPE  "\\"
# define url(s) s
#endif
# define INCODE() (yyextra->incode)
  
  static int check_tag (char * text, struct MyScanner * scan) {
    Tag * t = lookup_tag (text, scan->call, scan->ncall);
    if (t) {
      // link to another page
      output_s1 (scan, ESCAPE "<a href=");
      if (!strcmp(t->file, "stdlib")) {
	output_s1 (scan, "http://man7.org/linux/man-pages/man3/");
	output_s1 (scan, text);
	output_s1 (scan, ".3.html");
      }
      else if (!strcmp(t->file, "basilisk")) {
	output_s1 (scan, url ("/Basilisk%20C"));
	output_c1 (scan, '#');
	output_s1 (scan, t->line);
      }
      else {
	output_s1 (scan, url (t->file));
	output_c1 (scan, '#');
	char s1[256];
#if STANDALONE
	if (t->file[0] == '/')
	  strcat (strcpy (s1, baseurl), t->file);
	else
#endif
	strcpy (s1, t->file);
	output_s1 (scan, t->line);
      }
      output_s1 (scan, ESCAPE ">");
      output_s1 (scan, text);
      output_s1 (scan, ESCAPE "</a" ESCAPE ">");
    }
    return (t != NULL);
  }

  static char * append_c (char * tmp, int c) {
    int len = strlen(tmp);
    tmp = realloc (tmp, len + 2);
    tmp[len] = c; tmp[len+1] = '\0';
    return tmp;
  }

  #define nonspace(s) { while (strchr(" \t\v\n\f", *s)) s++; }
  #define space(s) { while (!strchr(" \t\v\n\f", *s)) s++; }
  static void comment(yyscan_t scanner);

  #define YY_INPUT(buf,result,max_size)			      \
  {							      \
    int i = 0;						      \
    char * s = buf;					      \
    while (i < max_size && *(yyextra->input)) {		      \
      *s++ = *(yyextra->input++); i++;			      \
    }							      \
    result = i;						      \
  }
%}

SYMID  [a-zA-Z0-9]
ID  [a-zA-Z0-9_]
SP  [ \t]
ES  (\\([\'\"\?\\abfnrtv]|[0-7]{1,3}|x[a-fA-F0-9]+))
WS  [ \t\v\n\f]

%%

"<code class=\"sourceCode c\">"  {
  yyextra->incode = 1; echo();
}

"</code>"                        {
  yyextra->incode = 0; echo();
}

^{SP}*"@def"{SP} {
  yyextra->incode = 0; echo();
}

^{SP}*"@"{SP}*$ {
  yyextra->incode = 1; echo();
}

\&quot;({ID}|[-/])+\.h\&quot; |
">"#{SP}*include{SP}+\".*\""<" |
^{SP}*#{SP}*include{SP}+\".*\" {
  if (!INCODE())
    REJECT;
  // add HTML links in include headers
  char c = '"', * s = strchr(yytext, c);
  if (!s)
    c = ';', s = strchr(yytext, c);
  *s++ = '\0';
  char c1 = '"', * s1 = strchr(s, c1);
  if (!s1)
    c1 = '&', s1 = strchr(s, c1);
  *s1++ = '\0';
  // look for header in tags
  char * header = NULL;
  int i = 0;
  for (i = 0; i < yyextra->nincl && !header; i++)
    if (strstr(yyextra->incl[i].id, s))
      header = yyextra->incl[i].id;
  echo();
  output_c (c);
  if (header) {
    output_s (ESCAPE "<a href=");
    output_s (url (header));
    output_s (ESCAPE ">");
  }
  output_s (s);
  if (header)
    output_s (ESCAPE "</a" ESCAPE ">");
  else
    fprintf (stderr, 
	     "codeblock: warning: %s: tag for \"%s\" not found\n",
	     yyextra->page, s);
  output_c (c1);
  output_s (s1);
}

{ID}+{SP}*\( {
  if (!INCODE() || yyextra->scope > 0)
    REJECT;
  // keyword  anchor (function definition)
  // this regexp needs to match that in src/include.lex (for tags)
  char * s = strdup(yytext), * s1 = s;
  while (!strchr(" \t\v\n\f(", *s1)) s1++;
  int c = *s1;
  *s1++ = '\0';
  char * tmp = malloc (sizeof(char)); *tmp = '\0';
  int p = 0, para = 1, c1;
  while (para > p && (c1 = input(yyscanner)) != EOF) {
    tmp = append_c (tmp, c1);
    if (c1 == '(') para++;
    else if (c1 == ')') para--;
  }
  if (c1 != ')') {
    output_s (s);
  }
  else {
    while ((c1 = input(yyscanner)) != EOF) {
      tmp = append_c (tmp, c1);
      if (c1 == '{' || c1 == ';')
	break;
      if (!strchr(" \t\v\n\f", c1))
	break;
    }
    if (c1 != '{') {
      output_s (s);
    }
    else {
      yyextra->scope++;
      Tag * t = lookup_decl (s);
      if (t == NULL || !strstr (yyextra->page, t->file)) {
	if (t != NULL)
	  fprintf (stderr, "file: %s page: %s\n", t->file, yyextra->page);
#if 1
	fprintf (stderr, 
		 "codeblock: warning: %s: tag for '%s' not found\n",
		 yyextra->page, s);
#endif
	output_s (s);
      }
      else {
	char buf[80];
	sprintf (buf, ESCAPE "<a id=%s" ESCAPE ">%s" ESCAPE "</a" ESCAPE ">",
		 t->id, s);
	//	fprintf (stderr, "%s: %s\n", yyextra->page, buf);
	output_s (buf);
      }
    }
  }
  output_c (c);
  output_s (s1);
  output_s (tmp);
  free (tmp);
  free (s);
}

\'.\' {
  echo(); // quoted character
}

"\&#39;"."\&#39;" {
  echo(); // quoted character in HTML
}

"{" {
  echo();
  yyextra->scope++;
}

"}" {
  echo();
  yyextra->scope--;
  if (yyextra->scope < 0) {
#if STANDALONE
    puts (yyextra->output);
    fprintf (stderr, "warning: %s: error: mismatched '}'\n", yyextra->page);
#endif
    yyextra->scope = 0;
  }
}

"<img src=\""[^\"]+\" |
"<a href=\""[^\"]+\" {
#if STANDALONE
  // URL
  yytext[strlen(yytext) - 1] = '\0';
  char * u = strchr (yytext, '"');
  *u = '\0';
  output_s (yytext);
  output_c ('"');
  output_s (url (u + 1));
  output_c ('"');
#else
  REJECT;
#endif
}

typedef {
  if (!INCODE())
    REJECT;
  echo();
  yyextra->intypedef = 1; yyextra->intypedefscope = yyextra->scope;
}

{ID}+{WS}*; {
  if (!INCODE())
    REJECT;
  // keyword  anchor (typedef)
  // this regexp needs to match that in src/include.lex (for tags)
  if (yyextra->intypedef && yyextra->intypedefscope == yyextra->scope) {
    char * s = yytext; space(s); *s-- = '\0';
    if (*s == ';')
      *s = '\0';
    yyextra->intypedef = 0;
    s = yytext;
    Tag * t = lookup_decl (s);
    if (t == NULL || !strstr (yyextra->page, t->file)) {
      if (t != NULL)
	fprintf (stderr, "file: %s page: %s\n", t->file, yyextra->page);
      fprintf (stderr, 
	       "codeblock: warning: %s: tag for '%s' not found\n",
	       yyextra->page, s);
      output_s (s);
    }
    else {
      char buf[80];
      sprintf (buf, "" ESCAPE "<a id=%s" ESCAPE ">%s" ESCAPE "</a" ESCAPE ">",
	       t->id, s);
      output_s (buf);
    }
    output_c (';');
  }
  else
    REJECT;  
}

scalar|vector|tensor {
#if STANDALONE
  output_s ("<span class=\"dt\">");
  echo();
  output_s ("</span>");
#else
  echo();
#endif
}

{ID}+ {
  if (!INCODE())
    REJECT;
  // keyword links
  if (!check_tag (yytext, yyextra))
    REJECT;
}

{SYMID}+ {
  if (!INCODE())
    REJECT;
  // Greek letters substitution
  Subst * sub = greek;
  while (sub->s) {
    if (!strcmp (sub->s, yytext)) {
      output_s (sub->sub);
      break;
    }
    sub++;
  }
  if (!sub->s)
    // not a greek letter
    REJECT;
}

"/*" {
  if (!INCODE())
    REJECT;
  echo();
  comment(yyscanner);
}

"//".* {
#if STANDALONE
  REJECT;
#else
  echo(); /* consume //-comment */
#endif
}

"//"[^<]* {
#if !STANDALONE
  REJECT;
#else
  if (!INCODE())
    REJECT;
  echo(); /* consume //-comment */
#endif
}

{ID}+                                 echo();
.                                     echo();
\n                                    echo();
({SP}?\"([^\"\\\n]|{ES})*\"{WS}*)+  { echo(); /* STRING_LITERAL */ }

%%

static int getput(yyscan_t scanner)
{
  int c = input(scanner);
  output_c1 (yyget_extra(scanner), c);
  return c;
}

static void comment(yyscan_t scanner)
{
  int c;
  while ((c = getput(scanner)) > 0) {
    if (c == '*') {
      while ((c = getput(scanner)) == '*')
	;
      if (c == '/')
	return;
      if (c == 0)
	break;
    }
  }
  fprintf (stderr, "codeblock: warning: %s: unterminated comment\n", 
	   yyget_extra(scanner)->page);
}

static void read_tagfile (struct MyScanner * scan)
{
  char files[2][200] = {"", ""};
  char * s = getenv ("BASILISK");
  if (s)
    strcpy (files[0], s);
  strcat (files[0], "/external.tags");
  strcpy (files[1], scan->page);
  strcat (files[1], ".tags");
  int i;
  for (i = 0; i < 2; i++) {
    char * s = files[i];
    FILE * fp = fopen (s, "r");
    if (fp == NULL)
      perror (s);
    else {
      Tag t;
      char type[10];
      while (fscanf (fp, "%s %s %s %s", type, t.id, t.file, t.line) == 4) {
	if (!strcmp (type, "decl")) {
	  scan->ndecl++;
	  scan->decl = realloc (scan->decl, scan->ndecl*sizeof(Tag));
	  scan->decl[scan->ndecl-1] = t;
	} else if (!strcmp (type, "call")) {
	  scan->ncall++;
	  scan->call = realloc (scan->call, scan->ncall*sizeof(Tag));
	  scan->call[scan->ncall-1] = t;
	}
	else {
	  scan->nincl++;
	  scan->incl = realloc (scan->incl, scan->nincl*sizeof(Tag));
	  scan->incl[scan->nincl-1] = t;
	}
      }
      fclose (fp);
    }
  }
}

char * codeblock (char * s, char * page)
{
  yyscan_t scanner;
  struct MyScanner sdata = {};
  sdata.input = s;
  sdata.len = strlen(s);
  sdata.output = malloc (sdata.len + 1);
  sdata.i = sdata.scope = sdata.intypedef = 0;
  sdata.output[0] = '\0';
  sdata.page = page;
#if STANDALONE
  sdata.incode = 0;
#else
  sdata.incode = 1;
#endif
  read_tagfile (&sdata);

  yylex_init_extra (&sdata, &scanner);
  FILE * fp = fopen ("/dev/null", "r");
  yyset_in (fp, scanner);
  yyset_out (stderr, scanner);
  yyset_debug (1, scanner);
  yylex (scanner);
  yylex_destroy (scanner);
  free (sdata.decl);
  free (sdata.call);
  free (sdata.incl);
  fclose (fp);

  //  fputs (sdata.output, stderr);
  
  return sdata.output;
}

#if STANDALONE || TEST
static void usage () {
#if STANDALONE
  fprintf (stderr,
	   "usage: codeblock BASEURL FILE.[ch] [EXT] < FILE.[ch].html\n");
#else // TEST
  fprintf (stderr,
	   "usage: codeblock BASEURL FILE.[ch] < FILE.markdown\n");
#endif
  exit (1);
}

int main (int argc, char * argv[])
{
  if (argc < 3)
    usage();
  char * name = argv[2];

#if STANDALONE
  baseurl = argv[1];
  if (argc >= 4)
    ext = argv[3];
#endif
  
#define BUF_SIZE 1024
  char buffer[BUF_SIZE];
  size_t contentSize = 1; // includes NULL
  /* Preallocate space.  We could just allocate one char here, 
     but that wouldn't be efficient. */
  char * content = malloc (sizeof(char)*BUF_SIZE);
  if (content == NULL) {
    perror ("Failed to allocate content");
    exit (1);
  }
  content[0] = '\0'; // make null-terminated
  while (fgets (buffer, BUF_SIZE, stdin)) {
    char * old = content;
    contentSize += strlen(buffer);
    content = realloc(content, contentSize);
    if (content == NULL) {
      perror("Failed to reallocate content");
      free (old);
      exit (2);
    }
    strcat (content, buffer);
  }
  if (ferror (stdin)) {
    free (content);
    perror ("Error reading from stdin.");
    exit (3);
  }

  char * output = codeblock (content, name);
  puts (output);
  free (output);
  free (content);
  return 0;
}
#endif // STANDALONE
