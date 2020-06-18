%option noyywrap
%{
  #include <unistd.h>
  #include <string.h>
  #include <ctype.h>
  #include <sys/stat.h>
  #include <sys/types.h>
  #include <sys/wait.h>
  #include <assert.h>

  enum { scalar, vector, tensor, bid,
	 struct_type = -1, other_type = -2 };

  typedef struct { int i; char * name; } Scalar;
  typedef struct { int x, y, face; char * name; } Vector;
  typedef struct { Vector x, y; char * name; } Tensor;

  int dimension = 2, bghosts = 0;
  
  int debug = 0, catch = 0, cadna = 0, nolineno = 0, events = 0, progress = 0;
  char dir[] = ".qccXXXXXX";

  char * autolink = NULL;
  int autolinks = 0, source = 0;
  
  int nvar = 0, nconst = 0, nevents = 0;
  int line;
  int scope, para, inforeach, foreachscope, foreachpara, 
    inforeach_boundary, inforeach_face, nmaybeconst = 0;
  int invardecl, vartype, varsymmetric, varface, varvertex, varmaybeconst;
  char * varconst, * type = NULL;
  int inval, invalpara, indef;
  int brack, inarray, inarraypara, inarrayargs;
  int infine;
  int inreturn;
  int inmalloc;

  #define EVMAX 100
  int inevent, eventscope, eventpara;
  char eventarray[EVMAX], * eventfile[EVMAX], * eventid[EVMAX];
  char * eventfunc[EVMAX];
  char * eventarray_elems[EVMAX];
  int nexpr[EVMAX], eventline[EVMAX], eventparent[EVMAX], eventchild[EVMAX];
  int eventlast[EVMAX];

  typedef struct {
    void * p;
    int n, size;
  } Stack;

  void stack_push (Stack * s, void * p) {
    s->n++;
    s->p = realloc (s->p, s->n*s->size);
    char * dest = ((char *)s->p) + (s->n - 1)*s->size;
    memcpy (dest, p, s->size);
  }

  void * stack_pop (Stack * s) {
    if (!s->n)
      return NULL;
    return ((char *)s->p) + --s->n*s->size;
  }

  void * stack_top (Stack * s) {
    if (!s->n)
      return NULL;
    return ((char *)s->p) + (s->n - 1)*s->size;
  }

  typedef struct {
    int scope, para, line, dim;
    FILE * fp;
  } foreachdim_t;

  Stack foreachdim_stack = {NULL, 0, sizeof(foreachdim_t)};

  int inattr, attrscope, attrline;
  FILE * attrfp;

  int inmap, mapscope, mapline;
  FILE * mapfp;

  int mallocpara;

  #define REDUCTMAX 10
  char foreachs[80], * fname;
  FILE * foreachfp;
  char reduction[REDUCTMAX][4], reductvar[REDUCTMAX][80];
  int nreduct;

  int inboundary;
  char boundaryvar[80], boundarydir[80];
  FILE * boundaryfp = NULL, * boundaryheader = NULL;
  int boundarycomponent, nboundary = 0, nsetboundary = 0;
  int boundaryindex[80], periodic[80];

  int infunction, infunctiondecl, functionscope, functionpara, inmain;
  int infunction_line;
  FILE * functionfp = NULL;
  char * return_type = NULL;
  FILE * dopen (const char * fname, const char * mode);

  int infunctionproto;

  FILE * tracefp = NULL;
  int intrace, traceon;
  char * tracefunc = NULL;
  
  int foreach_line;
  enum { face_x = 0, face_y, face_z, face_xyz };
  int foreach_face_xyz, foreach_face_start;

  char ** args = NULL, ** argss = NULL;
  int nargs = 0, inarg;
  int incadna = 0, incadnaargs = -1, incadnanargs = 0, incadnaarg[80];
  #define doubletype (cadna ? "double_st" : "double")
  #define floattype (cadna ? "float_st" : "float")
  
  typedef struct { 
    char * v, * constant, * block;
    int type, args, scope, automatic, symmetric, face, vertex, maybeconst;
    int wasmaybeconst;
    int i[9];
    char * conditional;
  } var_t;
  var_t * _varstack = NULL; int varstack = -1, varstackmax = 0;
  var_t * varpush (const char * s, int type, int scope, int maybeconst) {
    var_t * v = NULL;
    if (s[0] != '\0') {
      char * f = malloc (strlen (s) + 1);
      strcpy (f, s);
      char * q = f;
      int na = 0;
      while ((q = strchr (q, '['))) {
	*q++ = '\0'; na++;
      }
      varstack++;
      if (varstack >= varstackmax) {
	varstackmax += 100;
	_varstack = realloc (_varstack, varstackmax*sizeof (var_t));
      }
      _varstack[varstack] = (var_t) { f, NULL, NULL, type, na, scope, 
				      0, 0, 0, 0, maybeconst, 0, {-1} };
      v = &(_varstack[varstack]);
    }
    return v;
  }
  var_t ** foreachconst = NULL;

  char * makelist (const char * input, int type);

  char * automatic_list (int scope, int conditional) {
    char * list = NULL;
    int i, size = 1;
    for (i = varstack; i >= 0 && _varstack[i].scope > scope; i--) {
      var_t var = _varstack[i];
      if (var.automatic && !var.constant && var.type >= 0 &&
	  (conditional || !var.conditional)) {
	size += strlen (var.v) + 2;
	if (list == NULL) {
	  list = malloc (size);
	  strcpy (list, "{");
	}
	else {
	  list = realloc (list, size);
	  strcat (list, ",");
	}
	strcat (list, var.v);
      }
    }
    return list;
  }

  void delete_automatic (int scope) {
    char * list = automatic_list (scope, 0);
    if (list) {
      strcat (list, "}");
      char * slist = makelist (list, scalar);
      if (debug)
	fprintf (stderr, "%s:%d: deleting %s\n", fname, line, list);
      fprintf (yyout, " delete (%s); ", slist);
      free (slist);
      free (list);
    }
    int i;
    for (i = varstack; i >= 0 && _varstack[i].scope > scope; i--) {
      var_t var = _varstack[i];
      if (var.automatic) {
	if (var.conditional) {
	  if (debug)
	    fprintf (stderr, "%s:%d: deleting conditional %s\n", 
		     fname, line, var.v);
	  char list[80];
	  sprintf (list, "{%s}", var.v);
	  char * slist = makelist (list, scalar);
	  fprintf (yyout, " { if (!%s) delete (%s); } ", 
		   var.conditional, slist);
	  free (slist);
	}
	else if (var.type < 0) {
	  if (debug)
	    fprintf (stderr, "%s:%d: freeing array %s\n", 
		     fname, line, var.v);
	  fprintf (yyout, " free (%s);", var.v);
	}
      }
    }
  }

  void varpop() {
    delete_automatic (scope);
    while (varstack >= 0 && _varstack[varstack].scope > scope) {
      if (debug)
	fprintf (stderr, "%s:%d: '%s' out of scope\n", fname, line, 
		 _varstack[varstack].v);
      free (_varstack[varstack].v);
      free (_varstack[varstack].conditional);
      free (_varstack[varstack--].constant);
    }
  }

  void varwasmaybeconst() {
    int i;
    var_t * var = _varstack;
    for (i = 0; i <= varstack; i++, var++)
      if (var->wasmaybeconst) {
	var->maybeconst = 1;
	var->wasmaybeconst = 0;
      }
  }
  
  var_t * vartop() {
    return varstack >= 0 ? &_varstack[varstack] : NULL;
  }

  typedef struct {
    int scope, para;
    char name[80];
  } foreach_child_t;
  
  Stack foreach_child_stack = {NULL, 0, sizeof(foreach_child_t)};
  
  char * inbegin = NULL, ** blocks = NULL;
  Stack inblock_stack = {NULL, 0, sizeof(foreach_child_t)};
  
  int identifier (int c) {
    return ((c >= 'a' && c <= 'z') || 
	    (c >= 'A' && c <= 'Z') || 
	    (c >= '0' && c <= '9'));
  }

  int component (char * s, int c) {
    while (*s != '\0') {
      if (*s == '.') {
	s++;
	while (strchr(" \t\v\n\f", *s)) s++;
	if (*s == c)
	  return 1;
      }
      s++;
    }
    return 0;
  }
  
  var_t * varlookup (char * s, int len) {
    int i;
    for (i = varstack; i >= 0; i--)
      if (strlen(_varstack[i].v) == len && 
	  !strncmp(s, _varstack[i].v, len)) {
	if (_varstack[i].type < 0) // not a field
	  return NULL;
	if (s[len-1] == '.' && _varstack[i].type != vector)
	  return NULL;
	return &_varstack[i];
      }
    return NULL;
  }

  int get_vartype (char * s) {
    int len = strlen (s), i;
    for (i = varstack; i >= 0; i--)
      if (strlen(_varstack[i].v) == len && 
	  !strncmp(s, _varstack[i].v, len)) {
	if (s[len-1] == '.' && _varstack[i].type != vector)
	  return -100;
	return _varstack[i].type;
      }
    return -100;
  }

  int rotate (FILE * fin, FILE * fout, int n, int dimension);

  void writefile (FILE * fp, int nrotate, int dim,
		  int line1, const char * condition) {
    if (condition)
      fprintf (yyout, " if (%s) {", condition);
    if (line1 >= 0)
      fprintf (yyout, "\n#line %d\n", line1);
    rotate (fp, yyout, nrotate, dim);
    if (condition)
      fputs (" } ", yyout);
  }

  void endforeachdim() {
    foreachdim_t * dim;
    while ((dim = stack_top(&foreachdim_stack)) &&
	   dim->scope == scope && dim->para == para) {
      if (debug)
	fprintf (stderr, "%s:%d: popping foreach_dim\n", fname, line);
      FILE * fp = yyout;
      yyout = dim->fp;
      if (dim->scope > 0)
	fputc ('{', yyout);
      int i;
      for (i = 0; i < dim->dim; i++)
	writefile (fp, i, dim->dim, dim->line, NULL);
      fclose (fp);
      if (dim->scope > 0)
	fputc ('}', yyout);
      stack_pop (&foreachdim_stack);
    }
  }

  void write_face (FILE * fp, int i, int n) {
    fprintf (yyout, " { int %cg = -1; VARIABLES; ", 'i' + i);
    char s[30];
    sprintf (s, "is_face_%c()", 'x' + i);
    writefile (fp, n, dimension, foreach_line, s);
    fputs (" }} ", yyout);
  }
  
  void foreachbody() {
    if (inforeach_face) {
      // foreach_face()
      fputs ("foreach_face_generic()", yyout);
      FILE * fp = dopen ("_foreach_body.h", "r");
      if (foreach_face_xyz == face_xyz) {
	int i;
	for (i = 0; i < dimension; i++)
	  write_face (fp, (i + foreach_face_start) % dimension, i);
      }
      else if (foreach_face_xyz == face_x)
	write_face (fp, 0, 0);
      else if (foreach_face_xyz == face_y)
	write_face (fp, 1, 0);
      else
	write_face (fp, 2, 0);
      fputs (" end_foreach_face_generic()\n", yyout);
      fprintf (yyout, "#line %d\n", line);
      fclose (fp);
    }
    else {
      // foreach()
      FILE * fp = dopen ("_foreach.h", "r");
      int c;
      while ((c = fgetc (fp)) != EOF)
	fputc (c, yyout);
      fclose (fp);
      fp = dopen ("_foreach_body.h", "r");
      while ((c = fgetc (fp)) != EOF)
	fputc (c, yyout);
      fputs (" }", yyout);
      fclose (fp);
    }
    fprintf (yyout, " end_%s();", foreachs);
  }

  char * macro_from_var (const char * name) {
    char * macro = strdup (name), * s = macro;
    while (*s != '\0') {
      if (*s == '.')
	*s = '_';
      s++;
    }
    return macro;
  }

  void maybeconst_macro (var_t * var, const char * name,
			 const char * id) {
    int found = 0, i;
    for (i = 0; i < nmaybeconst && !found; i++)
      found = !strcmp (foreachconst[i]->v, var->v);
    if (!found) {
      foreachconst = realloc (foreachconst, 
			      ++nmaybeconst*sizeof (var_t *));
      foreachconst[nmaybeconst - 1] = var;
      if (debug)
	fprintf (stderr, "%s:%d: '%s' may be const\n", 
		 fname, line, var->v);
    }
    char * macro = macro_from_var (id);
    fprintf (yyout, "%s_%s(%s", name, macro, id);
    free (macro);
  }

  // combinations for each constant field
  void maybeconst_combinations (int line, void (* body)(void)) {
    if (nmaybeconst) {
      int n = 1 << nmaybeconst, bits;
      for (bits = 0; bits < n; bits++) {
	fputs ("\nif (", yyout);
	int i;
	for (i = 0; i < nmaybeconst; i++) {
	  fprintf (yyout, "%sis_constant(%s%s)",
		   (bits & (1 << i)) ? "" : "!",
		   foreachconst[i]->v,
		   foreachconst[i]->type == vector ? ".x" : "");
	  if (i < nmaybeconst - 1)
	    fputs (" && ", yyout);
	}
	fputs (") {\n", yyout);
	for (i = 0; i < nmaybeconst; i++)
	  if (foreachconst[i]->type == scalar) {
	    if (bits & (1 << i))
	      fprintf (yyout,
		       "const %s _const_%s = _constant[%s.i -_NVARMAX];\n"
		       "NOT_UNUSED(_const_%s);\n"
		       "#undef val_%s\n"
		       "#define val_%s(a,i,j,k) _const_%s\n"
		       "#undef fine_%s\n"
		       "#define fine_%s(a,i,j,k) _const_%s\n"
		       "#undef coarse_%s\n"
		       "#define coarse_%s(a,i,j,k) _const_%s\n",
		       doubletype,
		       foreachconst[i]->v, foreachconst[i]->v, 
		       foreachconst[i]->v, foreachconst[i]->v, 
		       foreachconst[i]->v, foreachconst[i]->v,
		       foreachconst[i]->v, foreachconst[i]->v, 
		       foreachconst[i]->v, foreachconst[i]->v, 
		       foreachconst[i]->v, foreachconst[i]->v);
	    else
	      fprintf (yyout, 
		       "#undef val_%s\n"
		       "#define val_%s(a,i,j,k) val(a,i,j,k)\n"
		       "#undef fine_%s\n"
		       "#define fine_%s(a,i,j,k) fine(a,i,j,k)\n"
		       "#undef coarse_%s\n"
		       "#define coarse_%s(a,i,j,k) coarse(a,i,j,k)\n",
		       foreachconst[i]->v, foreachconst[i]->v,
		       foreachconst[i]->v, foreachconst[i]->v,
		       foreachconst[i]->v, foreachconst[i]->v);
	  }
	  else if (foreachconst[i]->type == vector) {
	    int c, j;
	    if (bits & (1 << i)) {
	      fprintf (yyout, "const struct { %s x", doubletype);
	      for (c = 'y', j = 1; j < dimension; c++, j++)
		fprintf (yyout, ", %c", c);
	      fprintf (yyout, "; } _const_%s = {_constant[%s.x.i -_NVARMAX]",
		       foreachconst[i]->v, foreachconst[i]->v);
	      for (c = 'y', j = 1; j < dimension; c++, j++)
		fprintf (yyout, ", _constant[%s.%c.i - _NVARMAX]",
			 foreachconst[i]->v, c);
	      fprintf (yyout, "};\n"
		       "NOT_UNUSED(_const_%s);\n",
		       foreachconst[i]->v);
	    }
	    for (c = 'x', j = 0; j < dimension; c++, j++)
	      if (bits & (1 << i))
		fprintf (yyout,
			 "#undef val_%s_%c\n"
			 "#define val_%s_%c(a,i,j,k) _const_%s.%c\n"
			 "#undef fine_%s_%c\n"
			 "#define fine_%s_%c(a,i,j,k) _const_%s.%c\n"
			 "#undef coarse_%s_%c\n"
			 "#define coarse_%s_%c(a,i,j,k) _const_%s.%c\n",
			 foreachconst[i]->v, c,
			 foreachconst[i]->v, c, foreachconst[i]->v, c,
			 foreachconst[i]->v, c,
			 foreachconst[i]->v, c, foreachconst[i]->v, c,
			 foreachconst[i]->v, c,
			 foreachconst[i]->v, c, foreachconst[i]->v, c);
	      else
		fprintf (yyout, 
			 "#undef val_%s_%c\n"
			 "#define val_%s_%c(a,i,j,k) val(a,i,j,k)\n"
			 "#undef fine_%s_%c\n"
			 "#define fine_%s_%c(a,i,j,k) fine(a,i,j,k)\n"
			 "#undef coarse_%s_%c\n"
			 "#define coarse_%s_%c(a,i,j,k) coarse(a,i,j,k)\n",
			 foreachconst[i]->v, c, foreachconst[i]->v, c,
			 foreachconst[i]->v, c, foreachconst[i]->v, c,
			 foreachconst[i]->v, c, foreachconst[i]->v, c);
	  }
	fprintf (yyout, "#line %d\n", line);
	body();
	fputs (" }", yyout);
      }
    }
    else
      body();
  }

  void endforeach () {
    fclose (yyout);
    yyout = foreachfp;
    maybeconst_combinations (foreach_line, foreachbody);
    if (nreduct > 0) {
      int i;
      for (i = 0; i < nreduct; i++) {
	if (strcmp (reduction[i], "+"))
	  fprintf (yyout,
		   "OMP(omp critical) if (_%s %s %s) %s = _%s;\n",
		   reductvar[i], strcmp(reduction[i], "min") ? ">" : "<",
		   reductvar[i], reductvar[i], reductvar[i]);
	else
	  fprintf (yyout,
		   "OMP(omp critical) %s += _%s;\n",
		   reductvar[i], reductvar[i]);
	fprintf (yyout,
		 "mpi_all_reduce_double (%s, %s);\n",
		 reductvar[i], 
		 !strcmp(reduction[i], "min") ? "MPI_MIN" : 
		 !strcmp(reduction[i], "max") ? "MPI_MAX" : 
		 "MPI_SUM");
      }
      fputs ("\n#undef OMP_PARALLEL\n"
	     "#define OMP_PARALLEL() OMP(omp parallel)\n"
	     "}\n", yyout);
      fprintf (yyout, "#line %d\n", line);
    }
    fputs (" }", yyout);
    inforeach = inforeach_boundary = inforeach_face = nreduct = 0;
  }

  void endtrace() {
    fprintf (yyout, " end_trace(\"%s\", \"%s\", %d); ",
	     tracefunc, nolineno ? "" : fname, nolineno ? 0 : line);
  }

  void endevent() {
    if (intrace)
      endtrace();
    fputs (" return 0; } ", yyout);
    intrace = inevent = 0;
    nevents++;
  }

  void endattr() {
    inattr = 0;
    fclose (yyout); yyout = attrfp;
    FILE * fp = dopen ("_attribute.h", "r");
    FILE * out = dopen ("_attributes.h", "a");
    fprintf (out, "\n#line %d \"%s\"\n", attrline, fname);
    int c;
    while ((c = fgetc (fp)) != EOF) {
      if (c == '\n')
	fputc (c, yyout);
      fputc (c, out);
    }
    fclose (fp);
    fclose (out);
  }

  void endmap() {
    inmap = 0;
    fclose (yyout); yyout = mapfp;
    FILE * fp = dopen ("_map.h", "r");
    FILE * out = dopen ("_maps.h", "a");
    fprintf (out, "\n#line %d \"%s\"\n", mapline, fname);
    int c;
    while ((c = fgetc (fp)) != EOF) {
      if (c == '\n')
	fputc (c, yyout);
      fputc (c, out);
    }
    fclose (fp);
    fclose (out);
  }

  void maps (int line) {
    fputc ('\n', yyout);
    FILE * fp = dopen ("_maps.h", "r");
    int c;
    while ((c = fgetc (fp)) != EOF)
      fputc (c, yyout);
    fclose (fp);
    fprintf (yyout, "#line %d \"%s\"\n", line, fname);
  }

  void infunction_declarations (int line1) {
    if (!infunctiondecl) {
      if (debug)
	fprintf (stderr, "%s:%d: function declarations\n", fname, line1);
      char * name[3] = {"ig", "jg", "kg"};
      int i;
      for (i = 0; i < dimension; i++)
	fprintf (yyout, " int %s = 0; NOT_UNUSED(%s);", name[i], name[i]);
      fputs (" POINT_VARIABLES; ", yyout);
      maps (line1);
      infunctiondecl = 1;
      assert (functionfp == NULL);
      functionfp = yyout;
      infunction_line = line1;
      nmaybeconst = 0;
      yyout = dopen ("_infunction.h", "w");
    }
  }

  void function_body() {
    FILE * fp = dopen ("_infunction.h", "r");
    int c;
    assert (fp);
    while ((c = fgetc (fp)) != EOF)
      fputc (c, yyout);
    fclose (fp);
  }

  int endfunction() {
    infunction = 0;
    if (debug)
      fprintf (stderr, "%s:%d: outfunction %p\n", fname, line, functionfp);
    if (functionfp) {
      fclose (yyout);
      yyout = functionfp;
      functionfp = NULL;
      maybeconst_combinations (infunction_line, function_body);
      return 1;
    }
    return 0;
  }

  char * boundaryrep (char * name) {
    if (!inboundary)
      return name;
    if (!strcmp(name, boundaryvar)) {
      static char s[] = "_s";
      return s;
    }
    if (!strcmp (boundarydir, "left") || !strcmp (boundarydir, "right")) {
      char * s;
      if ((s = strstr (name, ".n")))
	strcpy (s, ".x");
      else if ((s = strstr (name, ".t")))
	strcpy (s, ".y");
      else if ((s = strstr (name, ".r")))
	strcpy (s, ".z");
    }
    else if (!strcmp (boundarydir, "top") || !strcmp (boundarydir, "bottom")) {
      char * s;
      if ((s = strstr (name, ".n")))
	strcpy (s, ".y");
      else if ((s = strstr (name, ".t")))
	strcpy (s, dimension > 2 ? ".z" : ".x");
      else if ((s = strstr (name, ".r")))
	strcpy (s, ".x");
    }
    else if (!strcmp (boundarydir, "front") || !strcmp (boundarydir, "back")) {
      char * s;
      if ((s = strstr (name, ".n")))
	strcpy (s, ".z");
      else if ((s = strstr (name, ".t")))
	strcpy (s, ".x");
      else if ((s = strstr (name, ".r")))
	strcpy (s, ".y");
    }
    return name;
  }

  int boundary_component (char * boundaryvar) {
    foreachdim_t * dim = stack_top(&foreachdim_stack);
    char * s1 = boundaryvar;
    int component = 0;
    while (*s1 != '\0') {
      if (*s1 == '.') {
	s1++;
	if (!dim) // replace .n, .t, .r with .x, .y, .z
	  *s1 = (*s1 == 'n' ? 'x' : *s1 == 't' ? 'y' : 'z');
	component = *s1;
      }
      else
	s1++;
    }
    return component;
  }

  char * opposite (const char * dir) {
    return (!strcmp(dir, "right") ? "left" :
	    !strcmp(dir, "left") ?  "right" :
	    !strcmp(dir, "top") ?   "bottom" :
	    !strcmp(dir, "bottom") ? "top" :
	    !strcmp(dir, "front") ?  "back" :
	    !strcmp(dir, "back") ? "front" :
	    "error");
  }
  
  void boundary_staggering (FILE * fp) {
    char index[3] = {'i','j','k'};
    char dir[3] = {'x','y','z'};
    int i;
    for (i = 0; i < dimension; i++)
      fprintf (fp,
	       " int %cg = neighbor.%c - point.%c; "
	       " if (%cg == 0) %cg = _attribute[_s.i].d.%c; "
	       " NOT_UNUSED(%cg);",
	       index[i], index[i], index[i],
	       index[i], index[i], dir[i],
	       index[i]);
  }

  char * makelist (const char * input, int type) {
    char * text = malloc (strlen(input) - 1);
    strncpy (text, &input[1], strlen (input) - 1);
    text[strlen (input) - 2] = '\0';
    char * s = strtok (text, " \t\v\n\f,");
    int listtype = 100;
    while (s) {
      char * dot = strchr (s, '.');
      var_t * var = varlookup (s, strlen(s) - (dot ? strlen(dot) : 0));
      if (!var) { // not a scalar or vector list, give up
	free (text);
	return NULL;
      }
      int vtype = var->type;
      while (dot) {
	vtype--;
	dot = strchr (dot+1, '.');
      }
      if (vtype < listtype)
	listtype = vtype;
      s = strtok (NULL, " \t\v\n\f,)");
    }
    if (listtype < type) { // incompatible list types
      free (text);
      return NULL;
    }      
    if (type >= 0)
      listtype = type;
    char * list = calloc (strlen("((scalar []){") + 1, sizeof (char));
    switch (listtype) {
    case scalar: strcpy (list, "((scalar []){"); break;
    case vector: strcpy (list, "((vector []){"); break;
    case tensor: strcpy (list, "((tensor []){"); break;
    default: assert (0); break;
    }
    strncpy (text, &input[1], strlen (input) - 1);
    text[strlen (input) - 2] = '\0';
    s = strtok (text, " \t\v\n\f,");
    while (s) {
      char * dot = strchr (s, '.');
      var_t * var = varlookup (s, strlen(s) - (dot ? strlen(dot) : 0));
      int vtype = var->type, si = 0;
      while (dot) {
	si += (dot[1] - 'x')*(vtype == 2 ? dimension : 1);
	vtype--;
	dot = strchr (dot+1, '.');
      }
      char member[80] = "";
      if (scope > 0 || var->i[0] < 0) { // dynamic allocation
	switch (listtype) {
	case scalar: {
	  switch (vtype) {
	  case 0: sprintf (member, "%s,", s); break;
	  case 1: {
	    int i, c;
	    for (i = 0, c = 'x'; i < dimension; i++, c++) {
	      char m[30];
	      sprintf (m, "%s.%c,", s, c);
	      strcat (member, m);
	    }
	    break;
	  }
	  case 2: {
	    int i, c;
	    for (i = 0, c = 'x'; i < dimension; i++, c++) {
	      int j, d;
	      for (j = 0, d = 'x'; j < dimension; j++, d++) {
		char m[30];
		sprintf (m, "%s.%c.%c,", s, c, d);
		strcat (member, m);
	      }
	    }
	    break;
	  }
	  default: assert (0);
	  }
	  break;
	}
	case vector: {
	  switch (vtype) {
	  case 1: {
	    sprintf (member, "{%s.x", s);
	    int i, c;
	    for (i = 1, c = 'y'; i < dimension; i++, c++) {
	      char m[30];
	      sprintf (m, ",%s.%c", s, c);
	      strcat (member, m);
	    }
	    strcat (member, "},");
	    break;
	  }
	  case 2: {
	    int i, c;
	    for (i = 0, c = 'x'; i < dimension; i++, c++) {
	      strcat (member, "{");
	      int j, d;
	      for (j = 0, d = 'x'; j < dimension; j++, d++) {
		char m[30];
		sprintf (m, "%s.%c.%c", s, c, d);
		strcat (member, m);
		if (j < dimension - 1)
		  strcat(member, ",");
	      }
	      strcat (member, "},");
	    }
	    break;
	  }
	  default: assert (0);
	  }
	  break;	  
	}
	case tensor: {
	  switch (vtype) {
	  case 2: {
	    int i, c;
	    strcat (member, "{");
	    for (i = 0, c = 'x'; i < dimension; i++, c++) {
	      strcat (member, "{");
	      int j, d;
	      for (j = 0, d = 'x'; j < dimension; j++, d++) {
		char m[30];
		sprintf (m, "%s.%c.%c", s, c, d);
		strcat (member, m);
		if (j < dimension - 1)
		  strcat(member, ",");
	      }
	      strcat (member, "}");
	      if (i < dimension - 1)
		strcat(member, ",");
	    }
	    strcat (member, "},");
	    break;
	  }
	  default: assert (0);
	  }
	  break;	  
	}
	default: assert(0);
	}
      }
      else { // static allocation
	int n = vtype - listtype, i;
	char * constant = var->constant ? "_NVARMAX + " : "";
	char coord[80] = "";
	int nd = 1;
	for (i = 0; i < n; i++)
	  nd *= dimension;
	for (i = 0; i < nd; i++) {
	  switch (listtype) {
	  case scalar:
	    sprintf (coord, "{%s%d},", constant, var->i[si+i]); break;
	  case vector: {
	    sprintf (coord, "{{%s%d}", constant, var->i[si+dimension*i]);
	    int j;
	    for (j = 1; j < dimension; j++) {
	      char s[80];
	      sprintf (s, ",{%s%d}", constant, var->i[si+dimension*i+j]);
	      strcat (coord, s);
	    }
	    strcat (coord, "},");
	    break;
	  }
	  case tensor: {
	    strcat (coord, "{");
	    int j, k, l = dimension*dimension*i;
	    for (j = 0; j < dimension; j++) {
	      strcat (coord, "{");
	      for (k = 0; k < dimension; k++) {
		char m[30];
		sprintf (m, "{%s%d}", constant, var->i[l++]);
		strcat (coord, m);
		if (k < dimension - 1)
		  strcat (coord, ",");
	      }
	      strcat (coord, "}");
	      if (j < dimension - 1)
		strcat (coord, ",");
	    }
	    strcat (coord, "},");
	    break;
	  }
	  default: assert (0);
	  }
	  strcat (member, coord);
	}
      }
      list = realloc (list, (strlen(list) + strlen(member) + 1)*sizeof(char));
      list = strcat (list, member);
      s = strtok (NULL, " \t\v\n\f,)");
    }
    free (text);
    char end[80]; // size should be OK at least up to dimension == 3
    switch (listtype) {
    case scalar: strcpy (end, "{-1}})"); break;
    case vector: {
      strcpy (end, "{{-1}");
      int i;
      for (i = 1; i < dimension; i++)
	strcat (end, ",{-1}");
      strcat (end, "}})");
      break;
    }
    case tensor: {
      strcpy (end, "{");
      int i, j;
      for (i = 0; i < dimension; i++) {
	strcat (end, "{{-1}");
	for (j = 1; j < dimension; j++)
	  strcat (end, ",{-1}");
	strcat (end, "}");
	if (i < dimension - 1)
	  strcat (end, ",");
      }
      strcat (end, "}})");
      break;
    }
    default: assert (0);
    }
    list = realloc (list, (strlen(list) + strlen(end) + 1)*sizeof(char));
    strcat (list, end);
    return list;
  }

  void new_field (var_t * var, char * block) {
    if (scope > 0) {
      // automatic variables
      fprintf (yyout, " new_%s%s%s%s(\"%s\"",
	       block ? "block_" : "",
	       var->constant ? "const_" : "",
	       var->symmetric ? "symmetric_" : 
	       var->face ? "face_" : 
	       var->vertex ? "vertex_" : 
	       "",
	       var->type == scalar ? "scalar" : 
	       var->type == vector ? "vector" : 
	       var->type == tensor ? "tensor" :
	       "internal_error",
	       var->v);
      if (block) {
	if (var->type == scalar)
	  fprintf (yyout, ", \"\"");
	fprintf (yyout, ", %s", block);
      }
      else if (var->constant) {
	if (var->type == scalar) {
	  fprintf (yyout, ", %d, %s", nconst, var->constant);
	  nconst++;
	}
	else if (var->type == vector) {
	  fprintf (yyout, ", %d, (%s [])%s", 
		   nconst, doubletype, var->constant);
	  nconst += dimension;
	}
	else
	  assert (0);
      }
      fputc (')', yyout);
    }
    else {
      // global constants and variables
      int * n = var->constant ? &nconst : &nvar;
      char * constant = var->constant ? "_NVARMAX + " : "";
      if (var->type == scalar) {
	fprintf (yyout, " {%s%d}", constant, (*n));
	var->i[0] = (*n)++;
      }
      else if (var->type == vector) {
	fputs (" {", yyout);
	int i;
	for (i = 0; i < dimension; i++) {
	  var->i[i] = (*n)++;
	  fprintf (yyout, "{%s%d}", constant, var->i[i]);
	  if (i < dimension - 1)
	    fputc (',', yyout);
	}
	fputc ('}', yyout);
      }
      else if (var->type == tensor) {
	int i, j, k = 0;
	for (i = 0; i < dimension; i++)
	  for (j = 0; j < dimension; j++) {
	    if (!var->symmetric || i <= j)
	      var->i[k] = (*n)++;
	    else
	      var->i[k] = var->i[j*dimension + i];
	    k++;
	  }
	fputs (" {", yyout);
	k = 0;
	for (i = 0; i < dimension; i++) {
	  fputc ('{', yyout);
	  for (j = 0; j < dimension; j++) {
	    fprintf (yyout, "{%s%d}", constant, var->i[k++]);
	    if (j < dimension - 1)
	      fputc (',', yyout);
	  }
	  fputc ('}', yyout);
	  if (i < dimension - 1)
	    fputc (',', yyout);
	}
	fputc ('}', yyout);
      }
      else
	assert (0);
    }
    if (!var->constant && var->maybeconst) {
      if (debug)
	fprintf (stderr, "%s:%d: '%s' cannot be a constant in this scope\n", 
		 fname, line, var->v);
      var->maybeconst = 0;
      var->wasmaybeconst = 1;
    }
  }

#define IDENTIFIER(c) (((c) >= 'A' && (c) <= 'Z') || \
                       ((c) >= 'a' && (c) <= 'z') || \
		       ((c) == '_'))
#define SEPARATOR(c)  (!IDENTIFIER(c) && c != '.')
  
  void cadna_echo (char * text) {
    if (cadna) {
      static char * rep[] = {"double", "float"};
      int i;
      for (i = 0; i < 2; i++) {
	char * s = strstr (text, rep[i]);
	if (s) {
	  char * end = s + strlen (rep[i]);
	  if ((s == text || SEPARATOR (*(s-1))) &&
	      (*end == '\0' || SEPARATOR (*end))) {
	    char c = *s;
	    *s = '\0';
	    cadna_echo (text);
	    *s = c;
	    fprintf (yyout, "%s_st", rep[i]);
	    cadna_echo (end);
	    return;
	  }
	}
      }
    }
    fputs (text, yyout);
  }

#define cadna_type(type) (!cadna ? type :			   \
			  !strcmp (type, "double") ? "double_st" : \
			  !strcmp (type, "float") ? "float_st" :   \
			  type)

  void declaration (char * var, char * text, int npointers) {
    var_t * v;
    if (vartype >= 0 && strchr (text, '[') && strchr (text, ']')) {
      // automatic
      *strchr (text, '[') = '\0';
      v = varpush (var, vartype, scope, 0);
      v->automatic = 1;
      v->symmetric = varsymmetric;
      v->face = varface;    
      v->vertex = varvertex;    
      if (varconst)
	v->constant = strdup (varconst);
      fputs (text, yyout);
      fputc ('=', yyout);
      new_field (v, NULL);
    }
    else {
      char * start = strchr (text, '['), * end = strchr (text, ']');
      if (start) {
	if (!end)
	  brack++;
	*start = '\0';
      }
      v = varpush (var, vartype, scope, varmaybeconst);
      v->symmetric = varsymmetric;
      v->face = varface;    
      v->vertex = varvertex;
      if (scope == 0 && v->type == bid)
	// global boundary id
	v->i[0] = 1;
      if (start) {
	if (cadna &&
	    scope > 0 && npointers == 0 && !end && type &&
	    (!strcmp (type, "double") ||
	     !strcmp (type, "float") ||
	     !strcmp (type, "coord") ||
	     !strcmp (type, "stats") ||
	     !strcmp (type, "norm"))) {
	  /* clang does not like variable-size arrays of "non-POD
	     types" e.g double_st or float_st of CADNA
	     we turn them into malloc'd arrays */
	  if (yytext[0] == ',')
	    fputc (',', yyout);
	  else
	    fputs (cadna_type (type), yyout);
	  fprintf (yyout, " *%s=(%s*)malloc(sizeof(%s)*(",
		   var, cadna_type (type), cadna_type (type));
	  inmalloc = brack;
	  v->automatic = 1;
	  if (debug)
	    fprintf (stderr, "%s:%d: malloc'd %s\n", fname, line, var);
	  return;
	}
	*start = '[';
      }
      cadna_echo (text);
    }
    if (debug)
      fprintf (stderr, "%s:%d: declaration: %s type: %d face: %d brack: %d "
	       "symmetric: %d\n",
	       fname, line, var, vartype, varface, brack, varsymmetric);
  }

  static int homogeneize (FILE * in, FILE * fp)
  {
    /* replace neumann/dirichlet(...) with neumann/dirichlet(0) */
    char * cond[2] = {"dirichlet", "neumann"};
    char str[1024];
    while (fgets (str, 1024, in)) {
      char * s = str;
      while (*s != '\0') {
	int i;
	for (i = 0; i < 2; i++)
	  if (!strncmp (s, cond[i], strlen (cond[i]))) {
	    fputs (cond[i], fp);
	    s += strlen (cond[i]);
	    while (strchr (" \t\v\n\f", *s))
	      s++;
	    if (*s == '(') {
	      fputs ("_homogeneous()", fp);
	      s++;
	      int para = 1;
	      while (para > 0 && *s != '\0') {
		if (*s == '(') para++;
		if (*s == ')') para--;
		s++;
	      }
	      if (para > 0)
		return 1;
	    }
	    break;
	  }
	fputc (*s++, fp);
      }
    }
    return 0;
  }
  
  void boundary_body (void) {
    FILE * fp = dopen ("_inboundary.h", "r");
    int c;
    while ((c = fgetc (fp)) != EOF)
      fputc (c, yyout);
    fclose (fp);
  }

  int yyerror(const char * s);

  void boundary_homogeneous_body (void) {
    FILE * fp = dopen ("_inboundary.h", "r");
    if (homogeneize (fp, yyout))
      yyerror ("expecting ')'");
    fclose (fp);
  }
  
#define nonspace(s) { while (strchr(" \t\v\n\f", *s)) s++; }
#define space(s) { while (!strchr(" \t\v\n\f", *s)) s++; }

#define YY_INPUT(buf,result,max_size)			      \
  {							      \
    int c = fgetc(yyin);				      \
    if (c == '\n') { line++; }				      \
    result = (c == EOF) ? YY_NULL : (buf[0] = c, 1);	      \
  }

  int getput(void);
  int comment(void);
%}

D			[0-9]
L			[a-zA-Z_]
H			[a-fA-F0-9]
E			[Ee][+-]?{D}+
FS			(f|F|l|L)
IS			(u|U|l|L)*

ID                      {L}({L}|{D})*
SP                      [ \t]
ES                      (\\([\'\"\?\\abfnrtv]|[0-7]{1,3}|x[a-fA-F0-9]+))
WS                      [ \t\v\n\f]
SCALAR                  {ID}({WS}*[.]{WS}*[xyzntr])*
CONSTANT                ({D}|{WS}|[<>+\-=*/^%&|!()~?:])+
TYPE                    [\*]*{WS}*{ID}({WS}*\[{WS}*{CONSTANT}?{WS}*\]|{WS}*\[)?

%%

"("                ECHO; para++;

")" {
  para--; if (para == invalpara) inval = 0;
  if (para < 0)
    return yyerror ("mismatched ')'");
  if (para == mallocpara - 1) {
    fprintf (yyout, ",__func__,__FILE__,%s", nolineno ? "0" : "__LINE__");
    mallocpara = 0;
  }
  if (inforeach == 1 && scope == foreachscope && para == foreachpara) {
    ECHO;
    fclose (yyout);
    yyout = foreachfp;
    if (nreduct > 0) {
      fputs ("\n#undef OMP_PARALLEL\n"
	     "#define OMP_PARALLEL()\n"
	     "OMP(omp parallel) {\n", yyout);
      int i;
      for (i = 0; i < nreduct; i++)
	fprintf (yyout, "%s _%s = %s; ",
		 doubletype, reductvar[i], reductvar[i]);
      fprintf (yyout, "\n#line %d\n", foreach_line);
    }
    yyout = dopen ("_foreach_body.h", "w");
    fputs ("{\n", yyout);
    maps (line);
    if (inforeach_face) {
      foreach_face_xyz = face_xyz;
      foreach_face_start = 0;
      FILE * fp = dopen ("_foreach.h", "r");
      int c, inbody = 0;
      while ((c = fgetc (fp)) != EOF) {
	if (c == '(')
	  inbody = 1;
	else if (c == ')')
	  break;
	else if (inbody) {
	  if (c == 'x' || c == 'y' || c == 'z') {
	    if (foreach_face_xyz == face_xyz)
	      foreach_face_xyz = face_x + c - 'x';
	  }
	  else if (c == ',') {
	    foreach_face_start = foreach_face_xyz;
	    foreach_face_xyz = face_xyz;
	    break;
	  }
	  else if (!strchr (" \t\v\n\f", c))
	    break;
	}
      }
      fclose (fp);
      if (foreach_face_xyz != face_xyz &&
	  ((dimension < 2 && foreach_face_xyz > face_x) ||
	   (dimension < 3 && foreach_face_xyz > face_y)))
	return yyerror ("face dimension is too high");
    }
    inforeach = 2;
  }
  else if (inevent > 0 && inevent < 4 && 
	   scope == eventscope && eventpara == para + 1) {
    if (!eventarray[nevents])
      fprintf (yyout, "); "
	       "  *ip = i; *tp = t; "
	       "  return ret; "
	       "} ");
    fprintf (yyout, 
	     "static int %s (const int i, const %s t, Event * _ev) { "
	     "trace (\"%s\", \"%s\", %d); ",
	     eventfunc[nevents], doubletype,
	     eventfunc[nevents], nolineno ? "" : fname, nolineno ? 0 : line);
    free (tracefunc);
    tracefunc = strdup (eventfunc[nevents]);
    intrace = 2; traceon = 0;
    assert (nevents < EVMAX);
    nexpr[nevents] = inevent;
    free (return_type);
    return_type = strdup ("int");
    inevent = 4;
  }
  else if (inarg == para + 1) {
    fputs ("})", yyout);
    inarg = 0;
  }
  else if (infine == para + 1) {
    int i;
    for (i = inarrayargs; i < 3; i++)
      fputs (",0", yyout);
    ECHO;
    infine = 0;
  }
  else
    ECHO;
  foreach_child_t * block = stack_top (&inblock_stack);
  if (block && block->scope == scope && block->para == para)
    fputc (';', yyout);
  if (incadna == para + 1) {
    incadna = 0;
    if (incadnanargs && incadnaargs == incadnaarg[incadnanargs-1] + 1)
      fputc (')', yyout);
  }
}

"{" {
  foreachdim_t * dim = stack_top(&foreachdim_stack);
  if ((!dim || dim->scope) || scope != 0 || infunctionproto)
    ECHO;
  if (infunction && functionpara == 1 && scope == functionscope)
    infunction_declarations (line - 1);
  if (inmain == 1 && scope == 0) {
    fputs (" _init_solver();", yyout);
    inmain = 2;
  }
  if (intrace == 1 && scope == 0) {
    fprintf (yyout, " trace (\"%s\", \"%s\", %d);",
	     tracefunc, nolineno ? "" : fname, nolineno ? 0 : line - 1);
    intrace = 2;
  }
  scope++;
}

"}" {
  scope--;
  if (scope < 0)
    return yyerror ("mismatched '}'");
  if (inmain == 2 && scope == 0) {
    fputs (" free_solver(); ",  yyout);
    inmain = 0;
  }
  varpop();
  varwasmaybeconst();
  if (infunction && scope <= functionscope)
    endfunction();
  foreachdim_t * dim = stack_top(&foreachdim_stack);
  if (dim && dim->scope == scope && dim->para == para) {
    if (scope != 0 || infunctionproto)
      ECHO;
    endforeachdim();
  }
  else if ((!inattr || scope != attrscope) &&
	   (!inmap || scope != mapscope)) {
    if (intrace && scope == 0) {
      endtrace ();
      intrace = 0;
    }
    ECHO;
  }
  foreach_child_t * i;
  while ((i = stack_top(&foreach_child_stack)) &&
	 i->scope == scope && i->para == para) {
    fprintf (yyout, " end_%s(); }", i->name);
    stack_pop (&foreach_child_stack);
  }
  while ((i = stack_top(&inblock_stack)) &&
	 i->scope == scope && i->para == para) {
    fprintf (yyout, " end_%s(); }", i->name);
    stack_pop (&inblock_stack);
  }
  if (inforeach && scope == foreachscope)
    endforeach ();
  else if (inevent && scope == eventscope)
    endevent();
  if (inattr && scope == attrscope)
    endattr ();
  if (inmap && scope == mapscope)
    endmap ();
  if (scope == 0)
    infunctionproto = 0;
}

(foreach_child|foreach_child_direction|foreach_neighbor|foreach_block) {
  if (indef)
    REJECT;
  foreach_child_t child = {scope, para};
  fputs (" { ", yyout);
  if (!strcmp(yytext, "foreach_block")) {
    if (inforeach) {
      fprintf (yyout, "foreach_block_inner");
      strcpy (child.name, "foreach_block_inner");
    }
    else {
      ECHO;
      strcpy (child.name, yytext);
    }      
  }
  else {
    ECHO;
    strcpy (child.name, yytext);
  }
  stack_push (&foreach_child_stack, &child);
}

break {
  if (foreach_child_stack.n)
    fprintf (yyout, "%s_break()",
	     ((foreach_child_t *)stack_top(&foreach_child_stack))->name);
  else
    ECHO;
}

foreach{ID}? {
  strcpy (foreachs, yytext);
  if (indef)
    REJECT;
  if (inforeach) {
    fprintf (stderr, "%s:%d: error: foreach* loops cannot be nested\n",
	     fname, line);
    return 1;    
  }
  fputs (" { ", yyout);
  inforeach = 1; foreachscope = scope; foreachpara = para;
  free (foreachconst); 
  foreachconst = NULL;
  nmaybeconst = 0;
  nreduct = 0;
  foreach_line = line;
  foreachfp = yyout;
  yyout = dopen ("_foreach.h", "w");
  inforeach_boundary = (!strncmp(foreachs, "foreach_boundary", 16));
  inforeach_face = (!strcmp(foreachs, "foreach_face"));
  ECHO;
}

end_foreach{ID}?{SP}*"()" {
  if (strncmp(&yytext[4], foreachs, strlen (foreachs))) {
    fprintf (stderr, 
	     "%s:%d: error: "
	     "%s() loop ended with %s\n", 
	     fname, line, foreachs, yytext);
    return 1;
  }
  if (indef)
    ECHO;
}

attribute{WS}+"{" {
  inattr = 1; attrscope = scope++; attrline = line;
  attrfp = yyout;
  yyout = dopen ("_attribute.h", "w");
}

(malloc|realloc|calloc|free|strdup) {
  fputc ('p', yyout);
  ECHO;
  mallocpara = para + 1;
}

map{WS}+"{" {
  inmap = 1; mapscope = scope++; mapline = line;
  mapfp = yyout;
  yyout = dopen ("_map.h", "w");
}

, {
  if (inarray && inarray == brack && inarraypara == para) {
    inarrayargs++;
    REJECT;
  }
  else if (infine == para) {
    inarrayargs++;
    REJECT;
  }
  else if (incadna == para && incadnanargs) {
    int i, found = 0;
    for (i = 0; i < incadnanargs && !found; i++) {
      if (incadnaargs == incadnaarg[i]) {
	fputs (",strp(", yyout);
	found = 1;
      }
      else if (incadnaargs == incadnaarg[i] + 1) {
	if (i < incadnanargs - 1 && incadnaargs == incadnaarg[i+1])
	  fputs ("),strp(", yyout);
	else
	  fputs ("),", yyout);
	found = 1;
      }
    }
    if (!found)
      ECHO;
    incadnaargs++;
  }
  else
    REJECT;
}

; {
  int insthg = 0;
  if (scope == 0)
    infunctionproto = 0;
  foreachdim_t * dim = stack_top(&foreachdim_stack);
  if (dim && dim->scope == scope && dim->para == para) {
    ECHO; insthg = 1;
    endforeachdim ();
  }
  foreach_child_t * i;
  while ((i = stack_top(&foreach_child_stack)) &&
	 i->scope == scope && i->para == para) {
    if (!insthg) {
      ECHO; insthg = 1;
    }
    fprintf (yyout, " end_%s(); }", i->name);
    stack_pop (&foreach_child_stack);
  }
  while ((i = stack_top(&inblock_stack)) &&
	 i->scope == scope && i->para == para) {
    if (!insthg) {
      ECHO; insthg = 1;
    }
    fprintf (yyout, " end_%s(); }", i->name);
    stack_pop (&inblock_stack);
  }
  if (infunction && scope == functionscope) {
    if (scope > 0 && !infunctiondecl) {
      fputs ("; ", yyout); insthg = 1;
      infunction_declarations (line);
    }
    else if (!infunctiondecl)
      endfunction();
  }
  if (inboundary) {
    ECHO;
    fclose (yyout);

    yyout = boundaryheader;
    maybeconst_combinations (line, boundary_body);
    fputs (" return 0.; } ", yyout);
    fprintf (yyout, 
	     "static %s _boundary%d_homogeneous "
	     "(Point point, Point neighbor, scalar _s, void * data) {",
	     doubletype, nboundary);
    boundary_staggering (yyout);
    fputs (" POINT_VARIABLES; ", yyout);
    maps (line - 1);
    maybeconst_combinations (line, boundary_homogeneous_body);
    fputs (" return 0.; }\n", yyout);

    yyout = boundaryfp;
    if (scope == 0) {
      /* file scope */
      fprintf (yyout, 
	       "static void _set_boundary%d (void) { ",
	       nboundary);
      boundaryindex[nsetboundary++] = nboundary;
    }
    /* function/file scope */
    fprintf (yyout,
	     "_attribute[%s.i].boundary[%s] = _boundary%d; "
	     "_attribute[%s.i].boundary_homogeneous[%s] = "
	     "_boundary%d_homogeneous;",
	     boundaryvar, boundarydir, nboundary,
	     boundaryvar, boundarydir, nboundary);
    if (scope == 0)
      /* file scope */
      fputs (" } ", yyout);
    periodic[nboundary++] = 0;
    inboundary = inforeach_boundary = inforeach = 0;
  }
  else if (inforeach && scope == foreachscope && para == foreachpara) {
    ECHO;
    endforeach ();
  }
  else if (inevent > 0 && inevent < 3 && para == eventpara) {
    if (eventarray[nevents])
      return yyerror ("cannot mix arrays and expressions");
    fprintf (yyout, "); "
	     "  *ip = i; *tp = t; "
	     "  return ret; "
	     "} "
	     "static int %s_expr%d (int * ip, %s * tp, Event * _ev) { "
	     "  int i = *ip; %s t = *tp; "
	     "  int ret = (", eventfunc[nevents], inevent++,
	     doubletype, doubletype);
  }
  else if (inevent == 4 && scope == eventscope && para == eventpara - 1) {
    ECHO;
    endevent ();
  }
  else if (inreturn) {
    fputs (";", yyout);
    delete_automatic (0);
    if (intrace)
      endtrace();
    fputs (return_type ? " return _ret; }" : " return; }", yyout);
    inreturn = 0;
  }
  else if (!insthg)
    ECHO;
  varpop();
  invardecl = varmaybeconst = 0;
}

{ID}{WS}*[=]{WS}*new{WS}+(symmetric|face|vertex){0,1}*{WS}*(scalar|vector|tensor)({WS}*\[[^\]]*\]){0,1} |
[=]{WS}*new{WS}+(symmetric|face|vertex){0,1}{WS}*(scalar|vector|tensor)({WS}*\[[^\]]*\]){0,1} {
  char * type = strchr (yytext, '=');
  type = strstr (type, "new"); space(type); nonspace(type);
  char * symmetric = strstr (type, "symmetric");
  if (symmetric) {
    space(type); nonspace(type);
  }
  char * face = strstr (type, "face");
  if (face) {
    space(type); nonspace(type);
  }
  char * vertex = strstr (type, "vertex");
  if (vertex) {
    space(type); nonspace(type);
  }
  char * block = strchr (type, '[');
  if (block) {
    char * s = type;
    while (s != block && !strchr(" \t\v\n\f", *s)) s++;
    if (strchr(" \t\v\n\f", *s)) *s = '\0';
    *block++ = '\0';
    *strchr (block, ']') = '\0';
  }
  int newvartype = (!strcmp (type, "scalar") ? scalar :
		    !strcmp (type, "vector") ? vector :
		    !strcmp (type, "tensor") ? tensor :
		    -1);
  var_t * var;
  if (yytext[0] == '=') {
    var = vartop();
    assert (var && var->scope == scope);
    if (var->type != newvartype)
      return yyerror ("type mismatch in new");
  }
  else {
    char * s = yytext; space (s); *s = '\0';
    var = varlookup (yytext, strlen(yytext));
    if (!var) {
      fprintf (stderr, "%s:%d: undeclared %s '%s'\n", fname, line, type, yytext);
      return 1;
    }
    if (var->type != newvartype)
      return yyerror ("type mismatch in new");
    fprintf (yyout, "%s ", yytext);
  }
  var->symmetric = (symmetric != NULL);
  var->face   = (face != NULL);
  var->vertex = (vertex != NULL);
  fputc ('=', yyout);
  new_field (var, block);
  if (debug)
    fprintf (stderr, "%s:%d: new %s%s: %s\n", 
	     fname, line, 
	     var->symmetric ? "symmetric " : 
	     var->face ? "face " :
	     var->vertex ? "vertex " :
	     "", 
	     type, var->v);
}

{ID}{WS}*[=]{WS}*automatic{WS}*[(][^)]*[)] |
[=]{WS}*automatic{WS}*[(][^)]*[)] {
  var_t * var;
  if (yytext[0] == '=') {
    var = vartop();
    assert (var && var->scope == scope);
  }
  else {
    char * s = yytext; space (s); *s++ = '\0';
    var = varlookup (yytext, strlen(yytext));
    if (!var) {
      fprintf (stderr, "%s:%d: '%s' undeclared", fname, line, yytext);
      return 1;
    }
    fprintf (yyout, "%s ", yytext);
    yytext = s;
  }
  var->automatic = 1;
  char * arg = strchr (yytext, '(');
  if (arg[1] == ')')
    fputc ('=', yyout);
  else {
    var->conditional = malloc (strlen(arg) + 5);
    sprintf (var->conditional, "%s%s.i", arg, 
	     var->type == vector ? ".x" : var->type == tensor ? ".x.x" : "");
    fprintf (yyout, "= %s ? %s :", var->conditional, arg);
  }
  new_field (var, NULL);
}

\({WS}*const{WS}*\)    varmaybeconst = 1;

(void|char|short|int|long|float|double|coord|stats|norm){WS}+{TYPE} {
  char * var = yytext;
  space (var);
  nonspace (var);
  int npointers = 0;
  while (*var == '*') var++, npointers++;
  nonspace (var);  if (para == 0) { /* declaration */
    varsymmetric = varface = varvertex = varmaybeconst = 0;
    varconst = NULL;
    vartype = other_type;
    free (type);
    char * s = yytext; space (s);
    type = strndup (yytext, s - yytext);
    declaration (var, yytext, npointers);
    invardecl = scope + 1;
  }
  else if (para == 1) { /* function prototype (no nested functions) */
    char * b = strchr (yytext, '[');
    if (b) {
      if (!strchr(yytext, ']'))
	brack++;
      *b = '\0';
    }
    varpush (var, other_type, scope + 1, 0);
    invardecl = varmaybeconst = 0;
    if (b)
      *b = '[';
    cadna_echo (yytext);
  }
  else
    cadna_echo (yytext);
}

struct{WS}+{ID}{WS}+{ID} {
  char * var = yytext;
  space (var); nonspace (var); space (var); nonspace (var);
  if (para == 0) /* declaration */
    varpush (var, struct_type, scope, 0);
  else if (para == 1) /* function prototype (no nested functions) */
    varpush (var, struct_type, scope + 1, 0);
  if (debug)
    fprintf (stderr, "%s:%d: %s\n", fname, line, yytext);
  ECHO;
}

symmetric{WS}+tensor{WS}+[a-zA-Z0-9_\[\]]+ |
face{WS}+vector{WS}+[a-zA-Z0-9_\[\]]+ |
vertex{WS}+scalar{WS}+[a-zA-Z0-9_\[\]]+ |
(scalar|vector|tensor){WS}+[a-zA-Z0-9_\[\]]+ |
bid{WS}+{ID} {
  varsymmetric = (strstr(yytext, "symmetric") == yytext);
  varface = (strstr(yytext, "face") == yytext);
  varvertex = (strstr(yytext, "vertex") == yytext);
  varconst = NULL;
  char * var = strstr(yytext,"scalar");
  vartype = scalar;
  if (!var) {
    var = strstr(yytext,"vector");
    vartype = vector;
    if (!var) {
      var = strstr(yytext,"tensor");
      vartype = tensor;
      if (!var) {
	var = strstr(yytext,"bid");
	vartype = bid;
      }
    }
  }
  char * text = var;
  if (vartype != bid)
    var = &var[7];
  else
    var = &var[4];
  nonspace (var);  
  if (*var == '[') {
    // scalar [..
    for (; *var != '\0'; var++)
      if (*var == '[') brack++;
      else if (*var == ']') brack--;
    fputs (text, yyout);
  }
  else if (para == 0) { /* declaration */
    free (type); type = NULL;
    declaration (var, text, 0);
    invardecl = scope + 1;
  }
  else if (para == 1) { /* function prototype (no nested functions) */
    fputs (text, yyout);
    if (debug)
      fprintf (stderr, "%s:%d: proto: %s\n", fname, line, var);
    varpush (var, vartype, scope + 1, varmaybeconst);
    invardecl = varmaybeconst = 0;
  }
  else
    fputs (text, yyout);
}

const{WS}+(symmetric{WS}+|face{WS}+|vertex{WS}+|{WS}*)(scalar|vector|tensor){WS}+{ID}{WS}*=[^;]+ {
  ECHO;
  char * s = strchr (yytext, '='); s--;
  while (strchr (" \t\v\n\f", *s)) s--; s++; *s = '\0';
  fprintf (stderr, "%s:%d: warning: did you mean '%s[]'?\n", 
	   fname, line, yytext);
}

const{WS}+(symmetric{WS}+|face{WS}+|vertex{WS}+|{WS}*)(scalar|vector|tensor){WS}+{ID}\[{WS}*\]{WS}*=[^;]+ {
  // const scalar a[] = 1.;
  char * var = strstr(yytext,"scalar");
  vartype = scalar;
  if (!var) {
    var = strstr(yytext,"vector");
    vartype = vector;
    if (!var) {
      var = strstr(yytext,"tensor");
      vartype = tensor;
    }
  }
  char * text = var;
  var = &var[7];
  nonspace (var);
  char * cst = strchr (var, ']'); cst++;
  if (*cst == '=')
    *cst++ = '\0';
  else {
    *cst++ = '\0';
    cst = strchr (cst, '=') + 1;
  }
  if (para != 0)
    return yyerror ("constant fields can only appear in declarations");
  if (debug)
    fprintf (stderr, "%s:%d: const '%s' '%s' '%s'\n",
	     fname, line, var, cst, text);
  varconst = cst;
  varsymmetric = varface = varvertex = 0;
  free (type); type = NULL;
  declaration (var, text, 0);
}

,{WS}*{TYPE} {
  if (para == 0 && brack == 0 && invardecl == scope + 1) {
    char * var = &yytext[1];
    nonspace (var);
    int npointers = 0;
    while (*var == '*') var++, npointers++;
    nonspace (var);
    declaration (var, yytext, npointers);
  }
  else
    REJECT;
}

return {
  char * list = automatic_list (0, 1);
  if (list || inevent || intrace) {
    // returning from a function: delete automatic fields before returning
    // note that this assumes that the function scope is always 1 
    // (i.e. no nested functions allowed).
    free (list);
    inreturn = 1;
    if (return_type)
      fprintf (yyout, "{ %s _ret = ", return_type);
    else
      fputs ("{ ", yyout);
  }
  else
    ECHO;
}

val{WS}*[(]    {
  inval = 1; invalpara = para++;
  assert (!infine && !inarray);
  infine = para;
  inarrayargs = 0;
  ECHO;
}

{SCALAR}{WS}*\[{WS}*. {
  /* v[... */
  var_t * var = NULL;
  if ((inforeach || infunction) && !inval) {
    char * s = yytext;
    while (!strchr(" \t\v\n\f[.", *s)) s++;
    if ((var = varlookup (yytext, s - yytext))) {
      s = yytext;
      while (*s != '[') s++;
      *s = '\0';
      if ((dimension < 2 && (component (yytext, 'y') ||
			     component (yytext, 't'))) ||
	  (dimension < 3 && (component (yytext, 'z')||
			     component (yytext, 'r')))) {
	if (debug)
	  fprintf (stderr, "%s:%d: the dimension of '%s' is too high\n",
		   fname, line, yytext);
	fprintf (yyout, "_val_higher_dimension(%s", yytext);
      }
      else if (var->constant)
	fprintf (yyout, "_val_constant(%s", boundaryrep(yytext));
      else if (var->maybeconst)
	maybeconst_macro (var, "val", boundaryrep(yytext));
      else
	fprintf (yyout, "val(%s", boundaryrep(yytext));
      if (yytext[yyleng-1] == ']')
	/* v[] */
	fputs (",0,0,0)", yyout);
      else {
	if (var->args > 0) {
	  /* v[...][... */
	  fputc ('[', yyout);
	  fputc (yytext[yyleng-1], yyout);
	  inarray = brack++; 
	  int j = var->args;
	  while (j) {
	    int c = getput();
	    if (c == EOF)
	      return yyerror ("unexpected EOF");
	    if (c == '[') brack++;
	    if (c == ']') {
	      brack--;
	      if (brack == inarray)
		j--;
	    }
	  }
	  int c = input();
	  if (c != '[') {
	    fprintf (stderr, "%s:%d: error: expecting '[' not '", 
		     fname, line);
	    fputc (c, stderr);
	    fputs ("'\n", stderr);
	    return 1;
	  }
	  fputc (',', yyout);
	  c = input();
	  if (c == ']')
	    /* v[...][] */
	    fputs ("0,0,0)", yyout);
	  else {
	    fputc (c, yyout);
	    inarray = ++brack;
	  }
	}
	else {
	  /* v[...] */
	  fputc (',', yyout);
	  unput (yytext[yyleng-1]);
	  inarray = ++brack;
	}
	inarraypara = para;
	inarrayargs = 1;
      }
    }
  }
  if (!var)
    REJECT;
}

[a-zA-Z_0-9]+[.ntr]*{WS}*\[{WS}*{ID}{WS}*\]{WS}*= {
  /* v[top] = ..., u.n[left] = ..., u.t[left] = ... */
  char * s = yytext;
  while (!strchr(" \t\v\n\f[.", *s)) s++;
  var_t * var;
  if ((var = varlookup (yytext, s - yytext))) {
    strcpy (boundaryvar, yytext);
    s = boundaryvar;
    while (!strchr(" \t\v\n\f[", *s)) s++;
    *s++ = '\0';
    nonspace (s);
    strcpy (boundarydir, s);
    s = boundarydir;
    while (!strchr(" \t\v\n\f]", *s)) s++;
    *s++ = '\0';
    var_t * dir = varlookup (boundarydir, strlen(boundarydir));
    if ((!dir || dir->type != bid) &&
	strcmp(boundarydir, "right") && strcmp(boundarydir, "left") &&
	(dimension < 2 || (strcmp(boundarydir, "top") &&
			   strcmp(boundarydir, "bottom"))) &&
	(dimension < 3 || (strcmp(boundarydir, "front") &&
			   strcmp(boundarydir, "back"))))
      REJECT;
    if (debug)
      fprintf (stderr, "%s:%d: boundarydir: %s yytext: %s\n",
               fname, line, boundarydir, yytext);
    boundarycomponent = boundary_component (boundaryvar);
    if (!var->face)
      boundarycomponent = 0;
    boundaryfp = yyout;

    yyout = boundaryheader;
    fprintf (yyout,
	     "#line %d \"%s\"\n"
	     "static %s _boundary%d"
             " (Point point, Point neighbor, scalar _s, void * data) {",
             line, fname, doubletype, nboundary);
    boundary_staggering (yyout);
    fputs (" POINT_VARIABLES; ", yyout);
    maps (line - 1);
    nmaybeconst = 0;
    inboundary = inforeach_boundary = 1;
    inforeach = 2;

    yyout = dopen ("_inboundary.h", "w");
    fputs ("return ", yyout);
  }
  else
    REJECT;
}

[a-zA-Z_0-9]+[.ntr]*{WS}*\[{WS}*(right|left|top|bottom|front|back){WS}*\]{WS}*={WS}*periodic{WS}*\({WS}*\){WS}*; {
  /* v[top] = periodic(); */
  char * s = yytext;
  while (!strchr(" \t\v\n\f[.", *s)) s++;
  var_t * var;
  if ((var = varlookup (yytext, s - yytext))) {
    strcpy (boundaryvar, yytext);
    s = boundaryvar;
    while (!strchr(" \t\v\n\f[", *s)) s++;
    *s++ = '\0';
    nonspace (s);
    strcpy (boundarydir, s);
    s = boundarydir;
    while (!strchr(" \t\v\n\f]", *s)) s++;
    *s++ = '\0';
    if ((dimension < 2 && (!strcmp(boundarydir, "top") ||
			   !strcmp(boundarydir, "bottom"))) ||
	(dimension < 3 && (!strcmp(boundarydir, "front") ||
			   !strcmp(boundarydir, "back"))))
      REJECT;
    if (debug)
      fprintf (stderr, "%s:%d: boundarydir: %s yytext: %s\n",
               fname, line, boundarydir, yytext);
    boundarycomponent = boundary_component (boundaryvar);
    if (!var->face)
      boundarycomponent = 0;

    if (scope == 0) {
      /* file scope */
      fprintf (yyout, 
	       "static void _set_boundary%d (void) { ",
	       nboundary);
      boundaryindex[nsetboundary++] = nboundary;
    }
    /* function/file scope */
    fprintf (yyout,
	     "_attribute[%s.i].boundary[%s] = "
	     "_attribute[%s.i].boundary_homogeneous[%s] = "
	     "_attribute[%s.i].boundary[%s] = "
	     "_attribute[%s.i].boundary_homogeneous[%s] = "
	     "periodic_bc;",
	     boundaryvar, boundarydir,
	     boundaryvar, boundarydir,
	     boundaryvar, opposite (boundarydir),
	     boundaryvar, opposite (boundarydir));
    if (scope == 0)
      /* file scope */
      fputs (" } ", yyout);
    periodic[nboundary++] = 1;

    /* keep newlines */
    s = yytext;
    while (*s != '\0')
      if (*s++ == '\n')
	fputc ('\n', yyout);
  }
  else
    REJECT;
}

(fine|coarse){WS}*\({WS}*{SCALAR}{WS}* {
  /* fine(v... */
  var_t * var = NULL;
  if ((inforeach || infunction) && !inval) {
    char * s = strchr (yytext, '('); s++;
    while (strchr(" \t\v\n\f", *s)) s++;
    char * id = s;
    while (!strchr(" \t\v\n\f.", *s)) s++;
    if ((var = varlookup (id, s - id))) {
      para++;
      s = id;
      while (!strchr(" \t\v\n\f", *s)) s++;
      *s = '\0';
      char * op = yytext[0] == 'f' ? "fine" : "coarse";
      if (debug)
	fprintf (stderr, "%s:%d: %s: %s\n", fname, line, op, id);
      if (var->maybeconst)
	maybeconst_macro (var, op, id);
      else
	fprintf (yyout, "%s(%s", op, boundaryrep(id));
      infine = para;
      inarrayargs = 0;
    }
  }
  if (!var)
    REJECT;
}

(allocated|allocated_child|neighbor|neighborp|aparent|child){WS}*\( {
  assert (!infine);
  infine = ++para;
  inarrayargs = 1;
  ECHO;
}
  
dirichlet{WS}*[(] {
  para++;
  if (inboundary && boundarycomponent == 'x')
    fputs (&yytext[9], yyout);
  else
    ECHO;
}

ghost {
  if (inforeach_boundary) {
    char * name[3] = {"(ig > 0 ? 1 : ig < 0 ? -1 : 0)",
		      "(ig > 0 ? 1 : ig < 0 ? -1 : 0),"
		      "(jg > 0 ? 1 : jg < 0 ? -1 : 0)",
		      "(ig > 0 ? 1 : ig < 0 ? -1 : 0),"
		      "(jg > 0 ? 1 : jg < 0 ? -1 : 0),"
		      "(kg > 0 ? 1 : kg < 0 ? -1 : 0)"};
    fputs (name[dimension - 1], yyout);
    if (infine || inarray)
      inarrayargs += dimension - 1;
  }
  else
    ECHO;
}

Point{WS}+point[^{ID}] {
  /* Point point */
  if (indef)
    REJECT;
  if (debug)
    fprintf (stderr, "%s:%d: infunction\n", fname, line);
  infunction = 1; infunctiondecl = 0;
  functionscope = scope; functionpara = para;
  if (yytext[yyleng-1] == ')') para--;
  ECHO;
}

for{WS}*[(]{WS}*(scalar|vector|tensor){WS}+{ID}{WS}+in{WS}+ {
  /* for (scalar .. in .. */
  char * s = strchr (&yytext[3], 'r'); s++;
  int vartype = s[-3] == 's' ? tensor : s[-3] == 't' ? vector : scalar;
  nonspace (s);
  char * id = s;
  while (!strchr (" \t\v\n\f", *s)) s++;
  *s = '\0';
  char * list = malloc (sizeof(char)); list[0] = '\0';
  int nl = 1, c;
  char last[] = ")";
  while ((c = input()) != EOF) {
    if (c == ')')
      break;
    else {
      list = realloc (list, (nl + 1)*sizeof(char));
      list[nl - 1] = c;
      list[nl++] = '\0';
    }
  }
  if (c != ')')
    return yyerror ("expecting ')'");
  static int i = 0;
  if (list[0] == '{') {
    char * oldlist = list;
    list = makelist (oldlist, vartype);
    free (oldlist);
    if (list == NULL)
      return yyerror ("invalid list");
  }
  else
    fprintf (yyout, "if (%s) ", list);
  if (debug)
    fprintf (stderr, "%s:%d: %s\n", fname, line, list);
  fprintf (yyout,
	   "for (%s %s = *%s, *_i%d = %s; ((scalar *)&%s)->i >= 0; %s = *++_i%d%s",
	   vartype == scalar ? "scalar" : 
	   vartype == vector ? "vector" : 
                               "tensor", 
	   id,
	   list, i, list, 
	   id, 
	   id, i, last);
  free (list);
  i++;
  varpush (id, vartype, scope + 1, 0);
}

for{WS}*[(][^)]+,[^)]+{WS}+in{WS}+[^)]+,[^)]+[)] {
  /* for (a,b in c,d) */
  char * id[10] = {NULL}, * list[10] = {NULL};
  int nid = 0, nlist = 0, inin = 0;
  char * s = strchr (yytext, '('); s++;
  s = strtok (s, " \t\v\n\f,");
  while (s) {
    if (!strcmp (s, "in"))
      inin = 1;
    else if (inin)
      list[nlist++] = s;
    else
      id[nid++] = s;
    s = strtok (NULL, " \t\v\n\f,)");
  }
  if (nlist != nid)
    return yyerror ("lists must be the same size");
  static int index = 0, i;
  char * sc = NULL, * vc = NULL, * tc = NULL;
  for (i = 0; i < nid; i++) {
    var_t * var = varlookup (id[i], strlen(id[i]));
    if (!var) {
      fprintf (stderr, 
	       "%s:%d: error: '%s' is not a scalar, vector or tensor\n", 
	       fname, line, id[i]);
      return 1;
    }
    fprintf (yyout, "%s * _i%d = %s; ", 
	     var->type == vector ? "vector" : 
	     var->type == scalar ? "scalar" :
	                           "tensor",
	     index + i, list[i]);
    if (!sc && var->type == scalar) sc = id[i];
    if (!vc && var->type == vector) vc = id[i];
    if (!tc && var->type == tensor) tc = id[i];
  }
  fprintf (yyout, "if (%s) for (%s = *%s", list[0], id[0], list[0]);
  for (i = 1; i < nid; i++)
    fprintf (yyout, ", %s = *%s", id[i], list[i]);
  fprintf (yyout, "; ((scalar *)&%s)->i >= 0; ", id[0]);
  fprintf (yyout, "%s = *++_i%d", id[0], index);
  for (i = 1; i < nid; i++)
    fprintf (yyout, ", %s = *++_i%d", id[i], index + i);
  fputc (')', yyout);
  index += nid;
}

[\n] {
  if (indef)
    fputc ('\\', yyout);
  ECHO;
}

"["        ECHO; brack++;
"]"        {
  if (inarray == brack) {
    int i;
    for (i = inarrayargs; i < 3; i++)
      fputs (",0", yyout);
    fputc (')', yyout);
    inarray = 0;
  }
  else if (inmalloc == brack) {
    fputs ("))", yyout);
    inmalloc = 0;
  }
  else
    ECHO;
  brack--;
  if (brack < 0)
    return yyerror ("mismatched ']'");  
}

event{WS}+{ID}{WS}*[(] {
  /* event (... */
  char ids[80], * id = ids;
  strncpy (id, yytext, 80);
  space (id); nonspace (id);
  char * s = id; while (!strchr(" \t\v\n\f(", *s)) s++; *s = '\0';
  eventparent[nevents] = eventchild[nevents] = -1;
  eventlast[nevents] = 0;
  eventid[nevents] = strdup (id);
  // check whether an event with the same name already exists
  // if it does, append a number (up to 9)
  char i = '0', found = 1;
  int len = strlen(id), lastfound = -1;
  id[len+2] = '\0';
  while (i <= '9' && found) {
    found = 0;
    int j;
    for (j = 0; j < nevents && !found; j++)
      if (!strcmp (id, eventfunc[j])) {
	lastfound = j;
	found = 1;
      }
    if (found) {
      id[len] = '_';
      id[len+1] = i++;
    }
  }
  if (lastfound >= 0) {
    eventparent[lastfound] = nevents;
    eventchild[nevents] = lastfound;
  }
  fprintf (yyout, 
	   "static int %s_expr%d (int * ip, %s * tp, Event * _ev) {"
	   "  int i = *ip; %s t = *tp;"
	   "  int ret = (", id, inevent++, doubletype, doubletype);
  eventscope = scope; eventpara = ++para;
  eventarray[nevents] = 0;
  eventfile[nevents] = strdup (fname);
  eventline[nevents] = line;
  eventfunc[nevents] = strdup (id);
}

[it]{WS}*={WS}*[{][^}]*[}] {
  if (inevent == 1) {
    eventarray[nevents] = yytext[0];
    yytext[yyleng-1] = '\0';
    eventarray_elems[nevents] = strdup (strchr (yytext, '{'));
    fputs ("1); "
	   "  *ip = i; *tp = t; "
	   "  return ret; "
	   "} ", yyout);
  }
  else
    REJECT;
}

,{WS}*last {
  if (inevent == 1)
    eventlast[nevents] = 1;
  else
    ECHO;
}

,{WS}*first {
  if (inevent == 1) {
    int i = nevents;
    while (i >= 0) {
      eventlast[i] = 0;
      i = eventchild[i];
    }
  }
  else
    ECHO;
}

end {
  if (inevent == 1)
    fputs ("1234567890", yyout);
  else
    ECHO;
}

trace {
  // function tracing
  if (scope == 0)
    traceon = 1;
  else
    REJECT;
}

stderr fputs ("qstderr()", yyout);
stdout fputs ("qstdout()", yyout);

foreach_dimension{WS}*[(]([1-3]|{WS})*[)] {
  int c;
  foreachdim_t dim;
  dim.dim = dimension;
  for (c = '1'; c <= '3'; c++)
    if (strchr(yytext, c))
      dim.dim = c - '1' + 1;
  if (debug)
    fprintf (stderr, "%s:%d: foreach_dimension (%d)\n", fname, line, dim.dim);
  dim.scope = scope; dim.para = para;
  dim.fp = yyout;
  char name[80];
  sprintf (name, "_dimension_%d.h", foreachdim_stack.n);
  yyout = dopen (name, "w+");
  dim.line = indef ? -1 : line;
  stack_push (&foreachdim_stack, &dim);
}

{ID}{WS}+{ID}{WS}*\( {
  if (scope == 0 && para == 0) {
    // function prototype
    para++;
    cadna_echo (yytext);
    char * s1 = yytext; space (s1); *s1++ = '\0';
    nonspace (s1);
    char * s2 = s1; while (!strchr (" \t\v\n\f(", *s2)) s2++; *s2 = '\0';
    free (return_type);
    return_type = strdup (cadna_type (yytext));
    infunctionproto = 1;
    if (debug)
      fprintf (stderr, "%s:%d: function '%s' returns '%s'\n", 
	       fname, line, s1, yytext);
    if (!strcmp (return_type, "void")) {
      free (return_type);
      return_type = NULL;
    }      
    if (!strcmp (s1, "main") && !strcmp (yytext, "int"))
      inmain = 1;
    else if (!strncmp (s1, "begin_", 6)) {
      if (debug)
	fprintf (stderr, "%s:%d: %s\n", fname, line, s1);
      free (inbegin);
      inbegin = strdup (s1 + 6);
    }
    else if (inbegin && !strncmp (s1, "end_", 4) &&
	     !strcmp (inbegin, s1 + 4)) {
      if (debug)
	fprintf (stderr, "%s:%d: %s\n", fname, line, s1);
      char ** b = blocks;
      int n = 0;
      while (b && *b)
	n++, b++;
      blocks = realloc (blocks, (n + 2)*sizeof (char *));
      blocks[n] = inbegin;
      blocks[n + 1] = NULL;
      inbegin = NULL;
    }
    // function tracing
    if (traceon) {
      intrace = 1; traceon = 0;
      free (tracefunc);
      tracefunc = strdup (s1);
    }
    else
      intrace = 0;
  }
  else
    REJECT;
}

,?{WS}*reduction{WS}*[(](min|max|\+):{ID}[)] {
  if (debug)
    fprintf (stderr, "%s:%d: '%s'\n", fname, line, yytext);
  if (yytext[0] == ',')
    yytext[0] = ' ';
  char * s = strchr (yytext, '('), * s1 = strchr (yytext, ':');
  *s1 = '\0'; s1++;
  assert (nreduct < REDUCTMAX);
  strcpy (reduction[nreduct], ++s);
  yytext[yyleng-1] = '\0';
  strcpy (reductvar[nreduct++], s1);
}

[{]({SP}*[a-zA-Z_0-9.]+{SP}*,{SP}*)*[a-zA-Z_0-9.]+{SP}*[}] {
  // {h, zb, ...}
  char * list = makelist (yytext, -1);
  if (list == NULL)
    ECHO;
  else {
    if (debug)
      fprintf (stderr, "%s:%d: %s\n", fname, line, list);
    fputs (list, yyout);
    free (list);
  }
}

{ID} {
  var_t * var;
  if (inforeach) {
    int i;
    for (i = 0; i < nreduct; i++)
      if (!strcmp (yytext, reductvar[i])) {
	fputc ('_', yyout);
	break;
      }
    cadna_echo (yytext);
  }
  else if (scope == 0 && para == 0 &&
	   (var = varlookup (yytext, strlen(yytext))) &&
	   var->constant) {
    // replace global scalar/vector constants
    if (var->type == scalar)
      fprintf (yyout, "{(_NVARMAX + %d)}", var->i[0]);
    else if (var->type == vector) {      
      fprintf (yyout, "{{_NVARMAX + %d}", var->i[0]);
      int i;
      for (i = 1; i < dimension; i++)
	fprintf (yyout, ",{_NVARMAX + %d}", var->i[i]);
      fputc ('}', yyout);
    }
    else
      assert (0);
  }
  else
    cadna_echo (yytext);
}

{SCALAR}[.][xyzi] ECHO;

{SCALAR}[.][a-wA-Z_0-9]+ {
  // scalar attributes
  char * dot1 = strchr (yytext, '.');
  var_t * var = varlookup (yytext, strlen(yytext) - strlen(dot1));
  if (var) {
    char * dot = dot1;
    while (dot1) {
      dot = dot1;
      dot1 = strchr (dot1 + 1, '.');
    }
    dot[0] = '\0'; dot++;
    fprintf (yyout, "_attribute[%s.i].%s", yytext, dot);
    if (debug)
      fprintf (stderr, "%s:%d: _attribute[%s.i].%s\n", fname, line, yytext, dot);
  }
  else
    ECHO;
}

{ID}{WS}+{ID}{WS}*[(]{WS}*struct{WS}+{ID}{WS}+{ID}{WS}*[)] {
  // function declaration with struct argument
  cadna_echo (yytext);
  char * s = yytext;
  space (s); *s++ = '\0'; nonspace (s);
  free (return_type);
  return_type = strdup (cadna_type (yytext));
  args = realloc (args, sizeof (char *)*++nargs);
  argss = realloc (argss, sizeof (char *)*nargs);
  char * s1 = s; while (!strchr(" \t\v\n\f(", *s1)) s1++;
  *s1++ = '\0'; s1 = strstr (s1, "struct"); space (s1); nonspace (s1);
  char * s2 = s1; space (s2); *s2++ = '\0';
  char * s3 = s2; space (s3); if (*s3 == '\0') s3[-1] = '\0';
  args[nargs-1] = strdup (s);
  argss[nargs-1] = strdup (s1);
  varpush (s2, struct_type, scope + 1, 0);

  if (!strncmp (s, "begin_", 6)) {
    if (debug)
      fprintf (stderr, "%s:%d: %s\n", fname, line, s);
    free (inbegin);
    inbegin = strdup (s + 6);
  }
  // function tracing
  if (traceon) {
    intrace = 1; traceon = 0;
    free (tracefunc);
    tracefunc = strdup (s);
  }
  else
    intrace = 0;
  if (debug)
    fprintf (stderr, 
	     "%s:%d: function '%s' with struct '%s' '%s' returns '%s'\n", 
	     fname, line, s, s1, s2, return_type);
}

{ID}{WS}*[(]{WS}*[)] {
  // function call with no 'args' assignment
  if (inarg)
    REJECT;
  char * s = yytext; while (!strchr(" \t\v\n\f(", *s)) s++;
  int len = s - yytext, i;
  for (i = 0; i < nargs && !inarg; i++)
    if (strlen(args[i]) == len && !strncmp (args[i], yytext, len)) {
      if (debug)
	fprintf (stderr, "%s:%d: no argument '%s'\n", fname, line, args[i]);
      fprintf (yyout, "%s ((struct %s){0})", args[i], argss[i]);
      inarg = 1;
    }
  if (!inarg)
    REJECT;
  inarg = 0;
}

{ID}{WS}*[(]{WS}*{ID}{WS}*[)] {
  // function call with a single 'args' assignment
  if (inarg)
    REJECT;
  char * s = yytext; while (!strchr(" \t\v\n\f(", *s)) s++;
  int len = s - yytext, i;
  for (i = 0; i < nargs && !inarg; i++) {
    if (strlen(args[i]) == len && !strncmp (args[i], yytext, len)) {
      char * s = strchr (yytext, '('); s++; nonspace (s);
      char * s1 = s; space (s1); *s1 = '\0';
      s1 = strchr (s, ')'); if (s1) *s1 = '\0';
      if (debug)
	fprintf (stderr, "%s:%d: single argument '%s'\n", fname, line, s);
      if (get_vartype (s) == struct_type)
	fprintf (yyout, "%s (%s)", args[i], s);
      else
	fprintf (yyout, "%s ((struct %s){%s})", args[i], argss[i], s);
      inarg = 1;
    }
  }
  if (!inarg)
    REJECT;
  inarg = 0;
}

{ID}{WS}*[(] {
  // function call with multiple 'args' assignment
  if (inarg)
    REJECT;
  char * s = yytext;
  while (!strchr(" \t\v\n\f(", *s)) s++;
  int len = s - yytext;
  char ** b = blocks, * func = yytext;
  while (b && *b) {
    if (strlen(*b) == len && !strncmp (*b, yytext, len)) {
      if (debug)
	fprintf (stderr, "%s:%d: block %s\n", fname, line, *b);
      func = malloc (strlen (yytext) + 7);
      strcpy (func, "begin_");
      strcat (func, yytext);
      len += 6;
      fputs ("{ ", yyout);
      break;
    }
    b++;
  }

  int i;
  for (i = 0; i < nargs && !inarg; i++)
    if (strlen(args[i]) == len && !strncmp (args[i], func, len)) {
      fputs (func, yyout);
      para++;
      inarg = para;
      fprintf (yyout, "(struct %s){", argss[i]);
    }
  
  if (func != yytext) {
    if (!inarg) {
      fputs (func, yyout); para++;
    }
    free (func);
    foreach_child_t block = {scope, para - 1};
    strcpy (block.name, *b);
    stack_push (&inblock_stack, &block);
  }
  else if (!inarg)
    REJECT;
}

[fs]?printf{WS}*[(] {
  if (cadna) {
    ECHO; para++;
    incadna = para;
    incadnanargs = 0;
  }
  else
    REJECT;
}

{ID}{WS}*=[^=] {
  // arguments of function call with 'args' assignment
  if (inarg && para == inarg) {
    fputc('.', yyout);
    ECHO;
  }
  else
    REJECT;
}

#{SP}+[0-9]+{SP}+\"[^\"]+\".*$ {
  /* line numbers */
  char * ln = yytext, * name, * quote;
  space (ln); nonspace(ln);
  name = ln; space (name); *name++ = '\0'; nonspace(name);
  line = atoi(ln);
  free (fname);
  name++; quote = strchr (name, '"'); *quote = '\0';
  fname = strdup (name);
  if (!indef)
    fprintf (yyout, "#line %d \"%s\"", line, fname);
  if (tracefp)
    fclose (tracefp);
  char * tracename = malloc (strlen(fname) + strlen(".trace") + 1);
  strcpy (tracename, fname);
  char * s = &tracename[strlen(tracename) - 1];
  while (*s != '.' && s != tracename) s--;
  if (s != tracename)
    strcpy (s, ".trace");
  else
    strcat (s, ".trace");
  tracefp = fopen ("tracename", "r");
  free (tracename);
}

^{SP}*@{SP}*def{SP}+(foreach|end_foreach|trace) |
^{SP}*@{SP}*def{SP}+ {
  // @def ...
  char * s = strstr (yytext, "def");
  space (s);
  fputs ("#define", yyout);
  fputs (s, yyout);
  indef = 1;
}

@ {
  if (indef) {
    fprintf (yyout, "\n#line %d\n", line - 1);
    indef = 0;
  }
  else
    ECHO;
}

^{SP}*#{SP}*pragma{SP}+autolink{SP}+.*$ {
  char * s = strstr (yytext, "autolink"); space (s);
  if (!autolink) {
    autolink = malloc (strlen(s) + 1);
    autolink[0] = '\0';
  }
  else
    autolink = realloc (autolink, strlen(s) + 1 + strlen (autolink));
  strcat (autolink, s);
  if (debug)
    fprintf (stderr, "%s:%d: %s\n", fname, line, autolink);
}
  
^{SP}*@{SP}*{ID}{SP}+(foreach|end_foreach|trace) |
^{SP}*@{SP}*{ID} {
  // @... foreach...
  // @... end_foreach...
  // @...
  yytext = strchr(yytext, '@'); yytext++;
  fprintf (yyout, "#%s", yytext);
  register int oldc = 0, c;
  char * text = malloc (100);
  int len = 0, maxlen = 100;
  for (;;) {
    while ((c = input()) != '\n' && c != EOF) {
      if (c == '(') para++;
      if (c == ')') para--;
      if (para < 0)
	return yyerror ("mismatched ')'");
      if (c == '{') scope++;
      if (c == '}') scope--;
      if (scope < 0)
	return yyerror ("mismatched '}'");
      oldc = c;    /* eat up text of preproc */
      if (len == maxlen - 1)
	maxlen += 100, text = realloc (text, maxlen);	
      text[len++] = c;
    }
    if (c == '\n') {
      if (len == maxlen - 1)
	maxlen += 100, text = realloc (text, maxlen);	
      text[len++] = c;
    }
    if (c == EOF || oldc != '\\')
      break;
  }
  text[len] = '\0';
  cadna_echo (text);
  free (text);
}

^{SP}*@.*" Pragma(" {
  yytext = strchr(yytext, '@'); yytext++;
  char * s = strstr (yytext, "Pragma("); *s++ = '\0';
  fprintf (yyout, "#%s", yytext);
  fputs ("_Pragma(", yyout);
  register int oldc = 0, c;
  for (;;) {
    while ((c = getput()) != '\n' && c != EOF)
      oldc = c;    /* eat up text of preproc */
    if (c == EOF || oldc != '\\')
      break;
  }
}

static{WS}+FILE{WS}*[*]{WS}*{ID}{WS}*= {
  if (scope == 0)
    REJECT;
  // static FILE * fp = ...
  ECHO;
  char * s = yytext;
  while (*s != '*') s++;
  s++;
  nonspace (s);
  char * id = s;
  space (s);
  *s = '\0';
  fprintf (yyout,
	   "NULL; if (!%s || i == 0) %s = pid() > 0 ? "
	   "fopen(\"/dev/null\", \"w\") : ", id, id);
}

{D}+{E}{FS}?		|
{D}*"."{D}+({E})?{FS}?	|
{D}+"."{D}*({E})?{FS}?	{
  if (cadna)
    /* Force (double_st) type casting of floating-point constants (for
       compatibility with CADNA) */
    fprintf (yyout, "(double_st)%s", yytext);
  else
    REJECT;
}

"/*"                                    { ECHO; if (comment()) return 1; }
"//".*                                  { ECHO; /* consume //-comment */ }

({SP}?\"([^\"\\\n]|{ES})*\"{WS}*)+	{
  /* STRING_LITERAL */
  if (incadna) {
    char * s = yytext;
    int nargs = 0;
    incadnanargs = 0;
    while (*s != '\0') {
      if (*s == '%') {
	char * f = s++;
	while (*s != '\0' && !strchr("diouxXeEfFgGaAcspn%", *s)) s++;

	switch (*s) {
	case 'e': case 'E': case 'f': case 'F': case 'g': case 'G':
	case 'a': case 'A':
	  fputs ("%s", yyout);
	  incadnaarg[incadnanargs++] = nargs++;
	  assert (incadnanargs < 80);
	  break;
	  
	default:
	  if (*s != 'm' && *s != '%')
	    nargs++;
	    
	  while (f != s)
	    fputc (*f++, yyout);
	  fputc (*s, yyout);
	}
      }
      else
	fputc (*s, yyout);
      s++;
    }
    incadnaargs = 0;
#if 0
    if (incadnanargs) {
      fprintf (stderr, "%s:%d: incadna: %d para: %d infine: %d %s:",
	       fname, line, incadna, para, infine, yytext);
      int i;
      for (i = 0; i < incadnanargs; i++)
	fprintf (stderr, " %d", incadnaarg[i]);
      fputc ('\n', stderr);
    }
#endif
  }
  else
    ECHO;
}

'.' { ECHO; /* character literal */ }

%%

int yyerror (const char * s)
{
  fprintf (stderr, "%s:%d: error: %s\n", fname, line, s);
  return 1;
}

int comment(void)
{
  int c;
  while ((c = getput())) {
    if (c == '*') {
      while ((c = getput()) == '*')
	;
      if (c == '/')
	return 0;
      if (c == 0)
	break;
    }
  }
  return yyerror ("unterminated comment");
}

int getput(void)
{
  int c = input();
  fputc (c, yyout);
  return c;
}

void stripname (char * path);
char * stripslash (char * path);
int includes (int argc, char ** argv, char ** out, 
	      char ** grid, int * default_grid, int * dimension, int * bg,
	      const char * dir);

int endfor (FILE * fin, FILE * fout)
{
  yyin = fin;
  yyout = fout;
  line = 1, scope = para = 0;
  inforeach = foreachscope = foreachpara = 
    inforeach_boundary = inforeach_face = 0;
  invardecl = 0;
  inval = invalpara = indef = 0;
  brack = inarray = infine = inmalloc = 0;
  inevent = inreturn = inattr = inmap = 0;
  mallocpara = 0;
  foreachdim_stack.n = 0;
  foreach_child_stack.n = 0;
  inblock_stack.n = 0;
  inboundary = nboundary = nsetboundary = 0;
  infunction = 0;
  infunctionproto = inmain = intrace = traceon = 0;
  if (tracefp)
    fclose (tracefp);
  tracefp = NULL;
  inarg = incadna = 0;
  boundaryheader = dopen ("_boundary.h", "w");
  int ret = yylex();
  if (!ret) {
    if (scope > 0)
      ret = yyerror ("mismatched '{'");
    else if (para > 0)
      ret = yyerror ("mismatched '('");
  }
  fclose (boundaryheader);
  return ret;
}

FILE * writepath (char * path, const char * mode)
{
  char * s = path;
  while ((s = strchr (s, '/'))) {
    *s = '\0';
    if (access (path, R_OK|W_OK|X_OK) && mkdir (path, 0700))
      return NULL;
    *s++ = '/';
  }
  return fopen (path, mode);
}

void cleanup (int status, const char * dir)
{
  if (!debug && dir) {
    char command[80] = "rm -r -f ";
    strcat (command, dir);
    int s = system (command); s = s;
  }
  exit (status);
}

FILE * dopen (const char * fname, const char * mode)
{
  char * out = malloc (strlen (dir) + strlen (fname) + 2);
  strcpy (out, dir); strcat (out, "/"); strcat (out, fname);
  FILE * fout = fopen (out, mode);
  free (out);
  return fout;
}

void write_event (int i, FILE * fout)
{
  char * func = eventfunc[i];
  fprintf (fout, "  event_register ((Event){ 0, %d, %s, {", nexpr[i], func);
  int j;
  for (j = 0; j < nexpr[i] - 1; j++)
    fprintf (fout, "%s_expr%d, ", func, j);
  fprintf (fout, "%s_expr%d}, ", func, j);
  if (eventarray[i] == 'i')
    fprintf (fout, "%s_array, ", func);
  else
    fprintf (fout, "((int *)0), ");
  if (eventarray[i] == 't')
    fprintf (fout, "%s_array,\n", func);
  else
    fprintf (fout, "((%s *)0),\n", doubletype);
  fprintf (fout, "    \"%s\", %d, \"%s\"});\n", nolineno ? "" : eventfile[i], 
	   nolineno ? 0 : eventline[i], eventid[i]);
}

static const char * typestr (var_t var) {
  return (var.type == scalar ? "scalar" :
	  var.type == vector ? "vector" :
	  var.type == tensor ? "tensor" : "");
}

void compdir (FILE * fin, FILE * fout, FILE * swigfp, 
	      char * swigname,
	      char * grid)
{
  if (endfor (fin, fout))
    cleanup (1, dir);
  fclose (fout);

  fout = dopen ("_boundarydecl.h", "w");
  int i;
  for (i = 0; i < nboundary; i++)
    if (!periodic[i])
      fprintf (fout, 
	       "static %s _boundary%d (Point point, Point neighbor,"
	       " scalar _s, void * data);\n"
	       "static %s _boundary%d_homogeneous (Point point,"
	       " Point neighbor, scalar _s, void * data);\n", 
	       doubletype, i, doubletype, i);
  fclose (fout);

  fout = dopen ("_grid.h", "w");
  /* new variables */
  fprintf (fout,
	   "size_t datasize = %d*sizeof (%s);\n",
	   nvar, doubletype);
  /* attributes */
  FILE * fp = dopen ("_attributes.h", "a");
  fputs ("\n"
	 "} _Attributes;\n"
	 "_Attributes * _attribute;\n",
	 fp);
  fclose (fp);
  /* event functions */
  for (i = 0; i < nevents; i++) {
    char * id = eventfunc[i];
    fprintf (fout, 
	     "static int %s (const int i, const %s t, Event * _ev);\n",
	     id, doubletype);
    int j;
    for (j = 0; j < nexpr[i]; j++)
      fprintf (fout,
	       "static int %s_expr%d (int * ip, %s * tp, Event * _ev);\n",
	       id, j, doubletype);
    if (eventarray[i])
      fprintf (fout, "static %s %s_array[] = %s,-1};\n", 
	       eventarray[i] == 'i' ? "int" : doubletype, id,
	       eventarray_elems[i]);
  }
  /* boundaries */
  for (i = 0; i < nsetboundary; i++)
    fprintf (fout, "static void _set_boundary%d (void);\n", 
	     boundaryindex[i]);
  /* _init_solver() */
  fputs ("void _init_solver (void) {\n"
	 "  void init_solver();\n"
	 "  init_solver();\n", fout);
  /* events */
  fprintf (fout,
	   "  Events = (Event *) pmalloc (sizeof (Event), __func__, __FILE__, %s);\n"
	   "  Events[0].last = 1;\n", nolineno ? "0" : "__LINE__");
  int last;
  for (last = 0; last <= 1; last++)
    for (i = 0; i < nevents; i++)
      if (eventchild[i] < 0 && eventlast[i] == last) {
	int j = i;
	while (j >= 0) {
	  write_event (j, fout);
	  j = eventparent[j];
	}
      }
  /* scalar attributes */
  fprintf (fout, "  _attribute = (_Attributes *) pcalloc (datasize/sizeof(%s), "
	   "sizeof (_Attributes), __func__, __FILE__, %s);\n",
	   doubletype, nolineno ? "0" : "__LINE__");
  /* list of all scalars */
  fprintf (fout, 
	   "  all = (scalar *) pmalloc (sizeof (scalar)*%d,__func__, __FILE__, %s);\n"
	   "  for (int i = 0; i < %d; i++)\n"
	   "    all[i].i = i;\n"
	   "  all[%d].i = -1;\n"
	   "  set_fpe();\n",
	   nvar + 1, nolineno ? "0" : "__LINE__", nvar, nvar);
  if (catch)
    fputs ("  catch_fpe();\n", fout);
  if (progress)
    fputs ("  last_events();\n", fout);
  fprintf (fout, "  %s_methods();\n", grid);
  for (i = varstack; i >= 0; i--) {
    var_t var = _varstack[i];
    if (var.i[0] >= 0) {
      if (var.constant) {
	// global constants
	if (var.type == scalar)
	  fprintf (fout, 
		   "  init_const_scalar ((scalar){_NVARMAX+%d}, \"%s\", %s);\n",
		   var.i[0], var.v, var.constant);
	else if (var.type == vector) {
	  int i;
	  fprintf (fout, "  init_const_vector ((vector){{_NVARMAX+%d}",
		   var.i[0]);
	  for (i = 1; i < dimension; i++)
	    fprintf (fout, ",{_NVARMAX+%d}", var.i[i]);
	  fprintf (fout, "}, \"%s\", (%s [])%s);\n",
		   var.v, doubletype, var.constant);
	}
	else
	  assert (0);
      }
      // global variables
      else if (var.type == scalar)
	fprintf (fout, "  init_scalar ((scalar){%d}, \"%s\");\n",
		 var.i[0], var.v);
      else if (var.type == vector) {
	fprintf (fout, "  init_%svector ((vector){{%d}",
		 var.face ? "face_" : 
		 var.vertex ? "vertex_" : 
		 "",
		 var.i[0]);
	int i;
	for (i = 1; i < dimension; i++)
	  fprintf (fout, ",{%d}", var.i[i]);
	fprintf (fout, "}, \"%s\");\n", var.v);
      }
      else if (var.type == tensor) {
	fprintf (fout, "  init_tensor ((tensor){");
	int i, j, k = 0;
	for (i = 0; i < dimension; i++) {
	  fprintf (fout, "{{%d}", var.i[k++]);
	  for (j = 1; j < dimension; j++)
	    fprintf (fout, ",{%d}", var.i[k++]);
	  fprintf (fout, "}");
	  if (i < dimension - 1)
	    fputc (',', fout);
	}
	fprintf (fout, "}, \"%s\");\n", var.v);
      }
      else if (var.type == bid)
	fprintf (fout, "  %s = new_bid();\n", var.v);
      else
	assert (0);
    }
  }
  for (i = 0; i < nsetboundary; i++)
    fprintf (fout, "  _set_boundary%d();\n", boundaryindex[i]);
  fputs ("}\n", fout);
  fclose (fout);

  if (source && autolinks && autolink)
    printf ("%s\n", autolink);
  
  /* SWIG interface */
  if (swigfp) {
    int i;
    fputs ("\n%{\n", swigfp);
    for (i = varstack; i >= 0; i--) {
      var_t var = _varstack[i];
      if (var.i[0] >= 0 && !var.constant)
	fprintf (swigfp, "  extern %s %s;\n", typestr (var), var.v);
    }
    fputs ("%}\n\n", swigfp);
    for (i = varstack; i >= 0; i--) {
      var_t var = _varstack[i];
      if (var.i[0] >= 0 && !var.constant)
	fprintf (swigfp, "extern %s %s;\n", typestr (var), var.v);
    }
    fputs ("\n%pythoncode %{\n", swigfp);
    for (i = varstack; i >= 0; i--) {
      var_t var = _varstack[i];
      if (var.i[0] >= 0 && !var.constant)
	fprintf (swigfp, "%s = %s(_%s.cvar.%s)\n",
		 var.v, typestr (var), swigname, var.v);
    }
    fputs ("%}\n", swigfp);
    fclose (swigfp);
  }
}

int main (int argc, char ** argv)
{
  char * cc = getenv ("CC99"), command[1000], command1[1000] = "";
  if (cc == NULL)
    strcpy (command, CC99);
  else
    strcpy (command, cc);
  char * file = NULL;
  int i, dep = 0, tags = 0, swig = 0;
  for (i = 1; i < argc; i++) {
    if (!strncmp (argv[i], "-grid=", 6))
      ;
    else if (!strcmp (argv[i], "-MD"))
      dep = 1;
    else if (!strcmp (argv[i], "-tags"))
      tags = 1;
    else if (!strcmp (argv[i], "-python"))
      swig = 1;
    else if (!strcmp (argv[i], "-debug"))
      debug = 1;
    else if (!strcmp (argv[i], "-events"))
      events = 1;
    else if (!strcmp (argv[i], "-catch"))
      catch = 1;
    else if (!strcmp (argv[i], "-source"))
      source = 1;
    else if (!strcmp (argv[i], "-autolink"))
      autolinks = 1;
    else if (!strcmp (argv[i], "-progress"))
      progress = 1;
    else if (!strcmp (argv[i], "-Wall")) {
      char * s = strchr (command, ' ');
      if (s) {
	char command1[1000];
	strcpy (command1, s);
	*(s+1) = '\0';
	strcat (command, argv[i]);
	strcat (command, command1);
      }
      else
	strcat (command, argv[i]);
    }
    else if (!strcmp (argv[i], "-cadna")) {
      cadna = 1;
      char * cc = getenv ("CADNACC");
      if (cc == NULL)
	strcpy (command, CADNACC);
      else
	strcpy (command, cc);
    }
    else if (!strncmp (argv[i], "-Ddimension=", 12))
      dimension = 1 + argv[i][12] - '1';
    else if (catch && !strncmp (argv[i], "-O", 2))
      ;
    else if (!strcmp (argv[i], "-nolineno")) {
      nolineno = 1;
      strcat (command1, " -D'assert(x)=((void)(x))'");
    }
    else if (!strcmp (argv[i], "-o")) {
      strcat (command1, " ");
      strcat (command1, argv[i++]);
      if (i < argc) {
	strcat (command1, " ");
	strcat (command1, argv[i]);
      }
    }
    else if (argv[i][0] != '-' && 
	     (tags || !strcmp (&argv[i][strlen(argv[i]) - 2], ".c"))) {
      if (file) {
	fprintf (stderr, "usage: qcc -grid=[GRID] [OPTIONS] FILE.c\n");
	return 1;
      }
      file = argv[i];
    }
    else if (!file) { 
      strcat (command, " ");
      strcat (command, argv[i]);
    }
    else {
      strcat (command1, " ");
      strcat (command1, argv[i]);
    }
  }
  char * openmp = strstr (command, "-fopenmp");
  if (openmp) {
    if (strstr (command, "-D_MPI")) {
      fprintf (stderr,
	       "qcc: warning: OpenMP cannot be used with MPI (yet): "
	       "switching it off\n");
      int i;
      for (i = 0; i < strlen("-fopenmp"); i++)
	openmp[i] = ' ';
    }
    else if (swig) {
      fprintf (stderr,
	       "qcc: warning: OpenMP cannot be used with Python (yet): "
	       "switching it off\n");
      int i;
      for (i = 0; i < strlen("-fopenmp"); i++)
	openmp[i] = ' ';
    }
  }
  int status;
  if (debug) {
    status = system ("rm -r -f .qcc");
    strcpy (dir, ".qcc");
    status = mkdir (dir, 0700);
  }
  else
    status = (mkdtemp (dir) == NULL);
  if (status) {
    perror (dir);
    return 1;
  }
  if (file) {
    char * out[100], * grid = NULL;
    int default_grid;
    includes (argc, argv, out, &grid, &default_grid, &dimension, &bghosts,
	      dep || tags ? NULL : dir);
    FILE * swigfp = NULL;
    char swigname[80] = "";
    if (swig) {
      strcpy (swigname, file);
      char * dot = strchr (swigname, '.');
      *dot = '\0'; strcat (swigname, ".i");
      swigfp = fopen (swigname, "a");
      if (!swigfp) {
	fprintf (stderr, "qcc: could not open '%s': ", swigname);
	return 1;
      }
      *dot = '\0';
    }
    if (!dep && !tags) {
      char * basename = strdup (file), * ext = basename;
      while (*ext != '\0' && *ext != '.') ext++;
      char * cpp = malloc (strlen(basename) + strlen("-cpp") + strlen(ext) + 1);
      if (*ext == '.') {
	*ext = '\0';
	strcpy (cpp, basename);
	strcat (cpp, "-cpp");
	*ext = '.';
	strcat (cpp, ext);
      }
      else {
	strcpy (cpp, basename);
	strcat (cpp, "-cpp");
      }
      FILE * fin = dopen (file, "r");
      if (!fin) {
	perror (file);
	cleanup (1, dir);
      }
      FILE * fp = dopen ("_attributes.h", "w");
      fputs ("typedef struct {\n", fp);
      fclose (fp);
      fp = dopen ("_maps.h", "w");
      fclose (fp);
      FILE * fout = dopen (cpp, "w");
      if (swig)
	fputs ("@include <Python.h>\n", fout);
      fputs ("@if _XOPEN_SOURCE < 700\n"
	     "  @undef _XOPEN_SOURCE\n"
	     "  @define _XOPEN_SOURCE 700\n"
	     "@endif\n"
	     "@if _GNU_SOURCE\n"
	     "@include <stdint.h>\n"
	     "@include <string.h>\n"
	     "@include <fenv.h>\n"
	     "@endif\n",
	     fout);
      if (catch)
	fputs ("#define TRASH 1\n"
	       "#define _CATCH last_point = point;\n", 
	       fout);
      else
	fputs ("#define _CATCH\n", fout);
      fprintf (fout, "#define dimension %d\n", dimension);
      if (bghosts)
	fprintf (fout, "#define BGHOSTS %d\n", bghosts);
      fputs ("#include \"common.h\"\n", fout);
      /* catch */
      if (catch)
	fputs ("void catch_fpe (void);\n", fout);
      /* grid */
      if (default_grid)
	fprintf (fout, "#include \"grid/%s.h\"\n", grid);
      /* declaration of boundary condition fonctions */
      fputs ("@include \"_boundarydecl.h\"\n", fout);
      char s[81];
      while (fgets (s, 81, fin)) {
	// replace '#include <' with '@include <'
	char * s1 = s; while (strchr (" \t", *s1)) s1++;
	if (*s1 == '#') {
	  char *s2 = s1 + 1; while (strchr (" \t", *s2)) s2++;
	  if (!strncmp (s2, "include", 7)) {
	    while (!strchr (" \t", *s2)) s2++;
	    while (strchr (" \t", *s2)) s2++;
	    if (*s2 == '<')
	      *s1 = '@';
	  }
	}
	fputs (s, fout);
      }
      if (swigfp)
	fputs ("#include \"python.h\"\n", fout);
      if (progress)
	fputs ("#include \"grid/progress.h\"\n", fout);
      fclose (fout);
      fclose (fin);
      fout = dopen (file, "w");
      if (!fout) {
	perror (file);
	cleanup (1, dir);
      }

      char preproc[1000], * cppcommand = getenv ("CPP99");
      strcpy (preproc, "cd ");
      strcat (preproc, dir);
      strcat (preproc, " && ");
      if (!cppcommand && strcmp (CPP99, ""))
	cppcommand = CPP99;
      if (!cppcommand) {
	if (source) {
	  /* remove -D_GNU_SOURCE flags from $CC99 */
	  char tmp[1000]; strcpy (tmp, command);
	  char * s = strtok (tmp, " \t");
	  while (s) {
	    if (strncmp (s, "-D_GNU_SOURCE", 13)) {
	      strcat (preproc, s);
	      strcat (preproc, " ");
	    }
	    s = strtok (NULL, " \t");
	  }
	}
	else
	  strcat (preproc, command);
	strcat (preproc, " -E");
      }
      else
	strcat (preproc, cppcommand);
      strcat (preproc, " -I. -I");
      strcat (preproc, LIBDIR);
      strcat (preproc, " ");
      if (events) {
	strcat (preproc, " -DDEBUG_EVENTS=1 -DBASILISK=\"\\\"");
	strcat (preproc, BASILISK);
	strcat (preproc, "\\\"\" ");
      }
      strcat (preproc, cpp);
      if (debug) {
	fprintf (stderr, "preproc: %s\n", preproc);
	strcat (preproc, " | tee _preproc.c");
      }

      fin = popen (preproc, "r");
      if (!fin) {
	perror (preproc);
	cleanup (1, dir);
      }

      compdir (fin, fout, swigfp, swigname, grid);
      int status = pclose (fin);
      if (status == -1 ||
	  (WIFSIGNALED (status) && (WTERMSIG (status) == SIGINT || 
				    WTERMSIG (status) == SIGQUIT)))
	cleanup (1, dir);

      fout = dopen ("_tmp", "w");

      fin = dopen (file, "r");
      char line[1024];
      int c, lineno = 1;
      char * fname = NULL;
      // includes _attributes.h
      while (fgets (line, 1024, fin)) {
	if (!strcmp (line, "#include \"_attributes.h\"\n"))
	  break;
	else {
	  if (!strncmp (line, "#line ", 6)) {
	    lineno = atoi (&line[6]) - 1;
	    char * s = &line[6]; space (s); nonspace(s);
	    free (fname);
	    fname = strdup (s);
	  }
	  fputs (line, fout);
	}
	lineno++;
      }
      fp = dopen ("_attributes.h", "r");
      while ((c = fgetc (fp)) != EOF)
	fputc (c, fout);
      fclose (fp);
      fprintf (fout, "#line %d %s\n", lineno, fname);
      free (fname);
      // includes _boundarydecl.h
      while (fgets (line, 1024, fin))
	if (!strcmp (line, "#include \"_boundarydecl.h\"\n"))
	  break;
	else
	  fputs (line, fout);
      fp = dopen ("_boundarydecl.h", "r");
      while ((c = fgetc (fp)) != EOF)
	fputc (c, fout);
      fclose (fp);
      // rest of the file
      while (fgets (line, 1024, fin))
	fputs (line, fout);
      fclose (fin);

      fin = dopen ("_boundary.h", "r");
      while ((c = fgetc (fin)) != EOF)
	fputc (c, fout);
      fclose (fin);

      fin = dopen ("_grid.h", "r");
      while ((c = fgetc (fin)) != EOF)
	fputc (c, fout);
      fclose (fin);

      fclose (fout);

      char src[80], dst[80];
      strcpy (src, dir); strcat (src, "/_tmp");
      if (source) {
	strcpy (dst, "_");
      }
      else {
	strcpy (dst, dir); strcat (dst, "/");
      }
      strcat (dst, file);
      rename (src, dst);

      strcat (command, " -I");
      strcat (command, LIBDIR);
      strcat (command, " ");
      strcat (command, dir);
      strcat (command, "/");
      strcat (command, file); 
      strcat (command, command1);
    }
  }
  else if (dep || tags) {
    fprintf (stderr, "usage: qcc -grid=[GRID] [OPTIONS] FILE.c\n");
    cleanup (1, dir);
  }
  else
    strcat (command, command1);
  /* compilation */
  if (!dep && !tags && !source) {
    if (autolinks && autolink)
      strcat (command, autolink);
    if (debug)
      fprintf (stderr, "command: %s\n", command);
    status = system (command);
    if (status == -1 ||
	(WIFSIGNALED (status) && (WTERMSIG (status) == SIGINT || 
				  WTERMSIG (status) == SIGQUIT)))
      cleanup (1, dir);
    cleanup (WEXITSTATUS (status), dir);
  }
  cleanup (0, dir);
  return 0;
}
