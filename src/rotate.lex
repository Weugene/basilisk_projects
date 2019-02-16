%option noyywrap
%option yylineno
%{
  static int n = 0, dimension = 2, para = 0, scope = 0;

  static int identifier (int c) {
    return ((c >= 'a' && c <= 'z') || 
	    (c >= 'A' && c <= 'Z') || 
	    (c >= '0' && c <= '9') || c == '_');
  }

  static void rotate_string (char * string) {
    char s[] = "123", * j = string;
    int i = 0, len = strlen (string), k;
    for (k = 0; k <= len + 1; k++, j++) {
      if (i < 3) {
	s[i++] = k <= len ? *j : '\0'; s[i] = '\0';
      }
      else {
	if ((s[0] == '.' || s[0] == '_') && 
	    s[1] >= 'x' && s[1] <= 'z' &&
	    !identifier(s[2]))
	  s[1] = 'x' + (s[1] + n - 'x') % dimension;
	fputc (s[0], yyout);
	s[0] = s[1]; s[1] = s[2]; 
	s[2] = s[1] != '\0' ? *j : '\0';
      }
    }
    fputs (s, yyout);
  }

%}

ID     [a-zA-Z0-9_]
SP     [ \t]
WS     [ \t\v\n\f]
ES     (\\([\'\"\?\\abfnrtv]|[0-7]{1,3}|x[a-fA-F0-9]+))

%%

\( { para++; scope = 0; ECHO; }
\) { para--; ECHO; }
\{ { scope++; ECHO; }
\} { scope--; ECHO; }

_attribute\[{ID}+.[ntr].i\] {
  // replace tangential and normal components
  char * s = strstr (yytext, ".n");
  if (s)
    s[1] = 'x';
  else if ((s = strstr (yytext, ".t")))
    s[1] = 'y';
  else if ((s = strstr (yytext, ".r")))
    s[1] = 'z';
  ECHO;
}

({ID}+\.)+x{WS}*,{WS}*({ID}+\.)+y |
({ID}+\.)+x{WS}*,{WS}*({ID}+\.)+y,{WS}*({ID}+\.)+z {
  // do not rotate lists of vector components
  char * id1 = strdup (yytext);
  char * s = strchr (id1, ','), * id2 = s + 1; 
  do *s-- = '\0'; while (strchr(" \t\v\n\f",*s));
  while (strchr(" \t\v\n\f",*id2)) id2++;
  s = strchr (id2, ',');
  char * id3 = s ? s + 1 : NULL;
  if (id3) {
    do *s-- = '\0'; while (strchr(" \t\v\n\f",*s));
    while (strchr(" \t\v\n\f",*id3)) id3++;
  }
  if (scope == 0 ||
      strlen(id1) != strlen(id2) || strncmp (id1, id2, strlen(id1) - 2) ||
      (id3 && (strlen(id1) != strlen(id3) ||
	       strncmp (id1, id3, strlen(id1) - 2))))
    rotate_string(yytext);
  else
    ECHO;
  free (id1);
}

{ID}*_[xyz] |
[.][xyz]  rotate_string (yytext);

{ID}*->[xyz] {
  int i = strlen(yytext) - 1;
  yytext[i] = 'x' + (yytext[i] + n - 'x') % dimension;
  ECHO;
}

left {
  char * s[3] = {"left", "bottom", "back"};
  fputs (s[n % dimension], yyout);
}
right {
  char * s[3] = {"right", "top", "front"};
  fputs (s[n % dimension], yyout);
}
top {
  char * s[3] = {"top", "front", "right"};
  fputs (s[n % dimension], yyout);
}
bottom {
  char * s[3] = {"bottom", "back", "left"};
  fputs (s[n % dimension], yyout);
}
front {
  char * s[3] = {"front", "right", "top"};
  fputs (s[n % dimension], yyout);
}
back {
  char * s[3] = {"back", "left", "bottom"};
  fputs (s[n % dimension], yyout);
}

{ID}+   ECHO;

val_{ID}*{WS}*\( |
(val|fine|coarse|allocated|neighbor|neighborp|aparent){WS}*\( {
  int para = 1;
  char * index[5];
  int len[5], i = 0, c = input(), j;
  for (j = 0; j < 5; j++) {
    index[j] = malloc (sizeof(char));
    index[j][0] = '\0';
    len[j] = 1;
  }
  while ((para > 1 || c != ')') && c != EOF) {
    if (c == '(') para++;
    else if (c == ')') para--;
    if (c == ',' && para == 1 && i < 5)
      i++;
    else {
      index[i] = realloc (index[i], sizeof(char)*++len[i]);
      index[i][len[i]-2] = c;
      index[i][len[i]-1] = '\0';
    }
    c = input();
  }
  int start = 0;
  if (!strncmp (yytext, "val", 3) || 
      !strncmp (yytext, "fine", 4) || 
      !strncmp (yytext, "coarse", 6)) {
    if (i == 3) {
      rotate_string (yytext);
      rotate_string (index[0]);
    }
    else {
      ECHO;
      fputs (index[0], yyout);
    }
    fputc (',', yyout);
    start = 1;
  }
  else
    ECHO;
  int ghost = 0;
  for (j = start; j <= i && !ghost; j++)
    if (strlen(index[j]) == 2 && 
	strchr ("ijk", index[j][0]) && index[j][1] == 'g')
      ghost = 1;
  if (!ghost && i == 2 + start) {
    for (j = 0; j < dimension; j++) {
      int k = (j + dimension - n) % dimension;
      rotate_string (index[k + start]);
      fputc (j < 2 ? ',' : ')', yyout);
    }
    for (j = dimension + 1; j <= 3; j++) {
      rotate_string (index[j - 1 + start]);
      fputc (j < 3 ? ',' : ')', yyout);
    }
  }
  else
    //  more than 3 indices or ghost indices: do not rotate
    for (j = start; j <= i; j++) {
      fputs (index[j], yyout);
      fputc (j < i ? ',' : ')', yyout);
    }
  for (j = 0; j < 5; j++)
    free (index[j]);
}

"//".*            { ECHO; /* consume //-comment */ }
.                   ECHO;
[\n]                ECHO;
({SP}?\"([^\"\\\n]|{ES})*\"{WS}*)+  ECHO; /* STRING_LITERAL */

%%

// rotate dimensions n times
int rotate (FILE * fin, FILE * fout, int nrotate, int dim)
{
  if (0) yyunput (0, NULL); // just prevents 'yyunput unused' compiler warning
  rewind (fin);
  yyin = fin;
  yyout = fout;
  yylineno = 1;
  n = nrotate;
  dimension = dim;
  para = scope = 0;
  return yylex();
}

#if TEST
int main (int argc, char * argv[])
{
  rotate (stdin, stdout, atoi (argv[1]), atoi (argv[2]));
  return 0;
}
#endif
