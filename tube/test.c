//The output of the following example program is approximately that of
//       "ls [a-c]*.c".

       #include <stdio.h>
       #include <stdlib.h>
       #include <wordexp.h>

       int
       main(int argc, char **argv)
       {
           wordexp_t p;
           char **w;
           int i;

           wordexp("dump-*", &p, 0);
           w = p.we_wordv;
           for (i = 0; i < p.we_wordc; i++)
               printf("%s\n", w[i]);
           wordfree(&p);
           exit(EXIT_SUCCESS);
       }
