#include "grid/cartesian1D.h"
#include "saint-venant.h"

double a, tmax;
scalar dh[];

int main(){
    X0 = -15.;
    L0 = 60.;
    N = 1024;
    G = 9.81;
    tmax = 50; // temps final
    a = 0.1; // amplitude de la vague
    run();
}

u.n[left] = neumann(0);
h[left] = neumann(0);

event init (i = 0)
{
    foreach(){
        zb[] =   (x>10)*(x-10)/(25);
        u.x[]=a*exp(-(x+12)*(x+12)) ; // pulse de vitesse pour créer la vague
        h[]=fmax(1+u.x[]-zb[],0);
    }
    boundary({zb,h,u});
}

event plot (t<tmax;t+=0.1) {
    printf("set title 'Vague sur une plage 1D ----- t= %.1lf '\n"
    "p[%g:%g][-0.1:1.5]  '-' u 1:($2+$4) t'free surface (m)' w l lt 3,"
    "'' u 1:3 t'velocity (m/s)' w l lt 4,"
    "'' u 1:4 t'topo' w l lt 1\n",t,X0+1,X0+L0);
    foreach()
        printf (" %g %g %g %g %g \n", x, h[], u.x[], zb[], t);
    printf ("e\n\n");
}


/**
## La compilation

Nous pouvons compiler à la main grâce au compilateur qcc qui a les
mêmes options que le compilateur gcc (+ quelques ajouts). Par exemple,
nous avons vu précédemment que l'option -events permet de garder une
trace de l'exécution des events dans la deuxième sortie du programme
(log). Compilez votre code grâce à la commande suivante dans votre
terminal :

~~~bash
qcc CM1.c -o CM1.x -lm -O2 -Wall
~~~

Si tout se passe bien, vous devriez pouvoir ensuite exécuter votre
code à l'aide de la commande suivante (en redirigeant la sortie
principale dans le fichiers "out" ):

~~~bash
./CM1.x > out
~~~

Vous pouvez voir le résultat de votre première simulation dans gnuplot
grâce à la commande :

~~~bash
load 'out'
~~~

## L'utilisation de Make

La commande "make code.tst" permet de compiler, d'executer et de
rediriger les fichiers de sortie en une seule commande. Copiez le
fichier Makefile se trouvant dans ~/basilisk/src/ dans votre
répertoire courant, puis testez make en executant la commande :

~~~bash
make CM1.tst
~~~

Si tout c'est bien passé, vous devriez maintenant avoir un dossier
"CM1" dans lequel vous allez trouver un fichier "out" et un fichier
"log", qui sont les deux sorties de votre programme. Vous pouvez faire
un "load 'out'" dans gnuplot pour vous assurez que le fichier "out"
est le même que précédemment. Pour se convaincre de l'utilité de la
commande make, relancez la commande :

~~~bash
make CM1.tst
~~~

Que vous a répondu le terminal ?

En effet, la commande make teste si le programme a changé depuis la
dernière fois qu'il a été executé. Ce qui est très utile en pratique :
on ne compte plus le nombre d'étudiants se plaignant d'un "bug" dans
leur programme alors qu'ils ne l'avaient simplement pas compilé et
qu'ils executaient une version antérieure de leur programme (qui,
elle, était buggée). Dans le doute : passez toujours par la commande
make  : si rien n'a changé, vous ne perdrez pas le temps d'une nouvelle
compilation/execution !.

## Lier les sorties avec votre programme

La commande make permet également de sortir des graphiques gnuplot grâce à la commande "make nom/plots". Ajoutez le code suivant à la fin de votre programme en commentaire et en remplaçant les "+" par des tildes :"~" :



Maintenant, testez la commande make /plots :

~~~bash
make CM1/plots
~~~

Vous pouvez admirer le résultat en ouvrant le fichier movie.gif dans
un explorateur internet de votre choix.  Voici le résultat chez moi :

~~~gnuplot Animation of the free surface.
reset
set xlabel 'X'
set ylabel 'Z'
set term gif animate
set output 'movie.gif'
load './out'
~~~
 */