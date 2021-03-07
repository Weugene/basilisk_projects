/* Conclusions:
 * fraction.h works if dimension>=2
 * the results worst to best:
 *  1. not set
 *  2. f.prolongation=refine_injection and 
 *	3. f.prolongation=fraction_refine; // probably refines if neighbors are changing!!! Do not use for smooth fields
 *	the result of comparison can be seen in adapt_comparison.png
 *	using of f.refine does not change the result. prolongation operator is more important.
 *	When the error is estimated for all grid cells they are marked to be either too fine, too coarse or just fine. All cells that are marked to be too coarse will be refined, also cells that require refinement in order to keep the levels at resolution boundaries differ a single level are refined. Coarsening is only done when all children of a single parent cell are marked to be too fine and coarsening those cells does not violate the mentioned requirement.

The techniques used for defing field values on the newly initialized grid cells are not always identical to the restriction and prolongation operation. For this purpose all fields defined on a tree grid also have a refinerefine and coarsencoarsen attributes to facilitate distiction between the methods used for prolongationprolongation and restrictionrestriction. The default for scalar fields is (as found in the common tree-grid header file):

...
s.refine = s.prolongation;
...
Furthermore, by default for a scalar field, the coarsening attribute is not even defined and the restriction operation is used instead when cells are coarsened.
 * */
//#include "grid/bitree.h"
#include "grid/tree.h"
#include "utils.h"
#include "fractions.h"
scalar f[]; 
int main(){
  FILE * fp0 = fopen("original","w");
  FILE * fp1 = fopen("adapted_no_set","w");
  FILE * fp2 = fopen("adapted_inj","w");
  FILE * fp3 = fopen("adapted_volume_frac","w");
  FILE * fp;

  X0=0;
  L0=1;

	for (int i=0; i<3; i++){
		init_grid(256);
		foreach() f[]= sin(((x)*M_PI*5)-11.5)* exp(-(sq((x-0.5)/0.15)));
		if (i==0){
			foreach() {
				if (fabs(y-0.501953) <1e-5){
					fprintf(fp0,"%g\t%g\n",x,f[]);
				}
			}
			fp = fp1;
		}
		if (i==1){
			f.prolongation=refine_injection;
			//f.refine = refine_injection;
			fp = fp2;
		}
		else if(i==2){
			f.prolongation=fraction_refine;
			//f.refine = fraction_refine;
			fp = fp3;
		}
		while(adapt_wavelet({f},(double[]){0.1},8).nc){
   	 		foreach()
      			f[]= sin(((x)*M_PI*5)-11.5)* exp(-(sq((x-0.5)/0.15)));
    		boundary({f});
		}
  		double xp=L0/512;
		while(xp<L0){
    		Point point = locate(xp, 0.5);
    		fprintf(fp,"%g\t%g\n",xp,f[]);
    		xp+= L0/256.;
 		}
	}
}


/* for gnuplot
 *
 * plot 'original' u 1:2 t 'original', 'adapted_no_set' u 1:2 t 'adapted no set', 'adapted_inj' u 1:2 t 'adapted inj', 'adapted_volume_frac' u 1:2 t 'adapted volume frac'
*/
