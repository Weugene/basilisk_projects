typedef struct {
    coord x, y, z;
} ttt;

#if dimension == 1
#define scalar_a_by_b(a, b) (a.x[]*b.x[])
#elif dimension == 2
#define scalar_a_by_b(a, b) (a.x[]*b.x[] + a.y[]*b.y[])
#else // dimension == 3
#define scalar_a_by_b(a, b) (a.x[]*b.x[] + a.y[]*b.y[] + a.z[]*b.z[])
#endif

int main(){
	init_grid(8);
	vector u[], n[], rm[];
	foreach(){

			rm.x[] = scalar_a_by_b(u, unityf);

	}
return 0;
}
