#ifndef NOT_ZERO
    #define NOT_ZERO 1.e-30
#endif

void MinMaxValues(scalar * list, double * arr_eps) {// for each scalar min and max
    double arr[10][2], small_val = 1e-10;
    int ilist = 0;
    for (scalar s in list) {
        double mina= HUGE, maxa= -HUGE;
        foreach( reduction(min:mina) reduction(max:maxa) ){
            if (fabs(s[]) < mina) mina = fabs(s[]);
            if (fabs(s[]) > maxa) maxa = fabs(s[]);
        }
#if _MPI
        MPI_Allreduce (MPI_IN_PLACE, &mina, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce (MPI_IN_PLACE, &maxa, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif
        arr[ilist][0] = mina;
        arr[ilist][1] = maxa;
        ilist++;
//        fprintf(stderr, "arr for i=%d", ilist);
    }
    int i = 0;
    for (scalar s in list){

#if EPS_MAXA == 1
        if (arr[i][1] > small_val){
            arr_eps[i] *=arr[i][1];
        }
#elif EPS_MAXA == 2
        if (arr[i][1] - arr[i][0] > small_val){
            arr_eps[i] *= arr[i][1] - arr[i][0];
        }
#else
        if (0.5*(arr[i][0] + arr[i][1]) > small_val){
            arr_eps[i] *= 0.5*(arr[i][0] + arr[i][1]);
        }
#endif
        else{
            arr_eps[i] *= 1;
        }
#ifdef DEBUG_MINMAXVALUES
        fprintf(stderr, "MinMaxValues: name=%s, min=%g, max=%g, eps=%g\n", s.name, arr[i][0], arr[i][1], arr_eps[i]);
#endif
        i++;
    }
}

int count_cells(double t, int i){
#if _MPI
    int rank, h_len;
    char hostname[MPI_MAX_PROCESSOR_NAME];
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Get_processor_name(hostname, &h_len);
    printf("i %d t %g hostname %s rank %d num cells %d total num cells %d compression rate %g\n", i, t, hostname, rank, grid->n, grid->tn, pow(2, dimension*grid->maxdepth)/grid->tn);
#else
    printf("i %d t %g total num cells %d\n", i, t, grid->tn);
#endif
    fflush(stdout);
    return grid->tn;
}
// statistical values inside cells with liquid
stats statsf_weugene (scalar f, scalar fs)
{
    double dvr, val, min = 1e100, max = -1e100, sum = 0., sum2 = 0., volume = 0.;
    foreach(reduction(+:sum) reduction(+:sum2) reduction(+:volume)
    reduction(max:max) reduction(min:min)){
        dvr = dv()*(1. - fs[]);
        val = f[]*(1. - fs[]);
        volume += dvr;
        sum    += f[]*dvr;
        sum2   += sq(f[])*dvr;
        if (val > max) max = val;
        if (val < min) min = val;
    }
	if (volume > 0.){
    	sum /= volume; sum2 /= volume;
	}
	sum2 -= sq(sum);
	fprintf(ferr, "***: %g %g\n", sum, sum2);
	stats s;
	s.min = min, s.max = max, s.sum = sum, s.volume = volume; //modified by Weugene
	s.stddev = sum2 > 0. ? sqrt(sum2) : 0.;
    return s;
}

// statistical values inside pure liquid
stats statsf_weugene2 (scalar f, scalar fs)
{
    double dvr, min = 1e100, max = -1e100, sum = 0., sum2 = 0., volume = 0.;
    foreach(reduction(+:sum) reduction(+:sum2) reduction(+:volume)
    reduction(max:max) reduction(min:min))
    if (fs[] == 0.) {
        dvr = dv()*(1. - fs[]);
        volume += dvr;
        sum    += dvr*f[];
        sum2   += dvr*sq(f[]);
        if (f[] > max) max = f[];
        if (f[] < min) min = f[];
    }
	if (volume > 0.){
    	sum /= volume; sum2 /= volume;
	}
	sum2 -= sq(sum);
    fprintf(ferr, "sum=%g\n", sum);
    stats s;
    s.min = min, s.max = max, s.sum = sum, s.volume = volume;
    return s;
}

norm normf_weugene (scalar f, scalar fs)
{
    double avg = 0., rms = 0., max = 0., volume = 0.;
    foreach(reduction(max:max) reduction(+:avg)
    reduction(+:rms) reduction(+:volume)){
        double dvr = dv()*(1. - fs[]);
        double v = fabs(f[])*(1. - fs[]);
        if (v > max) max = v;
        volume += dvr;
        avg    += dvr*v;
        rms    += dvr*sq(v);
    }
    norm n;
    n.avg = volume ? avg/volume : 0.;
    n.rms = volume ? sqrt(rms/volume) : 0.;
    n.max = max;
    n.volume = volume;
    return n;
}

double change_weugene (scalar s, scalar sn, scalar fs)
{
    double max = 0., ds;
    foreach(reduction(max:max)) {
        ds = fabs (s[] - sn[])*(1. - fs[]);
        if (ds > max)
            max = ds;
        sn[] = s[];
    }
    return max;
}
/**
 * Smoothing function
 * f insput scalar field
 * sf - output smoothed scalar field
 */
void filter_scalar(scalar f, scalar sf){
#if dimension <= 2
    foreach()
    sf[] = (4.*f[] +
            2.*(f[0,1] + f[0,-1] + f[1,0] + f[-1,0]) +
            f[-1,-1] + f[1,-1] + f[1,1] + f[-1,1])/16.;
#else // dimension == 3
    foreach()
        sf[] = (8.*f[] +
            4.*(f[-1] + f[1] + f[0,1] + f[0,-1] + f[0,0,1] + f[0,0,-1]) +
            2.*(f[-1,1] + f[-1,0,1] + f[-1,0,-1] + f[-1,-1] +
                f[0,1,1] + f[0,1,-1] + f[0,-1,1] + f[0,-1,-1] +
                f[1,1] + f[1,0,1] + f[1,-1] + f[1,0,-1]) +
            f[1,-1,1] + f[-1,1,1] + f[-1,1,-1] + f[1,1,1] +
            f[1,1,-1] + f[-1,-1,-1] + f[1,-1,-1] + f[-1,-1,1])/64.;
#endif
#if TREE
    sf.prolongation = refine_bilinear;
    boundary ({sf});
#endif
}

/**
Calculate scalar from face vector. compute viscosity in a cell
*/

void calc_scalar_from_face(const face vector vf, scalar vs){
    double vsum = 0;
    foreach() {
        vsum = 0;
        foreach_dimension() {
            vsum += vf.x[] + vf.x[1];
        }
        vs[] = vsum/(2.0*dimension);
    }
    boundary((scalar *){vs});
}

/**
A function to rescale normals so that they are unit vectors w.r.t. the
2-norm (by default, the 1-norm is adopted for efficiency purposes). */
coord normalize_coord(coord n){
    double nn = NOT_ZERO;
    foreach_dimension() nn += sq(n.x); // (sqrt(sq(nf.x[]) + sq(nf.y[])))
    nn = sqrt(nn);
    foreach_dimension() n.x /= nn;
    return n;
}

//coord normal (Point point, scalar c) {
//    coord n = mycs (point, c);
//    return normalize_coord(n);
//}
//
//coord normal_face (Point point, scalar f){
//    coord nf;
//    foreach_dimension() nf.x = 0.;
//    bool interface_plus = interfacial(point, f);
//    bool interface_minus = interfacial(neighborp(-1), f);
//    if (interface_minus || interface_plus) {
//        if (interface_plus) {
//            coord n = normal (point, f);
//            foreach_dimension() nf.x += n.x;
//        }
//        if (interface_minus) {
//            coord n = normal (neighborp(-1), f);
//            nf.x += n.x;
//            nf.y += n.y;
//            #if dimension > 2
//                nf.z += n.z;
//            #endif
//        }
////        double norm_nf = 0;
////        foreach_dimension() norm_nf += sq(nf.x); // (sqrt(sq(nf.x[]) + sq(nf.y[])))
////        norm_nf = sqrt(norm_nf);
////        foreach_dimension() nf.x /= norm_nf;
//        nf = normalize_coord(nf);
//    }
//    return nf;
//}
//
///**
// *
// * @param point
// * @param f - volume fraction of fluid $1$
// * @param fs - volume fraction of solid
// * @return  tau_w - is the tangential surface normal pointing along the surface, into the liquid.
// */
//coord tangential_wall_and_normal_CL_face (Point point, scalar f, scalar fs) {
//    coord n_f = normal_face (point, f);
//    coord n_fs = normal_face (point, fs);
//    double n_f_dot_n_fs = 0;
//    foreach_dimension() n_f_dot_n_fs += n_f.x*n_fs.x;
//    coord tau_w;
//    foreach_dimension() tau_w.x = n_f.x - n_f_dot_n_fs*n_fs.x;
//    tau_w = normalize_coord(tau_w);
//    return tau_w;
//}
//
//coord normal_face_correction (Point point, scalar f, scalar fs, double theta) {
//    coord n_fs = normal_face (point, fs);
//    coord tau_w = tangential_wall_and_normal_CL_face (point, f, fs);
//    coord n_cor;
//    foreach_dimension() n_cor.x = n_fs.x*cos(theta) - tau_w.x*sin(theta);
//    return n_cor;
//}
//
///**
//A function to compute 2-norm cell-centered/face-centered normals in every cell/face. */
//
//void compute_normal (scalar f, vector normal_vector) {
//    foreach() {
//        coord n = normal (point, f);
//        foreach_dimension()
//        normal_vector.x[] = n.x;
//    }
//    boundary((scalar*){normal_vector});
//}
//
//void compute_normal_face (scalar f, face vector normal_vector_face) {
//    foreach_face() {
//        coord n = normal_face (point, f);
//        foreach_dimension()
//        normal_vector_face.x[] = n.x;
//    }
//    boundary((scalar*){normal_vector_face});
//}
//
//void compute_tangential_wall_and_normal_CL_face (scalar f, scalar fs, face vector tau_w) {
//    foreach_face() {
//        coord tau = tangential_wall_and_normal_CL_face (point, f, fs);
//        foreach_dimension()
//        tau_w.x[] = tau.x;
//    }
//    boundary((scalar*){tau_w});
//}
//
///**
//A function to suppress glitches after an advection. */
//
//void magnet (scalar f, double error) {
//    foreach() {
//        f[] = clamp(f[], 0., 1.);
//        f[] = (f[] < error ? 0. : (f[] > 1. - error ? 1. : f[]));
//    }
//    boundary ({f});
//}
//
///**
//A function to compute in each point the divergence of a gradient based flux. */
//
//void my_laplacian (scalar f, scalar l, face vector D) {
//    boundary({f, D});
//    foreach() {
//        l[] = 0.;
//        foreach_dimension() l[] += (f[1] - f[0])*D.x[1] - (f[] - f[-1])*D.x[];
//        l[] /= sq(Delta);
//    }
//    boundary({l});
//}

/**
A function to mesure the length of the interface in the cell. Warning: the
length is normalised by the size of the cell. To get the real interface length
you have to multiplie it by the cell size $\Delta$. */

//double interface_length (Point point, scalar c)
//{
//    coord n = mycs (point, c);
//    double alpha = line_alpha (c[], n);
//    coord coord_centroid = {0, 0};
//    return line_length_center(n, alpha, &coord_centroid);
//}

#define betwf(val) ((val > 0) && (val < 1))

#define F_LIQ_EPS 1e-6
//void correct_f(scalar f, scalar f_corr){
//    double neighb, f1, f2;
//    foreach(){
//        f_corr[] = f[];
//    }
//    boundary((scalar *){f_corr});
////    // Correction step
//    foreach(){
//        if (interfacial (point, f)){
//            neighb = f[-1,1];
//            f1 = f[-1,0];
//            f2 = f[0,1];
//            if ( !betwf(f1) && !betwf(f2) && betwf(neighb) ){
//                f_corr[-1,0] = f[-1,0]*fabs(f[-1,0] - F_LIQ_EPS);
//                f_corr[0, 1] = f[0, 1]*fabs(f[0, 1] - F_LIQ_EPS);
//            }
//            neighb = f[1,1];
//            f1 = f[0,1];
//            f2 = f[1,0];
//            if ( !betwf(f1) && !betwf(f2) && betwf(neighb) ){
//                f_corr[0,1] = f[0,1]*fabs(f[0,1] - F_LIQ_EPS);
//                f_corr[1,0] = f[1,0]*fabs(f[1,0] - F_LIQ_EPS);
//            }
//            neighb = f[1,-1];
//            f1 = f[1,0];
//            f2 = f[0,-1];
//            if ( !betwf(f1) && !betwf(f2) && betwf(neighb) ){
//                f_corr[1,0] = f[1,0]*fabs(f[1,0] - F_LIQ_EPS);
//                f_corr[0,-1] = f[0,-1]*fabs(f[0,-1] - F_LIQ_EPS);
//            }
//            neighb = f[-1,-1];
//            f1 = f[0,-1];
//            f2 = f[-1,0];
//            if ( !betwf(f1) && !betwf(f2) && betwf(neighb) ){
//                f_corr[0,-1] = f[0,-1]*fabs(f[0,-1] - F_LIQ_EPS);
//                f_corr[-1,0] = f[-1,0]*fabs(f[-1,0] - F_LIQ_EPS);
//            }
//        }
//    }
//    boundary((scalar *){f_corr});
//}
//
//double average_neighbors(Point point, scalar c, int i, int j, int k){
//    double res = 0;
//    int kk=0;
//    for (int ii = -1; ii <= 1; ii++)
//        for (int jj = -1; jj <= 1; jj++)
//#if dimension>2
//            for (int kk = -1; kk <= 1; kk++)
//#endif
//                res += c[i + ii, j + jj, k + kk];
//    return res/pow(3, dimension);
//}



















//void correct_f(scalar f, scalar f_corr){
//    double avg1, avg2, neighb, f1, f2;
//    foreach(){
//        f_corr[] = f[];
//    }
//    boundary((scalar *){f_corr});
//    // Correction step
//    foreach(){
//        if (interfacial (point, f)){
//            neighb = f[-1,1];
//            f1 = f[-1,0];
//            f2 = f[0,1];
//            if ( !betwf(f1) && !betwf(f2) && betwf(neighb) ){
//                avg1 = average_neighbors(point, f, -1, 0, 0);
//                avg2 = average_neighbors(point, f, 0, 1, 0);
//                f_corr[-1,0] = (avg1/(avg1 + avg2));
//                f_corr[0,1] = (avg2/(avg1 + avg2));
//            }
//            neighb = f[1,1];
//            f1 = f[0,1];
//            f2 = f[1,0];
//            if ( !betwf(f1) && !betwf(f2) && betwf(neighb) ){
//                avg1 = average_neighbors(point, f, 0, 1, 0);
//                avg2 = average_neighbors(point, f, 1, 0, 0);
//                f_corr[0,1] = (avg1/(avg1 + avg2));
//                f_corr[1,0] = (avg2/(avg1 + avg2));
//            }
//            neighb = f[1,-1];
//            f1 = f[1,0];
//            f2 = f[0,-1];
//            if ( !betwf(f1) && !betwf(f2) && betwf(neighb) ){
//                avg1 = average_neighbors(point, f, 1, 0, 0);
//                avg2 = average_neighbors(point, f, 0, -1, 0);
//                f_corr[1,0] = (avg1/(avg1 + avg2));
//                f_corr[0,-1] = (avg2/(avg1 + avg2));
//            }
//            neighb = f[-1,-1];
//            f1 = f[0,-1];
//            f2 = f[-1,0];
//            if ( !betwf(f1) && !betwf(f2) && betwf(neighb) ){
//                avg1 = average_neighbors(point, f, 0, -1, 0);
//                avg2 = average_neighbors(point, f, -1, 0, 0);
//                f_corr[0,-1] = (avg1/(avg1 + avg2));
//                f_corr[-1,0] = (avg2/(avg1 + avg2));
//            }
//        }
//    }
//    boundary((scalar *){f_corr});
//}

//void correct_f(scalar f, scalar f_corr){
//    double avg1, avg2, neighb, f1, f2;
//    foreach(){
//        f_corr[] = f[];
//    }
//    boundary((scalar *){f_corr});
//    // Correction step
//    foreach(){
//        if (interfacial (point, f)){
//            neighb = f[-1,1];
//            f1 = f[-1,0];
//            f2 = f[0,1];
//            if ( !betwf(f1) && !betwf(f2) && betwf(neighb) ){
//                f_corr[-1,0] = fabs(f[-1,0] - F_LIQ_EPS);
//                f_corr[0,1] = fabs(f[0,1] - F_LIQ_EPS);
//            }
//            neighb = f[1,1];
//            f1 = f[0,1];
//            f2 = f[1,0];
//            if ( !betwf(f1) && !betwf(f2) && betwf(neighb) ){
//                f_corr[0,1] = fabs(f[0,1] - F_LIQ_EPS);
//                f_corr[1,0] = fabs(f[1,0] - F_LIQ_EPS);
//            }
//            neighb = f[1,-1];
//            f1 = f[1,0];
//            f2 = f[0,-1];
//            if ( !betwf(f1) && !betwf(f2) && betwf(neighb) ){
//                f_corr[1,0] = fabs(f[1,0] - F_LIQ_EPS);
//                f_corr[0,-1] = fabs(f[0,-1] - F_LIQ_EPS);
//            }
//            neighb = f[-1,-1];
//            f1 = f[0,-1];
//            f2 = f[-1,0];
//            if ( !betwf(f1) && !betwf(f2) && betwf(neighb) ){
//                f_corr[0,-1] = fabs(f[0,-1] - F_LIQ_EPS);
//                f_corr[-1,0] = fabs(f[-1,0] - F_LIQ_EPS);
//            }
//        }
//    }
//    boundary((scalar *){f_corr});
//}
