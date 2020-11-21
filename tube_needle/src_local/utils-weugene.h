#ifndef BASILISK_HEADER_24
#define BASILISK_HEADER_24
#line 1 "./../src_local/utils-weugene.h"
void MinMaxValues(scalar * list, double * arr_eps) {// for each scalar min and max
    double arr[10][2];
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
        arr_eps[i] *=arr[i][1];
#elif EPS_MAXA == 2
        arr_eps[i] *= arr[i][1] - arr[i][0];
#else
        arr_eps[i] *= 0.5*(arr[i][0] + arr[i][1]);
#endif
#ifdef DEBUG_MINMAXVALUES
        fprintf(stderr, "MinMaxValues: name=%s, min=%g, max=%g, eps=%g\n", s.name, arr[i][0], arr[i][1], arr_eps[i]);
#endif
        i++;
    }
}

int count_cells(double t, int i){
    int tnc = 0, nc = 0;
    foreach( reduction(+:tnc) )
        tnc++;
#if _MPI
    foreach()
        nc++;
    int rank, h_len;
    char hostname[MPI_MAX_PROCESSOR_NAME];
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Get_processor_name(hostname, &h_len);
    printf("i %d t %g hostname %s rank %d num cells %d total num cells %d\n", i, t, hostname, rank, nc, tnc);
#else
    printf("i %d t %g total num cells %d\n", i, t, tnc);
#endif
    fflush(stdout);
    return tnc;
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
        sum2   += dvr*sq(f[]);
        if (val > max) max = val;
        if (val < min) min = val;
//        fprintf(ferr, "val=%g\n", val);
    }
    fprintf(ferr, "sum=%g\n", sum);
    stats s;
    s.min = min, s.max = max, s.sum = sum, s.volume = volume;
    if (volume > 0.)
        sum2 -= sum*sum/volume;
    s.stddev = sum2 > 0. ? sqrt(sum2/volume) : 0.;
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
    fprintf(ferr, "sum=%g\n", sum);
    stats s;
    s.min = min, s.max = max, s.sum = sum, s.volume = volume;
    if (volume > 0.)
        sum2 -= sum*sum/volume;
    s.stddev = sum2 > 0. ? sqrt(sum2/volume) : 0.;
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

#endif
