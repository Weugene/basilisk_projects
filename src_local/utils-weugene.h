void MinMaxValues(scalar * list, double * arr_eps) {// for each scalar min and max
    double arr[10][2];
    int ilist = 0;
    for (scalar s in list) {
        double mina= HUGE, maxa= -HUGE;
        foreach( reduction(min:mina) reduction(max:maxa) ){
            if (fabs(s[]) < mina) mina = fabs(s[]);
            if (fabs(s[]) > maxa) maxa = fabs(s[]);
        }
        arr[ilist][0] = mina;
        arr[ilist][1] = maxa;
        ilist++;
//        fprintf(stderr, "arr for i=%d", ilist);
    }

    for (int i = 0; i < ilist; i++){
#if EPS_MAXA == 1
        arr_eps[i] *=arr[i][1];
#elif EPS_MAXA == 2
        arr_eps[i] *= arr[i][1] - arr[i][0];
#else
        arr_eps[i] *= 0.5*(arr[i][0] + arr[i][1]);
#endif
#ifdef DEBUG_MINMAXVALUES
        fprintf(stderr, "MinMaxValues: i=%d, min=%g, max=%g, eps=%g\n", i, arr[i][0], arr[i][1], arr_eps[i]);
#endif
    }
}

stats statsf_weugene (scalar f, scalar fs)
{
    double dvr, min = 1e100, max = -1e100, sum = 0., sum2 = 0., volume = 0.;
    foreach(reduction(+:sum) reduction(+:sum2) reduction(+:volume)
    reduction(max:max) reduction(min:min))
    if (fs[] < 1. && f[] != nodata) {
        dvr = dv()*(1. - fs[]);
        volume += dvr;
        sum    += dvr*f[];
        sum2   += dvr*sq(f[]);
        if (f[] > max) max = f[];
        if (f[] < min) min = f[];
    }
    stats s;
    s.min = min, s.max = max, s.sum = sum, s.volume = volume;
    if (volume > 0.)
        sum2 -= sum*sum/volume;
    s.stddev = sum2 > 0. ? sqrt(sum2/volume) : 0.;
    return s;
}

stats statsf_weugene2 (scalar f, scalar fs)
{
    double dvr, min = 1e100, max = -1e100, sum = 0., sum2 = 0., volume = 0.;
    foreach(reduction(+:sum) reduction(+:sum2) reduction(+:volume)
    reduction(max:max) reduction(min:min))
    if (fs[] == 0. && f[] != nodata) {
        dvr = dv()*(1. - fs[]);
        volume += dvr;
        sum    += dvr*f[];
        sum2   += dvr*sq(f[]);
        if (f[] > max) max = f[];
        if (f[] < min) min = f[];
    }
    stats s;
    s.min = min, s.max = max, s.sum = sum, s.volume = volume;
    if (volume > 0.)
        sum2 -= sum*sum/volume;
    s.stddev = sum2 > 0. ? sqrt(sum2/volume) : 0.;
    return s;
}

norm normf_weugene (scalar f, scalar fs)
{
    double dvr, avg = 0., rms = 0., max = 0., volume = 0.;
    foreach(reduction(max:max) reduction(+:avg)
    reduction(+:rms) reduction(+:volume))
    if (fs[] < 1. && f[] != nodata) {
        dvr = dv()*(1. - fs[]);
        double v = fabs(f[]);
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
    double max = 0.;
    foreach(reduction(max:max)) {
        if (fs[] < 1) {
            double ds = fabs (s[] - sn[]);
            if (ds > max)
                max = ds;
        }
        sn[] = s[];
    }
    return max;
}