int main (int argc, char * argv[])
{
  while (events (true))
    iter = inext, t = tnext;
}


event event1 (t = 4; t <= 12.23; t += 0.1)
  {
    fprintf (stderr, "event1 t = 4; t <= 12.23; t += 0.1:  t=%g i=%d\n", t, i);
  }

  event event2 (i = 5) fprintf (stderr, "event2 i = 5:  t=%g i=%d\n", t, i);
  event whaha (i++) fprintf (stderr, "whaha i++: t=%g i=%d\n", t, i);
  event event3 (i = 2; i *= 2) fprintf (stderr, "event3 i = 2; i *= 2: t=%g i=%d\n", t, i);
  event event4_0 (i++) fprintf (stderr, "event4_0 i++: t=%g i=%d\n", t, i);
  event event4_1 (i++) fprintf (stderr, "event4_1 i++: t=%g i=%d\n", t, i);
  event event4_2 (i++) fprintf (stderr, "event4_2 i++: t=%g i=%d\n", t, i);
  event init (i++) fprintf (stderr, "init i++: t=%g i=%d\n", t, i);

  event event5 (t += 1; t <= 20) fprintf (stderr, "event5 t += 1; t <= 20: t=%g i=%d\n", t, i);

  event event6 (t = end) fprintf (stderr, "event6 end:  t=%g i=%d\n", t, i);

  event event7 (t = 10; t <= 12; i++)
    fprintf (stderr, "event7 t = 10; t <= 12; i++:  t=%g i=%d\n", t, i);

  event event8 (i = end)
    fprintf (stderr, "event8 i = end:  t=%g i=%d\n", t, i);


//event adapt (i++){
//    adapt_wavelet ({cs, u.x, u.y}, (double[]){0.001, 0.01, 0.01}, 10);
//}


//event stop (t = 10) {
//    return 1;
//}

