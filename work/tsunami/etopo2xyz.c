#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <arpa/inet.h>
 
/* check that this matches ETOPO2v2c_i2_LSB.hdr */
#define NCOLS 10800
#define NROWS 5400
#define XLLCORNER -180.000000
#define YLLCORNER -90.000000
#define CELLSIZE 0.03333333333
#define NODATA_VALUE -32768
#define BYTEORDER LSBFIRST
#define NUMBERTYPE 2_BYTE_INTEGER
#define MIN_VALUE -10791
#define MAX_VALUE 8440
 
int main (int argc, char * argv[])
{
  double lat, lon;
  int16_t v;
  int i, j;
 
  for (j = 0; j < NROWS; j++) {
    lat = YLLCORNER + CELLSIZE*j;
    for (i = 0; i < NCOLS; i++) {
      lon = XLLCORNER + CELLSIZE*i;
      assert (fread (&v, sizeof (int16_t), 1, stdin));
      assert (v >= MIN_VALUE && v <= MAX_VALUE);
      printf ("%.8f %.8f %d\n", lon + CELLSIZE/2., - (lat + CELLSIZE/2.), v);
    }
    fprintf (stderr, "\rRow %d/%d              ", j + 1, NROWS);
  }
  fputc ('\n', stderr);
  return 0;
}
