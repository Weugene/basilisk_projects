#include "utils.h"

int main()
{
  init_grid (16);
  scalar a[];
  trash ({a});
  foreach()
    a[] = x;
  vector ga[];
  gradients ({a}, {ga});
}
