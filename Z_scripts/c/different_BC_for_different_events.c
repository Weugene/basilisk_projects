#include "run.h"

// The "Solver" code
scalar s[];

event event1 (i++) {
  boundary({s});
  foreach_boundary(left)
    assert (s[-1] == 2);
}

event event2 (i++) {
  boundary({s});
  foreach_boundary(left)
    assert (s[-1] == 1);
  return 1;
}
// The setup code:
int main() {
  run();
}
// Hooks for the solver events..
event event1 (i++) {
  s[left] = 2;
}

event event2 (i++) {
  s[left] = 1;
}
