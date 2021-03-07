#include "run.h"
event event1(i++,last){
    fprintf(ferr, "+event i++ last: i=%d\n",i);
}

int main() {
    dt=1;
    fprintf(ferr, "events will be performed the rarest one!!! here 2 events each and one i+=2 events but only the last one will be performed!!!");
    run();
}

event event1 (i++) {
    fprintf(ferr, "+event i++: i=%d\n",i);
}

event event1 (i+=2) {
    fprintf(ferr, "+event i+=2: i=%d\n",i);
}

event stop(i=10){

}
