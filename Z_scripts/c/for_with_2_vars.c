//it is imipossible to put another var, e.g. "j" into for loop parenthesis
void main(){
    double j=-0.1;
    for(int i=1; i<10; i++){
        j += 0.2;
        fprintf(stderr, "i=%d j=%g \n", i, j);
    }
}
