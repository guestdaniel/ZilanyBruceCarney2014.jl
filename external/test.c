#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

void test_change_input_vector(double *input) {
    /* Do stuff to input here */
    input[0] = 99.99;
}

void test_manipulate_cvector_with_julia_function(double *output, double *(*testfunc)(double *)) {
    /* Create vector of zeros here */
    double *input; 
    input = (double*)calloc((int) 10, sizeof(double));

    /* Call testfunc here */
    output = testfunc(input);
    free(input);
}