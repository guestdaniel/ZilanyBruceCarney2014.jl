#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

void test_change_input_vector(double *input) {
    /* Do stuff to input here */
    input[0] = 99.99;
}

void test_manipulate_cvector_with_julia_function(void(*testfunc)(double *)) {
    /* Create vector of zeros here */
    double *input; 
    input = (double*)calloc((int) 10, sizeof(double));

    /* Call testfunc here */
    testfunc(input);
    free(input);
}

void test_manipulate_julia_vector_with_julia_function(double *input, void(*testfunc)(double *)) {
    /* Call testfunc here */
    testfunc(input);
}