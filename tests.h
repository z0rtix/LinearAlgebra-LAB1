#ifndef TESTS_H
#define TESTS_H

#include "matrix.h"
#include "gauss.h"
#include "LU.h"


double *generateVector(int size, int seed);
double vectorNorm(double *v, int size);
double residual(Matrix &A, double *x, double *b);
double relativeError(double *x, double *xTrue, int size);

void experimentSingleSystem();
void experimentMultipleRHS();
void experimentHilbert();
void debugSmallExample();


#endif