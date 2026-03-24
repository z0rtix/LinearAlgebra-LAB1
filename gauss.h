#ifndef GAUSS_H
#define GAUSS_H

#include "matrix.h"

class GaussSolver {
    public:
        static double *solveNoPivot(Matrix matrix, double *b);
        static double *solvePartialPivot(Matrix matrix, double *b);
};

#endif