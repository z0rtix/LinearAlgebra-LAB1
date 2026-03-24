#ifndef LU_H
#define LU_H

#include "matrix.h"

class LUSolver {
    public:
        static void decompose(Matrix matrix, Matrix &L, Matrix &U);
        static double* forwardSubstitution(Matrix &L, double *b);
        static double* backwardSubstitution(Matrix &U, double *y);
        static double* solve(Matrix matrix, double *b);
};

#endif