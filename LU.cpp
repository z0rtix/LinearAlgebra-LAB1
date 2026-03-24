#include "LU.h"


void LUSolver::decompose(Matrix matrix, Matrix &L, Matrix &U) {
    int size = matrix.get_size();
    for (int row = 0; row < size; row++) {
        for (int col = row; col < size; col++) {
            double sum = 0;
            for (int k = 0; k < row; k++) {
                sum += L.get(row, k) * U.get(k, col);
            }
            double value = matrix.get(row, col) - sum;
            U.set(row, col, value);
        }

        for (int col = row; col < size; col++) {
            if (row == col) {
                L.set(row, row, 1.0);
            } 
            else {
                double sum = 0;
                for (int k = 0; k < row; k++) {
                    sum += L.get(col, k) * U.get(k, row);
                }
                double value = (matrix.get(col, row) - sum) / U.get(row, row);
                L.set(col, row, value);
            }
        }
    }
}

double *LUSolver::forwardSubstitution(Matrix &L, double *b) {
    int size = L.get_size();
    double *y = new double[size];
    for (int row = 0; row < size; row++) {
        double sum = 0;
        for (int col = 0; col < row; col++) {
            sum += L.get(row, col) * y[col];
        }
        y[row] = b[row] - sum;
    }
    return y;
}

double *LUSolver::backwardSubstitution(Matrix &U, double *y) {
    int size = U.get_size();
    double *x = new double[size];
    for (int row = size - 1; row >= 0; row--) {
        double sum = 0;
        for (int col = row + 1; col < size; col++) {
            sum += U.get(row, col) * x[col];
        }
        x[row] = (y[row] - sum) / U.get(row, row);
    }
    return x;
}

double *LUSolver::solve(Matrix matrix, double *b) {
    int size = matrix.get_size();

    Matrix L(size);
    Matrix U(size);

    decompose(matrix, L, U);
    double *y = forwardSubstitution(L, b);
    double *x = backwardSubstitution(U, y);
    delete[] y;

    return x;
}