#include "gauss.h"
#include <cmath>


double *GaussSolver::solveNoPivot(Matrix matrix, double *b) {
    int size = matrix.get_size();
    double *x = new double[size];

    for (int row = 0; row < size - 1; row++) {
        double pivot = matrix.get(row, row);
        for (int i = row + 1; i < size; i++) {
            double factor = matrix.get(i, row) / pivot;
            for (int col = row; col < size; col++) {
                double value = matrix.get(i, col) - factor * matrix.get(row, col);
                matrix.set(i, col, value);
            }
            b[i] = b[i] - factor * b[row];
        }
    }

    for (int row = size - 1; row >= 0; row--) {
        double sum = 0.0;
        for (int col = row + 1; col < size; col++) {
            sum += matrix.get(row, col) * x[col];
        }
        x[row] = (b[row] - sum) / matrix.get(row, row);
    }

    return x;
}

double *GaussSolver::solvePartialPivot(Matrix matrix, double *b) {
    int size = matrix.get_size();
    double *x = new double[size];

    for (int row = 0; row < size - 1; row++) {
        int maxRow = row;
        double maxValue = std::abs(matrix.get(row, row));

        for (int i = row + 1; i < size; i++) {
            double value = std::abs(matrix.get(i, row));
            if (value > maxValue) {
                maxValue = value;
                maxRow = i;
            }
        }

        if (maxRow != row) {
            for (int col = 0; col < size; col++) {
                double temp = matrix.get(row, col);
                matrix.set(row, col, matrix.get(maxRow, col));
                matrix.set(maxRow, col, temp);
            }
            double temp = b[row];
            b[row] = b[maxRow];
            b[maxRow] = temp;
        }

        double pivot = matrix.get(row, row);
        for (int i = row + 1; i < size; i++) {
            double factor = matrix.get(i, row) / pivot;
            for (int col = row; col < size; col++) {
                double value = matrix.get(i, col) - factor * matrix.get(row, col);
                matrix.set(i, col, value);
            }
            b[i] = b[i] - factor * b[row];
        }
    }

    for (int row = size - 1; row >= 0; row--) {
        double sum = 0.0;
        for (int col = row + 1; col < size; col++) {
            sum += matrix.get(row, col) * x[col];
        }
        x[row] = (b[row] - sum) / matrix.get(row, row);
    }

    return x;
}