#include <iostream>
#include <iomanip>
#include <random>
#include <cmath>
#include <ctime>

#include "tests.h"


double *generateVector(int size, int seed) {
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    double *b = new double[size];
    for (int i = 0; i < size; i++)
        b[i] = dist(gen);

    return b;
}

double vectorNorm(double *v, int size) {
    double sum = 0;
    for (int i = 0; i < size; i++)
        sum += v[i] * v[i];

    return std::sqrt(sum);
}

double residual(Matrix& A, double *x, double *b) {
    int n = A.get_size();
    double *Ax = new double[n];
    for (int i = 0; i < n; i++) {
        Ax[i] = 0;
        for (int j = 0; j < n; j++)
            Ax[i] += A.get(i, j)  *x[j];
    }

    for (int i = 0; i < n; i++)
        Ax[i] -= b[i];

    double r = vectorNorm(Ax, n);
    delete[] Ax;

    return r;
}

double relativeError(double *x, double *xTrue, int size) {
    double *diff = new double[size];
    for (int i = 0; i < size; i++)
        diff[i] = x[i] - xTrue[i];

    double err = vectorNorm(diff, size) / vectorNorm(xTrue, size);
    delete[] diff;

    return err;
}

void experimentSingleSystem() {
    std::cout << "\n=== СРАВНЕНИЕ ВРЕМЕНИ ДЛЯ ОДНОЙ СИСТЕМЫ ===\n\n";
    int sizes[] = {100, 200, 500, 1000};
    std::cout << "n\tGauss\t\tGaussPivot\tLU total\tLU decomp\tLU solve\n";
    for (int s = 0; s < 4; s++) {
        int n = sizes[s];
        Matrix A = Matrix::random(n, 42);
        double *b1 = generateVector(n, 1);
        double *b2 = generateVector(n, 1);
        double *b3 = generateVector(n, 1);

        clock_t start, end;

        start = clock();
        double *x1 = GaussSolver::solveNoPivot(A, b1);
        end = clock();
        double gaussTime = (double)(end - start) / CLOCKS_PER_SEC;
        
        start = clock();
        double *x2 = GaussSolver::solvePartialPivot(A, b2);
        end = clock();
        double gaussPivotTime = (double)(end - start) / CLOCKS_PER_SEC;
        
        Matrix L(n), U(n);
        start = clock();
        LUSolver::decompose(A, L, U);
        end = clock();
        double luDecomp = (double)(end - start) / CLOCKS_PER_SEC;
        
        start = clock();
        double *y = LUSolver::forwardSubstitution(L, b3);
        double *x3 = LUSolver::backwardSubstitution(U, y);
        end = clock();
        double luSolve = (double)(end - start) / CLOCKS_PER_SEC;
        
        double luTotal = luDecomp + luSolve;
        std::cout << n << "\t"
                  << gaussTime << "\t"
                  << gaussPivotTime << "\t"
                  << luTotal << "\t"
                  << luDecomp << "\t"
                  << luSolve << "\n";

        delete[] b1;
        delete[] b2;
        delete[] b3;
        delete[] x1;
        delete[] x2;
        delete[] x3;
        delete[] y;
    }
}

void experimentMultipleRHS() {
    std::cout << "\n=== МНОЖЕСТВЕННЫЕ ПРАВЫЕ ЧАСТИ ===\n\n";
    int n = 500;
    int ks[] = {1, 10, 100};
    Matrix A = Matrix::random(n, 42);
    std::cout << "k\tGaussPivot\tLU\n";
    for (int t = 0; t < 3; t++) {
        int k = ks[t];

        clock_t start, end;
        start = clock();
        for (int i = 0; i < k; i++) {
            double *b = generateVector(n, i);
            double *x = GaussSolver::solvePartialPivot(A, b);
            delete[] b;
            delete[] x;
        }
        end = clock();
        double gaussTime = (double)(end - start) / CLOCKS_PER_SEC;

        Matrix L(n), U(n);
        start = clock();
        LUSolver::decompose(A, L, U);
        for (int i = 0; i < k; i++) {
            double *b = generateVector(n, i);
            double *y = LUSolver::forwardSubstitution(L, b);
            double *x = LUSolver::backwardSubstitution(U, y);
            delete[] b;
            delete[] x;
            delete[] y;
        }
        end = clock();
        double luTime = (double)(end - start) / CLOCKS_PER_SEC;

        std::cout << k << "\t" << gaussTime << "\t" << luTime << "\n";
    }
}

void experimentHilbert() {
    std::cout << "\n=== ТОЧНОСТЬ НА МАТРИЦЕ ГИЛЬБЕРТА ===\n\n";
    int sizes[] = {5, 10, 15};
    std::cout << "n\tError(no pivot)\tError(pivot)\tResidual\n";
    for (int s = 0; s < 3; s++) {
        int n = sizes[s];
        Matrix H = Matrix::hilbert(n);
        double *xTrue = new double[n];
        for (int i = 0; i < n; i++)
            xTrue[i] = 1.0;
        double *b = new double[n];
        for (int i = 0; i < n; i++) {
            b[i] = 0;
            for (int j = 0; j < n; j++)
                b[i] += H.get(i, j) * xTrue[j];
        }
        double *b1 = new double[n];
        double *b2 = new double[n];
        for (int i = 0; i < n; i++) {
            b1[i] = b[i];
            b2[i] = b[i];
        }

        double *x1 = GaussSolver::solveNoPivot(H, b1);
        double *x2 = GaussSolver::solvePartialPivot(H, b2);
        double err1 = relativeError(x1, xTrue, n);
        double err2 = relativeError(x2, xTrue, n);
        double res = residual(H, x2, b);

        std::cout << n << "\t"
                  << err1 << "\t"
                  << err2 << "\t"
                  << res << "\n";

        delete[] xTrue;
        delete[] b;
        delete[] b1;
        delete[] b2;
        delete[] x1;
        delete[] x2;
    }
}

void debugSmallExample() {
    std::cout << "\n=== DEBUG ПРИМЕР ===\n\n";

    int n = 5;
    std::random_device rd;
    Matrix A = Matrix::random(n, rd());
    double* b = generateVector(n, 1);
    std::cout << "Матрица A:\n";
    A.print();

    std::cout << "\nВектор b:\n";
    for (int i = 0; i < n; i++)
        std::cout << b[i] << " ";
    std::cout << "\n";

    double* x = GaussSolver::solvePartialPivot(A, b);
    std::cout << "\nРешение x:\n";
    for (int i = 0; i < n; i++)
        std::cout << x[i] << " ";
    std::cout << "\n";

    double r = residual(A, x, b);
    std::cout << "\nНевязка: " << r << "\n";

    delete[] b;
    delete[] x;
}