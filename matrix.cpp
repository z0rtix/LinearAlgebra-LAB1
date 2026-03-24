#include "matrix.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <random>


Matrix::Matrix(int new_size) : size(new_size) {
    data = new double *[size];
    for (int i = 0; i < size; i++) {
        data[i] = new double[size] {0.0};
    }
}

Matrix::Matrix(const Matrix &other_matrix) : size(other_matrix.size) {
    data = new double *[size];
    for (int row = 0; row < size; row++) {
        data[row] = new double[size];
        for (int col = 0; col < size; col++) {
            data[row][col] = other_matrix.data[row][col];
        }
    }
}

Matrix::~Matrix() {
    for (int i = 0; i < size; i++) {
        delete[] data[i];
    }
    delete[] data;
}

Matrix Matrix::random(int size, int seed) {
    Matrix random_matrix(size);
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    for(int row = 0; row < size; row++)
        for(int col = 0; col < size; col++)
            random_matrix.data[row][col] = dist(gen);
    return random_matrix;
}

Matrix Matrix::hilbert(int size) {
    Matrix hilbert_matrix(size);
    for (int row = 0; row < size; row++)
        for (int col = 0; col < size; col++)
            hilbert_matrix.data[row][col] = 1.0 / (row + col + 1);
    return hilbert_matrix;
}

double Matrix::get(int row, int col) const {
    if (row < 0 || row >= size || col < 0 || col >= size) {
        return 0.0;
    }
    return data[row][col];
}

void Matrix::set(int row, int col, double value) {
    if (row < 0 || row >= size || col < 0 || col >= size) {
        return;
    }
    data[row][col] = value;
}

void Matrix::print() const {
    for (int row = 0; row < size; row++) {
        for(int col = 0; col < size; col++)
            std::cout << data[row][col] << " ";
        std::cout << std::endl;
    }
}

int Matrix::get_size() const {
    return size;
}