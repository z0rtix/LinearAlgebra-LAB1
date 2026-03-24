#ifndef MATRIX_H
#define MATRIX_H


class Matrix {
    private:
        double **data;
        int size;
    
    public:
        Matrix(int size);
        Matrix(const Matrix &other_matrix);
        ~Matrix();

        static Matrix random(int size, int seed);
        static Matrix hilbert(int size);

        double get(int row, int col) const;
        void set(int row, int col, double value);

        void print() const;
        int get_size() const;
};

#endif