#ifndef MATRIX_H_
#define MATRIX_H_

#include "common.h"
#include "util.h"

#include <string>
#include <cassert>

template <typename T>
struct Matrix {
    u64 width;
    u64 height;
    T *data; // Laid out sequentially, row by row.
};

template <typename T>
Matrix<T> alloc_matrix(i32 width, i32 height) {
    assert(width > 0 && height > 0);

    Matrix<T> m;
    m.width  = width;
    m.height = height;
    m.data = alloc_array<T>(width * height);
    return m;
}

template <typename T>
void free_matrix(Matrix<T> &matrix) {
    free_and_reset_array(matrix.data);
}

template <typename T>
void zero_matrix(Matrix<T> &matrix) {
    assert(matrix.data != NULL);
    memzero(matrix.data, matrix.width * matrix.height);
}

template <typename T>
inline T &at(Matrix<T> &matrix, i32 i, i32 j) {
    return matrix.data[i*matrix.width + j];
}

template <typename T>
inline T const &at(const Matrix<T> &matrix, i32 i, i32 j) {
    return matrix.data[i*matrix.width + j];
}

#endif // MATRIX_H_
