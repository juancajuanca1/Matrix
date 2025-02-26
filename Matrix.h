#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <stdexcept>
#include <initializer_list>
#include <iostream>
#include <utility>

// Templated Matrix class
template<typename T>
class Matrix {
public:
    // Constructors
    Matrix();  // Default (empty) matrix
    Matrix(int rows, int cols, const T& initial = T());
    Matrix(std::initializer_list<std::initializer_list<T>> init);
    
    // Element access with bounds checking
    T& operator()(int row, int col);
    const T& operator()(int row, int col) const;
    
    // Dimensions
    int numRows() const;
    int numCols() const;
    
    // Basic arithmetic operations
    Matrix<T> operator+(const Matrix<T>& other) const;
    Matrix<T> operator-(const Matrix<T>& other) const;
    Matrix<T> operator*(const Matrix<T>& other) const;   // Matrix multiplication
    Matrix<T> operator*(const T& scalar) const;            // Scalar multiplication
    Matrix<T> operator/(const T& scalar) const;            // Scalar division
    Matrix<T> operator-() const;                           // Unary minus

    // Element-wise (Hadamard) product
    Matrix<T> hadamard(const Matrix<T>& other) const;
    
    // Transpose of the matrix
    Matrix<T> transpose() const;
    
    // Determinant (only for square matrices)
    T determinant() const;
    
    // Inverse (only for square matrices; using Gaussian elimination)
    Matrix<T> inverse() const;
    
    // Compute the Reduced Row–Echelon Form (RREF)
    Matrix<T> rref() const;
    
    // Solve linear system Ax = b; b is a vector (size equal to number of rows)
    std::vector<T> solveSystem(const std::vector<T>& b) const;
    
    // LU decomposition (Doolittle method) returns a pair (L, U)
    std::pair<Matrix<T>, Matrix<T>> luDecomposition() const;
    
    // Equality operators
    bool operator==(const Matrix<T>& other) const;
    bool operator!=(const Matrix<T>& other) const;
    
    // Stream output for easy printing
    friend std::ostream& operator<<(std::ostream& os, const Matrix<T>& mat) {
        for (int i = 0; i < mat.numRows(); ++i) {
            for (int j = 0; j < mat.numCols(); ++j)
                os << mat(i, j) << " ";
            os << "\n";
        }
        return os;
    }
    
private:
    int rows_, cols_;
    std::vector<T> data_;
};

#include "Matrix.hpp"
#endif // MATRIX_H
