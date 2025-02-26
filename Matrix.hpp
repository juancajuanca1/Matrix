#ifndef MATRIX_HPP
#define MATRIX_HPP

#include "Matrix.h"
#include <cmath>
#include <limits>

// -----------------------
// Constructors
// -----------------------

// Default constructor: creates an empty matrix.
template<typename T>
Matrix<T>::Matrix() : rows_(0), cols_(0), data_() {}

// Constructor with dimensions and initial value.
template<typename T>
Matrix<T>::Matrix(int rows, int cols, const T& initial)
    : rows_(rows), cols_(cols), data_(rows * cols, initial) {
    if (rows < 0 || cols < 0)
        throw std::invalid_argument("Matrix dimensions must be non-negative");
}

// Constructor from nested initializer list.
// Example: Matrix<double> A = {{1, 2}, {3, 4}};
template<typename T>
Matrix<T>::Matrix(std::initializer_list<std::initializer_list<T>> init) {
    rows_ = init.size();
    if (rows_ == 0)
        cols_ = 0;
    else
        cols_ = init.begin()->size();
    data_.reserve(rows_ * cols_);
    for (auto row : init) {
        if (row.size() != static_cast<size_t>(cols_))
            throw std::invalid_argument("All rows must have the same number of columns");
        for (auto val : row)
            data_.push_back(val);
    }
}

// -----------------------
// Element Access and Dimensions
// -----------------------
template<typename T>
T& Matrix<T>::operator()(int row, int col) {
    if (row < 0 || row >= rows_ || col < 0 || col >= cols_)
        throw std::out_of_range("Matrix index out of range");
    return data_[row * cols_ + col];
}

template<typename T>
const T& Matrix<T>::operator()(int row, int col) const {
    if (row < 0 || row >= rows_ || col < 0 || col >= cols_)
        throw std::out_of_range("Matrix index out of range");
    return data_[row * cols_ + col];
}

template<typename T>
int Matrix<T>::numRows() const {
    return rows_;
}

template<typename T>
int Matrix<T>::numCols() const {
    return cols_;
}

// -----------------------
// Arithmetic Operators
// -----------------------
template<typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T>& other) const {
    if (rows_ != other.rows_ || cols_ != other.cols_)
        throw std::invalid_argument("Matrix dimensions must agree for addition");
    Matrix<T> result(rows_, cols_);
    for (int i = 0; i < rows_ * cols_; ++i)
        result.data_[i] = data_[i] + other.data_[i];
    return result;
}

template<typename T>
Matrix<T> Matrix<T>::operator-(const Matrix<T>& other) const {
    if (rows_ != other.rows_ || cols_ != other.cols_)
        throw std::invalid_argument("Matrix dimensions must agree for subtraction");
    Matrix<T> result(rows_, cols_);
    for (int i = 0; i < rows_ * cols_; ++i)
        result.data_[i] = data_[i] - other.data_[i];
    return result;
}

// Matrix multiplication (dot-product style)
template<typename T>
Matrix<T> Matrix<T>::operator*(const Matrix<T>& other) const {
    if (cols_ != other.rows_)
        throw std::invalid_argument("Matrix dimensions must agree for multiplication");
    Matrix<T> result(rows_, other.cols_, T());
    for (int i = 0; i < rows_; ++i) {
        for (int j = 0; j < other.cols_; ++j) {
            T sum = T();
            for (int k = 0; k < cols_; ++k)
                sum += (*this)(i, k) * other(k, j);
            result(i, j) = sum;
        }
    }
    return result;
}

// Scalar multiplication
template<typename T>
Matrix<T> Matrix<T>::operator*(const T& scalar) const {
    Matrix<T> result(rows_, cols_);
    for (int i = 0; i < rows_ * cols_; ++i)
        result.data_[i] = data_[i] * scalar;
    return result;
}

// Scalar division
template<typename T>
Matrix<T> Matrix<T>::operator/(const T& scalar) const {
    if (scalar == T())
        throw std::invalid_argument("Division by zero");
    Matrix<T> result(rows_, cols_);
    for (int i = 0; i < rows_ * cols_; ++i)
        result.data_[i] = data_[i] / scalar;
    return result;
}

// Unary minus
template<typename T>
Matrix<T> Matrix<T>::operator-() const {
    Matrix<T> result(rows_, cols_);
    for (int i = 0; i < rows_ * cols_; ++i)
        result.data_[i] = -data_[i];
    return result;
}

// -----------------------
// Element-wise (Hadamard) Product
// -----------------------
template<typename T>
Matrix<T> Matrix<T>::hadamard(const Matrix<T>& other) const {
    if (rows_ != other.rows_ || cols_ != other.cols_)
        throw std::invalid_argument("Matrix dimensions must agree for Hadamard product");
    Matrix<T> result(rows_, cols_);
    for (int i = 0; i < rows_ * cols_; ++i)
        result.data_[i] = data_[i] * other.data_[i];
    return result;
}

// -----------------------
// Transpose
// -----------------------
template<typename T>
Matrix<T> Matrix<T>::transpose() const {
    Matrix<T> result(cols_, rows_);
    for (int i = 0; i < rows_; ++i)
        for (int j = 0; j < cols_; ++j)
            result(j, i) = (*this)(i, j);
    return result;
}

// -----------------------
// Determinant (Gaussian elimination)
// -----------------------
template<typename T>
T Matrix<T>::determinant() const {
    if (rows_ != cols_)
        throw std::domain_error("Determinant is defined only for square matrices");
    if (rows_ == 0) return T();
    if (rows_ == 1) return (*this)(0, 0);
    if (rows_ == 2)
        return (*this)(0,0) * (*this)(1,1) - (*this)(0,1) * (*this)(1,0);
    
    Matrix<T> temp(*this);
    T det = T(1);
    const T eps = std::numeric_limits<T>::epsilon();
    
    for (int i = 0; i < rows_; ++i) {
        int pivot = i;
        while (pivot < rows_ && std::abs(temp(pivot, i)) < eps)
            ++pivot;
        if (pivot == rows_)
            return T(); // singular matrix
        if (pivot != i) {
            for (int j = 0; j < cols_; ++j)
                std::swap(temp(i, j), temp(pivot, j));
            det = -det;
        }
        det *= temp(i, i);
        for (int j = i + 1; j < rows_; ++j) {
            T factor = temp(j, i) / temp(i, i);
            for (int k = i; k < cols_; ++k)
                temp(j, k) -= factor * temp(i, k);
        }
    }
    return det;
}

// -----------------------
// Inverse (Gaussian elimination)
// -----------------------
template<typename T>
Matrix<T> Matrix<T>::inverse() const {
    if (rows_ != cols_)
        throw std::domain_error("Inverse is defined only for square matrices");
    int n = rows_;
    Matrix<T> aug(n, 2 * n);
    // Form augmented matrix [A | I]
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j)
            aug(i, j) = (*this)(i, j);
        for (int j = n; j < 2 * n; ++j)
            aug(i, j) = (i == (j - n)) ? T(1) : T(0);
    }
    const T eps = std::numeric_limits<T>::epsilon();
    // Forward elimination
    for (int i = 0; i < n; ++i) {
        int pivot = i;
        while (pivot < n && std::abs(aug(pivot, i)) < eps)
            ++pivot;
        if (pivot == n)
            throw std::domain_error("Matrix is singular and cannot be inverted");
        if (pivot != i) {
            for (int j = 0; j < 2 * n; ++j)
                std::swap(aug(i, j), aug(pivot, j));
        }
        T pivotVal = aug(i, i);
        for (int j = 0; j < 2 * n; ++j)
            aug(i, j) /= pivotVal;
        for (int r = 0; r < n; ++r) {
            if (r != i) {
                T factor = aug(r, i);
                for (int j = 0; j < 2 * n; ++j)
                    aug(r, j) -= factor * aug(i, j);
            }
        }
    }
    // Extract inverse from augmented matrix
    Matrix<T> inv(n, n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            inv(i, j) = aug(i, j + n);
    return inv;
}

// -----------------------
// Reduced Row–Echelon Form (RREF)
// -----------------------
template<typename T>
Matrix<T> Matrix<T>::rref() const {
    Matrix<T> result(*this);
    int lead = 0;
    const T eps = std::numeric_limits<T>::epsilon();
    for (int r = 0; r < rows_; ++r) {
        if (cols_ <= lead)
            break;
        int i = r;
        while (std::abs(result(i, lead)) < eps) {
            i++;
            if (i == rows_) {
                i = r;
                lead++;
                if (lead == cols_)
                    return result;
            }
        }
        // Swap rows i and r
        for (int k = 0; k < cols_; ++k)
            std::swap(result(i, k), result(r, k));
        T lv = result(r, lead);
        for (int k = 0; k < cols_; ++k)
            result(r, k) /= lv;
        for (int i = 0; i < rows_; ++i) {
            if (i != r) {
                T lv2 = result(i, lead);
                for (int k = 0; k < cols_; ++k)
                    result(i, k) -= lv2 * result(r, k);
            }
        }
        lead++;
    }
    return result;
}

// -----------------------
// Solve Linear System Ax = b using RREF
// -----------------------
template<typename T>
std::vector<T> Matrix<T>::solveSystem(const std::vector<T>& b) const {
    if (rows_ != static_cast<int>(b.size()))
        throw std::invalid_argument("Incompatible dimensions for solving linear system");
    // Augment this matrix with vector b
    Matrix<T> aug(rows_, cols_ + 1);
    for (int i = 0; i < rows_; ++i) {
        for (int j = 0; j < cols_; ++j)
            aug(i, j) = (*this)(i, j);
        aug(i, cols_) = b[i];
    }
    // Compute RREF of the augmented matrix
    Matrix<T> rrefMat = aug.rref();
    // Check for no solution (inconsistent system)
    for (int i = 0; i < rows_; ++i) {
        bool allZeros = true;
        for (int j = 0; j < cols_; ++j) {
            if (std::abs(rrefMat(i, j)) > std::numeric_limits<T>::epsilon())
                allZeros = false;
        }
        if (allZeros && std::abs(rrefMat(i, cols_)) > std::numeric_limits<T>::epsilon())
            throw std::domain_error("No solution exists for the system");
    }
    // Assuming a unique solution (square matrix and full rank)
    std::vector<T> solution(cols_);
    for (int i = 0; i < cols_; ++i)
        solution[i] = rrefMat(i, cols_);
    return solution;
}

// -----------------------
// LU Decomposition (Doolittle method)
// Returns a pair (L, U) such that A = L * U
// -----------------------
template<typename T>
std::pair<Matrix<T>, Matrix<T>> Matrix<T>::luDecomposition() const {
    if (rows_ != cols_)
        throw std::domain_error("LU decomposition requires a square matrix");
    int n = rows_;
    Matrix<T> L(n, n, T());
    Matrix<T> U(n, n, T());
    // Initialize L to identity matrix
    for (int i = 0; i < n; ++i)
        L(i, i) = T(1);
    for (int i = 0; i < n; ++i) {
        // Compute U
        for (int j = i; j < n; ++j) {
            T sum = T();
            for (int k = 0; k < i; ++k)
                sum += L(i, k) * U(k, j);
            U(i, j) = (*this)(i, j) - sum;
        }
        // Compute L
        for (int j = i + 1; j < n; ++j) {
            T sum = T();
            for (int k = 0; k < i; ++k)
                sum += L(j, k) * U(k, i);
            if (std::abs(U(i, i)) < std::numeric_limits<T>::epsilon())
                throw std::domain_error("Zero pivot encountered in LU decomposition");
            L(j, i) = ((*this)(j, i) - sum) / U(i, i);
        }
    }
    return std::make_pair(L, U);
}

// -----------------------
// Equality Operators
// -----------------------
template<typename T>
bool Matrix<T>::operator==(const Matrix<T>& other) const {
    if (rows_ != other.rows_ || cols_ != other.cols_)
        return false;
    for (int i = 0; i < rows_ * cols_; ++i)
        if (data_[i] != other.data_[i])
            return false;
    return true;
}

template<typename T>
bool Matrix<T>::operator!=(const Matrix<T>& other) const {
    return !(*this == other);
}

// -----------------------
// Explicit template instantiation for common types (optional)
// -----------------------
template class Matrix<int>;
template class Matrix<double>;

#endif // MATRIX_HPP
