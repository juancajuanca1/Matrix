Matrix Functions Overview

Constructors:
Matrix()
Default constructor; creates an empty matrix.

Matrix(int rows, int cols, const T& initial = T())
Creates a matrix with the specified number of rows and columns, with all elements initialized to a given value.

Matrix(std::initializer_list<std::initializer_list<T>> init)
Creates a matrix from a nested initializer list (e.g., Matrix<double> A = {{1, 2}, {3, 4}};).

Element Access:
T& operator()(int row, int col)
Access (and modify) the element at the specified row and column (with bounds checking).

const T& operator()(int row, int col) const
Access an element in a read-only context.

int numRows() const
Returns the number of rows in the matrix.

int numCols() const
Returns the number of columns in the matrix.

Arithmetic Operations:
Matrix<T> operator+(const Matrix<T>& other) const
Matrix addition.

Matrix<T> operator-(const Matrix<T>& other) const
Matrix subtraction.

Matrix<T> operator(const Matrix<T>& other) const*
*Matrix multiplication (dot product style).
For each element (i, j):





Additional Operations:
Matrix<T> hadamard(const Matrix<T>& other) const
*Element-wise (Hadamard) product.

Matrix<T> transpose() const
Returns the transpose of the matrix.

T determinant() const
Returns the determinant (only for square matrices).
Matrix<T> inverse() const
Returns the inverse (only for square, non-singular matrices).

Matrix<T> rref() const
Returns the matrix in Reduced Row–Echelon Form (RREF).

std::vector<T> solveSystem(const std::vector<T>& b) const
Solves the linear system Ax = b and returns the solution vector (assuming a unique solution).

std::pair<Matrix<T>, Matrix<T>> luDecomposition() const
Performs LU decomposition (Doolittle method) and returns a pair (L, U) such that A = L * U.

Equality & Output:
bool operator==(const Matrix<T>& other) const
Checks if two matrices are equal elementwise.

bool operator!=(const Matrix<T>& other) const
Checks if two matrices are not equal.

friend std::ostream& operator<<(std::ostream& os, const Matrix<T>& mat)
Allows printing the matrix to an output stream.
