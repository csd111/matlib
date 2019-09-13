# matlib

A C++/templated/header-only matrix computation library with column-major 
ordering, a few cache-aware implementations and Matlab style indexing and functions.

## Introduction
This was a small personal project I did back in 2015 as an exercise to improve
my C++ skills and be able to port my Matlab code to C++ easily. It was never really 
finished, and I found it again when going through my archives recently. I thought
 it would be nice to upload it in case it is useful to someone. Compared to my 
 original code, I just added a quick cmake, some tests, and put together this readme file. 

Please note it is very suboptimal as there are a lot of memory copies involved 
and no hardware-specific optimizations. The only places where I optimized the code 
to get a speed up are the cache-aware matrix transpose and matrix multiply. 
The matlib::Vector class is pretty much a re-implementation of std::vector with 
a few extra stuff. 

## Usage
``` cppp
#include "vector.h"
#include "matrix.h"
#include "matrix_op.h"

#include <iostream>

int main(int argc, const char * argv[]) {
    // Construct a vector of size 5 filled with value 0.2
    matlib::Vector <float> V(5, 0.2);
    // Construct a copy of this vector
    matlib::Vector <float> V2(V);
    // Change the value of the third element
    V2[3] = 15;
    // Clear the content of a vector
    V.clear();
    // push_back and pop_back
    V2.push_back(10);
    V2.pop_back();
    // insert and erase at index 2
    V2.insert(2, 10);
    V2.erase(2);

    // Construct a matrix of size (5, 4) filled with value 6
    matlib::Matrix <float> M(5, 4, 6.0);
    // Change the value at row 3 and column 2
    M(3, 2) = 15;
    // Perform some elementwise operations
    M.randn(6,8);
    matlib::Matrix<float> Unit(6, 8, 1.0);
    matlib::Matrix<float> Result = (M + Unit * 2 - Unit * 3 + 1) / M;

    // Test the cholesky decomposition of a matrix
    matlib::Matrix<double> M2; M2.rand(4,5);
    matlib::Matrix<double> M2T(M2); M2T.transpose();
    matlib::Matrix<double> M = M2.mul(M2T);
    matlib::Matrix<double> M_chol = matlib::cholesky(M);
    matlib::Matrix<double> M_chol_T(M_chol);
    M_chol_T.transpose();
    matlib::Matrix<double> M_res = M_chol.mul(M_chol_T) - M;
    std::cout << matlib::abs(M_res).sum() << std::endl; // should be close to 0
}
```

## Build

To build you don't need to pre-install anything, and the following commands 
should be enough: 
```bash
mkdir build && cd build
cmake ..
```

## Doc

--------------------------------------------------------------------------------

### Vector

#### Constructors

```Vector()```

```Vector(size_type n, const T& val = T())```

```Vector(const Vector& v)```

#### Size member function

```size_type size()```

#### Iterators

```iterator begin()```

```const_iterator begin()```

```iterator end()```

```const_iterator end()```

#### Overloaded operators

Use of i-1 to have the same index referencing as Matlab, with assertion to 
check to see if the index is within the bounds

```T& operator[](size_type i)```

```const T& operator[](size_type i)```

#### Adding elements

```void push_back(const T& val)```

```void insert(size_type idx, iterator it, size_type size)```

```void insert(size_type idx, value_type val)```

#### Removing elements

```void pop_back()```

```void erase(size_type start_idx, size_type end_idx)```

```void erase(size_type idx)```

Clear the whole vector

```void clear()```

#### Resizing

```void resize(size_type newSize, const value_type& value)```

```void resize(size_type newSize)```

#### Printing

```void print(std::ostream& out = std::cout)```

#### Random vectors

!!!! Complex types not permitted !!!!

Uniformly distributed

```void rand(size_type size)```

Normally distributed

```void randn(size_type size)```

#### Utilities

Sum of all the elements

```T sum()```

Mean value of all the elements

```T mean()```

--------------------------------------------------------------------------------

### Matrix

#### Constructors

```Matrix()```

```Matrix(size_type rows, size_type cols, const T& val = T())```

```Matrix(const Matrix& m)```

```Matrix(const Vector<T>& v)```

#### Size member functions

Total size

```size_type size()```

```size_type size_tot()```

Number of rows and columns

```size_type size_rows()```

```size_type size_cols()```

#### Overloaded operators

Use of i-1 to have the same index referencing as Matlab, with assertion to 
check to see if the index is within the bounds

```T& operator[](size_type i)```

```const T& operator[](size_type i)```

```T& operator()(size_type i, size_type j)```

```const T& operator()(size_type i, size_type j)```

Assignment operator

```Matrix& operator=(const Matrix&)```

Overloading of standard operators to act elementwise on elements (thus the 
multiplication operator acts as Hadamard product)

```Matrix operator+(const Matrix&)```

```Matrix operator+(const value_type)```

```Matrix operator-(const Matrix&)```

```Matrix operator-(const value_type)```

```Matrix operator*(const Matrix&)```

```Matrix operator*(const value_type)```

```Matrix operator/(const Matrix&)```

```Matrix operator/(const value_type)```

```Matrix operator^(const value_type)```

```void operator+=(const Matrix&)```

```void operator+=(const value_type)```

```void operator-=(const Matrix&)```

```void operator-=(const value_type)```

```void operator*=(const Matrix&)```

```void operator*=(const value_type)```

```void operator/=(const Matrix&)```

```void operator/=(const value_type)```

#### Iterators

```iterator begin()```

```const_iterator begin()```

```iterator end()```

```const_iterator end()```

#### Matrix creation

Identity matrix

```void eye(size_type rows)```

```void eye(size_type rows, size_type cols)```

Random matrices (complex type not permitted, int gives random high numbers)

Uniformly distributed

```void rand(size_type rows, size_type cols)```

Normally distributed

```void rand(size_type rows, size_type cols)```

Create a square diagonal matrix from a vector or Matrix (with one of its dim = 1)

```void diag(const Vector<T>& v)```

```void diag(const Matrix& m)```

Reverse operation, i.e. get a Matrix::column-vector from the main diagonal of the matrix

```void diag(const Matrix& m)```

#### Utilities

Clear a matrix

```void clear()```

Reshape a matrix

```void reshape(size_type newRows, size_type newCols)```

Prints its elements to the console

```void print()```

Write the matrix to a csv file - in a format readable by MATLAB

```void write_mat(std::ofstream& out)```

Read a comma-separated csv file (from MATLAB) and save it as a Matrix

```void read_mat(std::ifstream& in)```

Concatenation of matrices

```void horzcat(const Matrix& m)```

```void horzcat(const Vector<T>& v)```

```void vertcat(const Matrix& m)```

```void vertcat(const Vector<T>& v)```

#### Operations on elements

Flip the matrix, either up/down or left/right

```void flipud()```

```void fliplr()```

Transpose the matrix (cache aware algorithm)

```void transpose()```

Elementwise exponential

```void exp()```

Elementwise natural logarithm

```void log()```

Elementwise base 10 logarithm

```void log10()```

Elementwise square root

```void sqrt()```

Elementwise Square

```void square()```

Elementwise Inverse

```void eleminv()```

Absolute value

```void abs()```

Matrix multiplication (cache aware algorithm)

```Matrix mul(const Matrix&)```

Returns the sum of all the elements of the matrix

```T sum()```

Returns the mean value of all the elements of the matrix

```T mean()```

Returns the square root of the sum of squares (L2-norm)

```T norm()```

Trace of the matrix

```T trace()```

--------------------------------------------------------------------------------

### Matrix Operations

Functions working on matlib::Matrix objects

```Matrix<T> matlib::cholesky(const Matrix<T>& M)```

```Matrix<T> matlib::transpose(const Matrix<T>& M)```

```Matrix<T> matlib::exp(const Matrix<T>& M)```

```Matrix<T> matlib::log(const Matrix<T>& M)```

```Matrix<T> matlib::log10(const Matrix<T>& M)```

```Matrix<T> matlib::sqrt(const Matrix<T>& M)```

```Matrix<T> matlib::square(const Matrix<T>& M)```

```Matrix<T> matlib::eleminv(const Matrix<T>& M)```

```Matrix<T> matlib::abs(const Matrix<T>& M)```
