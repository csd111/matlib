//
//  matrix.h
//
//  Copyright (c) 2015 Clement DOIRE. All rights reserved.
//
// This program is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later 
// version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT 
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details
//
// You can obtain a copy of the GNU General Public License from
// http://www.gnu.org/copyleft/gpl.html or by writing to Free Software 
// Foundation, Inc.,675 Mass Ave, Cambridge, MA 02139, USA.

#ifndef MCL_MATRIX_H_
#define MCL_MATRIX_H_

#include <iostream>
#include <iomanip>
#include <memory>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <random>
#include <cassert>
#include <ctime>

#include "vector.h"

#define ROUND_UP_4(x) ((x + 3) & ~(3)) // rounds UP x to a multiple of 4

inline const int max_value(ptrdiff_t const val1, ptrdiff_t const val2, 
                           ptrdiff_t const val3, ptrdiff_t const leaf) {
    // to get the index of the maximum value between 4 choices
    return val1 > val2 ? (val1 > val3 ? (val1 > leaf ? 1 : 4) : 
           (val3 > leaf ? 3 : 4)) : (val2 > val3 ? (val2 > leaf ? 2 : 4) : 
           (val3 > leaf ? 3 : 4));
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */

/* Please note that dynamic allocation is deliberately not allowed with matrices
 as this would be either too time or memory consuming. If you really need to 
 'push back' some values, if you're reading a file for example, please use a 
 vector, or find a way to obtain the dimensions first. */

namespace matlib {

template <class T> class Matrix {
public:
    typedef T* iterator;
    typedef const T* const_iterator;
    typedef ptrdiff_t size_type;
    typedef T value_type;
    /* ---------------------------------------------------------------------- */
    // Constructors
    Matrix() {Create();};
    explicit Matrix(size_type rows, size_type cols, const T& val = T()) {
        Create(rows,cols,val);
    };
    // Copy constructor
    Matrix(const Matrix& m){Create(m.begin(),m.end(),m.size_rows(),m.size_cols());};
    // Copy constructor from vector
    Matrix(const Vector<T>& v){Create(v.begin(),v.end(),v.size(),1);};
    // Destructor
    ~Matrix() {Uncreate();};
    /* ---------------------------------------------------------------------- */
    // Size member functions
    size_type size() const {return avail - data;};
    size_type size_tot() const {return limit - data;};
    size_type size_rows() const {return mRows;};
    size_type size_cols() const {return mCols;};
    /* ---------------------------------------------------------------------- */
    // Overloaded operators -- use of i-1 to have the same index referencing as 
    // Matlab, with assertion to check to see if the index is within the bounds
    T& operator[](size_type i) { 
        assert(static_cast<size_t>(i-1) < static_cast<size_t>(size_tot())); 
        return data[i-1]; 
    };
    const T& operator[](size_type i) const {
        assert(static_cast<size_t>(i-1) < static_cast<size_t>(size_tot())); 
        return data[i-1];
    };
    T& operator()(size_type i, size_type j) {
        assert(static_cast<size_t>((j-1)*mRows+i-1) < 
               static_cast<size_t>(size_tot())); 
        return data[(j-1)*mRows+i-1]; 
    };
    const T& operator()(size_type i, size_type j) const {
        assert(static_cast<size_t>((j-1)*mRows+i-1) < 
               static_cast<size_t>(size_tot())); 
        return data[(j-1)*mRows+i-1];
    };
    // Assignement operator
    Matrix& operator=(const Matrix&);
    /* ---------------------------------------------------------------------- */
    // Begin iterator
    iterator begin() {return data;};
    const_iterator begin() const {return data;};
    // End iterator
    iterator end() {return avail;};
    const_iterator end() const {return avail;};
    /* ---------------------------------------------------------------------- */
    // Clear a matrix
    inline void clear() {Uncreate();};
    // Reshape the matrix
    void reshape(size_type newRows, size_type newCols);
    // Print the elements to the console
    void print();
    // Write the matrix to a csv file - in a format readable by MATLAB
    void write_mat(std::ofstream& out);
    // Read a comma-separated csv file (from MATLAB) and save it in memory
    void read_mat(std::ifstream& in);
    // Create a matrix with 1s on the main diagonal and 0s elsewhere
    void eye(size_type rows, size_type cols);
    void eye(size_type rows);
    // Create random matrices -- !!! Complex types not permitted !!!
    //                        -- !!! Integer types give random high numbers !!!
    void rand(size_type rows, size_type cols);  //uniformly distributed
    void randn(size_type rows, size_type cols); //normally distributed
    //Concatenation of matrices
    void horzcat(const Matrix& m);
    void horzcat(const Vector<T>& v);
    void vertcat(const Matrix& m);
    void vertcat(const Vector<T>& v);
    // Create a square diagonal matrix from a vector/Matrix::vector + 
    // reverse operation, i.e. Matrix::column-vector from main diagonal of matrix
    void diag(const Vector<T>& v);
    void diag(const Matrix& m); //if m is a row or col vector, this->square diagonal matrix
                                //if m is a matrix, this->col vector from the main diagonal
    /* --------------------------Operations on elements---------------------- */
    // Flip the matrix, either up/down or left/right
    void flipud();
    void fliplr();
    // Transpose the matrix (cache aware algorithm)
    void transpose();
    // Overloading of standard operators to act elementwise on elements
    Matrix operator+(const Matrix&);
    Matrix operator+(const value_type);
    Matrix operator-(const Matrix&);
    Matrix operator-(const value_type);
    Matrix operator*(const Matrix&);    // hadamard product
    Matrix operator*(const value_type);
    Matrix operator/(const Matrix&);    // hadamard division
    Matrix operator/(const value_type);
    Matrix operator^(const value_type); // elementwise power
    void operator+=(const Matrix&);
    void operator+=(const value_type);
    void operator-=(const Matrix&);
    void operator-=(const value_type);
    void operator*=(const Matrix&);
    void operator*=(const value_type);
    void operator/=(const Matrix&);
    void operator/=(const value_type);
    // Elementwise exponential
    void exp();
    // Elementwise natural logarithm
    void log();
    // Elementwise base 10 logarithm
    void log10();
    // Elementwise square root
    void sqrt();
    // Elementwise Square
    void square();
    // Elementwise Inverse
    void eleminv();
    // Absolute value
    void abs();
    // Matrix multiplication (cache aware algorithm)
    Matrix mul(const Matrix&);
    // Returns the sum of all the elements of the matrix
    T sum();
    // Returns the mean value of all the elements of the matrix
    T mean();
    // Returns the square root of the sum of squares (L2-norm)
    T norm();
    // Trace of the matrix
    T trace();
    
private:
    iterator data;   // First element in the matrix
    iterator avail;  // one past the last element in the matrix
    iterator limit;  // avail rounded to the nearest multiple of 4
    
    size_type mRows; // number of rows
    size_type mCols; // number of cols
    
    // Facilities for memory allocation
    std::allocator<T> alloc; // object to handle memory allocation
    // Allocate and initialize the underlying array
    void Create() ;
    void Create(size_type,size_type, const T&);
    void Create(const_iterator, const_iterator,size_type,size_type);
    // Destroy the elements in the array and free the memory
    void Uncreate();
};


/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* ------------------- Definition of the member functions ------------------- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */


/* -------------------------------- Create ---------------------------------- */
template <class T>
void Matrix<T>::Create() {
    data = avail = limit = 0;
    mRows = mCols = 0;
}

template <class T>
void Matrix<T>::Create(size_type rows, size_type cols, const T& val) {
    data = alloc.allocate(ROUND_UP_4(rows*cols));
    avail = data + rows*cols;
    limit = data + ROUND_UP_4(rows*cols);
    std::uninitialized_fill(data, limit, val);
    mRows = rows;
    mCols = cols;
}

template <class T>
void Matrix<T>::Create(const_iterator i, const_iterator j, 
                       size_type rows, size_type cols) {
    // Check for consistency of the dimensions, abort if problem
    assert((j-i) == (rows*cols)); 
    data = alloc.allocate(ROUND_UP_4(j - i));
    avail = std::uninitialized_copy(i, j, data);
    limit = data + ROUND_UP_4(j-i);
    std::uninitialized_fill(avail, limit , T());
    mRows = rows;
    mCols = cols;
}

/* ------------------------------- Uncreate --------------------------------- */
template <class T>
void Matrix<T>::Uncreate() {
    if (data) {
        // destroy (in reverse order) the elements that were constructed
        iterator it = limit;
        while (it != data)
            alloc.destroy(--it);
        // return all the space that was allocated
        alloc.deallocate(data, limit - data);
    }
    // reset pointers to indicate that the matrix is empty again
    data = avail = limit = 0;
    mRows = mCols = 0;
}

/* ------------------------------ Assignement ------------------------------- */
template <class T>
Matrix<T>& Matrix<T>::operator=(const Matrix& rhs) {
    // check for self-assignment
    if (&rhs != this) {
        // free the array in the left-hand side
        Uncreate();
        // copy elements from the right-hand to the left-hand side
        Create(rhs.begin(), rhs.end(), rhs.size_rows(), rhs.size_cols());
    }
    return *this;
}

/* ------------------------------- reshape ---------------------------------- */
template <class T>
void Matrix<T>::reshape(size_type newRows, size_type newCols) {
    // check for consistency of the dimensions, abort if problem
    assert((newRows*newCols) == (mRows*mCols)); 
    mRows = newRows;
    mCols = newCols;
}

/* -------------------------------- print ----------------------------------- */
template <class T>
void Matrix<T>::print() {
    for (size_type i = 0;i<mRows;++i) {
        for (size_type j = 0;j<mCols;++j) {
            std::cout << data[(j*mRows) + i] << std::setw(12); //not sure about the use of setw here
        }
        std::cout << '\n';
    }
}

/* ------------------------------- write_mat -------------------------------- */
template <class T>
void Matrix<T>::write_mat(std::ofstream& out) {
    for (size_type i = 0;i<mRows;++i) {
        for (size_type j = 0;j<mCols-1;++j) {
            out << data[(j*mRows) + i] <<',';
        }
        out << data[((mCols-1)*mRows) + i] << '\n';
    }
    
}

/* ------------------------------- read_mat --------------------------------- */
template <class T>
void Matrix<T>::read_mat(std::ifstream& in) {
    // First create a dynamically allocated vector
    Vector<T> vec;
    // Create the variables containing the nb of rows and elements
    size_type rows = 0;
    size_type idx = 0;
    // Now read
    std::string line;
    if (in.is_open()) {
        while (std::getline(in,line)) {
            ++rows;
            // Convert the line (string format) to a stream
            std::stringstream iss(line); 
            while (iss.good()) { //while there are characters to read
                std::string substr;
                // Get each of the elements separated by a comma
                getline( iss, substr, ',' ); 
                std::stringstream ss(substr);
                // Convert the value to the correct type
                T mVal;
                ss >> mVal; 
                // Push back in the vector
                vec.push_back(mVal); 
                ++idx;
            }
        }
    }
    // Now convert the vector into a Matrix with correct dimensions
    Uncreate();
    // Copy elements from the right-hand side to the left-hand side
    Create(vec.begin(), vec.end(), idx/rows, rows);
    vec.clear();
    // Allocate new space
    size_type totSize = limit - data;
    size_type avSize = avail - data;
    iterator new_data = alloc.allocate(totSize);
    // Start the recursive cache oblivious transpose operation
    recursiveCOtranspose(data,new_data,mRows,mCols);
    std::uninitialized_fill(new_data + avSize, new_data + totSize, T());
    size_type newRows = mCols;
    size_type newCols = mRows;
    // Return the old space
    Uncreate();
    // Reset pointers to point to the newly allocated space
    data = new_data;
    avail = new_data + avSize;
    limit = new_data + totSize;
    mRows = newRows;
    mCols = newCols;
}

/* ---------------------------------- eye ----------------------------------- */
template<class T>
void Matrix<T>::eye(size_type rows, size_type cols) {
    //allocate new space
    size_type newSize = rows*cols;
    iterator new_data = alloc.allocate(ROUND_UP_4(newSize));
    //fill the space with zeros
    iterator new_avail = new_data + newSize;
    iterator new_limit = new_data + ROUND_UP_4(newSize);
    std::uninitialized_fill(new_data, new_limit, T());
    //get the diagonal equal to 1s
    for (size_type i=0;i<rows;++i) {
        new_data[i*rows+i] = T(1.0);
    }
    // return the old space
    Uncreate();
    // reset pointers to point to the newly allocated space
    data = new_data;
    avail = new_avail;
    limit = new_limit;
    mRows = rows;
    mCols = cols;
}

template<class T>
void Matrix<T>::eye(size_type rows) {
    eye(rows, rows);
}

/* -------------------------------- rand ------------------------------------ */
//!!!! Complex types not permitted !!!!
template <class T>
void Matrix<T>::rand(size_type rows, size_type cols){
    size_type newSize = ROUND_UP_4(rows*cols);
    //allocate new space
    iterator new_data = alloc.allocate(newSize);
    //initialize the random number generator
    std::default_random_engine generator(time(NULL));
    std::uniform_real_distribution<T> distribution(T(0.0),T(1.0));
    for (size_type i=0;i<newSize;++i) {
        new_data[i] = distribution(generator);
    }
    // return the old space
    Uncreate();
    // reset pointers to point to the newly allocated space
    data = new_data;
    avail = new_data + rows*cols;
    limit = new_data + newSize;
    mRows = rows;
    mCols = cols;
}

/* -------------------------------- randn ----------------------------------- */
//!!!! Complex types not permitted !!!!
template <class T>
void Matrix<T>::randn(size_type rows, size_type cols){
    size_type newSize = ROUND_UP_4(rows*cols);
    //allocate new space
    iterator new_data = alloc.allocate(newSize);
    //initialize the random number generator
    std::default_random_engine generator(time(NULL));
    std::normal_distribution<T> distribution(T(0.0),T(1.0));
    for (size_type i=0;i<newSize;++i) {
        new_data[i] = distribution(generator);
    }
    // return the old space
    Uncreate();
    // reset pointers to point to the newly allocated space
    data = new_data;
    avail = new_data + rows*cols;
    limit = new_data + newSize;
    mRows = rows;
    mCols = cols;
}

/* -------------------------------- horzcat --------------------------------- */
template<class T>
void Matrix<T>::horzcat(const Matrix& m){
    // check compatibility
    assert(mRows == m.size_rows());
    // allocate new space
    size_type cols = mCols + m.size_cols();
    size_type UpSize = ROUND_UP_4(mRows*cols);
    iterator new_data = alloc.allocate(UpSize);
    // copy elements
    iterator tmp_avail = std::uninitialized_copy(data, avail, new_data);
    iterator new_avail = std::uninitialized_copy(m.begin(),m.end(),tmp_avail);
    iterator new_limit = std::uninitialized_fill_n(new_avail, UpSize - mRows*cols, T());
    // return the old space
    Uncreate();
    // reset pointers to point to the newly allocated space
    data = new_data;
    avail = new_avail;
    limit = new_limit;
    mRows = m.size_rows();
    mCols = cols;
}

template<class T>
void Matrix<T>::horzcat(const Vector<T>& v){
    // check compatibility
    assert(mRows == v.size());
    // allocate new space
    size_type cols = mCols+1;
    size_type UpSize = ROUND_UP_4(mRows*cols);
    iterator new_data = alloc.allocate(UpSize);
    // copy elements
    iterator tmp_avail = std::uninitialized_copy(data, avail, new_data);
    iterator new_avail = std::uninitialized_copy(v.begin(),v.end(),tmp_avail);
    iterator new_limit = std::uninitialized_fill_n(new_avail, 
                                                   UpSize - mRows*cols, 
                                                   T());
    // return the old space
    Uncreate();
    // reset pointers to point to the newly allocated space
    data = new_data;
    avail = new_avail;
    limit = new_limit;
    mRows = v.size();
    mCols = cols;
}

/* ------------------------------- vertcat ---------------------------------- */
template<class T>
void Matrix<T>::vertcat(const Matrix& m) {
    // check compatibility
    assert(mCols == m.size_cols());
    // get useful stuff from the other matrix
    size_type rows = m.size_rows();
    size_type cols = m.size_cols();
    const_iterator data2 = m.begin();
    // allocate new space
    size_type new_rows = rows + mRows;
    size_type UpSize = ROUND_UP_4(mCols*new_rows);
    iterator new_data = alloc.allocate(UpSize);
    // copy elements
    for (size_type i=0;i<mCols;++i) {
        std::uninitialized_copy(data + mRows*i, 
                                data + mRows*(i+1), 
                                new_data + new_rows*i);
        std::uninitialized_copy(data2 + rows*i, 
                                data2 + rows*(i+1), 
                                new_data + new_rows*i + mRows);
    }
    // return the old space
    Uncreate();
    // reset pointers to point to the newly allocated space
    data = new_data;
    avail = new_data + cols*new_rows;
    limit = std::uninitialized_fill_n(avail, UpSize - cols*new_rows, T());
    mRows = new_rows;
    mCols = cols;
}

template<class T>
void Matrix<T>::vertcat(const Vector<T>& v) {
    // check compatibility
    assert(mCols == v.size());
    // get useful stuff from the other matrix
    size_type cols = mCols;
    const_iterator data2 = v.begin();
    // allocate new space
    size_type new_rows = 1 + mRows;
    size_type UpSize = ROUND_UP_4(mCols*new_rows);
    iterator new_data = alloc.allocate(UpSize);
    // copy elements
    for (size_type i=0;i<mCols;++i) {
        std::uninitialized_copy(data + mRows*i, 
                                data + mRows*(i+1), 
                                new_data + new_rows*i);
        std::uninitialized_copy(data2 + i, 
                                data2 + i + 1, 
                                new_data + new_rows*i + mRows);
    }
    // return the old space
    Uncreate();
    // reset pointers to point to the newly allocated space
    data = new_data;
    avail = new_data + cols*new_rows;
    limit = std::uninitialized_fill_n(avail, UpSize - cols*new_rows, T());
    mRows = new_rows;
    mCols = cols;
}

/* -------------------------------- diag ------------------------------------ */
//square diagonal matrix from vector
template<class T>
void Matrix<T>::diag(const Vector<T>& v) {
    size_type rows = v.size();
    //allocate new space
    size_type newSize = rows*rows;
    iterator new_data = alloc.allocate(ROUND_UP_4(newSize));
    //fill the space with zeros
    iterator new_avail = new_data + newSize;
    iterator new_limit = new_data + ROUND_UP_4(newSize);
    std::uninitialized_fill(new_data, new_limit, T());
    //get the diagonal equal to the values of the vector
    for (size_type i=0;i<rows;++i) {
        new_data[i*rows+i] = v[i+1];
    }
    // return the old space
    Uncreate();
    // reset pointers to point to the newly allocated space
    data = new_data;
    avail = new_avail;
    limit = new_limit;
    mRows = rows;
    mCols = rows;
}

template<class T>
void Matrix<T>::diag(const Matrix<T>& m) {
    if (m.size_cols()==1 || m.size_rows()==1) { 
        //square diagonal matrix from column/row matrix
        size_type rows = m.size();
        //allocate new space
        size_type newSize = rows*rows;
        iterator new_data = alloc.allocate(ROUND_UP_4(newSize));
        //fill the space with zeros
        iterator new_avail = new_data + newSize;
        iterator new_limit = new_data + ROUND_UP_4(newSize);
        std::uninitialized_fill(new_data, new_limit, T());
        //get the diagonal equal to the values of the vector
        for (size_type i=0;i<rows;++i) {
            new_data[i*rows+i] = m[i+1];
        }
        // return the old space
        Uncreate();
        // reset pointers to point to the newly allocated space
        data = new_data;
        avail = new_avail;
        limit = new_limit;
        mRows = rows;
        mCols = rows;
    }
    else { 
        // vector from main diagonal of matrix 
        // (!! vector is actually a column matrix here !!)
        size_type rows = std::min(m.size_rows(),m.size_cols());
        //allocate new space
        iterator new_data = alloc.allocate(ROUND_UP_4(rows));
        //fill the space with zeros
        iterator new_avail = new_data + rows;
        iterator new_limit = new_data + ROUND_UP_4(rows);
        std::uninitialized_fill(new_data, new_limit, T());
        //get the vector equal to the values of the main diagonal
        for (size_type i=0;i<rows;++i) {
            new_data[i] = m(i+1,i+1);
        }
        // return the old space
        Uncreate();
        // reset pointers to point to the newly allocated space
        data = new_data;
        avail = new_avail;
        limit = new_limit;
        mRows = rows;
        mCols = 1;
    }
}


/* -------------------------------------------------------------------------- */
/* ----------------------------Operations on elements------------------------ */
/* -------------------------------------------------------------------------- */


/* ------------------------------- flipud ----------------------------------- */
template<class T>
void Matrix<T>::flipud() {
    //allocate new space
    size_type newSize = limit - data;
    iterator new_data = alloc.allocate(newSize);
    iterator new_avail = new_data + (mRows*mCols);
    iterator new_limit = new_data + newSize;
    std::uninitialized_fill(new_data, new_limit, T());
    //copy elements backward
    for (size_type i=0;i<mCols;++i) {
        std::reverse_copy(data + mRows*i, data + mRows*(i+1), new_data + mRows*i);
    }
    // return the old space
    size_type cols = mCols;
    size_type rows = mRows;
    Uncreate();
    // reset pointers to point to the newly allocated space
    data = new_data;
    avail = new_avail;
    limit = new_limit;
    mRows = rows;
    mCols = cols;
}

/* ------------------------------- fliplr ----------------------------------- */
template<class T>
void Matrix<T>::fliplr() {
    //allocate new space
    size_type newSize = limit - data;
    iterator new_data = alloc.allocate(newSize);
    iterator new_limit = new_data + newSize;
    std::uninitialized_fill(new_data, new_limit, T());
    //copy elements
    size_type j=0;
    for (size_type i=mCols-1;i>=0;--i,++j) {
        std::uninitialized_copy(data + mRows*i, data + mRows*(i+1), new_data + mRows*j);
    }
    // return the old space
    size_type cols = mCols;
    size_type rows = mRows;
    Uncreate();
    // reset pointers to point to the newly allocated space
    data = new_data;
    avail = new_data + cols*rows;
    limit = new_limit;
    mRows = rows;
    mCols = cols;
}


/* ------------------------------ transpose --------------------------------- */
// Implementation of a cache aware transpose operation. 
// Adjust 'leaf' size to tailor performances
template<class T> 
void recursiveCOtranspose(const T* input, T* output, ptrdiff_t const rows, 
                          ptrdiff_t const cols, ptrdiff_t const r1 = 0, 
                          ptrdiff_t const c1 = 0, ptrdiff_t r2 = 0, 
                          ptrdiff_t c2 = 0, ptrdiff_t const leaf = 32)
{
    if (!c2) { c2 = cols - c1; }
    if (!r2) { r2 = rows - r1; }
    ptrdiff_t const diff_r = r2 - r1, diff_c = c2 - c1;
    // Check which size of the current block is the biggest
    if (diff_r >= diff_c && diff_r > leaf) { // Divide into 2 sub-problems
        recursiveCOtranspose(input, output, rows, cols, r1, c1, (r1 + r2) / 2, c2);
        recursiveCOtranspose(input, output, rows, cols, (r1 + r2) / 2, c1, r2, c2);
    }
    else if (diff_c > leaf) { // Divide into 2 sub-problems
        recursiveCOtranspose(input, output, rows, cols, r1, c1, r2, (c1 + c2) / 2);
        recursiveCOtranspose(input, output, rows, cols, r1, (c1 + c2) / 2, r2, c2);
    }
    else // Minimum 'leaf' size reached (base case), perform iterative transpose
    {
        for (ptrdiff_t i1 = r1, i2 = (i1 * cols); i1 < r2; ++i1, i2 += cols) {
            for (ptrdiff_t j1 = c1, j2 = (j1 * rows); j1 < c2; ++j1, j2 += rows) {
                output[i2 + j1] = input[j2 + i1];
            }
        }
    }
}

template <class T>
void Matrix<T>::transpose() {
    // allocate new space
    size_type totSize = limit - data;
    size_type avSize = avail - data;
    iterator new_data = alloc.allocate(totSize);
    //start the recursive cache oblivious transpose operation
    recursiveCOtranspose(data,new_data,mRows,mCols);
    std::uninitialized_fill(new_data + avSize, new_data + totSize, T());
    size_type newRows = mCols;
    size_type newCols = mRows;
    // return the old space
    Uncreate();
    // reset pointers to point to the newly allocated space
    data = new_data;
    avail = new_data + avSize;
    limit = new_data + totSize;
    mRows = newRows;
    mCols = newCols;
}

/* ------------------------------- addition --------------------------------- */

template <class T>
Matrix<T> Matrix<T>::operator+(const Matrix& rhs) {
    assert(mRows == rhs.size_rows() && mCols == rhs.size_cols());
    // create output matrix
    Matrix<T> output(mRows,mCols);
    iterator outdat = output.begin();
    const_iterator rhsdata = rhs.begin();
    for (size_type i = 0; i < limit-data; ++i) {
        outdat[i] = data[i] + rhsdata[i];
    }
    return output;
}

template <class T>
Matrix<T> Matrix<T>::operator+(const value_type rhs) {
    // create output matrix
    Matrix<T> output(mRows,mCols);
    iterator outdat = output.begin();
    for (size_type i = 0; i < limit-data; ++i) {
        outdat[i] = data[i] + rhs;
    }
    return output;
}

template <class T>
void Matrix<T>::operator+=(const Matrix& rhs) {
    assert(mRows == rhs.size_rows() && mCols == rhs.size_cols());
    const_iterator rhsdata = rhs.begin();
    for (size_type i = 0; i < limit-data; ++i) {
        data[i] += rhsdata[i];
    }
    return;
}

template <class T>
void Matrix<T>::operator+=(const value_type rhs) {
    for (size_type i = 0; i < limit-data; ++i) {
        data[i] += rhs;
    }
    return;
}

/* ----------------------------- subtraction -------------------------------- */

template <class T>
Matrix<T> Matrix<T>::operator-(const Matrix& rhs) {
    assert(mRows == rhs.size_rows() && mCols == rhs.size_cols());
    // create output matrix
    Matrix<T> output(mRows,mCols);
    iterator outdat = output.begin();
    const_iterator rhsdata = rhs.begin();
    for (size_type i = 0; i < limit-data; ++i) {
        outdat[i] = data[i] - rhsdata[i];
    }
    return output;
}

template <class T>
Matrix<T> Matrix<T>::operator-(const value_type rhs) {
    // create output matrix
    Matrix<T> output(mRows,mCols);
    iterator outdat = output.begin();
    for (size_type i = 0; i < limit-data; ++i) {
        outdat[i] = data[i] - rhs;
    }
    return output;
}

template <class T>
void Matrix<T>::operator-=(const Matrix& rhs) {
    assert(mRows == rhs.size_rows() && mCols == rhs.size_cols());
    const_iterator rhsdata = rhs.begin();
    for (size_type i = 0; i < limit-data; ++i) {
        data[i] -= rhsdata[i];
    }
    return;
}

template <class T>
void Matrix<T>::operator-=(const value_type rhs) {
    for (size_type i = 0; i < limit-data; ++i) {
        data[i] -= rhs;
    }
    return;
}

/* --------------------------- Hadamard product ----------------------------- */

template <class T>
Matrix<T> Matrix<T>::operator*(const Matrix& rhs) {
    assert(mRows == rhs.size_rows() && mCols == rhs.size_cols());
    // create output matrix
    Matrix<T> output(mRows,mCols);
    iterator outdat = output.begin();
    const_iterator rhsdata = rhs.begin();
    for (size_type i = 0; i < limit-data; ++i) {
        outdat[i] = data[i] * rhsdata[i];
    }
    return output;
}

template <class T>
Matrix<T> Matrix<T>::operator*(const value_type rhs) {
    // create output matrix
    Matrix<T> output(mRows,mCols);
    iterator outdat = output.begin();
    for (size_type i = 0; i < limit-data; ++i) {
        outdat[i] = data[i] * rhs;
    }
    return output;
}

template <class T>
void Matrix<T>::operator*=(const Matrix& rhs) {
    assert(mRows == rhs.size_rows() && mCols == rhs.size_cols());
    const_iterator rhsdata = rhs.begin();
    for (size_type i = 0; i < limit-data; ++i) {
        data[i] *= rhsdata[i];
    }
    return;
}

template <class T>
void Matrix<T>::operator*=(const value_type rhs) {
    for (size_type i = 0; i < limit-data; ++i) {
        data[i] *= rhs;
    }
    return;
}

/* ------------------------- Elementwise divide ----------------------------- */

template <class T>
Matrix<T> Matrix<T>::operator/(const Matrix& rhs) {
    assert(mRows == rhs.size_rows() && mCols == rhs.size_cols());
    // create output matrix
    Matrix<T> output(mRows,mCols);
    iterator outdat = output.begin();
    const_iterator rhsdata = rhs.begin();
    for (size_type i = 0; i < limit-data; ++i) {
        outdat[i] = data[i] / rhsdata[i];
    }
    return output;
}

template <class T>
Matrix<T> Matrix<T>::operator/(const value_type rhs) {
    assert(rhs != 0);
    // create output matrix
    Matrix<T> output(mRows,mCols);
    iterator outdat = output.begin();
    for (size_type i = 0; i < limit-data; ++i) {
        outdat[i] = data[i] / rhs;
    }
    return output;
}

template <class T>
void Matrix<T>::operator/=(const Matrix& rhs) {
    assert(mRows == rhs.size_rows() && mCols == rhs.size_cols());
    const_iterator rhsdata = rhs.begin();
    for (size_type i = 0; i < limit-data; ++i) {
        data[i] /= rhsdata[i];
    }
    return;
}

template <class T>
void Matrix<T>::operator/=(const value_type rhs) {
    assert(rhs != 0);
    for (size_type i = 0; i < limit-data; ++i) {
        data[i] /= rhs;
    }
    return;
}


/* ------------------------- Elementwise power ------------------------------ */

template <class T>
Matrix<T> Matrix<T>::operator^(const value_type rhs) {
    // create output matrix
    Matrix<T> output(mRows,mCols);
    iterator outdat = output.begin();
    for (size_type i = 0; i < limit-data; ++i) {
        outdat[i] = pow(data[i],rhs);
    }
    return output;
}

/* ------------------------ Elementwise square ------------------------------ */

template <class T>
void Matrix<T>::square() {
    for (size_type i = 0; i < limit-data; ++i) {
        data[i] = data[i] * data[i];
    }
}

/* ------------------------ Elementwise inverse ----------------------------- */

template <class T>
void Matrix<T>::eleminv() {
    for (size_type i = 0; i < limit-data; ++i) {
        data[i] = T(1.0) / data[i];
    }
}

/* ---------------------- Elementwise exponential --------------------------- */

template <class T>
void Matrix<T>::exp() {
    for (size_type i = 0; i < limit-data; ++i) {
        data[i] = exp(data[i]);
    }
}

/* ------------------- Elementwise natural logarithm ------------------------ */

template <class T>
void Matrix<T>::log() {
    for (size_type i = 0; i < limit-data; ++i) {
        data[i] = log(data[i]);
    }
}

/* ------------------- Elementwise base 10 logarithm ------------------------ */

template <class T>
void Matrix<T>::log10() {
    for (size_type i = 0; i < limit-data; ++i) {
        data[i] = log10(data[i]);
    }
}

/* --------------------- Elementwise square root ---------------------------- */

template <class T>
void Matrix<T>::sqrt() {
    for (size_type i = 0; i < limit-data; ++i) {
        data[i] = sqrt(data[i]);
    }
}

/* ------------------------- Absolute value --------------------------------- */

template <class T>
void Matrix<T>::abs() {
    for (size_type i = 0; i < limit-data; ++i) {
        data[i] = std::abs(data[i]);
    }
}

/* -------------------- Multiplication (dot product) ------------------------ */
// Implementation of a cache aware matrix multiplication operation. 
// Adjust 'leaf' size to tailor performances
template<class T> 
void recursiveCOmultiply(const T* left, const T* right, T* output, 
                         ptrdiff_t const lRows, ptrdiff_t const lCols, 
                         ptrdiff_t const rCols, ptrdiff_t const rl1 = 0, 
                         ptrdiff_t rl2 = 0, ptrdiff_t const cl1 = 0, 
                         ptrdiff_t cl2 = 0, ptrdiff_t const rr1 = 0, 
                         ptrdiff_t rr2 = 0, ptrdiff_t const cr1 = 0, 
                         ptrdiff_t cr2 = 0, ptrdiff_t const leaf = 12)
{
    // Initialization stuff
    if (!cl2) { cl2 = lCols - cl1; }
    if (!rl2) { rl2 = lRows - rl1; }
    if (!cr2) { cr2 = rCols - cr1; }
    if (!rr2) { rr2 = lCols - rr1; }
    // Compute the current block limits
    ptrdiff_t const lDiff_r = rl2 - rl1, lDiff_c = cl2 - cl1, rDiff_c = cr2 - cr1;
    // Check which size of the current block is the biggest
    const int which = max_value(lDiff_r,lDiff_c,rDiff_c,leaf);
    switch (which) { // We are doing A*B with A(n,m) and B(m,p)
        case 1: // n is largest
            recursiveCOmultiply(left, right, output, lRows, lCols, rCols,
                                rl1, (rl1 + rl2) / 2, cl1, cl2, rr1, rr2, cr1, 
                                cr2);
            recursiveCOmultiply(left, right, output, lRows, lCols, rCols,
                                (rl1 + rl2) / 2, rl2, cl1, cl2, rr1, rr2, cr1, 
                                cr2);
            break;
        case 2: // m is largest
            recursiveCOmultiply(left, right, output, lRows, lCols, rCols,
                                rl1, rl2, cl1, (cl1 + cl2) / 2, rr1, 
                                (rr1 + rr2) / 2, cr1, cr2);
            recursiveCOmultiply(left, right, output, lRows, lCols, rCols, 
                                rl1, rl2, (cl1 + cl2) / 2, cl2, (rr1 + rr2) / 2,
                                rr2, cr1, cr2);
            break;
        case 3: // p is largest
            recursiveCOmultiply(left, right, output, lRows, lCols, rCols, 
                                rl1, rl2, cl1, cl2, rr1, rr2, cr1, 
                                (cr1 + cr2) / 2);
            recursiveCOmultiply(left, right, output, lRows, lCols, rCols, 
                                rl1, rl2, cl1, cl2, rr1, rr2, (cr1 + cr2) / 2,
                                cr2);
            break;
        case 4: // max(n,m,p) < leaf (base case)
            for (ptrdiff_t i=rl1; i<rl2; ++i) {
                for (ptrdiff_t j=cr1; j< cr2;++j) {
                    for (ptrdiff_t k=0; k<lDiff_c ;++k) {
                        output[j*lRows+i] += left[(k+cl1)*lRows+i] * 
                            right[j*lCols+k+rr1];
                    }
                }
            }
            break;
        default: 
            // in case we're in none of the cases above - shouldn't happen in practice
            break;
    }
}

template <class T>
Matrix<T> Matrix<T>::mul(const Matrix& rhs) {
    // mCols of left hand side must be equal to mRows of right hand side !!
    assert(mCols == rhs.size_rows()); 
    // create output matrix
    size_type newCols = rhs.size_cols();
    Matrix<T> output(mRows,newCols);
    iterator outdat = output.begin();
    const_iterator rhsdata = rhs.begin();
    recursiveCOmultiply(data,rhsdata,outdat,mRows,mCols,newCols);
    return output;
}

/* -------------------------------- sum ------------------------------------- */

template <class T>
T Matrix<T>::sum() {
    T Sum = 0;
    for (size_type i = 0; i < avail-data; ++i) {
        Sum += data[i];
    }
    return Sum;
}

/* ------------------------------- mean ------------------------------------- */

template <class T>
T Matrix<T>::mean() {
    T Sum = 0;
    for (size_type i = 0; i < avail-data; ++i) {
        Sum += data[i];
    }
    return Sum / (avail-data);
}

/* -------------------------------- norm ------------------------------------ */

template <class T>
T Matrix<T>::norm() {
    T Sum = 0;
    for (size_type i = 0; i < avail-data; ++i) {
        Sum += data[i]*data[i];
    }
    return sqrt(Sum);
}

/* ------------------------------- trace ------------------------------------ */

template <class T>
T Matrix<T>::trace() {
    assert(mRows == mCols);  // limited to square matrices only !!
    T output = 0;
    for (size_type i = 0; i < mRows; ++i) {
        output += data[i+i*mRows];
    }
    return output;
}

}

#endif /* MCL_MATRIX_H_ */