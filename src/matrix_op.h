//
//  matrix_op.h
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

#ifndef MCL_MATRIX_OP_H_
#define MCL_MATRIX_OP_H_

#include "matrix.h"

namespace matlib{

/* ------------------------------- Cholesky --------------------------------- */

template <class T>
Matrix<T> cholesky(const Matrix<T>& M) {
    // The matrix needs to be square (actually it needs to be Hermitian and 
    // positive definite but at least use the method on a square matrix) !!
    assert(M.size_cols() == M.size_rows()); 
    // create output matrix
    ptrdiff_t size = M.size_rows();
    Matrix<T> output(size,size);
    T* outdat = output.begin();
    const T* indat = M.begin();
    T tmp;
    for (ptrdiff_t j= 0; j < size; ++j) {
        tmp = 0.0;
        for (ptrdiff_t k = 0; k < j ; ++k) {
            tmp += outdat[k*size + j]*outdat[k*size + j];
        }
        outdat[j*size + j] = sqrt(indat[j*size + j] - tmp);
        for (ptrdiff_t i = j + 1 ;i < size; ++i) {
            tmp = 0.0;
            for (ptrdiff_t k2 = 0; k2 < j ; ++k2) {
                tmp += outdat[k2*size + i]*outdat[k2*size + j];
            }
            outdat[j*size + i] = (indat[j*size + i] - tmp)/outdat[j*size + j];
        }
    }
    return output;
}

/* ------------------------------ Transpose --------------------------------- */

template <class T>
Matrix<T> transpose(const Matrix<T>& M) {
    // create output matrix
    ptrdiff_t rows = M.size_cols();
    ptrdiff_t cols = M.size_rows();
    Matrix<T> output(rows,cols);
    T* outdat = output.begin();
    const T* indat = M.begin();
    //start the recursive cache oblivious transpose operation
    recursiveCOtranspose(indat,outdat,cols,rows);
    return output;
}

/* ----------------------------- Exponential -------------------------------- */

template <class T>
Matrix<T> exp(const Matrix<T>& M) {
    // create output matrix
    Matrix<T> output(M);
    output.exp();
    return output;
}

/* ------------------------------- Logarithm -------------------------------- */

template <class T>
Matrix<T> log(const Matrix<T>& M) {
    // create output matrix
    Matrix<T> output(M);
    output.log();
    return output;
}

/* --------------------------- Base-10 Logarithm ---------------------------- */

template <class T>
Matrix<T> log10(const Matrix<T>& M) {
    // create output matrix
    Matrix<T> output(M);
    output.log10();
    return output;
}

/* ----------------------------- Square-root -------------------------------- */

template <class T>
Matrix<T> sqrt(const Matrix<T>& M) {
    // create output matrix
    Matrix<T> output(M);
    output.sqrt();
    return output;
}

/* -------------------------------- Square ---------------------------------- */

template <class T>
Matrix<T> square(const Matrix<T>& M) {
    // create output matrix
    Matrix<T> output(M);
    output.square();
    return output;
}

/* -------------------------- Elementwise Inverse --------------------------- */

template <class T>
Matrix<T> eleminv(const Matrix<T>& M) {
    // create output matrix
    Matrix<T> output(M);
    output.eleminv();
    return output;
}

/* ----------------------------- Absolute Value ----------------------------- */

template <class T>
Matrix<T> abs(const Matrix<T>& M) {
    // create output matrix
    Matrix<T> output(M);
    output.abs();
    return output;   
}

}

#endif /* MCL_MATRIX_OP_H_ */
