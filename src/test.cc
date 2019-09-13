#include <iostream>
#include <gtest/gtest.h>
#include <complex>

#include "vector.h"
#include "matrix.h"
#include "matrix_op.h"

/* ---------------------------------------------------------------------------*/

TEST(MATLIB, VectorBasics) {
    // Test constructor by value, copy constructor, and empty constructor
    matlib::Vector <float> V(5, 0.2);
    matlib::Vector <float> V1(V);
    matlib::Vector <float> V2;
    // Test the assignement operator
    V2 = V;
    ASSERT_EQ(V2.size(), V.size());
    ASSERT_EQ(V.size(), V1.size());
    ASSERT_EQ(V[2], V1[2]);
    ASSERT_EQ(V.sum(), V2.sum());
    V2[3] = 15;
    ASSERT_NE(V2[3], V[3]);
    // Test the destructor
    V.clear();
    ASSERT_EQ(V.begin(), V.end());
    ASSERT_EQ(V.size(), 0);
    // Test push_back and pop_back
    V2.push_back(10);
    ASSERT_EQ(V2[V2.size()], 10);
    V2.pop_back();
    ASSERT_NE(V2[V2.size()], 10);
    // Test insert and erase
    V2.insert(2, 10);
    ASSERT_EQ(V2[2], 10);
    V2.erase(2);
    ASSERT_NE(V2[2], 10);
    // Test resize
    V2.resize(8);
    ASSERT_EQ(V2.size(), 8);
    ASSERT_EQ(V2[V2.size()], 0);
    // Test rand and randn
    V2.rand(250);
    V1.randn(250);
    ASSERT_LT(std::abs(V2.mean() - 0.5), 0.1);
    ASSERT_LT(V1.mean(), 0.1); // Note : this sometimes fails (rarely)
    // Test sum and mean
    matlib::Vector <int> Vint(3, 3);
    Vint[1] = 2;
    Vint[3] = 4;
    ASSERT_EQ(Vint.sum(), 9);
    ASSERT_EQ(Vint.mean(), 3);
}

/* ---------------------------------------------------------------------------*/

TEST(MATLIB, MatrixBasics) {
    // Test the constructor by value, the copy constructor, the empty constructor
    matlib::Matrix <double> M(5, 4, 5.0);
    matlib::Matrix <double> M1(M);
    matlib::Matrix <double> M2;
    // Test the assignement operator
    M2 = M;
    ASSERT_EQ(M2.size(), M.size());
    ASSERT_EQ(M.size(), M1.size());
    ASSERT_EQ(M(2, 3), M1(2, 3));
    ASSERT_EQ(M.sum(), M2.sum());
    M2(3, 2) = 15;
    ASSERT_NE(M2(3, 2), M(3, 2));
    // Test the destructor
    M.clear();
    ASSERT_EQ(M.begin(), M.end());
    ASSERT_EQ(M.size(), 0);
    // Test the size member functions
    ASSERT_EQ(M1.size(), 20);
    ASSERT_EQ(M1.size_tot(), 20);
    ASSERT_EQ(M1.size_rows(), 5);
    ASSERT_EQ(M1.size_cols(), 4);
    // Test the overloading operator for accessing one particular element, 
    // in a vector-fashion (column-major ordering)
    M1[1]=2.0;M1[3]=2.0;M1[4]=4.0;M1[7]=9.0;M1[12]=0.0;
    ASSERT_EQ(M1(1,1), 2.0);
    ASSERT_EQ(M1(3,1), 2.0);
    ASSERT_EQ(M1(4,1), 4.0);
    ASSERT_EQ(M1(2,2), 9.0);
    ASSERT_EQ(M1(2,3), 0.0);
    // Test the elementwise operations
    M.randn(6,8);
    matlib::Matrix<double> Unit(6, 8, 1.0);
    matlib::Matrix<double> Result = (M + Unit * 2 - Unit * 3 + 1) / M;
    ASSERT_EQ(Result.sum(), 48);
    Result.clear();
    Unit.clear();
    // Test matrix multiplication
    matlib::Matrix<double> lhs(2, 2, 2.0);
    matlib::Matrix<double> rhs(lhs);
    matlib::Matrix<double> expected_res(lhs);
    lhs(1, 1) = 3; lhs(2, 1) = 1; lhs(1, 2) = 2; lhs(2, 2) = 0;
    expected_res(1, 1) = 10; expected_res(2, 1) = 2; 
    expected_res(1, 2) = 10; expected_res(2, 2) = 2;
    Result = lhs.mul(rhs);
    ASSERT_EQ((Result - expected_res).sum(), 0);
    expected_res.clear();
    // Test identity matrix
    lhs.eye(9);
    rhs.rand(9, 9);
    Result = lhs.mul(rhs);
    ASSERT_EQ((Result - rhs).sum(), 0);
    // Check we can create a matrix of std::complex
    matlib::Matrix<std::complex<float>> Mat_c(3, 3, std::complex<float>(1, 3));
}

/* ---------------------------------------------------------------------------*/

TEST(MATLIB, MatrixOperations) {
    // Test the cholesky decomposition works
    matlib::Matrix<double> M2;
    M2.rand(4,5);
    matlib::Matrix<double> M2T(M2);
    M2T.transpose();
    matlib::Matrix<double> M = M2.mul(M2T);

    matlib::Matrix<double> M_chol = matlib::cholesky(M);
    matlib::Matrix<double> M_chol_T(M_chol);
    M_chol_T.transpose();
    
    matlib::Matrix<double> M_res = M_chol.mul(M_chol_T) - M;
    ASSERT_LT(matlib::abs(M_res).sum(), 1e-12);
}