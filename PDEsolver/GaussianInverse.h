/*!
    \file
    The header file for functor Gaussian elimination matrix inversion class
    \brief The header file for functor Gaussian elimination matrix inversion class
    \author Myles Dillon and Kyle Cusack
*/


#ifndef GAUSSIANINVERSE_H
#define GAUSSIANINVERSE_H

#include <vector>
#include "Math_Vector.h"

/*!
    Functor class to find the inverse of a matrix via Gaussian Elimination
    \brief Gaussian elimination matrix inverter class
*/
class GaussianInverse
{
    public:
        /*!
            Functor used to solve for the inverse of a square matrix using Gaussian elimination.
            \pre Matrix A has at least one row and one column,
                 size of b == height of A,
                 T != 0 (not equal w/ literal 0) defined,
                 T = T (assignment) defined,
                 T = 0 (assignment with literal 0) defined,
                 T = 1 (assignment with literal 1) defined,
                 T * T (multiplication) defined,
                 T + T (addition) defined,
                 T - T (subtraction) defined,
                 1 / T (division with literal 1) defined
            \post see return
            \return The matrix inverse of "A" is returned
            \throws std::invalid_argument if "A" is not a square matrix
        */
        template<typename T>
        denseMatrix<T> operator () (const absMatrix<T>& A) const;

    private:
        /*!
            Function used to find the index of a matrix row that should become the pivot row
            \pre A.getHeight() and A.getWidth() >= 1,
                rowNumber < A.getHeight(),
                columnNumber < A.getWidth(),
                T usable as parameter in std::abs(),
                T = T (assignment) defined,
                T > T (greater than) defined
            \post see return
            \return The row number from "rowNumber" to A.getHeight()-1 in matrix A, that has the greatest absolute value in column "columnNumber" of matrix A.
            \throws std::invalid_argument if any of the first three pre-conditions are not met.
                
        */
        template <typename T>
        unsigned getPivotRowNum(const absMatrix<T>& A, const unsigned rowNumber, const unsigned columnNumber) const;
    
};



#include "GaussianInverse.hpp"

#endif