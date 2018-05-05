/*!
    \file
    This is the header file for the Gaussian Elimination Library
    \brief Header file for the Guassian Elimination class

*/


#ifndef GAUSSIANSOLVER_H
#define GAUSSIANSOLVER_H

#include "Matrix.h"
#include "Math_Vector.h"

/*!
    A class with functions that can be used to solve a system of equations of the form Ax=b 
    where a A is coeficient matrix, "x" is the system solution, and "b" contains the values
    that the equations are equal to.

    \date 3/26/18
    \author Myles Dillon
    \brief Guassian Elimination class

*/
class GaussianSolver
{
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
        unsigned getPivotRowNum(const Matrix<T>& A, const unsigned rowNumber, const unsigned columnNumber) const;
    
    public:
        /*!
            Functor used to solve for vector x in the system Ax=b using Gaussian elimination.
            \pre Matrix A has at least one row and one column,
                 size of b == height of A,
                 T != 0 (not equal w/ literal 0) defined,
                 T = T (assignment) defined,
                 T = 0 (assignment with literal 0) defined,
                 T += T defined,
                 T * T (multiplication) defined,
                 T + T (addition) defined,
                 T - T (subtraction) defined,
                 1 / T (division with literal 1) defined
            \post If the system defined by Ax=b is consistent the "x" (solution) vector is returned.
            If the system is inconsistent "NO SOLUTION" is inserted into cout and an empty "x" vector is returned.
            \return The solution vector "x" in Ax=b is returned such that for each index "i" in returned vector,
            returnedVector[i] == xi. Note the solutions with infinitly many solutions will have their
            free variables set to zero. "x" will be an empty vector if the system is inconsistent.
            \throws std::invalid_argument if A does not have at least one row and at least one column, or if b.getSize() != A.getHeight().
        */
        template <typename T>
        Math_Vector<T> operator () (const Matrix<T>& A, const Math_Vector<T>& b);


};

#include "GaussianSolver.hpp"

#endif