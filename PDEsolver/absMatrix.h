/*!
    \file
    This is the header file for the abstract base matrix library
    \brief Header file for the matrix abstract base class
    \author Myles Dillon and Kyle Cusack
*/

#ifndef ABSMATRIX_H
#define ABSMATRIX_H

#include <string>
#include <iostream>
#include "Math_Vector.h"
using namespace std;

/*!
    abstract base class for matrices.
    \brief abstract base class for matrices
*/
template <typename T>
class absMatrix
{
    public:

        /*!
            right hand accessor
            \pre varies
            \post see return
            \return the element at row "r" and column "c" in the matrix is returned
        */
        virtual T operator () (const unsigned r, const unsigned c) const = 0;
        
        /*!
            matrix-vector multiplication
            \pre varies
            \post see return
            \return The product of the calling object and vect is returned
        */
        virtual Math_Vector<T> operator * (const Math_Vector<T>& vect) const = 0;

        /*!
            a functions for getting the columns of a matrix as a 2-d vector
            \pre varies
            \post see return
            \return The calling matrix is returned in the form of a 2-d vector with each sub-vector
                containing the values of the columns from the matrix. Note matrix(i,j) == returnVector[j][i] for each i and j in matrix.
        */
        virtual Math_Vector<Math_Vector<T>> getColumns() const = 0;

        /*!
            A function for returning a matrix as a vector of matrix rows
           \pre varies
           \post see return
           \return The calling matrix is returned in the form of a 2-d vector with each sub-vector
                containing the values of the rows from the matrix. Note matrix(i,j) == returnVector[i][j] for each i and j in matrix.
        */
        virtual Math_Vector<Math_Vector<T>> getRows() const = 0;

        /*!
            A function to determine if a matrix is diagonally dominate
            \pre T = 0 (assignment with literal 0) defined,
                T += T (addition/assignment) defined,
                T -= T (subtraction/assignment) defined,
                T <= T (less than equal to) defined,
                std::abs(T) defined
            \post see return
            \return True if the calling object is diagonally dominate, false otherwise
        */
        bool diagonallyDominate () const;

        /*!
            Accessor for matrix width
            \pre None
            \post see return
            \return The number of columns in the matrix is returned
        */
        inline unsigned getWidth() const {return this->m_width;}

        /*!
            Accessor for matrix height
            \pre None
            \post see return
            \return the number of rows in the matrix
        */
        inline unsigned getHeight() const {return this->m_height;}


        /*!
            Insertion operator for printing matrices
            \pre none
            \post matrix is inserted in os row-by-row from column 0 to column width-1
            \return A reference to "os" for chaining
        */
        template <typename U>
        friend ostream& operator << (ostream& os, const absMatrix<U>& matrix);

    protected:
        /*!
            Number of columns in the matrix
        */
        unsigned m_width;
        /*!
            Number of rows in the matrix
        */
        unsigned m_height;

};

#include "absMatrix.hpp"

#endif