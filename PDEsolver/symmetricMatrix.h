/*!
    \file
    This is the header file for the symmetric matrix class
    \brief Header file for the symmetric matrix class
    \author Myles Dillon and Kyle Cusack
*/


#ifndef SYMMETRICMATRIX_H
#define SYMMETRICMATRIX_H

#include <string>
#include <iostream>
#include "absMatrix.h"
#include "denseMatrix.h"
#include "Math_Vector.h"
using namespace std;


/*!
    A templated symmetric matrix class. Smmetric matrix data is stored in a triangular 2-d vector (m_data)
    representing the bottom triangle of the symmetric matrix. symmetricMatrix(i,j) == symmetricMatrix(j,i) for
    all "i" and "j" columns in the matrix.

    \brief symmetric matrix template class
*/
template <typename T>
class symmetricMatrix : public absMatrix<T>
{
    public:
        

        /*!
            Default constructor
            \pre none
            \post matrix width == 0, matrix height == 0, m_data is resized to 0
        */
        symmetricMatrix ();


        /*!
            constructor for creating a square matrix
            \pre size >= 1
            \post m_width == size, m_height == size, m_data resized to size, 
                    m_data[i] resized to i+1 for all rows "i" [0,size) in the matrix
            \throws std::invalid_argument if size < 1
        */
        symmetricMatrix(const unsigned size);


        /*!
            right hand side accessor for matrix elements
            \pre r < m_height,
                c < m_width,
                T = T (assignment) defined
            \post see return
            \return The value of the matrix at row "r" and column "c" is returned
            \throws std::invalid_argument if r >= m_height or c >= m_width
        */
        virtual T operator () (const unsigned r, const unsigned c) const;


        /*!
            matrix addition operator
            \pre rhs must have same dimensions as the calling object,
                T = T (assignment) defined,
                T += T (addition/assignment) defined
            \post see return
            \return The matrix sum of the calling object and rhs is returned
            \throws std::invalid_argument if rhs does not have the same dimensions as the calling object
        */
        denseMatrix<T> operator + (const absMatrix<T>& rhs) const;

        /*!
            matrix subtraction operator
            \pre rhs dimensions must match calling object dimensions,
                T = T (assignment) defined,
                T - T (subtration) defined
            \post see return
            \return The matrix difference between the calling object and rhs is returned
            \throws std::invalid_argument if rhs dimensions do not match the calling object's dimensions
        */
        denseMatrix<T> operator - (const absMatrix<T>& rhs) const;

        /*!
            matrix multiplication operator
            \pre m_width == rhs height,
                T = T (assignment) defined,
                T + T (addition) defined,
                T * T (multiplication) defined
            \post see return
            \return The matrix product of the calling object and rhs is returned
            \throws std::invalid_argument if m_width != rhs height
        */
        denseMatrix<T> operator * (const absMatrix<T>& rhs) const;

        /*!
            matrix-scalar multiplication operator
            \pre T = T (assignment) defined,
                T * T (multiplication) defined
            \post see return
            \return A copy of the calling object is returned with each entry multiplied by scalar
        
        */
        symmetricMatrix<T> operator * (const T& scalar) const;


        /*!
            matrix-vector multiplication operator
            \pre m_width == vect size,
                T = T (assignment) defined,
                T * T (mutliplication) defined,
                T + T (addition) defined
            \post see return
            \return The product of the calling object and vect is returned
            \throws std::invalid argument if m_width != vect size
        */
        virtual Math_Vector<T> operator * (const Math_Vector<T>& vect) const;


        /*!
            A function for returning a matrix as a vector of matrix rows
           \pre  T = T (assignment) defined
           \post see return
           \return The calling matrix is returned in the form of a 2-d vector with each sub-vector
                containing the values of the rows from the matrix. Note matrix(i,j) == returnVector[i][j] for each i and j in matrix.
        */
        virtual Math_Vector<Math_Vector<T>> getRows() const;

        /*!
            A function for returning a matrix as a vector of matrix columns
           \pre  T = T (assignment) defined
           \post see return
           \return The calling matrix is returned in the form of a 2-d vector with each sub-vector
                containing the values of the columns from the matrix. Note matrix(i,j) == returnVector[j][i] for each i and j in the matrix.
        */
        virtual Math_Vector<Math_Vector<T>> getColumns() const;


        /*!
            Left hand accessor for symmetric matrices
            \pre c < r,
                c < matrix width,
                r < matrix height
            \post see return
            \return A non-const reference to index (row = r, column = c) of the matrix
            \throws std::invalid_argument if any of the preconditions are violated
        */
        T& getRef(const unsigned r, const unsigned c);


        /*!
            Extraction operator for reading in a matrix from a file
           \pre in >> U defined
           \post matrix.getHeight()*matrix.getWidth() values are read in from "in" and inserted
                into matrix row by row starting from column 0 to column width-1
           \return reference to the ifstream for chaining
        */
        template <typename U>
        friend istream& operator >> (ifstream& in, symmetricMatrix<U>& matrix);
        


    private:
        Math_Vector<Math_Vector<T>> m_data;
};

#include "symmetricMatrix.hpp"

#endif