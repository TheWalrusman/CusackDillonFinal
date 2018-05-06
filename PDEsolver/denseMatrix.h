/*!
    \file
    This is the header file for the denseMatrix Library
    \brief Header file for the denseMatrix<T> class
    \author Myles Dillon and Kyle Cusack
*/

#ifndef DENSEMATRIX_H
#define DENSEMATRIX_H

#include <iostream>
#include <fstream>
#include <vector>
#include "Math_Vector.h"
using namespace std;


/*!
    A templated denseMatrix class. This class stores data in a 2-d Math_Vector with each sub-Math_Vector
    representing a row in the denseMatrix. This means that denseMatrix[x][y] is the value
    contained within row x and column y of the denseMatrix.

    \date 3/26/18
    \author Myles Dillon and Kyle Cusack
    \brief Template dense matrix class
*/
template <class T>
class denseMatrix : public absMatrix<T>
{
    private:
        Math_Vector<Math_Vector<T>> m_data;     //data of the denseMatrix

    public:
        /*!
            Default contructor
            \pre None
            \post m_width == 0, m_height == 0
        */
        denseMatrix();

        /*!
            Copy contructor
            \pre T = T (assignment) defined (due to Math_Vector<T> assignment call)
            \post m_width == other.m_width, m_height == other.m_height, m_data is assigned the value of other.m_data
        */
        denseMatrix(const denseMatrix<T>& other);

        /*!
            A constructor for declaring a sqaure denseMatrix of sizeXsize dimensions
            \pre None
            \post m_width == size, m_height == size, m_data is resized to "size", m_data[i] is resized to "size" for each index "i" in m_data
        */
        denseMatrix(const unsigned size);




        /*!
            A constructor that allows for a custom height and width
            \pre none
            \post m_height == height, m_width == width
        */
        denseMatrix(const unsigned height, const unsigned width);





        /*!
            A constructor for declaring a vect.getSize() by 1 matrix. This create a matrix from a math_vector
            \pre vect.m_size != 0,
                 T = T (assignment) defined
            \post A a vect.getSize() by 1 matrix is created. m_heght == vect.getSize(), m_widht == 1, m_data is resized to vect.getSize(), m_data[i] is resized to 1 for all index "i" in m_data, m_data[i][0] == vect[i] for all index "i" in vect    
            \throws std::invalid_argument if vect.getSize() == 0 (vect is empty)    
        */
        denseMatrix(const Math_Vector<T>& vect);



        /*!
            A function for returning a matrix as a vector of matrix columns
            \pre  T = T (assignment) defined, T = 0 (assigment with literal 0) defined
            \post see return
            \return The calling matrix is returned in the form of a 2-d vector with each sub-vector
                containing the values of the columns from the matrix. Note matrix(i,j) == returnVector[j][i] for each i and j in matrix.
        */
        virtual Math_Vector<Math_Vector<T>> getColumns() const;


        /*!
            A function for returning a matrix as a vector of matrix rows
            \pre  T = T (assignment) defined, T = 0 (assignment with literal 0) defined
            \post see return
            \return The calling matrix is returned in the form of a 2-d vector with each sub-vector
                containing the values of the rows from the matrix. Note matrix(i,j) == returnVector[i][j] for each i and j in matrix.
        */
        virtual Math_Vector<Math_Vector<T>> getRows() const;


        /*!
            The right hand side overload for operator [] to access denseMatrix rows (sub-arrays of m_data). Note that a second opertor [] can be used on what is returned by this operator to access a specific element of the denseMatrix row using Math_Vector<T>::operator [].
            \pre i < m_height
            \post see return
            \return vector m_data[i] (row "i" of the denseMatrix) is returned
            \throws std::invalid_argument if i >= m_height
        */
        const Math_Vector<T>& operator [] (const unsigned i) const;

        /*!
            The left hand side overload for operator [] to modify denseMatrix rows. Note that a second opertor [] can be used on what is returned by this operator to access a specific element of the denseMatrix row using Math_Vector<T>::operator [].
            \pre i < m_height
            \post see return
            \return vector m_data[i] (row "i" of the denseMatrix) is returned
            \throws std::invalid_argument if i >= m_height
        */
        Math_Vector<T>& operator [] (const unsigned i);

      

        /*!
            Matrix multiplication operator
            \pre matrix width == rhs height, 
                T = T (assignment) defined,
                 T + T (addition) defined,
                 T * T (multiplication) defined
            \post The matrix product of the calling object and rhs is returned
            \return The matrix product of the calling object and rhs is retruned
            \throws std::invalid_argument if c.o. width != rhs height
        */
        denseMatrix<T> operator * (const absMatrix<T>& rhs) const;



        /*!
            denseMatrix vector multiplication operator.
            \pre size of vect must == m_width, 
                size of vect != 0,
                T = T (assignment) defined,
                T * T (multiplication) defined,
                T + T (addition) defined
            \post see return
            \return vect is converted to a single column denseMatrix and the denseMatrix product of this denseMatrix and the calling denseMatrix is returned.
            \throws std::invalid_argument if vect size != m_width or vect size == 0
        */
        virtual Math_Vector<T> operator * (const Math_Vector<T>& vect) const;

        /*!
            denseMatrix scalar multiplication operator.
            \pre T = T (assignment) defined,
                T * T (multiplication) defined
            \post see return
            \return The scalar product of "scalar" and the calling object is returned.
            
        */
        denseMatrix<T> operator * (const T& scalar) const;



        /*!
            denseMatrix addition operator.
            \pre rhs has the same dimensions as the calling object (implies m_width == rhs.m_width and m_height == rhs.m_height),
                rhs has at least one column and one row, 
                T = T (assignment) defined,
                T + T (addition) defined
            \post see return
            \return The matrix sum of the calling object and rhs is return
            \throws std::invalid_argument if rhs dimensions are not the same as calling object dimensions.
        */
        denseMatrix<T> operator + (const absMatrix<T>& rhs) const;



        /*!
            matrix subtraction operator.
            \pre rhs has the same dimensions as the calling object (implies m_width == rhs.m_width and m_height == rhs.m_height),
                T = T (assignment) defined,
                T + T (addition) defined,
                T * T (multiplication) defined
            \post see return
            \return The result of the denseMatrix subraction between the calling object minus rhs is returned.
            \throws std::invalid_argument if rhs dimensions are not the same as calling object dimensions or rhs does not have at least one row and one column.
        */
        denseMatrix<T> operator - (const absMatrix<T>& rhs) const;


        /*!
            right hand accessor for matrix elements
           \pre r < m_width, 
                c < m_height,
           \post see return
           \return returns the value of the matrix at (r,c)
           \throws std::invalid_argument if r >= m_width or c >= m_height
        */
        virtual T operator () (const unsigned r, const unsigned c) const;


        /*!
            assignment operator for 2-d vectors
            \pre rhs.getSize() >= 1,
                rhs[0] >= 1,
                T = T (assignment) defined
            \post m_data == rhs
            \throws std::invalid_argument if rhs.getSize() < 1 or rhs[0].getSize() < 1
        */
        void operator = (const Math_Vector<Math_Vector<T>>& rhs);


        /*!
            A function for caring out matrix row interchanges on a matrix.
            \pre rowIndex1 and rowIndex2 < m_height,
                T = T (assignment) defined since this function calls std::swap(Math_Vector<T>,Math_Vector<T>)
            \post std::swap is called with m_data[rowIndex1] and m_data[rowIndex2] so that the contents of two rows in the matrix are swapped entry for entry
            \throws std::invalid_argument if rowIndex1 or rowIndex2 are >= m_height
        */
        void rowInterchange(const unsigned rowIndex1, const unsigned rowIndex2);


        /*!
            A function for matrix rowAddition.
            \pre T * T (multiplication) defined, T + T (addition) defined, T = T (assignment) defined, rowIndex1 and rowIndex2 < m_height, rowIndex != rowIndex2
            \post m_data[rowIndex1] == m_data[rowIndex1] + (m_data[rowIndex2] * scalar)
            \throws std::invalid_argument if rowIndex1 and rowIndex2 do not meet their pre-conditions
        */
        void rowAddition(const unsigned rowIndex1, const unsigned rowIndex2, const T& scalar = 1);        


        /*!
            A function for multiplying a row of a matrix by a scalar value
            \pre T = T (assignment) defined, T * T (multiplication) defined, rowIndex < m_height
            \post m_data[rowIndex] == (m_data[rowIndex] * scalar) (the row specified by rowIndex is scalar multiplied)
            \throws std::invalid_argument if rowIndex >= m_height
        */
        void rowMult(const unsigned rowIndex, const T& scalar);


        /*!
            Assingment operator for dense matrices
            \pre none
            \post m_width == rhs.m_width, m_height == rhs.m_height, m_data == rhs.m_data
        */
        denseMatrix<T>& operator = (const denseMatrix<T>& rhs);


        /*!
            Extraction operator for reading in files
            \pre see Math_Vector operator >>
            \post The values from the input file assigned to denseMatrix row by row in order of increasing column number
            \return A reference to "in" for chaining
        */
        template <typename U>
        friend ifstream& operator >> (ifstream& in, denseMatrix<U>& matrix);

};

#include "denseMatrix.hpp"

#endif