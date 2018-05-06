/*!
    \file
    This is the header file for the Gauss-Seidel functor class
    \brief Header file for the  Gauss-Seidel functor class
    \Author Myles Dillon and Kyle Cusack
*/

#ifndef GAUSSSEIDEL_H
#define GAUSSSEIDEL_H

#include "Math_Vector.h"
#include "symmetricMatrix.h"
#include "GaussianInverse.h"
#include <iostream>
#include <math.h>
using namespace std;

/*!
    A functor class for solving a system of equations using the method of steepest descent
    \brief steepest descent functor class
*/
class GaussSeidel
{
    private:
        double m_tollerance;            //The value the the approximation error vector's magnitude has to drop below to stop the iterations
        unsigned m_finalIterations;     //how many iterations the last calculation took to get below the tolerance
        double m_finalResidualMag;      //The final approximation error vector magnitude from the last calculation
        unsigned m_maxIterations;       //the maximum number of iterations a calculation can take

    public:
        /*!
            default constructor
            \pre none
            \post sets m_tollerance to 0.0001 and m_maxItertions to 10000
        */
        GaussSeidel():m_tollerance(0.00001), m_maxIterations(10000) {}

        /*!
            mutator for tollerance
            \pre toll >= 0
            \post m_tollerance == toll
            \throws std::invalid_argument if toll < 0
        */
        void setTollerance (const double toll);

        /*!
            mutator for m_maxIterations
            \pre none
            \post m_maxIterations = max
        */
        void setMaxIters (const unsigned max) {m_maxIterations = max; return;}


        /*!
            A function that report how many iterations the last calculation took.
            \pre none
            \post m_finalIterations, the number of iterations it took to complete the last calculation, is returned
        */
        inline unsigned getFinalIterations() {return m_finalIterations;}

        /*!
            A function that reports the magnitude of the approximation error vector on the last iteration of the last calculation
            \pre none
            \post m_finalResidualMat is returned
        */
        inline double getFinalResidualMag() {return m_finalResidualMag;}


        /*!
            Functor to solve Ax=b using steepest decent
            \pre A must be diagonally dominate,
                b.getSize() == A.getHeight(),
                T = 0 (assignment with literal 0) defined,
                T = T (assignment) defined,
                T - T (subtraction) defined,
                T * T (multiplication) defined,
                T + T (addition) defined,
                T / T (division) defined
            \post The solution to the system Ax=b is approximated using the method of steepest descent.
                The iterations of the method are capped at m_maxIterations. Iterating will stop if the magnitude of the 
                error is <= m_tollerance.
            \return The approximate solution to Ax=b
            \throws std::invalid_argument if A is not diagonally dominate or b.getSize() != A.getHeight()
        */
        template <typename T>
        Math_Vector<T> operator () (const symmetricMatrix<T>& A, const Math_Vector<T>& b);

};

#include "GaussSeidel.hpp"

#endif
