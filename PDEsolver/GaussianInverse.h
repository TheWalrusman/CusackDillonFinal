#ifndef GAUSSIANINVERSE_H
#define GAUSSIANINVERSE_H

#include <vector>
#include "Math_Vector.h"

class GaussianInverse
{
    public:
        template<typename T>
        denseMatrix<T> operator () (const absMatrix<T>& A) const;

    private:
        template <typename T>
        unsigned getPivotRowNum(const absMatrix<T>& A, const unsigned rowNumber, const unsigned columnNumber) const;
    
};



#include "GaussianInverse.hpp"

#endif