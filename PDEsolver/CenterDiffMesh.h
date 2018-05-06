#ifndef CENTERDIFFMESH_H
#define CENTERDIFFMESH_H

#include "symmetricMatrix.h"

template< class T, class U, class V, T xLower(const U,const V), T xUpper(const U, const V), T yLower(const U, const V), T yUpper(const U, const V)>
class CenterDiffMesh
{
    private:
        unsigned m_subDivisions;
    public:
        void setSubdivisions(const unsigned numSubs);

        void operator()(symmetricMatrix<T>& A, Math_Vector<T>& B);
};


#include "CenterDiffMesh.hpp"
#endif // !CENTERDIFFMESH_H