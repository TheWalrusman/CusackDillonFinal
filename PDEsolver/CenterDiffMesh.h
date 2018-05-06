/*!
\file
This is the header file for the CenterDiffMesh Library
\brief Header file for the CenterDiffMesh<T, U, V, xLower, xUpper, yLower, yUpper> class

*/

#ifndef CENTERDIFFMESH_H
#define CENTERDIFFMESH_H

#include "symmetricMatrix.h"
/*!
	A templated CenterDiffMesh class. This functor class will take a symmetric matrix and a math_vector and will fill them out
	accordingly with the defined Mesh Density(m_subDivisions).

\date 5/06/18
\author Myles Dillon & Kyle Cusack
\brief Template dense matrix class
*/
template< class T, class U, class V, T xLower(const U,const V), T xUpper(const U, const V), T yLower(const U, const V), T yUpper(const U, const V)>
class CenterDiffMesh
{
    private:
        unsigned m_subDivisions;	//Mesh density to be applied
    public:
		/*!
		setter
		\pre None
		\post set m_subDivisons to numSubs
		*/
        void setSubdivisions(const unsigned numSubs);
		/*!
		Functor that creates the mesh for the bounded problem
		\pre symmetric matrix and math_vector of any size
		\post will have mutated the given matrix and math_vector according to the mesh Density
		*/
        void operator()(symmetricMatrix<T>& A, Math_Vector<T>& B);
};


#include "CenterDiffMesh.hpp"
#endif // !CENTERDIFFMESH_H