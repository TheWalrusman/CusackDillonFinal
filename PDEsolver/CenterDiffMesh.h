#pragma once
#ifndef CENTERDIFFMESH_H
#define CENTERDIFFMESH_H

template< class T, class U, class V, T xLower(const U,const V), T xUpper(const U, const V), T yLower(const U, const V), T yUpper(const U, const V)>
class CenterDiffMesh
{
private:
	unsigned m_stepsize;
public:
	void setSubdivisions(const unsigned numSubs);

	void operator()(denseMatrix<T>& A, Math_Vector<T>& B);






};




#include "CenterDiffMesh.hpp"
#endif // !CENTERDIFFMESH_H

