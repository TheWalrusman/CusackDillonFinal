#pragma once
#ifndef CENTERDIFFMESH_HPP
#define CENTERDIFFMESH_HPP


template< class T, class U, class V, T xLower(const U, const V), T xUpper(const U, const V), T yLower(const U, const V), T yUpper(const U, const V)>
void CenterDiffMesh<T, U, V, xLower, xUpper, yLower, yUpper>::setSubdivisions(const unsigned numSubs)
{
	if (numSubs < 3)
		throw std::invalid_argument("In CenterDiffMesh::setSubdivisions, numSubs must be >= 3");
	m_subDivisions = numSubs;
}


template< class T, class U, class V, T xLower(const U, const V), T xUpper(const U, const V), T yLower(const U, const V), T yUpper(const U, const V)>
void CenterDiffMesh<T, U, V, xLower, xUpper, yLower, yUpper>::operator()(symmetricMatrix<T>& A, Math_Vector<T>& B)
{


	double h = 1.0 / (m_subDivisions);
	double onefourth = 1.0 / 4.0;

	double x = 1;

	double y = 1;
	double xp;

	double xm;

	double yp;

	double ym;

	double bSum;

	bool firstOnBound;
	bool secOnBound;
	bool thirdOnBound;
	bool fourthOnBound;

	symmetricMatrix<double> correctA((m_subDivisions - 1)*(m_subDivisions - 1));
	A = correctA;
	B.resize((m_subDivisions - 1)*(m_subDivisions - 1));
	Math_Vector<double> nb;

	nb.resize((m_subDivisions - 1)*(m_subDivisions - 1));




	for (unsigned i = 0; i < (m_subDivisions - 1)*(m_subDivisions - 1); i++)
	{

		//===== GENERATE THE PARAM VALUES USED IN EQUATION ======
		xp = x + 1;

		xm = x - 1;

		yp = y + 1;

		ym = y - 1;

		//====================================================== 


		//============ ASSUME NO TERMS ARE ON BOUNDARY AT FIRST =====
		firstOnBound = false;
		secOnBound = false;
		thirdOnBound = false;
		fourthOnBound = false;
		//======================================================


		//================= GENERATE THE NUMBER THAT NEEDS TO GO INTO THE B VECTOR =========
		// We do not need to look at the upper and right location because of the banded and symmetric
		//nature of the matrix we only need to look at the lower and left locations to fill out the matrix
		bSum = 0;

		//see first term
		if (xm == 0)
		{

			bSum += xLower(xm*h, y*h);
			firstOnBound = true;
		}
		if (xm == m_subDivisions)
		{

			bSum += xUpper(xm*h, y*h);
			firstOnBound = true;
		}
		//see second term
		if (xp == 0)
		{

			bSum += xLower(xp*h, y*h);
			secOnBound = true;
		}
		if (xp == m_subDivisions)
		{

			bSum += xUpper(xp*h, y*h);
			secOnBound = true;
		}

		//see third term
		if (ym == 0)
		{

			bSum += yLower(x*h, ym*h);
			thirdOnBound = true;
		}
		if (ym == m_subDivisions)
		{

			bSum += yUpper(x*h, ym*h);
			thirdOnBound = true;
		}
		//see fourth term
		if (yp == 0)
		{

			bSum += yLower(x*h, yp*h);
			fourthOnBound = true;
		}
		if (yp == m_subDivisions)
		{

			bSum += yUpper(x*h, yp*h);
			fourthOnBound = true;
		}

		B[i] = bSum;

		//===================== END GENERATE B VALUE ===========================


		//=================== FILL IN MATRIX =======================
		//mapping equation -----> col = ( ((y/h)-1) * (m_subDivisions-1) ) + (x/h) - 1 
		A.getRef(i, i) = 1; //diagonal is always 1

		unsigned col;

		if (!firstOnBound)
		{
			//xm   y
			col = static_cast<unsigned>((((y)-1) * (m_subDivisions - 1)) + (xm)-1);
			A.getRef(i, col) = -onefourth;

		}

		if (!thirdOnBound)
		{
			//x   ym
			col = static_cast<unsigned>((((ym)-1) * (m_subDivisions - 1)) + (x)-1);
			A.getRef(i, col) = -onefourth;

		}

		//=============================================================



		//======================== UPDATE X AND Y ============================

		x += 1;

		if (x == m_subDivisions)
		{
			x = 1;
			y += 1;

		}
		//========================================================================

	}

	cout << secOnBound << fourthOnBound << endl;



	B = B * onefourth;



}










#endif // !CENTERDIFFMESH_HPP