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

	//unsigned m_subDivisions = atoi(argv[1]);
	double h = 1.0 / (m_subDivisions);
	double onefourth = 1.0 / 4.0;
	//double stepsize = 1.0 / (m_subDivisions);
	//double newx = 1;
	double x = 1;
	//double newy = 1;
	double y = 1;
	double xp;
	//double nxp;
	double xm;
	//double nxm;
	double yp;
	//double nyp;
	double ym;
	//double nym;
	double bSum;
	//double nbSum;
	bool firstOnBound;
	bool secOnBound;
	bool thirdOnBound;
	bool fourthOnBound;
	//Math_Vector<double> b;
	symmetricMatrix<double> correctA((m_subDivisions - 1)*(m_subDivisions - 1));
	A = correctA;
	B.resize((m_subDivisions - 1)*(m_subDivisions - 1));
	Math_Vector<double> nb;
	//denseMatrix<double> nA((m_subDivisions - 1)*(m_subDivisions - 1));
	nb.resize((m_subDivisions - 1)*(m_subDivisions - 1));




	for (unsigned i = 0; i < (m_subDivisions - 1)*(m_subDivisions - 1); i++)
	{

		//===== GENERATE THE PARAM VALUES USED IN EQUATION ======
		xp = x + 1;
		//nxp = newx + 1;
		xm = x - 1;
		//nxm = newx - 1;
		yp = y + 1;
		//nyp = newy + 1;
		ym = y - 1;
		//nym = newy - 1;
		//====================================================== 


		//============ ASSUME NO TERMS ARE ON BOUNDARY AT FIRST =====
		firstOnBound = false;
		secOnBound = false;
		thirdOnBound = false;
		fourthOnBound = false;
		//======================================================


		//================= GENERATE THE NUMBER THAT NEEDS TO GO INTO THE B VECTOR =========
		bSum = 0;
		//nbSum = 0;
		//see third and forth term
		//if (x == 0)
		//{
		//	//nbSum += xLower(newx*h, nym*h);
		//	//nbSum += xLower(newx*h, nyp*h);
		//	bSum += xLower(x*h, ym*h);
		//	bSum += xLower(x*h, yp*h);
		//	thirdOnBound = true;
		//	fourthOnBound = true;
		//}
		//if (x == m_subDivisions)
		//{
		//	//nbSum += xUpper(newx*h, nym*h);
		//	//nbSum += xUpper(newx*h, nyp*h);
		//	bSum += xUpper(x*h, ym*h);
		//	bSum += xUpper(x*h, yp*h);
		//	thirdOnBound = true;
		//	fourthOnBound = true;
		//}
		//see first term
		if (xm == 0)
		{
			//nbSum += xLower(nxm*h, newy*h);
			bSum += xLower(xm*h, y*h);
			firstOnBound = true;
		}
		if (xm == m_subDivisions)
		{
			//nbSum += xUpper(nxm*h, newy*h);
			bSum += xUpper(xm*h, y*h);
			firstOnBound = true;
		}
		//see second term
		if (xp == 0)
		{
			//nbSum += xLower(nxp*h, newy*h);
			bSum += xLower(xp*h, y*h);
			secOnBound = true;
		}
		if (xp == m_subDivisions)
		{
			//nbSum += xUpper(nxp*h, newy*h);
			bSum += xUpper(xp*h, y*h);
			secOnBound = true;
		}
		//see first and second term
		//if (y == 0)
		//{
		//	//nbSum += yLower(nxm*h, newy*h);
		//	//nbSum += yLower(nxp*h, newy*h);
		//	bSum += yLower(xm*h, y*h);
		//	bSum += yLower(xp*h, y*h);
		//	firstOnBound = true;
		//	secOnBound = true;
		//}
		//if (y == m_subDivisions)
		//{
		//	//nbSum += yUpper(nxm*h, newy*h);
		//	//nbSum += yUpper(nxp*h, newy*h);
		//	bSum += yUpper(xm*h, y*h);
		//	bSum += yUpper(xp*h, y*h);
		//	firstOnBound = true;
		//	secOnBound = true;
		//}
		//see third term
		if (ym == 0)
		{
			//nbSum += yLower(newx*h, nym*h);
			bSum += yLower(x*h, ym*h);
			thirdOnBound = true;
		}
		if (ym == m_subDivisions)
		{
			//nbSum += yUpper(newx*h, nym*h);
			bSum += yUpper(x*h, ym*h);
			thirdOnBound = true;
		}
		//see fourth term
		if (yp == 0)
		{
			//nbSum += yLower(newx*h, nyp*h);
			bSum += yLower(x*h, yp*h);
			fourthOnBound = true;
		}
		if (yp == m_subDivisions)
		{
			//nbSum += yUpper(newx*h, nyp*h);
			bSum += yUpper(x*h, yp*h);
			fourthOnBound = true;
		}

		B[i] = bSum;
		//nb[i] = nbSum;
		//===================== END GENERATE B VALUE ===========================


		//=================== FILL IN MATRIX =======================
		//mapping equation -----> col = ( ((y/h)-1) * (m_subDivisions-1) ) + (x/h) - 1 
		A.getRef(i,i) = 1; //diagonal is always 1
					 //nA[i][i] = 1;
		unsigned col;
		//unsigned ncol;
		if (!firstOnBound)
		{
			//xm   y
			col = static_cast<unsigned>((((y)-1) * (m_subDivisions - 1)) + (xm)-1);
			A.getRef(i,col) = -onefourth;
			//ncol = static_cast<unsigned>(((((newy*h) / h) - 1) * (m_subDivisions - 1)) + ((nxm*h) / h) - 1);
			//nA[i][ncol] = -h;
		}
		//if (!secOnBound)
		//{
		//	//xp   y
		//	col = static_cast<unsigned>((((y ) - 1) * (m_subDivisions - 1)) + (xp ) - 1);
		//	A[i][col] = -onefourth;
		//	//ncol = static_cast<unsigned>(((((newy*h) / h) - 1) * (m_subDivisions - 1)) + ((nxp*h) / h) - 1);
		//	//nA[i][ncol] = -h;
		//}
		if (!thirdOnBound)
		{
			//x   ym
			col = static_cast<unsigned>((((ym)-1) * (m_subDivisions - 1)) + (x)-1);
			A.getRef(i, col) = -onefourth;
			//ncol = static_cast<unsigned>(((((nym*h) / h) - 1) * (m_subDivisions - 1)) + ((newx*h) / h) - 1);
			//nA[i][ncol] = -h;
		}
		//if (!fourthOnBound)
		//{
		//	//x   yp
		//	col = static_cast<unsigned>((((yp ) - 1) * (m_subDivisions - 1)) + (x ) - 1);
		//	A[i][col] = -onefourth;
		//	//ncol = static_cast<unsigned>(((((nyp*h) / h) - 1) * (m_subDivisions - 1)) + ((newx*h) / h) - 1);
		//	//nA[i][ncol] = -h;
		//}
		//=============================================================



		//======================== UPDATE X AND Y ============================
		//update "x" and "y" as needed
		cout << "(x,y) = " << "(" << x << "," << y << ")" << "     i = " << i << endl;
		//cout << "(newx,newy) = " << "(" << newx * h << "," << newy * h << ")" << "     i = " << i << endl;
		x += 1;
		//newx += 1;
		if (x == m_subDivisions)
		{
			x = 1;
			y += 1;
			//newx = 1;
			//newy += 1;
		}
		//========================================================================

	}

	cout << secOnBound << fourthOnBound << endl;
	cout << "A = " << endl;
	cout << A << endl;
	cout << endl;
	//cout << "nA = " << endl;
	//cout << nA << endl;


	B = B * onefourth;
	cout << "B = " << endl;
	cout << B << endl;
	//nb = nb * h;
	//cout << "nb = " << endl;
	//cout << nb << endl;


}










#endif // !CENTERDIFFMESH_HPP