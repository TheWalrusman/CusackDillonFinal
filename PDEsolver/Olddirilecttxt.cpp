#pragma once
#ifndef CENTERDIFFMESH_HPP
#define CENTERDIFFMESH_HPP


template< class T, class U, class V, T xLower(const U, const V), T xUpper(const U, const V), T yLower(const U, const V), T yUpper(const U, const V)>
void CenterDiffMesh<T, U, V, xLower, xUpper, yLower, yUpper>::setSubdivisions(const unsigned numSubs)
{
	m_stepsize = numSubs;
}


template< class T, class U, class V, T xLower(const U, const V), T xUpper(const U, const V), T yLower(const U, const V), T yUpper(const U, const V)>
void CenterDiffMesh<T, U, V, xLower, xUpper, yLower, yUpper>::operator()(denseMatrix<T>& A, Math_Vector<T>& B)
{

	//unsigned m_stepsize = atoi(argv[1]);
	double h = 1.0 / (m_stepsize);
	//double stepsize = 1.0 / (m_stepsize);
	double newx = 1;
	double x = h;
	double newy = 1;
	double y = h;
	double xp;
	double nxp;
	double xm;
	double nxm;
	double yp;
	double nyp;
	double ym;
	double nym;
	double bSum;
	double nbSum;
	bool firstOnBound;
	bool secOnBound;
	bool thirdOnBound;
	bool fourthOnBound;
	//Math_Vector<double> b;
	denseMatrix<double> correctA((m_stepsize - 1)*(m_stepsize - 1));
	A = correctA;
	B.resize((m_stepsize - 1)*(m_stepsize - 1));
	Math_Vector<double> nb;
	denseMatrix<double> nA((m_stepsize - 1)*(m_stepsize - 1));
	nb.resize((m_stepsize - 1)*(m_stepsize - 1));




	for (unsigned i = 0; i < (m_stepsize - 1)*(m_stepsize - 1); i++)
	{

		//===== GENERATE THE PARAM VALUES USED IN EQUATION ======
		xp = x + h;
		nxp = newx + 1;
		xm = x - h;
		nxm = newx - 1;
		yp = y + h;
		nyp = newy + 1;
		ym = y - h;
		nym = newy - 1;
		//====================================================== 


		//============ ASSUME NO TERMS ARE ON BOUNDARY AT FIRST =====
		firstOnBound = false;
		secOnBound = false;
		thirdOnBound = false;
		fourthOnBound = false;
		//======================================================


		//================= GENERATE THE NUMBER THAT NEEDS TO GO INTO THE B VECTOR =========
		bSum = 0;
		nbSum = 0;
		//see third and forth term
		if (newx == 0)
		{
			nbSum += xLower(newx*h, nym*h);
			nbSum += xLower(newx*h, nyp*h);
			bSum += xLower(x, ym);
			bSum += xLower(x, yp);
			thirdOnBound = true;
			fourthOnBound = true;
		}
		if (newx == m_stepsize)
		{
			nbSum += xUpper(newx*h, nym*h);
			nbSum += xUpper(newx*h, nyp*h);
			bSum += xUpper(x, ym);
			bSum += xUpper(x, yp);
			thirdOnBound = true;
			fourthOnBound = true;
		}
		//see first term
		if (nxm == 0)
		{
			nbSum += xLower(nxm*h, newy*h);
			bSum += xLower(xm, y);
			firstOnBound = true;
		}
		if (nxm == m_stepsize)
		{
			nbSum += xUpper(nxm*h, newy*h);
			bSum += xUpper(xm, y);
			firstOnBound = true;
		}
		//see second term
		if (nxp == 0)
		{
			nbSum += xLower(nxp*h, newy*h);
			bSum += xLower(xp, y);
			secOnBound = true;
		}
		if (nxp == m_stepsize)
		{
			nbSum += xUpper(nxp*h, newy*h);
			bSum += xUpper(xp, y);
			secOnBound = true;
		}
		//see first and second term
		if (newy == 0)
		{
			nbSum += yLower(nxm*h, newy*h);
			nbSum += yLower(nxp*h, newy*h);
			bSum += yLower(xm, y);
			bSum += yLower(xp, y);
			firstOnBound = true;
			secOnBound = true;
		}
		if (newy == m_stepsize)
		{
			nbSum += yUpper(nxm*h, newy*h);
			nbSum += yUpper(nxp*h, newy*h);
			bSum += yUpper(xm, y);
			bSum += yUpper(xp, y);
			firstOnBound = true;
			secOnBound = true;
		}
		//see third term
		if (nym == 0)
		{
			nbSum += yLower(newx*h, nym*h);
			bSum += yLower(x, ym);
			thirdOnBound = true;
		}
		if (nym == m_stepsize)
		{
			nbSum += yUpper(newx*h, nym*h);
			bSum += yUpper(x, ym);
			thirdOnBound = true;
		}
		//see fourth term
		if (nyp == 0)
		{
			nbSum += yLower(newx*h, nyp*h);
			bSum += yLower(x, yp);
			fourthOnBound = true;
		}
		if (nyp == m_stepsize)
		{
			nbSum += yUpper(newx*h, nyp*h);
			bSum += yUpper(x, yp);
			fourthOnBound = true;
		}

		B[i] = bSum;
		nb[i] = nbSum;
		//===================== END GENERATE B VALUE ===========================


		//=================== FILL IN MATRIX =======================
		//mapping equation -----> col = ( ((y/h)-1) * (m_stepsize-1) ) + (x/h) - 1 
		A[i][i] = 1; //diagonal is always 1
		nA[i][i] = 1;
		unsigned col;
		unsigned ncol;
		if (!firstOnBound)
		{
			//xm   y
			col = static_cast<unsigned>((((y / h) - 1) * (m_stepsize - 1)) + (xm / h) - 1);
			A[i][col] = -h;
			ncol = static_cast<unsigned>(((((newy*h) / h) - 1) * (m_stepsize - 1)) + ((nxm*h) / h) - 1);
			nA[i][ncol] = -h;
		}
		if (!secOnBound)
		{
			//xp   y
			col = static_cast<unsigned>((((y / h) - 1) * (m_stepsize - 1)) + (xp / h) - 1);
			A[i][col] = -h;
			ncol = static_cast<unsigned>(((((newy*h) / h) - 1) * (m_stepsize - 1)) + ((nxp*h) / h) - 1);
			nA[i][ncol] = -h;
		}
		if (!thirdOnBound)
		{
			//x   ym
			col = static_cast<unsigned>((((ym / h) - 1) * (m_stepsize - 1)) + (x / h) - 1);
			A[i][col] = -h;
			ncol = static_cast<unsigned>(((((nym*h) / h) - 1) * (m_stepsize - 1)) + ((newx*h) / h) - 1);
			nA[i][ncol] = -h;
		}
		if (!fourthOnBound)
		{
			//x   yp
			col = static_cast<unsigned>((((yp / h) - 1) * (m_stepsize - 1)) + (x / h) - 1);
			A[i][col] = -h;
			ncol = static_cast<unsigned>(((((nyp*h) / h) - 1) * (m_stepsize - 1)) + ((newx*h) / h) - 1);
			nA[i][ncol] = -h;
		}
		//=============================================================



		//======================== UPDATE X AND Y ============================
		//update "x" and "y" as needed
		cout << "(x,y) = " << "(" << x << "," << y << ")" << "     i = " << i << endl;
		cout << "(newx,newy) = " << "(" << newx * h << "," << newy * h << ")" << "     i = " << i << endl;
		x += h;
		newx += 1;
		if (newx == m_stepsize)
		{
			x = h;
			y += h;
			newx = 1;
			newy += 1;
		}
		//========================================================================

	}

	cout << endl;
	cout << "A = " << endl;
	cout << A << endl;
	cout << endl;
	cout << "nA = " << endl;
	cout << nA << endl;


	B = B * h;
	cout << "B = " << endl;
	cout << B << endl;
	nb = nb * h;
	cout << "nb = " << endl;
	cout << nb << endl;


}










#endif // !CENTERDIFFMESH_H