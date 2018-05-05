#include <iostream>
#include "Math_Vector.h"
#include "denseMatrix.h"
#include "BoundFuncts.hpp"
#include "CenterDiffMesh.h"
using namespace std;

double xLower(const double x, const double y)
{
	return 0 * x*y;
}

double xUpper(const double x, const double y)
{
	return 0 * x*y;//
}

double yLower(const double x, const double y)
{
	return (1 - (4 * (x - 0.5)*(x - 0.5))) + 0 * y;
}

double yUpper(const double x, const double y)
{
	return 0 * x*y;
}


int main(int argc, char* argv[])
{
	cout << argc <<argv[0]<< endl;
	denseMatrix<double> A;
	Math_Vector<double> b;
	CenterDiffMesh<double, double, double, xLower, xUpper, yLower, yUpper> dirichlet;
	dirichlet.setSubdivisions(atoi(argv[1]));
	dirichlet(A, b);



//    if(argc){}
//
//    unsigned n = atoi(argv[1]);
//    double h = 1.0/(n);
//	//double stepsize = 1.0 / (n);
//	double newx = 1;
//    double x = h;
//	double newy = 1;
//    double y = h;
//    double xp;
//	double nxp;
//    double xm;
//	double nxm;
//    double yp;
//	double nyp;
//    double ym;
//	double nym;
//    double bSum;
//	double nbSum;
//    bool firstOnBound;
//    bool secOnBound;
//    bool thirdOnBound;
//    bool fourthOnBound;
//    Math_Vector<double> b;
//    denseMatrix<double> A((n-1)*(n-1));
//    b.resize((n-1)*(n-1));
//	Math_Vector<double> nb;
//	denseMatrix<double> nA((n - 1)*(n - 1));
//	nb.resize((n - 1)*(n - 1));
//
//    //"i" is the row of the symmetric matrix
//    for(unsigned i = 0; i < (n-1)*(n-1); i++)
//    {
//
//        //===== GENERATE THE PARAM VALUES USED IN EQUATION ======
//        xp = x+h;
//		nxp = newx + 1;
//        xm = x-h;
//		nxm = newx - 1;
//        yp = y+h;
//		nyp = newy + 1;
//        ym = y-h;
//		nym = newy - 1;
//        //====================================================== 
//
//
//        //============ ASSUME NO TERMS ARE ON BOUNDARY AT FIRST =====
//        firstOnBound = false;
//        secOnBound = false;
//        thirdOnBound = false;
//        fourthOnBound = false;
//        //======================================================
//
//
//        //================= GENERATE THE NUMBER THAT NEEDS TO GO INTO THE B VECTOR =========
//        bSum = 0;
//		nbSum = 0;
//        //see third and forth term
//        if(newx == 0)
//        {
//			nbSum += xLower(newx*h, nym*h);
//			nbSum += xLower(newx*h, nyp*h);
//            bSum += xLower(x,ym);
//            bSum += xLower(x,yp);
//            thirdOnBound = true;
//            fourthOnBound = true;
//        }
//        if(newx == n)
//        {
//			nbSum += xUpper(newx*h, nym*h);
//			nbSum += xUpper(newx*h, nyp*h);
//            bSum += xUpper(x,ym);
//            bSum += xUpper(x,yp);
//            thirdOnBound = true;
//            fourthOnBound = true;
//        }
//        //see first term
//        if(nxm == 0)
//        {
//			nbSum += xLower(nxm*h, newy*h);
//            bSum += xLower(xm,y);
//            firstOnBound = true;
//        }
//        if(nxm == n)
//        {
//			nbSum += xUpper(nxm*h, newy*h);
//            bSum += xUpper(xm,y);
//            firstOnBound = true;
//        }
//        //see second term
//        if(nxp == 0)
//        {
//			nbSum += xLower(nxp*h, newy*h);
//            bSum += xLower(xp,y);
//            secOnBound = true;
//        }
//        if(nxp == n)
//        {
//			nbSum += xUpper(nxp*h, newy*h);
//            bSum += xUpper(xp,y);
//            secOnBound = true;
//        }
//        //see first and second term
//        if(newy == 0)
//        {
//			nbSum += yLower(nxm*h, newy*h);
//			nbSum += yLower(nxp*h, newy*h);
//            bSum += yLower(xm,y);
//            bSum += yLower(xp,y);
//            firstOnBound = true;
//            secOnBound = true;
//        }
//        if(newy == n)
//        {
//			nbSum += yUpper(nxm*h, newy*h);
//			nbSum += yUpper(nxp*h, newy*h);
//            bSum += yUpper(xm,y);
//            bSum += yUpper(xp,y);
//            firstOnBound = true;
//            secOnBound = true;
//        }
//        //see third term
//        if(nym == 0)
//        {
//			nbSum += yLower(newx*h, nym*h);
//            bSum += yLower(x,ym);
//            thirdOnBound = true;
//        }
//        if(nym == n)
//        {
//			nbSum += yUpper(newx*h, nym*h);
//            bSum += yUpper(x,ym);
//            thirdOnBound = true;
//        }
//        //see fourth term
//        if(nyp == 0)
//        {
//			nbSum += yLower(newx*h, nyp*h);
//            bSum += yLower(x,yp);
//            fourthOnBound = true;
//        }
//        if(nyp == n)
//        {
//			nbSum += yUpper(newx*h, nyp*h);
//            bSum += yUpper(x,yp);
//            fourthOnBound = true;
//        }
//
//        b[i] = bSum;
//		nb[i] = nbSum;
//        //===================== END GENERATE B VALUE ===========================
//
//      
//        //=================== FILL IN MATRIX =======================
//        //mapping equation -----> col = ( ((y/h)-1) * (n-1) ) + (x/h) - 1 
//        A[i][i] = 1; //diagonal is always 1
//		nA[i][i] = 1;
//        unsigned col;
//		unsigned ncol;
//        if(!firstOnBound)
//        {
//            //xm   y
//            col = static_cast<unsigned>( ( ((y/h)-1) * (n-1) ) + (xm/h) - 1 );
//            A[i][col] = -h;
//			ncol = static_cast<unsigned>(((((newy*h) / h) - 1) * (n - 1)) + ((nxm*h) / h) - 1);
//			nA[i][ncol] = -h;
//        }
//        if(!secOnBound)
//        {
//            //xp   y
//            col = static_cast<unsigned>( ( ((y/h)-1) * (n-1) ) + (xp/h) - 1 );
//            A[i][col] = -h;
//			ncol = static_cast<unsigned>(((((newy*h) / h) - 1) * (n - 1)) + ((nxp*h) / h) - 1);
//			nA[i][ncol] = -h;
//        }
//        if(!thirdOnBound)
//        {
//            //x   ym
//            col = static_cast<unsigned>( ( ((ym/h)-1) * (n-1) ) + (x/h) - 1 );
//            A[i][col] = -h;
//			ncol = static_cast<unsigned>(((((nym*h) / h) - 1) * (n - 1)) + ((newx*h) / h) - 1);
//			nA[i][ncol] = -h;
//        }
//        if(!fourthOnBound)
//        {
//            //x   yp
//            col = static_cast<unsigned>( ( ((yp/h)-1) * (n-1) ) + (x/h) - 1 );
//            A[i][col] = -h;
//			ncol = static_cast<unsigned>(((((nyp*h) / h) - 1) * (n - 1)) + ((newx*h) / h) - 1);
//			nA[i][ncol] = -h;
//        }
//        //=============================================================
//
//
//
//        //======================== UPDATE X AND Y ============================
//        //update "x" and "y" as needed
//        cout << "(x,y) = " << "(" << x << "," << y << ")" << "     i = " << i << endl;
//		cout << "(newx,newy) = " << "(" << newx*h << "," << newy*h << ")" << "     i = " << i << endl;
//        x += h;
//		newx += 1;
//        if(newx ==n)
//        {
//            x = h;
//            y += h;
//			newx = 1;
//			newy += 1;
//        }
//        //========================================================================
//
//    }
//
//    cout << endl;
//    cout << "A = " << endl;
//    cout << A << endl;
//	cout << endl;
//	cout << "nA = " << endl;
//	cout << nA << endl;
//
//
//    b = b * h;
//    cout << "b = " << endl;
//    cout << b << endl;
//	nb = nb * h;
//	cout << "nb = " << endl;
//	cout << nb << endl;
//
	return 0;
}