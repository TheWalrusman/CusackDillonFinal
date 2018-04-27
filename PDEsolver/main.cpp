#include <iostream>
#include "Math_Vector.h"
#include "denseMatrix.h"
using namespace std;

double xLower (const double x, const double y)
{
    return 0*x*y;
}

double xUpper (const double x, const double y)
{
    return 0*x*y;//
}

double yLower (const double x, const double y)
{
    return ( 1 - (4*(x-0.5)*(x-0.5)) ) + 0*y; 
}

double yUpper (const double x, const double y)
{
    return 0*x*y;
}


int main(int argc, char* argv[])
{
    if(argc){}

    unsigned n = atoi(argv[1]);
    double h = (1.0/n);
    double x = h;
    double y = h;
    double xp;
    double xm;
    double yp;
    double ym;
    double bSum;
    bool firstOnBound;
    bool secOnBound;
    bool thirdOnBound;
    bool fourthOnBound;
    Math_Vector<double> b;
    denseMatrix<double> A((n-1)*(n-1));
    b.resize((n-1)*(n-1));

    //"i" is the row of the symmetric matrix
    for(unsigned i = 0; i < (n-1)*(n-1); i++)
    {

        //===== GENERATE THE PARAM VALUES USED IN EQUATION ======
        xp = x+h;
        xm = x-h;
        yp = y+h;
        ym = y-h;
        //====================================================== 


        //============ ASSUME NO TERMS ARE ON BOUNDARY AT FIRST =====
        firstOnBound = false;
        secOnBound = false;
        thirdOnBound = false;
        fourthOnBound = false;
        //======================================================


        //================= GENERATE THE NUMBER THAT NEEDS TO GO INTO THE B VECTOR =========
        bSum = 0;
        //see third and forth term
        if(x == 0)
        {
            bSum += xLower(x,ym);
            bSum += xLower(x,yp);
            thirdOnBound = true;
            fourthOnBound = true;
        }
        if(x == 1)
        {
            bSum += xUpper(x,ym);
            bSum += xUpper(x,yp);
            thirdOnBound = true;
            fourthOnBound = true;
        }
        //see first term
        if(xm == 0)
        {
            bSum += xLower(xm,y);
            firstOnBound = true;
        }
        if(xm == 1)
        {
            bSum += xUpper(xm,y);
            firstOnBound = true;
        }
        //see second term
        if(xp == 0)
        {
            bSum += xLower(xp,y);
            secOnBound = true;
        }
        if(xp == 1)
        {
            bSum += xUpper(xp,y);
            secOnBound = true;
        }
        //see first and second term
        if(y == 0)
        {
            bSum += yLower(xm,y);
            bSum += yLower(xp,y);
            firstOnBound = true;
            secOnBound = true;
        }
        if(y == 1)
        {
            bSum += yUpper(xm,y);
            bSum += yUpper(xp,y);
            firstOnBound = true;
            secOnBound = true;
        }
        //see third term
        if(ym == 0)
        {
            bSum += yLower(x,ym);
            thirdOnBound = true;
        }
        if(ym == 1)
        {
            bSum += yUpper(x,ym);
            thirdOnBound = true;
        }
        //see fourth term
        if(yp == 0)
        {
            bSum += yLower(x,yp);
            fourthOnBound = true;
        }
        if(yp == 1)
        {
            bSum += yUpper(x,yp);
            fourthOnBound = true;
        }

        b[i] = bSum;
        //===================== END GENERATE B VALUE ===========================

      
        //=================== FILL IN MATRIX =======================
        //mapping equation -----> col = ( ((y/h)-1) * (n-1) ) + (x/h) - 1 
        A[i][i] = 1; //diagonal is always 1
        unsigned col;
        if(!firstOnBound)
        {
            //xm   y
            col = static_cast<unsigned>( ( ((y/h)-1) * (n-1) ) + (xm/h) - 1 );
            A[i][col] = -h;
        }
        if(!secOnBound)
        {
            //xp   y
            col = static_cast<unsigned>( ( ((y/h)-1) * (n-1) ) + (xp/h) - 1 );
            A[i][col] = -h;
        }
        if(!thirdOnBound)
        {
            //x   ym
            col = static_cast<unsigned>( ( ((ym/h)-1) * (n-1) ) + (x/h) - 1 );
            A[i][col] = -h;
        }
        if(!fourthOnBound)
        {
            //x   yp
            col = static_cast<unsigned>( ( ((yp/h)-1) * (n-1) ) + (x/h) - 1 );
            A[i][col] = -h;
        }
        //=============================================================



        //======================== UPDATE X AND Y ============================
        //update "x" and "y" as needed
        cout << "(x,y) = " << "(" << x << "," << y << ")" << "     i = " << i << endl;
        x += h;
        if(x >=0.999999)
        {
            x = h;
            y += h;
        }
        //========================================================================

    }

    cout << endl;
    cout << "A = " << endl;
    cout << A << endl;


    b = b * h;
    cout << "b = " << endl;
    cout << b << endl;

    return 0;
}