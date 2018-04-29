#include <iostream>
#include <fstream>
#include "absMatrix.h"
#include "denseMatrix.h"
#include "Math_Vector.h"
#include "steepestDescent.h"
#include "symmetricMatrix.h"
#include "GaussSeidel.h"
#include "GaussianInverse.h"
using namespace std;


int main()
{
    GaussSeidel gsSolver;
    //symmetricMatrix<double> tester(4);
    Math_Vector<double> blah(4);
    GaussianInverse gauInverse;
    denseMatrix<double> inverseTest(3);
    ifstream in;
    denseMatrix<double> invTest2(3);
    denseMatrix<double> invTest3(3);
    symmetricMatrix<double> A(10);
    Math_Vector<double> b(10);

    in.open("inverseTestCases.txt");
        in >> inverseTest;
        in >> invTest2;
        in >> invTest3;
    in.close();

    in.open("gaussSeidelTest.txt");
        in >> A;
        in >> b;
    in.close();


    cout << inverseTest << endl;
    cout << invTest2 << endl;
    cout << invTest3 << endl;

    cout << inverseTest * gauInverse(inverseTest) << endl;
    cout << invTest2 * gauInverse(invTest2) << endl;
    cout << invTest3 * gauInverse(invTest3) << endl;

    cout << endl << endl << endl;
    cout << "A = " << endl;
    cout << A << endl;
    cout << "B = " << endl;
    cout << b << endl;
    cout << "x (approximation) = " << endl;
    cout << gsSolver(A,b) << endl;


    return 0;
}