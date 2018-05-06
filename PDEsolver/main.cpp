#include <iostream>
#include <fstream>
#include <ctime>
#include <iomanip>
#include "absMatrix.h"
#include "denseMatrix.h"
#include "Math_Vector.h"
#include "steepestDescent.h"
#include "symmetricMatrix.h"
#include "GaussSeidel.h"
#include "GaussianInverse.h"
#include "CenterDiffMesh.h"
using namespace std;

double xLower(const double x, const double y)
{
	return 0 * x*y;
}

double xUpper(const double x, const double y)
{
	return 0 * x*y;
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
    /*
    ifstream in;


    GaussSeidel gsSolver;
    GaussianInverse gauInverse;

    denseMatrix<double> inverseTest(3);
    denseMatrix<double> invTest2(3);
    denseMatrix<double> invTest3(3);

    symmetricMatrix<double> A(10);
    Math_Vector<double> b(10);


	cout << argc <<argv[0]<< endl;
	
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
    */
	denseMatrix<double> AA;
	Math_Vector<double> bb;
	CenterDiffMesh<double, double, double, xLower, xUpper, yLower, yUpper> dirichlet;
	dirichlet.setSubdivisions(atoi(argv[1]));
	dirichlet(AA, bb);
    cout << setprecision(10);

    if(argc != 2)
        throw std::invalid_argument("Incorrect number of command line args");

    //1) get A and b dimensions from user
    int dimensions = atoi(argv[1]);
    cout << endl << "dimensions = " << dimensions << endl;

    //2) generate random A
    srand((unsigned int)time(NULL));
    symmetricMatrix<double> A(dimensions);
    for(unsigned i = 0; i < A.getHeight(); i++) //for each row of A
    {
        for(unsigned j = 0; j <= i; j++) //for each column up to the diagonal of the row
        {
            if(j == i)
            {
                double number = static_cast<double>((rand() % 10000) + 123456.0);
                double decimal = (static_cast<double>(rand()%100)/100.0f);
                cout << number << endl;
                cout << decimal << endl;
                cout << (number + decimal) << endl;
                A.getRef(i,j) = ((rand() % 10000) + 123456) + static_cast<double>(rand()%1000000)/1000000; //force A to be diagonally dominate
            }
            else
            {
                A.getRef(i,j) = ((rand() % 10000)) + static_cast<double>(rand()%1000000)/1000000;
            }
        }
    }
    
    //3) generate random b
    Math_Vector<double> b(dimensions);
    for(int i = 0; i < b.getSize(); i++)
        b[i] = ((rand() % 10000)) + static_cast<double>(rand()%1000000)/1000000;
    
    //4) calculate x approximation
    GaussSeidel gsSolver;
    Math_Vector<double> Xapprox(dimensions);
    Xapprox = gsSolver(A,b);

    //5) output A
    cout << endl << "A = " << endl;
    cout << A << endl;

    //6) output b
    cout << "b = " << endl;
    cout << b << endl;

    //7) output (A * Xapprox)   (This should b almost equal or exactly equal to b)
    cout << "(A * Xapprox) = " << endl;
    cout << (A*Xapprox) << endl;

    //8) output (b - (A * Xapprox))    (This should be close to or equal to the zero vector)
    cout << "(b - (A*Xapprox)) = " << endl;
    cout << (b - (A*Xapprox)) << endl;

    return 0;
}