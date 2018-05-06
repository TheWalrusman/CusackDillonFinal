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
#include "boundaryFunctions.h"
using namespace std;


int main(int argc, char* argv[])
{

    if(argc != 2)
        throw std::invalid_argument("The must be at exactly two command line arguments");


    int numOfSubDivisions = atoi(argv[1]);
    symmetricMatrix<double> A;
    Math_Vector<double> b;
    CenterDiffMesh<double,double,double, xLower, xUpper, yLower, yUpper> dirichlet;
    steepestDescent steepDescent;
    GaussSeidel gausSied;

    dirichlet.setSubdivisions(numOfSubDivisions);
    dirichlet(A,b);

    cout << "=== Steepest Descent ===" << endl;
    cout << "x = " << endl;
    cout << steepDescent(A,b) << endl;
    cout << "b = " << endl;
    cout << b << endl;
    cout << "A*x = " << endl;
    cout << A*steepDescent(A,b) << endl;

    cout << "=== Gauss Seidel ===" << endl;
    cout << "x = " << endl;
    cout << gausSied(A,b) << endl;
    cout << "b = " << endl;
    cout << b << endl;
    cout << "A*x = " << endl;
    cout << A*gausSied(A,b) << endl;


    // cout << "A = " << endl;
    // cout << A << endl;
    // cout << "b = " << endl;
    // cout << b << endl;

    /*
    //===================== TESTING GAUSS-SEIDEL WITH A RANDOM MATRIX A AND VECTOR B =========================
    cout << setprecision(10);

    if(argc != 2)
        throw std::invalid_argument("Incorrect number of command line args");

    //1) get A and b dimensions from user
    int dimensions = atoi(argv[1]);
    cout << endl << "dimensions = " << dimensions << endl;

    //2) generate random A
    srand(time(NULL));
    symmetricMatrix<double> A(dimensions);
    for(unsigned i = 0; i < A.getHeight(); i++) //for each row of A
    {
        for(unsigned j = 0; j <= i; j++) //for each column up to the diagonal of the row
        {
            if(j == i)
            {
                double number = static_cast<double>((rand() % 10000) + 123456.0f);
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
    //===========================================================================================
    */


    return 0;
}