template <typename T>
Math_Vector<T> GaussSeidel::operator () (const symmetricMatrix<T>& A, const Math_Vector<T>& b) const
{
    Math_Vector<T> Xapprox(b.getSize());    //the approximation of vector x in Ax=b
    denseMatrix<T> L (A.getHeight());       //lower matrix part of matrix A
    denseMatrix<T> U (A.getHeight());       //upper matrix part of matrix A
    denseMatrix<T> negU;                    //the negative version of U
    denseMatrix<T> invL;                    //the matrix inverse of matrix L
    GaussianInverse inverse;                //find the inverse of L using gaussian elimination

    //generate lower matrix part of matrix A
    for(unsigned i = 0; i < A.getHeight(); i++)
    {
        for(unsigned j = 0; j <= i; j++)
            L[i][j] = A(i,j);
    }

    //generate upper matrix part of matrix A
    for(unsigned i = 0; i < A.getHeight(); i++)
    {
        for(unsigned j = i+1; j < A.getWidth(); j++)
            U[i][j] = A(i,j);
    }

    //make the zero vector the intial guess
    for(unsigned i = 0; i < static_cast<unsigned>(Xapprox.getSize()); i++)
        Xapprox[i] = 0;


    // cout << U << endl;
    // cout << L << endl;
    // cout << b << endl;

    negU = U * -1;
    invL = inverse(L).getRows();

    // cout << negU << endl;
    // cout << invL << endl;


    for(unsigned i = 0; i < 10; i++)
    {
        Xapprox = invL * ((negU*Xapprox) + b);
    }


    return Xapprox;
}


void GaussSeidel::setTollerance (const double toll)
{
    if(toll < 0) //if toll is negative
        throw std::invalid_argument("tollerance cannot be negative");

    m_tollerance = toll;

    return;
}