template <typename T>
Math_Vector<T> GaussSeidel::operator () (const symmetricMatrix<T>& A, const Math_Vector<T>& b)
{
	using namespace std;
    if(!A.diagonallyDominate())
        throw std::invalid_argument("In GaussSeidel::operator(), matrix A is not diagonally dominate");
    if(static_cast<unsigned>(b.getSize()) != A.getHeight())
        throw std::invalid_argument("In steepestDescent<T>::operator() length of vector b must == height of matrix A");


    Math_Vector<T> Xapprox(b.getSize());    //the approximation of vector x in Ax=b
    denseMatrix<T> L (A.getHeight());       //lower matrix part of matrix A
    denseMatrix<T> U (A.getHeight());       //upper matrix part of matrix A
    denseMatrix<T> negU;                    //the negative version of U
    denseMatrix<T> invL;                    //the matrix inverse of matrix L
    GaussianInverse inverse;                //find the inverse of L using gaussian elimination
    Math_Vector<T> r(b.getSize());          //approximation error
    unsigned iterations = 0;                //number of iterations that have occured

    //generate lower matrix part of matrix A
	L = L * 0; //set all entries to zero to clear out junk
    for(unsigned i = 0; i < A.getHeight(); i++)
    {
		for (unsigned j = 0; j <= i; j++)
		{
			L[i][j] = A(i, j);
		}
    }
	
    //generate upper matrix part of matrix A
    U = U * 0; //set all entries to zero to clear out junk
    for(unsigned i = 0; i < A.getHeight(); i++)
    {
        for(unsigned j = i+1; j < A.getWidth(); j++)
            U[i][j] = A(i,j);
    }
	
    //make the zero vector the intial guess
    for(unsigned i = 0; i < static_cast<unsigned>(Xapprox.getSize()); i++)
        Xapprox[i] = 0;

    negU = U * -1;
    invL = inverse(L).getRows();

    do
    {
        Xapprox = invL * ((negU*Xapprox) + b); //calculate new approximation
        r = b - (A*Xapprox); //calculate approximation error vector
        iterations++;       
        m_finalIterations = iterations;
        m_finalResidualMag = sqrt(r.innerProduct(r));
    } while(iterations < m_maxIterations && sqrt(r.innerProduct(r)) > m_tollerance); //while not over max number of iterations and magnitude of remainder is not less than or equal to tolerance

    return Xapprox;
}


void GaussSeidel::setTollerance (const double toll)
{
    if(toll < 0) //if toll is negative
        throw std::invalid_argument("tolerance cannot be negative");

    m_tollerance = toll;

    return;
}