template <typename T>
Math_Vector<T> steepestDescent::operator () (const symmetricMatrix<T>& A, const Math_Vector<T>& b) const
{
    if(!A.diagonallyDominate())
        throw std::invalid_argument("In steepestDescent<T>::operator() parameter A must be diagonally dominate");
    if(static_cast<unsigned>(b.getSize()) != A.getHeight())
        throw std::invalid_argument("In steepestDescent<T>::operator() length of vector b must == height of matrix A");

    Math_Vector<T> xApprox(A.getHeight());    //the aproximate value of x in Ax=b
    Math_Vector<T> r(b.getSize());           //the error in the approximation
    T a; 
    unsigned iterations = 0; //number of steepest decent iterations

    //guess that solution is the zero vector first
    for(unsigned i = 0; i < static_cast<unsigned>(xApprox.getSize()); i++)
    {
        xApprox[i] = 0; //generate dat zero vector
    }

    //try to get a better solution using steepest decent
    do
    {
        r = b - (A*xApprox);                            //error
        a = r.innerProduct(r) / (A*r).innerProduct(r);  //magical "a" from proof
        xApprox = xApprox + (r*a);                      //the new approximation generated from ancient magic passed on by many mathematicians for centuries
        iterations++;
    } while(iterations < m_maxIterations && sqrt(r.innerProduct(r)) > m_tollerance); //while not over max number of iterations and magnitude of remainder is not less than or equal to tolerance
      

    return xApprox;
}


void steepestDescent::setTollerance (const double toll)
{
    if(toll < 0) //if toll is negative
        throw std::invalid_argument("tollerance cannot be negative");

    m_tollerance = toll;

    return;
}
