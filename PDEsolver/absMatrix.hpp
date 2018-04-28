template <typename U>
ostream& operator << (ostream& os, const absMatrix<U>& matrix)
{
    for(unsigned i = 0; i < matrix.getHeight(); i++)
    {
        for(unsigned j = 0; j < matrix.getWidth(); j++)
        {
            os << matrix(i,j) << " ";
        }

        os << endl;
    }

    return os;
}


template <typename T>
bool absMatrix<T>::diagonallyDominate () const
{
    bool diagDominate = true;
    T rowSum; //sum of the absolute values of each row except the diagonal element

    //for each row in matrix
    for(unsigned i = 0; i < this->getHeight(); i++)
    {
        rowSum = 0;
        //for each column in row i
        for(unsigned j = 0; j < this->getWidth(); j++)
        {
            rowSum += abs((*this)(i,j)); //sum the absolute values of the elements in the row
        }
        rowSum -= abs((*this)(i,i)); //subtract out the diagonal element

        if(abs((*this)(i,i)) <= rowSum) //if row is not diagonally dominate
        {
            diagDominate = false;
            break;
        }
    }

    return diagDominate;
}