template <typename T>
denseMatrix<T>::denseMatrix()
{
    this->m_width = 0;
    this->m_height = 0;
}


template <typename T>
denseMatrix<T>::denseMatrix(const denseMatrix<T>& other)
{
    this->m_width = other.m_width;
    this->m_height = other.m_height;
    m_data = other.m_data;
}


template <typename T>
denseMatrix<T>::denseMatrix(const unsigned size)
{
    this->m_width = size;
    this->m_height = size;
    m_data.resize(size);
    for(int i = 0; i < m_data.getSize(); i++) //resize each row
    {
        m_data[i].resize(size);
    }
}

template <typename T>
denseMatrix<T>::denseMatrix(const Math_Vector<T>& vect)
{
    if(vect.getSize() == 0)
        throw std::invalid_argument("   -In Matrix(const Math_Vector<T>& vect), size of vect must be >= 1");

    this->m_height = vect.getSize();
    this->m_width = 1;
    m_data.resize(vect.getSize());
    for(int i = 0; i < m_data.getSize(); i++) //copy vector entries into the first and only column of matrix
    {
        m_data[i].resize(1);
        m_data[i][0] = vect[i];
    }
}


template <typename T>
Math_Vector<Math_Vector<T>> denseMatrix<T>::getColumns() const
{
    Math_Vector<Math_Vector<T>> matColumns(this->m_width); //vector of vector to hold the columns

    for(int i = 0; i < matColumns.getSize(); i++)   //for each column in c.o.
    {
        matColumns[i].resize(this->m_height);             //each column vector is equal to matrix height in length
        for(unsigned j = 0; j < this->m_height; j++)      //for each row in c.o.
        {
            matColumns[i][j] = m_data[j][i];        //fill in matColumns from left to right and top to bottom
        }
    }

    return matColumns;
}


template <typename T>
Math_Vector<Math_Vector<T>> denseMatrix<T>::getRows() const
{
    return m_data;
}


template <typename T>
const Math_Vector<T>& denseMatrix<T>::operator [] (const unsigned i) const
{
    if(i >= this->m_height)
        throw std::invalid_argument("   -i in denseMatrix::operator[] must be < the height of the matrix");
    return m_data[i];
}


template <typename T>
Math_Vector<T>& denseMatrix<T>::operator [] (const unsigned i)
{
    if(i >= this->m_height)
        throw std::invalid_argument("   -i in denseMatrix::operator[] must be < the height of the matrix");
    return m_data[i];
}


template <typename T>
Math_Vector<T> denseMatrix<T>::operator * (const Math_Vector<T>& vect) const
{
    if(static_cast<unsigned>(vect.getSize()) != this->m_width)
        throw std::invalid_argument("   -In denseMatrix<T>::operator * (const Math_Vector<T>& vect), size of vect must match calling objects height");
    if(vect.getSize() == 0)
        throw std::invalid_argument("   -In denseMatrix<T>::operator * (const Math_Vector<T>& vect), size of vect must be at least 1");

    denseMatrix<T> matVect (vect);   //explicit conversion
    denseMatrix<T> product((*this) * matVect);
    Math_Vector<T> vProd(this->m_height);
    for(unsigned i = 0; i < this->m_height; i++)
    {
        vProd[i] = product[i][0];
    }

    return vProd; //matrix mult  
}


template <typename T>
denseMatrix<T> denseMatrix<T>::operator * (const T& scalar) const
{
    denseMatrix<T> scaProduct(*this);

    for(unsigned i = 0; i < this->m_height; i++)
        scaProduct[i] = m_data[i] * scalar;
    
    return scaProduct;
}


template <typename T>
denseMatrix<T> denseMatrix<T>::operator + (const absMatrix<T>& rhs) const
{
    if(this->getWidth() != rhs.getWidth() || this->getHeight() != rhs.getHeight())
    throw std::invalid_argument("-In denseMatrix<T>::operator +, rhs must have same dimensions as calling object");

    denseMatrix<T> sum(this->getHeight()); //the sum of the two matrices
    sum = rhs.getRows();

    for(unsigned i = 0; i < this->getWidth(); i++) //for each column
    {
        for(unsigned j = 0; j < this->getHeight(); j++) //for each row
        {
            sum[i][j] += (*this)(i,j); //add elements
        }
    }

    return sum;
}


template <typename T>
denseMatrix<T> denseMatrix<T>::operator - (const absMatrix<T>& rhs) const
{
    if(this->getWidth() != rhs.getWidth() || this->getHeight() != rhs.getHeight())
    throw std::invalid_argument("-In denseMatrix<T>::operator -, rhs must have same dimensions as calling object");

    denseMatrix<T> diff(this->getHeight()); //the sum of the two matrices
    diff = rhs.getRows();

    for(unsigned i = 0; i < this->getWidth(); i++) //each column
    {
        for(unsigned j = 0; j < this->getHeight(); j++) //each row
        {
            diff[i][j] = (*this)(i,j) - diff[i][j];
        }
    }

    return diff;
}


template <typename T>
denseMatrix<T> denseMatrix<T>::operator * (const absMatrix<T>& rhs) const
{
    if(this->getWidth() != rhs.getHeight())
    throw std::invalid_argument("-In denseMatrix matrix multiplication, rhs must have the same number of rows as there are columns in the calling object");

    denseMatrix<T> product(this->getWidth());
    Math_Vector<Math_Vector<T>> coRows;
    Math_Vector<Math_Vector<T>> rhsCols;
    coRows = this->getRows();
    rhsCols = rhs.getColumns();

    for(int i = 0; i < coRows.getSize(); i++) //for each row in c.o.
    {
        for(int j = 0; j < rhsCols.getSize(); j++) //for each col in rhs
        {
            product[i][j] = +(coRows[i]^rhsCols[j]); //multiply and sum each row with each column
        }
    }

    return product;
}



template <typename T>
T denseMatrix<T>::operator () (const unsigned r, const unsigned c) const
{
    if(r >= this->m_height)
        throw std::invalid_argument("-r in denseMatrix<T>::operator () must be < height of the matrix");
    if(c >= this->m_width)
        throw std::invalid_argument("-c in denseMatrix<T>::operator () must be < width of matrix");
    
    return m_data[r][c];
}


template <typename T>
void denseMatrix<T>::operator = (const Math_Vector<Math_Vector<T>>& rhs)
{
    if(rhs.getSize() < 1)
        throw std::invalid_argument("In denseMatrix<T> operator = , Math_Vector<Math_Vector<T>> rhs must have at least one entry");
    if(rhs[0].getSize() < 1)
        throw std::invalid_argument("In denseMatrix<T> operator = , Math_Vector<Math_Vector<T>> rhs must have at least one entry");
    
    this->m_height = rhs.getSize();
    this->m_width = rhs[0].getSize();
    m_data = rhs;
    return;
}



template <typename T>
ifstream& operator >> (ifstream& in, denseMatrix<T>& matrix)
{
    in >> matrix.m_data; //read in the vectors making up the matrix
    return in;
}


