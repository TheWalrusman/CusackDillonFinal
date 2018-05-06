template <typename T>
symmetricMatrix<T>::symmetricMatrix()
{
    this->m_width = 0;
    this->m_height = 0;
    m_data.resize(0);
}


template <typename T>
symmetricMatrix<T>::symmetricMatrix(const unsigned size)
{
    if(size < 1)
        throw std::invalid_argument("-size must be >= 1 in symmetricMatrix(const unsigned size)");

    this->m_width = size;
    this->m_height = size;
    m_data.resize(size);
    for(unsigned i = 0; i < size; i++)
        m_data[i].resize(i+1);
}


template <typename T>
T symmetricMatrix<T>::operator () (const unsigned r, const unsigned c) const
{
    if(r >= this->getHeight())
        throw std::invalid_argument("-r must be < matrix width in diagonalMatrix<T>::operator()");
    if(c >= this->getWidth())
        throw std::invalid_argument("-c must be < matrix height in diagonalMatrix<T>::operator()");

    T retVal;

    if(c <= r) //if column number is less than row number
    {
        retVal = m_data[r][c];
    }
    else
    {
        retVal = m_data[c][r];
    }

    return retVal;

}



template <typename T>
denseMatrix<T> symmetricMatrix<T>::operator + (const absMatrix<T>& rhs) const
{
    if(this->getWidth() != rhs.getWidth() || this->getHeight() != rhs.getHeight())
        throw std::invalid_argument("-In symmetricMatrix<T>::operator +, rhs must have same dimensions as calling object");

    denseMatrix<T> sum(this->getHeight()); //the sum of the two matrices
    sum = rhs.getRows();

    for(unsigned i = 0; i < this->getHeight(); i++) //for each row in c.o.
    {
        for(unsigned j = 0; j <= i; j++) //for some columns c.o.
        {
            if(i == j) //if diagonal entry
            {
                sum[i][j] += (*this)(i,j);
            }
            else
            {
                sum[i][j] += (*this)(i,j);
                sum[j][i] += (*this)(i,j);
            }
        }
    }

    return sum;
}


template <typename T>
denseMatrix<T> symmetricMatrix<T>::operator - (const absMatrix<T>& rhs) const
{
    if(this->getWidth() != rhs.getWidth() || this->getHeight() != rhs.getHeight())
        throw std::invalid_argument("-In symmetricMatrix<T>::operator -, rhs must have same dimensions as calling object");

    denseMatrix<T> diff(this->getHeight()); //the sum of the two matrices
    diff = rhs.getRows();

    for(unsigned i = 0; i < this->getHeight(); i++) //for each row in c.o.
    {
        for(unsigned j = 0; j <= i; j++) //for some columns c.o.
        {
            if(i == j) //if diagonal entry
            {
                diff[i][j] = (*this)(i,j) - rhs(i,j);
            }
            else
            {
                diff[i][j] = (*this)(i,j) - rhs(i,j);
                diff[j][i] = (*this)(i,j) - rhs(j,i);
            }
        }
    }

    return diff;

}


template <typename T>
denseMatrix<T> symmetricMatrix<T>::operator * (const absMatrix<T>& rhs) const
{
    if(this->getWidth() != rhs.getHeight())
        throw std::invalid_argument("-In symmetricMatrix matrix multiplication, rhs must have the same number of rows as there are columns in the calling object");

    denseMatrix<T> product(this->getWidth());
    Math_Vector<Math_Vector<T>> coRows; //calling object rows
    Math_Vector<Math_Vector<T>> rhsCols;    //rhs columns
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
symmetricMatrix<T> symmetricMatrix<T>::operator * (const T& scalar) const
{
    symmetricMatrix<T> product(this->getHeight());

    //for each row
    for(unsigned i = 0; i < this->getHeight(); i++)
    {
        //for some columns
        for(unsigned j = 0; j <= i; j++)
        {
            product.m_data[i][j] = m_data[i][j] * scalar;
        }
    }

    return product;
}


template <typename T>
Math_Vector<T> symmetricMatrix<T>::operator * (const Math_Vector<T>& vect) const
{
    if(this->getWidth() != static_cast<unsigned>(vect.getSize()))
        throw std::invalid_argument("-In symmetricMatrix vector multiplication, vect must have the same number of rows as there are columns in the calling object");

    Math_Vector<T> product(vect.getSize());
    Math_Vector<Math_Vector<T>> rows;

    rows = this->getRows();

    //for each row in the matrix
    for(unsigned i = 0; i < static_cast<unsigned>(rows.getSize()); i++)
    {
        product[i] = +(rows[i]^vect);
    }

    return product;
}


template <typename T>
Math_Vector<Math_Vector<T>> symmetricMatrix<T>::getRows() const
{
    Math_Vector<Math_Vector<T>> rows;

    rows.resize(this->getHeight());

    //for each row in the matrix
    for(unsigned i = 0; i < this->getHeight(); i++)
    {
        rows[i].resize(this->getWidth());
        //for each column
        for(unsigned j = 0; j < this->getWidth(); j++)
        {
            rows[i][j] = (*this)(i,j);
        }
    }

    return rows;
}


template <typename T>
Math_Vector<Math_Vector<T>> symmetricMatrix<T>::getColumns() const
{
    Math_Vector<Math_Vector<T>> cols;

    cols.resize(this->getWidth());

    //for each col in the matrix
    for(unsigned i = 0; i < this->getWidth(); i++)
    {
        cols[i].resize(this->getHeight());
        //for each row
        for(unsigned j = 0; j < this->getWidth(); j++)
        {
            cols[i][j] = (*this)(j,i);
        }
    }

    return cols;
}


template <typename T>
T& symmetricMatrix<T>::getRef(const unsigned r, const unsigned c)
{
    if(c > r)
        throw std::invalid_argument("-In symmetricMatrix<T>::getRef(), c must be <= r");
    if(c >= this->m_width)
        throw std::invalid_argument("-In symmetricMatrix<T>::getRef(), c must be < matrix width");
    if(r >= this->m_height)
        throw std::invalid_argument("-In symmetricMatrix<T>::getRef(), r must be < m_height");
    
    return m_data[r][c];
}





template <typename U>
istream& operator >> (ifstream& in, symmetricMatrix<U>& matrix)
{
    U junk;
    for(unsigned i = 0; i < matrix.getHeight(); i++) //for each row
    {
        for(unsigned j = 0; j < matrix.getWidth(); j++) //for each column
        {
            if(j <= i) //if in lower part of matrix
            {
                in >> matrix.m_data[i][j];
            }
            else
            {
                in >> junk;   
            }
        }
    }

    return in;
}


