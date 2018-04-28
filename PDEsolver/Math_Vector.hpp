/*
    Author: Myles Dillon
    Date: 2/24/18
    Purpose: Implementation file for the Math_Vector class
*/

template <typename T>
Math_Vector<T>::Math_Vector()
{
    m_dataPtr = NULL;
    m_size = 0;
}



template <typename T>
Math_Vector<T>::Math_Vector(const unsigned s)
{
    m_dataPtr = new T[s];
    m_size = s;
}



template <typename T>
Math_Vector<T>::Math_Vector(const Math_Vector<T>& Math_VectorToCopy)
{
    m_size = Math_VectorToCopy.getSize();
    m_dataPtr = new T[m_size];
    copy(Math_VectorToCopy); //copy Math_Vector contents
}


template <typename T>
Math_Vector<T>::Math_Vector(Math_Vector<T>&& other)
{
    m_dataPtr = other.m_dataPtr;
    m_size = other.m_size;
    other.m_dataPtr = NULL;
}


template <typename T>
Math_Vector<T>::~Math_Vector()
{
    delete[] m_dataPtr;
}


template<typename T>
ostream& operator << (ostream& os, const Math_Vector<T>& arr)
{
    for(int i = 0; i < arr.getSize(); i++) //print Math_Vector elements seperated by white space
    {
        os << arr[i] << endl;
    }
    
    return os;
}



template <typename T>
void Math_Vector<T>::resize(const unsigned s)
{
    clear();
    m_dataPtr = new T[s];
    m_size = s;

    return;
}




template <typename T>
void Math_Vector<T>::clear()
{
    delete[] m_dataPtr;
    m_dataPtr = NULL;
    m_size = 0;

    return;
}




template <typename T>
Math_Vector<T>& Math_Vector<T>::operator = (const Math_Vector<T>& rhs)
{
    if(rhs.m_dataPtr != m_dataPtr)
    {
        delete[] m_dataPtr;
        m_size = rhs.getSize();

        if(m_size == 0)
        {
            m_dataPtr = NULL;
        }
        else
        {
            m_dataPtr = new T[m_size];
            copy(rhs);
        }
    }

    return *this;

}




template <typename T>
Math_Vector<T>& Math_Vector<T>::operator = (Math_Vector<T>&& rhs)
{
    std::swap(m_dataPtr,rhs.m_dataPtr);
    std::swap(m_size,rhs.m_size);
    return *this;
}




template <typename T>
void Math_Vector<T>::copy(const Math_Vector<T>& arr)
{
    if(m_size >= arr.getSize())
    {
        for(int i = 0; i < arr.getSize(); i++)
        {
            m_dataPtr[i] = arr[i];
        }
    }
    else
    {
        throw std::invalid_argument("  -Calling Math_Vector must be the same m_size or bigger than the Math_Vector it is trying to copy");          
    }

    return;
}



template <typename T>
const T& Math_Vector<T>::operator [] (const unsigned i) const 
{
    if(i >= static_cast<unsigned>(m_size))
        throw std::invalid_argument("  - i cannot be greater than or equal to Math_Vector size");

    return m_dataPtr[i];
}




template <typename T>
T& Math_Vector<T>::operator [] (const unsigned i)
{
    //if(i >= static_cast<unsigned>(m_size))
        //throw std::invalid_argument("  - i cannot be greater than or equal to Math_Vector size");

    return m_dataPtr[i];
}






template <typename T>
Math_Vector<T> Math_Vector<T>::operator - (const Math_Vector<T>& rhs) const
{
    if(m_size != rhs.getSize())
        throw std::invalid_argument("Math_Vector sizes don't match in Math_Vector<T>::operator = ");

    Math_Vector<T> diffMath_Vector(m_size);     
        
    for(int i = 0; i < rhs.getSize(); i++)      //subtract each element
    {
        diffMath_Vector[i] = m_dataPtr[i] - rhs[i];
    }

    return std::move(diffMath_Vector);
}


template <typename T>
Math_Vector<T> Math_Vector<T>::operator + (const Math_Vector<T>& rhs) const
{
    if(rhs.getSize() != m_size) //if sizes of two vector don't match
        throw std::invalid_argument("   -Vectors must be of equal length for vector addition to work");

    //Math_Vector<T> coCopy(*this);   //copy of calling object
    return ((*this) - (rhs * -1));   //multiply rhs by -1 and use vect sub to create addition
}



template <typename T>
Math_Vector<T> Math_Vector<T>::operator * (const T& scalar) const
{
    Math_Vector<T> scaledArr(m_size);
    for(int i = 0; i < m_size; i++) //multiply each element by the scalar
    {
        scaledArr[i] = m_dataPtr[i] * scalar;
    }
    return std::move(scaledArr);
}



template <typename T>
ifstream& operator >> (ifstream& in, Math_Vector<T>& vect)
{
    for(int i = 0; i < vect.m_size; i++) //read in m_size number of vals from input stream
    {
        in >> vect[i];
    }

    return in;
}


template <typename T>
Math_Vector<T> Math_Vector<T>::operator ^ (const Math_Vector<T>& rhs) const
{
    if(m_size != rhs.getSize())
        throw std::invalid_argument("   -Vectors must be of the same length for component-wise multiplication (operator ^) to function correctly");

    Math_Vector<T> product(m_size);                 //product of the two vectors

    for(int i = 0; i < product.getSize(); i++)      //for each index of product
    {
        product[i] = m_dataPtr[i] * rhs[i];         //store product of i'th components
    }

    return std::move(product);
}


template <typename T>
T Math_Vector<T>::operator + () const
{
    if(m_size == 0)
        throw std::invalid_argument("   -Vector sum operator (operator +) requires as least one component to find the sum");

    T sum = 0;
    for(int i = 0; i < m_size; i++) //get sum of components
    {
        sum += m_dataPtr[i]; 
    }
    return sum;
}



template <typename T>
T Math_Vector<T>::innerProduct(const Math_Vector<T>& other) const
{
    if(m_size != other.getSize())
        throw std::invalid_argument("Vector sizes must match in Math_Vector<T>::innerProduct");

    return +((*this)^other);
}