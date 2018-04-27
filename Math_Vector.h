/*!
    \file
    This is the header file for the Math_Vector Library
    \brief Header file for the Math_Vector<T> class

*/


#ifndef Math_Vector_H
#define Math_Vector_H

#include <iostream>
using namespace std;

/*!
    This is a templated Math_Vector class. The class accesses data through an
    Math_Vector pointer (m_dataPtr) that is of type T. The Math_Vector can be resized
    and cleared as needed by the user.

    \date 2/24/18
    \author Myles Dillon
    \brief Templated Math_Vector class

*/
template <class T>
class Math_Vector
{
    private:
        T* m_dataPtr;   //pointer to an Math_Vector
        int m_size;    //the size (AKA dimension) of the Math_Vector that m_dataPtr points to

    public:
        
        /*!
            Default Math_Vector Constructor
            \pre None
            \post m_dataPtr points to NULL, m_size == 0
        */
        Math_Vector();


        /*!
            Constructor the sets Math_Vector m_size to a specific size
            \pre None
            \post m_size == s 
                  m_dataPtr = new T[s]
        */
        Math_Vector(const unsigned s);
        

        /*!
            Copy Constructor
            \pre Since Math_Vector<T>::copy() is called, T = T must be defined
            \post m_size = Math_VectorToCopy.getSize()
                  m_dataPtr = new T[Math_VectorToCopy.getSize()]
                  The contents of *(Math_VectorToCopy.m_dataPtr) are copied into *m_dataPtr using Math_Vector<T>::copy
            \throws See Math_Vector<T>::copy(const Math_Vector<T>& arr)
        */
        Math_Vector(const Math_Vector<T>& Math_VectorToCopy);

        /*!
            Move contructor
            \pre none
            \post m_dataPtr is assigned other.m_dataPtr,
                  m_size is assigned other.m_size,
                  other.m_dataPtr is set to null
        */
        Math_Vector(Math_Vector<T>&& other);

        /*!
            Math_Vector Destructor
            \pre None
            \post The Math_Vector that m_dataPtr points to is delete via delete[] m_dataPtr
        */
        ~Math_Vector();

        
        /*!
            The right-hand overload of operator [] used to access m_dataPtr elements
            \pre i < m_size
            \post See return
            \return A const reference to m_dataPtr[i] if no exception is thrown, else a const reference to m_dataPtr[0] 
            \throws std::invalide_argument if i >= m_size
        */
        const T& operator [] (const unsigned i) const;
        

        /*!
            The left-hand overload of operator [] used to change m_dataPtr elements.
            \pre i < m_size
            \post See return
            \return A reference to m_dataPtr[i] if no exception is thrown, else a reference to m_dataPtr[0] 
            \throws std::invalide_argument if i >= m_size
        */
        T& operator [] (const unsigned i);


        /*!
            An insertion operator overload.
            \pre operator << defined for T
            \post Each element of 'arr' is inserted into 'os' followed by a whitespace.
                  If 'arr'.getSize() > 0 an endl is also inserted into 'os'.
            \return An ostream& to allow for operator << chaining
            \throws See const T& operator [] (const unsigned i) const
        */
        template<typename U>
        friend ostream& operator << (ostream& os, const Math_Vector<U>& arr);


        /*!
            Assignment operator.
            \pre T = T defined since Math_Vector<T>::copy is called
            \post m_dataPtr is deallocated and then reallocated via m_dataPtr = new T[m_size] if rhs.getSize() != 0
                    or set to NULL if rhs.getSize() == 0. 
                    m_size == rhs.getSize().
                    Contents of rhs is copied in this Math_Vector via a call to Math_Vector<T>::copy 
            \return A reference to the calling object for chaining
            \throws See Math_Vector<T>::copy(const Math_Vector<T>& arr)
        */
        Math_Vector<T>& operator = (const Math_Vector<T>& rhs);
        

        /*!
            Move Assignment operator
            \pre None
            \post std::swap called on m_dataPtr and rhs.m_dataPtr,
                  std::swap called on m_size and rhs.m_size
            \returns A reference to the calling object
        */
        Math_Vector<T>& operator = (Math_Vector<T>&& rhs);


        /*!
            Accessor for m_size.
            \pre None
            \post See return
            \return The value of m_size
        */
        inline int getSize() const {return m_size;}


        /*!
            Clears the contents of an Math_Vector
            \pre None
            \post Contents of m_dataPtr is delete via delete[] m_dataPtr,
                  m_dataPtr == NULL,
                  m_size == 0.
        */
        void clear();


        /*!
            Clears and resizes a Math_Vector
            \pre None
            \post Math_Vector<T> clear() is called,
                  m_dataPtr = new T[s],
                  m_size == s
        */
        void resize(const unsigned s);


        /*!
            A function that allows the user to copy the contents of one Math_Vector into
            the calling Math_Vector without delete what m_dataPtr points to first.
            \pre m_size >= arr.getSize(),
                 T = T (assignment) defined
            \post From i = 0 to i = (arr.getSize() - 1) m_dataPtr[i] == arr[i] is true
            \throws std::invalid_argument if m_size < arr.getSize()
        */
        void copy(const Math_Vector<T>& arr);



        /*!
            Vector subtraction operator
            \pre rhs must be the same size as the calling object, 
                 T - T (subtraction) defined
            \post See return
            \return A new Math_Vector that is the same size as the calling object is returned
            that has a value of m_dataPtr[i] - rhs[i] for each index i
            \throws std::invalid_argument if rhs size != calling object size
        */
        Math_Vector<T> operator - (const Math_Vector<T>& rhs) const;




        /*!
            Vector addition operator
            \pre rhs.m_size == m_size,
                 T * -1 (multiplication with literal -1) defined,
                 T - T (subraction) defined
            \post see return
            \return A Math_Vector that is the same size as the calling object such that
            for each in "i" in the returned vector, returnedVector[i] == m_data[i] + rhs.m_data[i];
            \throws std::invalid_argument if rhs.m_size != m_size    
        */
        Math_Vector<T> operator + (const Math_Vector<T>& rhs) const;

        



        /*!
            Scalar multiplication for a Math_Vector.
            \pre T * T (multiplication) defined
            \post See return
            \return An Math_Vector of the same size as the calling object is returned such that for
            each index i, returnedMath_Vector[i] == callingObject[i] * scalar
        */
        Math_Vector<T> operator * (const T& scalar) const;


        /*!
            Extraction operator for vectors
            \pre ifstream >> T (Extraction) defined
            \post reads in vect.m_size values from "in" and assigns the values to vect from vect[i] to vect[m_size-1]
            \return reference to "in" for chaining
        */
        template <typename U>
        friend ifstream& operator >> (ifstream& in, Math_Vector<U>& vect);


        /*!
            component-wise multiplication operator for vectors
            \pre m_size == rhs.m_size, 
                T * T (multiplication) defined,
                T = T (assignment) defined
            \post See return
            \return A Math_Vector of the same size as the calling object is returned such that
                for each index i, returnedMath_Vector[i] == m_dataPtr[i] * rhs.m_data[i]
        */
        Math_Vector<T> operator ^ (const Math_Vector<T>& rhs) const;


        /*!
            operator that returns the sum of the component values of a vector
            \pre m_size > 0,
                T += T (addition/assignment) defined
            \post see return
            \return The sum of all the values in the vector
            \throws std::invalid_argument if m_size == 0
        */
        T operator + () const;


        /*!
            function for calculating the inner product of two vectors
            \pre m_size == other size,
                T * T (mutliplication) defined,
                T + T (addtion) defined
            \post see return
            \return The inner product of the calling object and other is calculated and returned
            \throws std::invalid_argument if m_size != other size
        */
        T innerProduct(const Math_Vector<T>& other) const;


};

#include "Math_Vector.hpp"

#endif 