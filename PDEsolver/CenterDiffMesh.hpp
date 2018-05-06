#ifndef CENTERDIFFMESH_HPP
#define CENTERDIFFMESH_HPP


template< class T, class U, class V, T xLower(const U, const V), T xUpper(const U, const V), T yLower(const U, const V), T yUpper(const U, const V)>
void CenterDiffMesh<T, U, V, xLower, xUpper, yLower, yUpper>::setSubdivisions(const unsigned numSubs)
{
    if(numSubs < 3)
        throw std::invalid_argument("In CenterDiffMesh::setSubdivisions, numSubs must be >= 3");

    m_subDivisions = numSubs;
}


template< class T, class U, class V, T xLower(const U, const V), T xUpper(const U, const V), T yLower(const U, const V), T yUpper(const U, const V)>
void CenterDiffMesh<T, U, V, xLower, xUpper, yLower, yUpper>::operator()(symmetricMatrix<T>& A, Math_Vector<T>& B)
{
    symmetricMatrix<double> correctA((m_subDivisions - 1)*(m_subDivisions - 1));
	double h = 1.0 / (m_subDivisions);
	double x = h;  //the x val of the point you are trying to approximate
	double y = h;   //the y val of the point you are trying to apporixmate
	double xp;      //== x+h
	double xm;      //== x-h
	double yp;      //== y+h
	double ym;      //== y-h
    double bSum;        //value to put in vector b
	bool firstOnBound;  //true if first bracketed is on bound
	bool secOnBound;    // dido but for second
	bool thirdOnBound;  // dido but for third
    bool fourthOnBound; // dido but for fourth
    bool iIsOnRightEdge;  //true if "i" is trying to calculate a point that is on the right edge of the interior of the bound 
    unsigned interiorColumn = 1;  //a value that is equal to x/h at any given point, used to map into A
    unsigned interiorRow = 1;       //a value that is equalt to y/h at any given point, used to map into A

    //for debugging ...
    //unsigned pointsThatHityLower = 0;

    A = correctA;
    B.resize((m_subDivisions - 1)*(m_subDivisions - 1));


	for (unsigned i = 0; i < (m_subDivisions - 1)*(m_subDivisions - 1); i++) //for each interior point withing the bounds
	{

		//===== GENERATE THE PARAM VALUES USED IN EQUATION ======
		xp = x + h;
		xm = x - h;
		yp = y + h;
		ym = y - h;
		//====================================================== 


		//============ ASSUME NO TERMS ARE ON BOUNDARY AT FIRST =====
		firstOnBound = false;
		secOnBound = false;
		thirdOnBound = false;
        fourthOnBound = false;
        iIsOnRightEdge = false; //assume current point not on rigth edge
		//======================================================


		//================= GENERATE THE NUMBER THAT NEEDS TO GO INTO THE B VECTOR =========
        //===================================== CHECK FOR SPECIAL VALUES OF I ==================================================
        if(i == 0) //bottom left corner
        {
            firstOnBound = true;
            thirdOnBound = true;
        }
        if(i == (m_subDivisions-2)) //bottom right corner
        {
            secOnBound = true;
            thirdOnBound = true;
            iIsOnRightEdge = true;
        }
        if(i == ( ( ((m_subDivisions-1)*(m_subDivisions-1)) -1) - (m_subDivisions-2) ) ) //top left corner
        {
            firstOnBound = true;
            fourthOnBound = true;
        }
        if(i == ( ((m_subDivisions-1)*(m_subDivisions-1)) - 1) ) //top right corner
        {
            secOnBound = true;
            fourthOnBound = true;
            iIsOnRightEdge = true;
        }
        if(i > 0 && i < (m_subDivisions-2))  //bottom edge
        {
            thirdOnBound = true;
        }
        if(i > ( ( ((m_subDivisions-1)*(m_subDivisions-1)) -1) - (m_subDivisions-2) ) && i < ( ((m_subDivisions-1)*(m_subDivisions-1)) -1) ) //top edge
        {
            fourthOnBound = true;
        }
        if(i != 0) //this must also be here to prevent a bug when i == 0
        {
            if(i % (m_subDivisions - 1) == 0) //left edge
            {
                firstOnBound = true;
            }
        }
        if(i > (m_subDivisions-2)) //since "i" is unsigned this if statement must be here to prevent a bug
        {
            if( ( (i-(m_subDivisions-2)) % (m_subDivisions-1)) == 0) //right edge
            {
                secOnBound = true;

                iIsOnRightEdge = true;
            }
        }
        
        //================ END CHECK FOR SPECIAL I VALUES ======================================================


        //================ FILL IN THE B VECTOR =========================================
        bSum = 0;
        if(firstOnBound)
        {
            bSum += xLower(xm,y);
        }
        if(secOnBound)
        {
            bSum += xUpper(xp,y);
        }
        if(thirdOnBound)
        {
            bSum += yLower(x,ym);
            //for debugging...
            //pointsThatHityLower++;
        }
        if(fourthOnBound)
        {
            bSum += yUpper(x,yp);
        }
        B[i] = bSum;
        //====================== END FILL IN THE B VECTOR
		//===================== END GENERATE B VALUE ===========================



        //=================== FILL IN MATRIX =======================
        //mapping equation -----> col = ( ((y/h)-1) * (n-1) ) + (x/h) - 1 
        A.getRef(i,i) = 1; //diagonal is always 1
        unsigned col;
        if(!firstOnBound)
        {
            //xm   y
            col = static_cast<unsigned>( ( (interiorRow-1) * (m_subDivisions-1) ) + (interiorColumn-1) - 1 );            
            if(col <= i)
                A.getRef(i,col) = -h;
        }
        if(!secOnBound)
        {
            //xp   y
            col = static_cast<unsigned>( ( (interiorRow-1) * (m_subDivisions-1) ) + (interiorColumn+1) - 1 );
            if(col <= i)
                A.getRef(i,col) = -h;
        }
        if(!thirdOnBound)
        {
            //x   ym
            col = static_cast<unsigned>( ( ((interiorRow-1)-1) * (m_subDivisions-1) ) + (interiorColumn) - 1 );
            if(col <= i)
                A.getRef(i,col) = -h;
        }
        if(!fourthOnBound)
        {
            //x   yp
            col = static_cast<unsigned>( ( ((interiorRow+1)-1) * (m_subDivisions-1) ) + (interiorColumn) - 1 );
            if(col <= i)
                A.getRef(i,col) = -h;
        }
        //=============================================================

        

		//======================== UPDATE X AND Y ============================
		//cout << "(x,y) = " << "(" << x << "," << y << ")" << "     i = " << i << endl;
        //cout << "(interiorRow, interiorCol) = " << "(" << interiorRow << ", " << interiorColumn << ")" << endl;
        x += h;
        interiorColumn++;
        //if "i" is a right edge interior point
		if (iIsOnRightEdge)
		{
            x = h;
            y += h;
            interiorColumn = 1;  //allows mapping via ints to A
            interiorRow++;          //allows mapping via ints to A
		}
		//========================================================================

    }
    

    //for debugging ...
    //cout << "pointsThatHityLower = " << pointsThatHityLower << " ... " << (m_subDivisions - 1) << endl;

    return;

}


#endif // !CENTERDIFFMESH_H