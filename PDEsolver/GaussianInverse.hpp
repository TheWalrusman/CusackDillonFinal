template <typename T>
denseMatrix<T> GaussianInverse::operator () (const absMatrix<T>& A) const
{
    denseMatrix<T> Ainv(A.getHeight());     //the inverse of matrix A
    denseMatrix<T> AI(A.getHeight(),A.getWidth()*2); //The augmented matrix [A|I]
    denseMatrix<T> I(A.getHeight());                    //The identity matrix
    unsigned cRow = 0;                   //the row you are trying to find a pivot for
    unsigned cColumn = 0;                //the column you are trying to find a pivot for
    unsigned pivotRowNum;                //the row number where the next pivot is before interchange
    vector<unsigned> pivColPos;          //keep track of where the pivots are
    vector<unsigned> pivRowPos;          //keep track of where the pivots are
    vector<T> pivValue;                  //the value of the pivots before RREF
    vector<T> pivValueRecips;            //recipricol values of the pivot value used to get from REF to RREF
    bool notInvertable = false;         //true if matrix A is not invertable


    //build the identity matrix
    for(unsigned i = 0; i < I.getHeight(); i++) //for each row
    {
        for(unsigned j = 0; j < I.getWidth(); j++) //for each column
        {
            if(i == j)          //if entry is on the diagonal set it equal to 1
                I[i][j] = 1;
            else                //all non-diagonal entries are 0
                I[i][j] = 0;
        }
    }

    //create the augmented matrix
    for(unsigned i = 0; i < AI.getHeight(); i++) //for each row in AI
    {
        for(unsigned j = 0; j < AI.getWidth(); j++) //for each column in Ai
        {
            if(j < A.getWidth())   //the front columns of the augmented matrix are set to the values of matrix A
                AI[i][j] = A(i,j);
            else                    //the back columns are set to the values of the identity matrix
                AI[i][j] = I[i][j-A.getWidth()];
        }
    }

    //cout << AI << endl;
    
    //================================ USE GAUSSIAN ELIMINATION TO GET TO REF ====================
    //while a pivot has not been found for each row and column count has not
    //  run over the edge of the matix while tring to find a pivot for a row
    while(cColumn < AI.getWidth() && cRow < AI.getHeight())
    {
        pivotRowNum = getPivotRowNum(AI,cRow,cColumn); //which row in column number cColumn at or below cRow has the greatest absolute value?

        if(AI[pivotRowNum][cColumn] != 0)            //pivot values cannot be zero
        {
            AI.rowInterchange(pivotRowNum,cRow);        //move the pivot value to the current row

            //for each row below row cRow in the augmented matric
            for(unsigned k = cRow + 1; k < AI.getHeight(); k++)
            {
                //make the leading entry of each row 0 if it is not already 0
                //  via the use of row addition
                if(AI[k][cColumn] != 0)
                {
                    AI.rowAddition(k,cRow,-(AI[k][cColumn]/AI[cRow][cColumn]));
                }
            }

            pivColPos.push_back(cColumn);           //remember the pivot location for later use
            pivRowPos.push_back(cRow);              //remember the pivot location for later use 
            pivValue.push_back(AI[cRow][cColumn]);  //remember the pivot location for later use

            cRow++;     //a pivot was found for the row, so increment
        }

        cColumn++;      //always increment column even if no pivot was found
    }
    //============================= END GAUSSIAN ELIMINATION (REF ACHIEVED) =========================

    //cout << AI << endl;


    //check that each column in the "A" part of AI has a pivot because if it does
    //  not then matrix "A" is not invertable
    for(unsigned i = 0; i < A.getHeight(); i++)
    {
        if(i != pivColPos[i])
            notInvertable = true;
    }

    
    //if matrix A is invertable
    if(!notInvertable)
    {
        //===============================CONVERT TO REDUCED ROW ECHELON FORM (RREF) ====================
        //for each row in pivRowPos
        for(unsigned i = 0; i < pivValue.size(); i++)
        {
            pivValueRecips.push_back(1/pivValue[i]);        //get recipricols of pivot values
            //multiply each row that contains a pivot by the recipricol of the pivot in 
            //  order to convert all pivot values to 1 ...
            AI.rowMult(pivRowPos[i],pivValueRecips[i]);   
        }    
        //*************EACH PIVOT VALUE IS NOW EQUAL TO 1****************************
        //eliminate non-zero numbers above the pivots .....
        //for each pivot position
        for(unsigned i = 0; i < pivRowPos.size(); i++)
        {
            //for each row above the pivot row
            for(int j = (pivRowPos[i]-1); j >= 0; j--)
            {
                if(AI[j][pivColPos[i]] != 0)         //eliminate only non-zero numbers
                    AI.rowAddition(j,pivRowPos[i],-AI[j][pivColPos[i]]);  //eliminate non-zero's above 1's via Rj = Rj + *negative above value* * Ri
            }
            
        }
        //======================================================================================
    }
    else //if matrix A is not invertables
    {
        throw std::invalid_argument("In GaussianInverse::operator(), matrix A is not invertable");
    }

    //cout << AI << endl;

    //copy the right matrix in AI into Ainv because it is the inverse of matrix A
    for(unsigned i = 0; i < A.getHeight(); i++)
    {
        for(unsigned j = A.getWidth(); j < AI.getWidth(); j++)
            Ainv[i][j-A.getWidth()] = AI[i][j];
    }


    // cout << "Ainv = " << endl;
    // cout << Ainv << endl;


    return Ainv;
}


template <typename T>
unsigned GaussianInverse::getPivotRowNum(const absMatrix<T>& A, const unsigned rowNumber, const unsigned columnNumber) const
{
    if(A.getHeight() < 1 || A.getWidth() < 1)
        throw std::invalid_argument("     -In call to GaussianInverse::getPivotRowNum from GaussianSolver::operator(), matirx A must have at least one column and at least one row");
    if(rowNumber >= A.getHeight())
        throw std::invalid_argument("     -In call to GaussianInverse::getPivotRowNum from GaussianSolver::operator(), rowNumber must be less than the matrix height");
    if(columnNumber >= A.getWidth())
        throw std::invalid_argument("     -In call to GaussianInverse::getPivotRowNum from GaussianSolver::operator(), columnNumber must be less than matrix width");


    T maxAbsVal;            //the max absolute value found in a column
    unsigned rowNumOfMax;   //the row number in "A" where the maxAbsVal was found

    maxAbsVal = abs(A(rowNumber, columnNumber));    //assume pivot is already in right place
    rowNumOfMax = rowNumber;
    for(unsigned i = rowNumber; i < A.getHeight() ; i++)    //for each value in column #columnNumber of A, look at each value at row #rowNumber and below to find the biggest
    {
        if(abs(A(i, columnNumber)) > maxAbsVal)             //if the abs of the entry is greater than the current max abs
        {
            maxAbsVal = abs(A(i, columnNumber));            //update the max
            rowNumOfMax = i;                                //assume row is the pivot until found otherwise
        }
    }

    return rowNumOfMax;
}