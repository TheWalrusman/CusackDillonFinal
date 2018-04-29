template <typename T>
Math_Vector<T> GaussianSolver::operator () (const Matrix<T>& A, const Math_Vector<T>& b)
{
    if(A.getHeight() == 0 || A.getWidth() == 0)
        throw std::invalid_argument("   -Matrix A in GaussianSolver::operator() must have at least one column and one row");
    if(static_cast<unsigned>(b.getSize()) != A.getHeight())
        throw std::invalid_argument("   -Size of vector b must match height of matrix A in GaussianSolver::operator()");


   
    Matrix<T> augMatrix(A|b);                   //create augmented matrix
    unsigned cRow = 0;                          //the row you are trying to find a pivot for
    unsigned cColumn = 0;                       //the column you are trying to find a pivot for
    unsigned pivotRowNum;                       //the row number where the next pivot is before interchange
    vector<unsigned> pivColPos;              //keep track of where the pivots are
    vector<unsigned> pivRowPos;              //keep track of where the pivots are
    vector<T> pivValue;                      //the value of the pivots before RREF
    vector<T> pivValueRecips;                //recipricol values of the pivot value used to get from REF to RREF
    Math_Vector<T> x;                        //the solution to A*x=b
    bool solutionExists = true;

    //================================ USE GAUSSIAN ELIMINATION TO GET TO REF ====================
    //while a pivot has not been found for each row and column count has not
    //  run over the edge of the matix while tring to find a pivot for a row
    while(cColumn < augMatrix.getWidth() && cRow < augMatrix.getHeight())
    {
        pivotRowNum = getPivotRowNum(augMatrix,cRow,cColumn); //which row in column number cColumn at or below cRow has the greatest absolute value?

        if(augMatrix[pivotRowNum][cColumn] != 0)            //pivot values cannot be zero
        {
            augMatrix.rowInterchange(pivotRowNum,cRow);        //move the pivot value to the current row

            //for each row below row cRow in the augmented matric
            for(unsigned k = cRow + 1; k < augMatrix.getHeight(); k++)
            {
                //make the leading entry of each row 0 if it is not already 0
                //  via the use of row addition
                if(augMatrix[k][cColumn] != 0)
                {
                    augMatrix.rowAddition(k,cRow,-(augMatrix[k][cColumn]/augMatrix[cRow][cColumn]));
                }
            }

            pivColPos.push_back(cColumn);
            pivRowPos.push_back(cRow);
            pivValue.push_back(augMatrix[cRow][cColumn]);

            cRow++;     //a pivot was found for the row, so increment
        }

        cColumn++;      //always increment column even if no pivot was found
    }
    //============================= END GAUSSIAN ELIMINATION (REF ACHIEVED) =========================

    
    //============================== CHECK FOR NO SOLUTION BEFORE GETTING RREF ==================================
    //check for no solution before getting RREF
    for(unsigned i = 0; i < pivColPos.size(); i++)
    {
        if(pivColPos[i] == (augMatrix.getWidth() - 1))
        {   
            cout << "NO SOLUTION" << endl;
            solutionExists = false;
            break;
        }
    }
    //==========================================================================================

    
    //only continue if a there is a solution to the system
    if(solutionExists)
    {
        //===============================CONVERT TO REDUCED ROW ECHELON FORM (RREF) ====================
        //for each row in pivRowPos
        for(unsigned i = 0; i < pivValue.size(); i++)
        {
            pivValueRecips.push_back(1/pivValue[i]);        //get recipricols of pivot values
            //multiply each row that contains a pivot by the recipricol of the pivot in 
            //  order to convert all pivot values to 1 ...
            augMatrix.rowMult(pivRowPos[i],pivValueRecips[i]);   
        }    
        //*************EACH PIVOT VALUE IS NOW EQUAL TO 1****************************
        //eliminate non-zero numbers above the pivots .....
        //for each pivot position
        for(unsigned i = 0; i < pivRowPos.size(); i++)
        {
            //for each row above the pivot row
            for(int j = (pivRowPos[i]-1); j >= 0; j--)
            {
                if(augMatrix[j][pivColPos[i]] != 0)         //eliminate only non-zero numbers
                    augMatrix.rowAddition(j,pivRowPos[i],-augMatrix[j][pivColPos[i]]);  //eliminate non-zero's above 1's via Rj = Rj + *negative above value* * Ri
            }
            
        }
        //======================================================================================
        


        //============================= SYSTEM IS CONSITENT. CHECK FOR FREE VARIABLES =========================
        //check for free variables
        if(pivColPos.size() == (augMatrix.getWidth()-1)) //if no free variables exist (no piv columns == no columns in A)
        {
            //=============== BACK SUBSTITUTION ==============================================
            x.resize(augMatrix.getWidth()-1);   //prep the x vector to hold the solutions
            for(int i = (pivColPos.size()-1); i >= 0; i--) //for each pivot position
            {
                T sumOfKnowns = 0;
                //for each number in A to the left of each pivot value in row i
                for(unsigned j = (pivColPos[i] + 1); j < (augMatrix.getWidth()-1); j++)
                {
                    sumOfKnowns += x[j] * augMatrix[pivRowPos[i]][j];
                }

                x[i] = augMatrix[pivRowPos[i]][augMatrix.getWidth()-1] - sumOfKnowns;
            }
            //==================================================================================
        }
        else //otherwise, there are free variables
        {
            //========================== FIGURE OUT WHICH VARIABLES ARE FREE VARIABLES ===============
            vector<unsigned> freeColumnNum; //stores which columns are free
            //for each column of the augmented matrix except the last column
            for(unsigned i = 0; i < (augMatrix.getWidth()-1); i++)
            {
                //=======figure out which columns are not pivot columns========
                unsigned wCounter = 0;  //temp while loop counter
                bool isPivot = false;   //temp bool used in while loop. true if column i is a pivot column
                //iterate through pivColPos to see if column "i" is a pivot column
                while(wCounter < pivColPos.size() && !isPivot)
                {
                    //if column "i" is a pivot
                    if(pivColPos[wCounter] == i)
                    {
                        isPivot = true;
                    }

                    wCounter++; //check next column
                }

                if(!isPivot) //if column "i" was not a pivot, column "i" is free column
                {
                    freeColumnNum.push_back(i);
                }
                
            }
            //====================================================================================



            x.resize(augMatrix.getWidth()-1);   //prep the x vector to hold the solutions
            //================================= SET FREE VARIBALE VALUES TO ZERO ==========================
            for(unsigned i = 0; i < freeColumnNum.size(); i++)  //for each free variable
            {
                x[freeColumnNum[i]] = 0; //set each free value variable to 0
            }
            //=============================================================================
            

            //================================== REGULAR BACK SUBSTITUTION =======================
            //this for loop conducts regular back substitution
            for(int i = (pivColPos.size()-1); i >= 0; i--) //for each pivot position
            {
                T sumOfKnowns = 0;
                //for each number in A to the left of each pivot value in row i
                for(unsigned j = (pivColPos[i] + 1); j < (augMatrix.getWidth()-1); j++)
                {
                    sumOfKnowns += x[j] * augMatrix[pivRowPos[i]][j];
                }

                x[pivColPos[i]] = augMatrix[pivRowPos[i]][augMatrix.getWidth()-1] - sumOfKnowns;
            }
            //======================================================================================
        }
    }    

    return x;
}



template <typename T>
unsigned GaussianSolver::getPivotRowNum(const Matrix<T>& A, const unsigned rowNumber, const unsigned columnNumber) const
{
    if(A.getHeight() < 1 || A.getWidth() < 1)
        throw std::invalid_argument("     -In call to GaussianSolver::getPivotRowNum from GaussianSolver::operator(), matirx A must have at least one column and at least one row");
    if(rowNumber >= A.getHeight())
        throw std::invalid_argument("     -In call to GaussianSolver::getPivotRowNum from GaussianSolver::operator(), rowNumber must be less than the matrix height");
    if(columnNumber >= A.getWidth())
        throw std::invalid_argument("     -In call to GaussianSolver::getPivotRowNum from GaussianSolver::operator(), columnNumber must be less than matrix width");
    

    T maxAbsVal;            //the max absolute value found in a column
    unsigned rowNumOfMax;   //the row number in "A" where the maxAbsVal was found

    maxAbsVal = abs(A[rowNumber][columnNumber]);    //assume pivot is already in right place
    rowNumOfMax = rowNumber;
    for(unsigned i = rowNumber; i < A.getHeight() ; i++)    //for each value in column #columnNumber of A, look at each value at row #rowNumber and below to find the biggest
    {
        if(abs(A[i][columnNumber]) > maxAbsVal)             //if the abs of the entry is greater than the current max abs
        {
            maxAbsVal = abs(A[i][columnNumber]);            //update the max
            rowNumOfMax = i;                                //assume row is the pivot until found otherwise
        }
    }
    
    return rowNumOfMax;
}

