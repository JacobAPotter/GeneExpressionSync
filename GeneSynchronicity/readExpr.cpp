/****************************************************************************
*
*	readExpr.cpp:	Code for reading in expression values and storing in matrix 
*                    
*                       Sharlee Climer
*                       March 2019
*
*	Edited by:
*	Jacob Potter
*	March 2019
*	Modified so that 'main' is now a function which returns the normalized data
*
****************************************************************************/
  


#include "readExpr.h"

float** ReadAndNormalize(char * file , int numGenes, int numInd)
{
  
  timer t;
  t.start("Timer started.");

  FILE *input;

  if ((input = fopen(file, "r")) == NULL)
    fatal("Input file could not be opened.\n");
  
  char strng[STRNGLENGTH];  // string for holding temporary values

  std::cout << "\n***Important: Assumed columns represent individuals and rows represent genes/proteins in input file.***\n" << std::endl;

  std::cout << numInd << " individuals and " << numGenes << " genes." << std::endl;

  std::cout << "Assumed " << NUMHEADCOLS << " header column(s) and " << NUMHEADROWS << " header row(s) in '" << file << "'." << std::endl;
  
  std::cout << "**Missing values must be indicated by NN or NA.**\n" << std::endl;
   
  if ((numInd < 1) || (numGenes < 1))
	  fatal("Too few individuals or genes");

  // allocate data array memory
  float **data; // hold expression values: MIN_EXP - 1 will indicate missing value 
  if ((data = new float*[numGenes]) == NULL)
    fatal("memory not allocated");
  for (int i = 0; i < numGenes; i++)
    if ((data[i] = new float[numInd]) == NULL)
      fatal("memory not allocated");

  // initialize data values to MIN_EXP-1, indicating missing values
  for (int i = 0; i < numGenes; i++)
    for (int j = 0; j < numInd; j++)
      data[i][j] = MIN_EXP - 1;

  // allocate memory for number of individuals without missing data (per gene)
  int *numIndNoMiss; 
  if ((numIndNoMiss = new int[numGenes]) == NULL)
    fatal("memory not allocated");

  // initialize values to all individuals
  for (int i = 0; i < numGenes; i++)
    numIndNoMiss[i] = numInd;

  // read in data
  float val; // value read in
  long int numMissing = 0; // count number of missing values
  float maxExpr = MIN_EXP; // initialize maximum expression value
  float minExpr = MAX_EXP; // initialize minimum expression value

  // read in header rows   
  for (int i = 0; i < NUMHEADROWS; i++)
    for (int j = 0; j < NUMHEADCOLS + numInd; j++)
      fscanf(input, "%s", strng); 
      
  // read in data rows
  for (int i = 0; i < numGenes; i++) {
    for (int j = 0; j < NUMHEADCOLS; j++) 
      fscanf(input, "%s", strng); // read in header columns

    for (int j = 0; j < numInd; j++) {
      if(feof(input))
	      fatal("Input file is missing data");

      fscanf(input, "%s", strng);
      val = atof(strng);

      // check that strng holds a number
      int ascii = (int)strng[0]; // ascii value of first character

      if((ascii < 48) || (ascii > 57)) { // not a digit
	      if((ascii != 45) && (ascii != 46)) { // not a '-' or '.'
	  
	        if(ascii == 78) { // starts with 'N' for missing data
	          if ((strncmp(strng, "NA", 2) == 0) || (strncmp(strng, "NN", 2) == 0)){
	            numMissing++; // add to count
	            val = MIN_EXP-1; // set missing value 
	            numIndNoMiss[i]--; // one less indiv. with no missing data for this gene
	          }
                     
	          else {
	          std::cout << strng << std::endl;
	          fatal("Invalid expression value in input file");
	          }
	        }
	
	        else {
	          std::cout << strng << std::endl;
	          fatal("Invalid expression value in input file");
	        }
	      }
      }

      if((ascii >= 48) && (ascii <= 57)) { // is a digit
	      if ((val > MAX_EXP) || (val < MIN_EXP)) {
	        std::cout << strng << std::endl;
	        fatal ("Invalid expression value in input data.  \n\tCheck MAX_EXP and MIN_EXP in header file.");
	      }
      }

      if((ascii == 45) || (ascii == 46)) // is a '-' or '.'
	      if ((val > MAX_EXP) || (val < MIN_EXP)) {
	        std::cout << strng << std::endl;
	        fatal ("Invalid expression value in input data.  \n\tCheck MAX_EXP and MIN_EXP in header file.");
	    }
      
      // update min and max, if appropriate
      if(val > maxExpr)
	      maxExpr = val;

      if(val < minExpr)
        if (val > MIN_EXP - 1 + TOL) // don't use missing value symbol
	        minExpr = val;

      // record expression in data matrix (with missing values encoded as MIN_EXP - 1)
      data[i][j] = val;
    } 
  } // end of read in data rows loop

  // check for end of file
  if (!feof(input)){
    fscanf(input, "%s", strng);
    if (!feof(input))
      fatal("Unread data in input file");
  }

  fclose(input); 

  std::cout << "\nExpression values range from " << minExpr << " to " << maxExpr << ".\n" ;
  std::cout << numMissing << " missing values in input file.\n"; 

  // print data
  /*if(VERBOSE) {
    std::cout << "\nInput Data:\n";
    for (int i = 0; i < numGenes; i++) {
      for (int j = 0; j < numInd; j++)
	      std::cout << data[i][j] << " ";
      std::cout << std::endl;
    }
  }
*/

  // normalize values to range between -1 and +1 using global feature scaling
  float range = maxExpr - minExpr; // range of input values
  
  for (int i = 0; i < numGenes; i++)
    for (int j = 0; j < numInd; j++)
      if (data[i][j] > MIN_EXP - 1 + TOL) // don't normalize missing data flags
        data[i][j] = ((data[i][j] * 2 - minExpr * 2) / range) - 1; // scale value

  // print normalized values
/*  if (VERBOSE) {
    std::cout << "\nNormalized values:\n";
    for (int i = 0; i < numGenes; i++) {
      for (int j = 0; j < numInd; j++)
	      std::cout << data[i][j] << " ";
      std::cout << std::endl;
    }
  }

*/

return data;


}

