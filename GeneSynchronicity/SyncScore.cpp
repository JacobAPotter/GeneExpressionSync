/*
SyncScore Project
Jacob Potter
3/29/2019

This program will read in a file of gene expression data and score the synchronicity of the data by choosing an archetype and calculating the deviation from it.
It does this by sorting an individuals genes from lowest expression to highest, applying this order to all other individuals and determining how well that individual represents the others. The order which best fits the data is chosen by looking at the discrete derivatives across rows. Ideally the rows would increase the same way the archetype does. Scoring is calculated by how consistently these derivatives are increasing. 

quick sort algorithm adapted from:
https://www.geeksforgeeks.org/cpp-program-for-quicksort/
*/

#include "readExpr.cpp"
#include "Operations.cpp"

void quickSortAndStore(vector<float> * , vector<int> * , int , int );
int partition (vector<float>* , vector<int> *, int , int );
void quickSortColumns(vector<int> * , vector<vector<float> *>* , int , int );
int partition (vector<int>* , vector<vector<float> *>*, int , int );


using namespace std;

int main(int argc, char** argv)
{
  //read in the data 
   if(argc != 4)
     cout << "Invalid argument list" << endl;

   int numGenes = atoi(argv[2]);
   int numInd = atoi(argv[3]);

//function from readExpr.cpp
   float **arr = ReadAndNormalize(argv[1], numGenes, numInd);




   //convert data to vector
   vector<vector<float> *> *data = new vector<vector<float>* >();

   for(int x = 0 ; x < numGenes; x++)
   {
	  vector<float> * v = new vector<float>();
	
	  for(int y = 0 ; y < numInd; y++)
	     v->push_back(arr[x][y]);

	  data->push_back(v);
   }


//benchmark the programs runtime
timer t;
t.start("SyncScore timer started");

/* For each individual, order the rows so that genes are ordered
 from least expressed to most expressed. (the derivative of this column will always be positive)   Store the order and apply it to every other individual. 
Then see if the derivatives of other individuals (columns) are also mostly positive.
The ordering which gives the most consistently increasing derivatives will be chosen.
The individual which produced this ‘best’ ordering can be thought of as the archetype */



vector<int> * order = new vector<int>();

for(int i = 0 ; i < numInd; i++)
   order->push_back(i);


//worst case for totalNegative is numInd * numGenes
int bestTotalNegative = numInd * numGenes;
int indexOfBest = -1;

for(int i = 0; i < data->size(); i++)
{
   //order the genes based on the individual data->at(i) 
   //so that data->at(i)’s genes go from least expressed to most 

//this is a modofied version of quicksort which will store the ordering in a vector
   quickSortAndStore(data->at(i), order, 0 , data->at(i)->size()-1);

   //apply this ordering to every other individual
   for(int x =0 ; x < data->size(); x++)
   {
      if(x != i)
     {
        for(int y =0 ; y < data->at(x)->size(); y++)
        {
 	      float temp = data->at(x)->at(y);
	   
 	      data->at(x)->at(y) = data->at(x)->at(order->at(y));
	      data->at(x)->at(order->at(y)) = temp;
        }
     }
   }


   //get discrete derivative of every column
   vector<vector<float> * > * derivative = Derivative(data);
   

   //count the number of times the genes in each column are not increasing in expression.
   //this will be used to choose which ordering is best. 
   //(least negative numbers in derivative is best)

   
   int totalNegativeCount =0;

  for(int x = 0 ; x < derivative->size(); x++)	
     for(int y = 0; y < derivative->at(x)->size(); y++)
       	if(derivative->at(x)->at(y) < 0)
	 totalNegativeCount++;

   if(totalNegativeCount < bestTotalNegative)
   {
     bestTotalNegative = totalNegativeCount;
     indexOfBest = i;
   }

//clean up garbage vector
   for(int d = 0 ; d < derivative->size(); d++)
      delete derivative->at(d);


   delete derivative;
}


//now that we have the best ordering, go ahead and sort one last time  
quickSortAndStore(data->at(indexOfBest), order, 0, data->at(indexOfBest)->size()-1);

//apply this ordering to every other individual
for(int x =0 ; x < data->size(); x++)
{
    if(x != indexOfBest)
    for(int y =0 ; y < data->at(x)->size(); y++)
      {
 	      float temp = data->at(x)->at(y);

 	      data->at(x)->at(y) = data->at(x)->at(order->at(y));
	      data->at(x)->at(order->at(y)) = temp;
     }
}


//sort the individuals based on how consistantly their derivatives are increasing.
//individuals whose expression is least like the archetype should be last.
//start by storing the number of negative values in the derivative for each column

vector<vector<float> *>* derivative = Derivative(data);

vector<int> * negativesInEachColumn = new vector<int>();

for(int x= 0 ; x < data->size(); x++)
{
 int negatives = 0;

  for(int y = 0; y < data->at(x)->size(); y++)
    if(derivative->at(x)->at(y) < 0)
       negatives++;
 
  negativesInEachColumn->push_back(negatives);
}

//now order the columns based on the number of negatives in their derivatives
 quickSortColumns(negativesInEachColumn, data,0, negativesInEachColumn->size()-1);

//Data is now sorted by the relationship of each columns derivatives with the derivative of the archetype.
//To score the data we will look at how similar the derivatives of each column are and how the 
//gene expression rows increase from top to bottom


//Sum the number of negative values in the derivatives and divide by the total number of data points
float negativeSum = 0;

for(int i=0 ; i< negativesInEachColumn ->size(); i++)
    negativeSum += negativesInEachColumn->at(i);

float score1 = 1.0f - (negativeSum/(numGenes * numInd));


//sum up all the rows and count how many times the sums increase from top to bottom
//the fewer times it increases the better the synchronicity 
vector<float> * rowSums = new vector<float>();

for(int x =0 ; x < data->size(); x++)
 {
     float rowSum =0;
  
    for(int y =0 ; y < data->at(x)->size(); y++)
       rowSum += data->at(x)->at(y);

    rowSums->push_back(rowSum);
}
 

//count the times row sums increase 
float rowNegatives = 0;

for(int i = 1; i < rowSums->size(); i++)
   if( rowSums->at(i) < rowSums->at(i-1))
      rowNegatives += 1;

float score2 = 1.0f - (rowNegatives/rowSums->size());

//the two scores are equally weighted
float syncScore = (score1 + score2) / 2;

t.stop("Timer Stopped.");
cout << endl << "Program runtime: " << t << " seconds." << endl;


std::cout << endl << "Column Synchronicity: " << score1 << endl << "Row Synchronicity: " << score2 << endl;
std::cout << endl<< "Final Synchronicity Score: " << syncScore << endl << endl;
}

//modified version of quick sort which will keep track of the ordering
void quickSortAndStore(vector<float>* vec, vector<int>*order, int low, int high)
{
    if (low < high)
    {
        int pi = partition(vec, order, low, high);
  
        quickSortAndStore(vec,order, low, pi - 1);
        quickSortAndStore(vec, order, pi + 1, high);
    }
}


int partition (vector<float>* vec, vector<int> *order, int low, int high)
{
    int pivot = vec->at(high);	// pivot
    int i = (low - 1);  // Index of smaller element
  
    for (int j = low; j <= high- 1; j++)
    {
       
        if (vec->at(j) <= pivot)
        {
            i++;	// increment index of smaller element
          
  	    float temp = vec->at(j);
	    int tempIndex = order->at(j);
	    vec->at(j) = vec->at(i);
  	    order->at(j) = order->at(i);
	    vec->at(i) = temp;
	    order->at(i) = tempIndex;
        }
    }

    float temp = vec->at(high);
    int tempIndex = order->at(high);
    vec->at(high) = vec->at(i + 1);
    order->at(high) = order->at(i + 1);
    vec->at(i + 1) = temp;
    order->at(i + 1) = tempIndex;

    return (i + 1);
}
 
//another modified quick sort which will sort the vector ‘data’ based on the ordering of ‘neg’
void quickSortColumns(vector<int>* neg, vector<vector<float> *>* data, int low, int high)
{
    if (low < high)
    {
        int pi = partition(neg, data, low, high);
  
        quickSortColumns(neg, data, low, pi - 1);
        quickSortColumns(neg, data, pi + 1, high);
    }

}

//overload partition to sort data based on another vector (neg)’s order
int partition (vector<int>* neg, vector<vector<float> *>* data, int low, int high)
{
    int pivot = neg->at(high);	// pivot
    int i = (low - 1);  // Index of smaller element
  
    for (int j = low; j <= high- 1; j++)
    {
        if (neg->at(j) <= pivot)
        {
            i++;	// increment index of smaller element
          
  	    float temp = neg->at(j);
	    vector<float> * tempVec = data->at(j);
	    neg->at(j) = neg->at(i);
  	    data->at(j) = data->at(i);
	    neg->at(i) = temp;
	    data->at(i) = tempVec;
        }
    }

    float temp = neg->at(high);
    vector<float> * tempVec = data->at(high);
    neg->at(high) = neg->at(i + 1);
    data->at(high) = data->at(i + 1);
    neg->at(i + 1) = temp;
    data->at(i + 1) = tempVec;

return (i + 1);
}
