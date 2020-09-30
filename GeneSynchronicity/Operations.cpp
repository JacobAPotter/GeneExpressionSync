#ifndef MATH
#define MATH
#include <math.h>
#endif

#include "Operations.h"


vector< vector<float> * > * Derivative(vector< vector<float> *>* data)
{
   vector< vector<float> *>* deriv = new vector<vector<float> *>; 

    for(int x = 0; x < data->size(); x++)
      { 
         vector<float> *v = new vector<float>();

         for(int y =0; y < data->at(x)->size() - 1 ; y++)
		v->push_back(data->at(x)->at(y+1) - data->at(x)->at(y));

         //padding for last entry, cant take dd
         v->push_back(0.0f);

         deriv->push_back(v);
	}

  return deriv;
}

   float StandardDeviation(vector<float>* vec)
{
    float avg = Mean(vec);  
    float variation = 0;
   for(int i = 0 ; i < vec->size(); i++)
       variation += pow(abs(vec->at(i) - avg), 2);
   variation /= (vec->size() - 1);

   return sqrt(variation);
}

float Mean(vector<float> *vec)
{
  float sum =0;

    for(int i = 0; i < vec->size(); i++)
       sum += vec->at(i);

  return sum/vec->size();
}


/*
vector<vector<bool> *>* BinaryThreshold(vector<vector<float> *>* data, float threshold)
{
   vector<vector<bool> *>* thresh;

   for(int x =0  ;x < data->size(); x++)
     {
       vector<bool> * v = new vector<bool>();

      for(int y =0; y < data->size(); y++)
      {
        if(data->at(x)->at(y) < threshold)
           v->push_back(false);
        else 
           v->push_back(true);
      }

      thresh->push_back(v);
    }
  
return thresh;
}

vector<vector<float> *>* InvertSign(vector<vector<float>* >* data)
{
   vector< vector<float> *>* inv; 

    for(int x = 0; x < data->size(); x++)
      { 
         vector<float> *v = new vector<float>();

         for(int y =0; y < data->at(x)->size(); y++)
	       	 v->push_back(-data->at(x)->at(y));

         inv->push_back(v);
	}

return inv;
}

vector<vector<bool> *>* HitOrMiss(vector<vector<bool> *>* data)
{
  vector<vector<bool>> str1 {
  vector<bool>{0, 0, 0, 1, 0, 0, 0 },
  vector<bool>{0, 0, 1, 1, 1, 0, 0 },
  vector<bool>{0, 1, 1, 1, 1, 1, 0 },
  vector<bool>{0, 0, 1, 1, 1, 0, 0 },
  vector<bool>{0, 0, 0, 1, 0, 0, 0 } };

  vector<vector<bool>> str2 {
  vector<bool>{0, 0, 0, 1, 0, 0, 0 },
  vector<bool>{0, 0, 1, 1, 1, 0, 0 },
  vector<bool>{0, 1, 1, 1, 1, 1, 0 },
  vector<bool>{0, 0, 1, 1, 1, 0, 0 },
  vector<bool>{0, 0, 0, 1, 0, 0, 0 } };

  vector<vector<bool> *>* A;

  //Erode data structuring element 1
  A = Erode(data,str1);

  //Erode the complement of image A ( Ac) with structuring element 2

  vector<vector<bool> *>* compliment;
   
  for(int x = 0; x < data->size(); x++)
      { 
         vector<bool> *v = new vector<bool>();

         for(int y =0; y < data->at(x)->size(); y++)
	       	 v->push_back(!A->at(x)->at(y));

         compliment->push_back(v);
      }

  vector<vector<bool> *>* B = Erode(compliment,str2);

//AND results from step 1 and step 2
  for(int x =0 ;x < data->size(); x++)
     for(int y =0 ; y < data->at(x)->size(); y++)
       compliment->at(x)->at(y) = A->at(x)->at(y) && B->at(x)->at(y);


 for(int x = 0; x < data->size(); x++)
    {
       delete A->at(x);
       delete B->at(x);
    }

 delete A;
 delete B;

  return compliment;
}

vector<vector<bool> *>* Erode(vector<vector<bool>*>* data, vector<vector<bool>> str)
{
 
   vector<vector<bool> *>* erd;

   for(int r =0 ; r < data->size(); r++)
      {

         vector<bool> * v = new vector<bool>();

          for(int c = 0; c < data->at(r)->size(); c++)
          {
           
        bool setPixel = false;

        for(int x = 0; x < str.size(); x++)
        for(int y = 0; y < str[x].size(); y++)
         {
           if(r + x - (str.size()/2) >= 0 &&
              r + x - (str.size()/2) < data->size() &&
              c + y - (str.size()/2) >= 0 &&
              c + y - (str.size()/2) < data->at(r)->size())              
           if(data->at(x)->at(y) > 0 &&
           data->at(r+x-(data->size()/2))->at(c+y-(data->size()/2)) == 0)
              setPixel = false;
         }

           v->push_back(setPixel);
          }

          erd->push_back(v);
      }
}
*/

