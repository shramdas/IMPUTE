#include "CalcGL.h"

void VectorNormalize(float * vector)
{
    float * total;
    for (int i=0; i<states; i++)  { total += vector[i]}
    for (int i=0; i<states; i++)  { vector[i] /= total}
}

void MatrixReplace(char * vector, int i){
    for (j=0; j<states; j++)   {matrix[i][j]=(double)vector[j]}
}


void RunLeftHmm(char * observed, char ** haplotypes, double ** freqs)
   {
   // Initialize likelihoods at first position ( used in MarkovModel::WalkLeft())
  // for (int i = 0; i < states; i++)
  //    matrix[0][i] = 1.;
   InitializeFirstVector(matrix[0])
 
   // Scan along chromosome
   for (int i = 0; i < markers - 1; i++)
      {
      if (observed[i])
         Condition(matrix[i], haplotypes, i, observed[i], E[i], freqs[observed[i]][i]);
         Transpose(matrix[i], matrix[i+1], R[i]);
      }
 
   if (observed[markers - 1])
      Condition(matrix[markers - 1], haplotypes, markers - 1, observed[markers - 1], E[markers - 1], freqs[observed[markers - 1]][markers - 1]);
   }


void RunRightHmmCombine(char * major, char * observed, char ** haplotypes, double ** freqs)
   {
   float * swap;
   float * vector = new float [states];
   float * extra = new float [states]; 
   for (int i = 0; i < states; i++)
      vector[i] = 1.;
   for (int i = markers - 1; i > 0; i--)
      {
      for (int j = 0; j < states; j++)
         extra[j] = vector[j] * matrix[i][j];
         VectorNormalize(extra)
         MatrixReplace(extra,i)
      if (observed[i])
         Condition(vector, haplotypes, i, observed[i], E[i], freqs[observed[i]][i]);
         Transpose(vector, extra, R[i - 1]);
         swap = vector; vector = extra; extra = swap;
       }
   if (observed[0])
      Condition(vector, haplotypes, 0, observed[0], E[0], freqs[observed[0]][0]);
      VectorNormalize(vector)
      MatrixReplace(extra,0)
   delete [] vector;
   delete [] extra;
   }








