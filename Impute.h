
using namespace std
#include <stdio.h>      
#include <stdlib.h>     
#include <time.h>
#include <math.h> 

void ImputationGenotype(float *matrix1, float *matrix2, float *matrix3){
	srand(time(NULL))
	// Most Likely Genotype in one iteration, if two genotypes with equal probability, then randomly choose one
	for (int i=0; i<marker; i++, matrix1++, matrix2++, matrix3++, ){
		int randnumber = rand()%2;
		if (randnumber == 0){
			if (*matrix1 >= *matrix2){
				if (*matrix1 >= *matrix3)
					GenotypeSampling[0][i]++;
			}
			else if(*matrix2 >= *matrix1){
				if (*matrix2 >= *matrix3)
					GenotypeSampling[1][i]++;
			}
			else if(*matrix3 >= *matrix1){
				if (*matrix3 >= *matrix2)
					GenotypeSampling[2][i]++;
			}
		}
		else if (randnumber == 1){
			if(*matrix2 >= *matrix1){
				if (*matrix2 >= *matrix3)
					GenotypeSampling[1][i]++;
			}
			else if(*matrix3 >= *matrix1){
				if (*matrix3 >= *matrix2)
					GenotypeSampling[2][i]++;
			}
			else if (*matrix1 >= *matrix2){
				if (*matrix1 >= *matrix3)
					GenotypeSampling[0][i]++;
			}
		}
		else {
			if(*matrix3 >= *matrix1){
				if (*matrix3 >= *matrix2)
					GenotypeSampling[2][i]++;
			}
			else if (*matrix1 >= *matrix2){
				if (*matrix1 >= *matrix3)
					GenotypeSampling[0][i]++;
			}
			else if(*matrix2 >= *matrix1){
				if (*matrix2 >= *matrix3)
					GenotypeSampling[1][i]++;
			}
		}
	}
}

void ImputationQuality(int iteration){
	int n = iteration;
	srand(time(NULL))
	for (int i=0; i<marker; i++){
		// Genotype score for certain marker
		GenotypeScore[i] = (2*posterior[2][i]+posterior[1][i])/n;
		// Most Likely Genotype after n iterations, if two genotypes with equal probability, then randomly choose one
		int randnumber = rand()%2;
		if (randnumber == 0){
			if (GenotypeSampling[0][i] >= GenotypeSampling[1][i]){
				if (GenotypeSampling[0][i] >= GenotypeSampling[2][i])
					MLGenotype[i] = 0;
			}
			else if(GenotypeSampling[1][i] >= GenotypeSampling[0][i]){
				if (GenotypeSampling[1][i] >= GenotypeSampling[2][i])
					MLGenotype[i] = 1;
			}
			else if(GenotypeSampling[2][i] >= GenotypeSampling[0][i]){
				if (GenotypeSampling[2][i] >= GenotypeSampling[1][i])
					MLGenotype[i] = 2;
			}
		}
		else if (randnumber == 1){
			if(GenotypeSampling[1][i] >= GenotypeSampling[0][i]){
				if (GenotypeSampling[1][i] >= GenotypeSampling[2][i])
					MLGenotype[i] = 1;
			}
			else if(GenotypeSampling[2][i] >= GenotypeSampling[0][i]){
				if (GenotypeSampling[2][i] >= GenotypeSampling[1][i])
					MLGenotype[i] = 2;
			}
			else if (GenotypeSampling[0][i] >= GenotypeSampling[1][i]){
				if (GenotypeSampling[0][i] >= GenotypeSampling[2][i])
					MLGenotype[i] = 0;
			}
		}
		else{
			if(GenotypeSampling[2][i] >= GenotypeSampling[0][i]){
				if (GenotypeSampling[2][i] >= GenotypeSampling[1][i])
					MLGenotype[i] = 2;
			}
			else if (GenotypeSampling[0][i] >= GenotypeSampling[1][i]){
				if (GenotypeSampling[0][i] >= GenotypeSampling[2][i])
					MLGenotype[i] = 0;
			}
			else if(GenotypeSampling[1][i] >= GenotypeSampling[0][i]){
				if (GenotypeSampling[1][i] >= GenotypeSampling[2][i])
					MLGenotype[i] = 1;
			}
		}
		// Genotype Quality score for certain marker
		GenotypeQualityScore[i] = GenotypeSampling[MLGenotype[i]][i]/n;
	}
	double sum, mean, variance;
	sum = 0.0;
	for (int i = 0; i < marker; ++i) {
		sum += GenotypeScore[i];
	}
	mean = sum/marker;
	for (int i = 0; i < marker; ++i) {
		variance += pow((GenotypeScore[i]-mean),2);
	}
	variance = variance/marker;
	R_square_expected = variance/(mean*(1-mean/2));
}

