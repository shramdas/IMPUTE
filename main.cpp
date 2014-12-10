#include "LoadVCF.h"
#include "CalcGL.h"
#include "HMM.h"
#include "Impute.h"


//marker/haplotype index starting from 0

//declared globally
char ** samplematrix;
char ** referencematrix;
double * E; // error vector
double ** matrix;
int states, markers;

int main(int argc, char ** argv){
  
  if(argc < 4){cout << "command format: ./impute sample.vcf reference.vcf output.vcf" << endl;}
  
  char *samplevcf = argv[1];
  char *referencevcf = argv[2];
  string outputvcf = argv[3];
    
    int number_markers = countMarkers(argv[1]);
    markers = number_markers;
    int num_sample_haplotypes = countHaplotypes(argv[1]);
    
    //load sample vcf
    //allocate memory
    samplematrix = AllocateCharMatrix(num_sample_haplotypes, number_markers);
    //load genotypes
    loadStudy(samplematrix, num_sample_haplotypes, number_markers, samplevcf);
    
    int number_reference_markers = countMarkers(argv[2]);
    int num_reference_haplotypes = countHaplotypes(referencevcf);
    states = num_reference_haplotypes;
    
    referencematrix = AllocateCharMatrix(number_reference_markers, num_reference_haplotypes);
    loadReference(referencematrix, num_reference_haplotypes, number_reference_markers, referencevcf);
      
    //allocate memory for HMM
    // Define double array to save GenotypeLikelyhood and Frequences
<<<<<<< HEAD
    matrix= AllocateDoubleMatrix(states, markers);
    double ** freqs= AllocateDoubleMatrix(5, markers);
    InitialFreqs(freqs, num_reference_haplotypes); // Initialize freqs
    E = new [markers]; // Initialize error vector
    for (int i=0; i < markers; i++) {
        E[i] = 0.001;
    }
=======
    double ** matrix= AllocateDoubleMatrix(num_reference_haplotypes,number_markers);
    double ** freqs= AllocateDoubleMatrix(5, number_markers);
    // Define Most Likely Genotype, imputation quality
    double ** GenotypeSampling =  AllocateDoubleMatrix(3, number_markers);
    double * GenotypeScore = AllocateMatrix(number_markers);
    int * MLGenotype = AllocateMatrix(number_markers);
    double * GenotypeQualityScore = AllocateMatrix(number_markers);
    double R_square_expected = 0.0;
>>>>>>> beb0b45a45792e696cb988bd4b4b9a749f10e227
    
    int numIterations=20; // number of iterations of HMM
    double ** probmatrix = AllocateCharMatrix(num_sample_haplotypes, number_markers)
    
    for (int iter=0; iter<numIterations; iter++) {
        
<<<<<<< HEAD
        for (int i=0; i < num_sample_haplotypes; i++) {
            
            //walk forward/backward to produce the prob matrix as "matrix"
            RunLeftHmm(samplematrix[i], referencematrix, freqs )
            RunRightHmmCombine(samplematrix[i],referencematrix, freqs)
            
            // Impute with prob matrix "matrix"
            
        }
=======
        
        // Impute
        ImputationGenotype(freqs[1], freqs[2], freqs[3]);
>>>>>>> beb0b45a45792e696cb988bd4b4b9a749f10e227
    }
    ImputationQuality(numIterations);

    //clean memory
    FreeDoubleMatrix(matrix);
    FreeDoubleMatrix(freqs);
    delete [] E;
    FreeDoubleMatrix(samplematrix);
    FreeCharMatrix(referencematrix);
    FreeDoubleMatrix(probmatrix);
    
  return 1;
}

