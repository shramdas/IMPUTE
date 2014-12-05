#include <stdio.h>
#include <limits.h>
#include <math.h>
#include <vector>

// HMM Functions

// Function to initialize the first vector of the Genotype Likelihood matrix
void InitialFirstVector(double * Sstart){
    int Nr = Sstart->size();
    for (int i = 0; i < Nr; i++) {
        Sstart[i] = 1.0;
    }
    return;
}

// Function to multiple a given Genotype Likelihood vector by the transition matrix
void  Transpose(double * Sfrom, double * Sto, double &theta)
{
    // For each study haplotype: StudyHap
    // Calculate the current state probabilities transitioned from a previous position with probability vector Sfrom
    // State space is consisted of N_r reference haplotypes
    // Sfrom, Sto is a vector of size N_r, Sfrom[i] is the probability of being in the state of reference haplotype S_i, i= 1, ..., N_r, on previous marker position
    // theta is transition rate, i.e., recombination rate, assume it is the same for all markers
    int Nr = Sfrom->size(); // total number of reference haplotypes
    double p_sum; // sum of vector Sfrom
    
    if (theta == 0) {
        for (int i=0; i < Nr; i++) {
            Sto[i] = Sfrom[i];
        }
    }
    
    else{
        p_sum = 0.0;
        for (int j=0; j < Nr; j++) {
            p_sum += Sfrom[j];
        }
        p_sum *= theta / (double)N_r;
        
        double q = 1.0 - theta;
        
        if (p_sum < 1e-10) {
            p_sum *= 1e15;
            q *= 1e15;
        }
        
        for(int i = 0; i < Nr; i++){
            Sto[i] = q * Sfrom[i] + p_sum;
        }
    }
    
    return;
    
}

// Function to calculate the conditional likelihood given the observed hyplotypes
void Condition(double * GV, char ** haplotypes, int position, char observed, double &epsilon, double &freq){
    
    if (observed == 0) {
        return;
    }
    
    int Nr = GV->size();
    char *study_allel = haplotypes[i][position];
    
    double prand = epsilon * freq;
    double pmatch = (1.0 - epsilon) + prand;
    
    for (int i=0; i < Nr; i++) {
        if (haplotypes[i][position] == observed) {
            GV[i] *= pmatch;
        }
        else GV[i] *= prand;
    }
    
    return;
    
}

