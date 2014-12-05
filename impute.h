#include "VcfFileReader.h"
#include "StringBasics.h"
#include "StringHash.h"
#include "MemoryAllocators.h"

#include <stdio.h>
#include <limits.h>
#include <math.h>

#include <iostream>     // std::cout, std::endl
#include <iomanip>

using namespace std

// Load Genotypes
#include <vector>
#include <map>
#include "libVcfVcfFile.h"
using namespace libVcf;
//std::unordered_map<std::string, int> SampleIdx;
//std::unordered_map<std::string, int> MarkerIdx;
void Initialize(char** genotypes, const string & filename, map<string, int> &SampleIdx, map<string,int>& MarkerIdx) 
{
		VcfFile* pVcf = new VcfFile;
		pVcf->bSiteOnly = false;
		pVcf->bParseGenotypes = false;
		pVcf->bParseDosages = false;
		pVcf->bParseValues = true;
		pVcf->openForRead(filename.c_str());
		int nSamples = pVcf->getSampleCount();
		int person = 0;
		for (int i = 0; i < nSamples; i++) {
			SampleIdx[pVcf->vpVcfInds[i]->sIndID + "." + pVcf->vpVcfInds[i]->sIndID] = person;
			person++;
		}
		int markerindex = 0;
		VcfMarker* pMarker = new VcfMarker;
		string markerName;
		while (pVcf->iterateMarker()) {//for each marker
			pMarker = pVcf->getLastMarker();
			markerName.printf("%s:%d", pMarker->sChrom.c_str(), pMarker->nPos);
			MarkerIdx[markerName] = markerindex;
			int AFidx = pMarker->asInfoKeys.Find("AF");
			int PLidx = pMarker->asFormatKeys.Find("PL");
			int GLflag = 0;
			if (PLidx < 0) {
				PLidx = pMarker->asFormatKeys.Find("GL");
				if (PLidx >= 0) GLflag = 1;
			}
			int formatLength = pMarker->asFormatKeys.Length();
			int idx11 = 0, idx12 = 1, idx22 = 2;
			string phred;
			int genoindex = markerindex * 3;
			for (int i = 0; i < nSamples; i++)//for each individual
			{
				phred.ReplaceTokens(pMarker->asSampleValues[PLidx + i*formatLength], ",");
				int phred11 = GLflag ? static_cast<int>(-10. * phred[idx11].AsDouble()) : phred[idx11].AsInteger();
				int phred12 = GLflag ? static_cast<int>(-10. * phred[idx12].AsDouble()) : phred[idx12].AsInteger();
				int phred22 = GLflag ? static_cast<int>(-10. * phred[idx22].AsDouble()) : phred[idx22].AsInteger();
				if ((phred11 < 0) || (phred11 < 0) || (phred12 < 0)) {
					printf(stderr,"Negative PL or Positive GL observed");
				}
				if (phred11 > maxPhred) phred11 = maxPhred;
				if (phred12 > maxPhred) phred12 = maxPhred;
				if (phred22 > maxPhred) phred22 = maxPhred;
				genotypes[i][genoindex] = phred11;
				genotypes[i][genoindex + 1] = phred12;
				genotypes[i][genoindex + 2] = phred22;
			}
			++markerindex;
		}
		delete pVcf;
		delete pMarker;
}



int CountStudyHap(){}
int CountRefHap(){}
int CountMarkers(){}

// Run HMM 
void RunLeftHMM(int StudyHap, double ** Lmat){
    
    InitialFirstVector(StudyHap, Lmat[0]);
    if (GetStudyHap(StudyHap, 0) != NULL) {
        CondGL(StudyHap[0], Lmat[0], epsilon);
    }
    
    for (int j=1; j<Ms; j++) {
        if (GetStudyHap(StudyHap, j) != NULL) {
            Transpose(StudyHap, Lmat[j-1], Lmat[j], theta);
            CondGL(StudyHap, j, Lmat[j], epsilon, freq);
        }
    }

}

RunRightHMM(int StudyHap, double ** Rmat){}

ConbineHMM(double ** Lmat, double ** Rmat, double ** V){}


// HMM Functions

void InitialFirstVector(double * Sstart){
    int Nr = Sstart->size();
    for (int i = 0; i < Nr; i++) {
        Sstart[i] = 1.0;
    }
    return;
}

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
Condition(matrix[i], haplotypes, i, observed[i], E[i], freqs[observed[i]][i])
void Condition(double * GV, char observe, double &epsilon, double &freq){

    int Nr = GV->size();
    char *study_allel = GetStudyHap(StudyHap, position);

    double prand = epsilon * freq;
    double pmatch = (1.0 - epsilon) + prand;
    
    for (int i=0; i<Nr; i++) {
        char *ref_allel = GetRefHap(RefHap[i], position);
        if (strcmp(study_allel, ref_allel) == 0) {
            GV[i] *= pmatch;
        }
        else GV[i] *= prand;
    }
    return;
}



// Imputation

void Impute(char * StudyHap, double ** V){}


