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
void Initialize((char** genotypes, const string & filename, map<string, int> &SampleIdx, map<string,int>& MarkerIdx) 
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
void RunLeftHMM(char * StudyHap, double ** Lmat){
    
    InitialFirstVector(StudyHap, Lmat[0]);
    for (int j=1; j<Ms; j++) {
        Transpose(StudyHap[j], Lmat[j-1], Lmat[j], theta);
        CondGL(StudyHap[j], glvec, epsilon);
        for (int i=0; i<Nr; i++) {
            Lmat[j][i] = Lmat[j][i] * glvec[i];
        }
    }

}

RunRightHMM(char * StudyHap, double ** Rmat){}

ConbineHMM(double ** Lmat, double ** Rmat, double ** V){}


// HMM Functions

void InitialFirstVector(int StudyHap, double * Sstart){
    int Nr = Sstart->size();
    for (int i = 0; i < Nr; i++) {
        Sstart[i] = 1.0/Nr;
    }
    return;
}

void  Transpose(int StudyHap, double * Sfrom, double * Sto, double &theta)
{
// For each study haplotype: StudyHap
// Calculate the current state probabilities transitioned from a previous position with probability vector Sfrom
// State space is consisted of N_r reference haplotypes
// Sfrom, Sto is a vector of size N_r, Sfrom[i] is the probability of being in the state of reference haplotype S_i, i= 1, ..., N_r, on previous marker position
// theta is transition rate, i.e., recombination rate, assume it is the same for all markers
    int Nr = Sfrom->size(); // total number of reference haplotypes
    double p_sum; // sum of vector Sfrom

	for(int i = 0; i < Nr; i++){
        p_sum = 0.0;
        for (int j=0; j < Nr; j++) {
                p_sum += Sfrom[j];
        }
        
        Sto[i] = (1.0 - theta) * Sfrom[i] + theat / (double)N_r * p_sum;
	}
    
    return;
}

void CondGL(int StudyHap, double * GV, double &espilon){

    int Nr = GV->size();
    char *study_allel = GetStudyHap(StudyHap);

    for (int i=0; i<Nr; i++) {
        char *ref_allel = GetRefHap(RefHap[i]);
        if (strcmp(study_allel, ref_allel) == 0) {
            GV[i] = 1 - epsilon;
        }
        else GV[i] = epsilon;
    }
    return;
}



// Imputation

void Impute(char * StudyHap, double ** V){}


