#include "VcfFileReader.h"
#include "StringBasics.h"
#include "StringHash.h"
#include "MemoryAllocators.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include "string.h"
#include <boost/algorithm/string.hpp>
#include <sstream>
#include <typeinfo>

using namespace std;


//function that counts the number of markers in a VCF. Requirement: The last line ends with a newline character.
int countMarkers(const char* filename){
	int num_markers = 0;
	ifstream infile;
	infile.open(filename);
	string s;
	string line;
	while(getline(infile, line)){
		if(line.substr(0,1) == "#"){
			continue;
		}
		else{
			++num_markers;
		}
	}
	infile.close();
	return(num_markers);
}

//counts the number of haplotypes in a vcf file (#genotypes*2)
int countHaplotypes(const char* filename){
	int num_haplotypes = 0;
	ifstream infile;
	infile.open(filename);
	vector<string> s;
	string line;
	while(getline(infile, line)){
		if(line.substr(0,6) == "#CHROM"){
			boost::split(s, line, boost::is_any_of(" "));
			return(2*(s.size() - 9));
			break;
		}
	}
}

char getSampleHaplotype(const int haploNum, const int numMarker){
	return samplematrix[haploNum][numMarker];
}

char getReferenceHaplotype(const int haploNum, const int numMarker){
	return referencematrix[numMarker][haploNum];
}

String genMarkerStr(VcfRecord& record){
	static String marker;
	cout << record.getChromStr() << "\n";
	marker = record.getChromStr();
	marker += ":";
	marker += record.get1BasedPosition();
	marker += ':';
	marker += record.getRefStr();
	marker += ":";
	marker += record.getAltStr();
	return(marker);
}

vector<char> getAlleles(const string genotype){
	vector<char> alleles;
	if(genotype == "0|0"){
		alleles.push_back('0');
		alleles.push_back('0');
	}
	else if(genotype == "0|1"){
		alleles.push_back('0');
		alleles.push_back('1');
	}
	else if(genotype == "1|0"){
		alleles.push_back('1');
		alleles.push_back('0');
	}
	else if(genotype == "1|1"){
		alleles.push_back('1');
		alleles.push_back('1');
	}
	else if(genotype == ".|."){
		alleles.push_back('.');
		alleles.push_back('.');
	}
	else if(genotype == "0|2"){
		alleles.push_back('0');
		alleles.push_back('.');
	}
	else if(genotype == "2|0"){
		alleles.push_back('.');
		alleles.push_back('0');
	}
	else if(genotype == "1|2"){
		alleles.push_back('1');
		alleles.push_back('.');
	}
	else if(genotype == "2|1"){
		alleles.push_back('.');
		alleles.push_back('1');
	}
	else if(genotype == "2|2"){
		alleles.push_back('2');
		alleles.push_back('2');
	}
	else{
		alleles.push_back('.');
		alleles.push_back('.');
	}
	return(alleles);
}

void loadStudy(char ** samplematrix, int numhaplotypes, int nummarkers, const char *samplevcf){
	ifstream infile;
	infile.open(samplevcf);
	int found_header = 0;
	string line;
	int markernum = 0;
    
	vector<string> sample_ids;
	while(getline(infile, line)){
		if(line.substr(0,6) == "#CHROM"){
			found_header = 1;
			boost::split(sample_ids, line, boost::is_space());
			sample_ids.erase(sample_ids.begin(), sample_ids.begin()+9);
		}
		else{		// marker line
			if(found_header){
				int haplotypenum = 0;
				vector<string> genotypes;
				boost::split(genotypes, line, boost::is_space());
				genotypes.erase(genotypes.begin(),genotypes.begin()+9);
				vector<string> genotypefield;
				
				for(int i =0 ;i < genotypes.size(); i++ ){
					boost::split(genotypefield,genotypes[i],boost::is_any_of(":"));
					vector<char> alleles = getAlleles(genotypefield[0]);
                    
					samplematrix[haplotypenum][markernum] = alleles[0];
					samplematrix[haplotypenum+1][markernum] = alleles[1];
					haplotypenum += 2;
                    
				}
				markernum++;
			}
		}
	}
}

void loadReference(char** referencematrix, const int numhaplotypes, const int num_reference_markers, const char *referencevcf){
	ifstream infile;
	infile.open(referencevcf);
	int found_header = 0;
	string line;
	int markernum = 0;
    
	vector<string>	sample_ids;
	while(getline(infile, line)){
		if(line.substr(0,6) == "#CHROM"){
			found_header = 1;
			boost::split(sample_ids, line, boost::is_space());
			sample_ids.erase(sample_ids.begin(), sample_ids.begin()+9);
		}
		else{
			if(found_header){
				cout << line <<"\n";
				int haplotypenum = 0;
				vector<string> genotypes;
				boost::split(genotypes, line,boost::is_space());
				genotypes.erase(genotypes.begin(), genotypes.begin()+9);
				vector<string> genotypefield;
				
				for(int i = 0; i < genotypes.size(); i++){
					boost::split(genotypefield, genotypes[i], boost::is_any_of(":"));
					vector<char> alleles = getAlleles(genotypefield[0]);
					
					referencematrix[markernum][haplotypenum] = alleles[0];
					referencematrix[markernum][haplotypenum+1] = alleles[1];
					haplotypenum += 2;
                    
				}
				markernum++;
			}
		}
	}
}
