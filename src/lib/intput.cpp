/*
 * intput.cpp
 *
 *  Created on: 27 Mar, 2015
 *      Author: Nok
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <set>
#include <vector>
#include <cmath>
#include <map>
#include <cstdlib>
#include <cstring>
#include "../header/input.h"

void readGeneExpression(const char* filename, GeneExpression* geneExpression,
		char delim, map<string, int>* geneSymbolToId, vector<string>* sampleIdToName){
	ifstream inFile;
	inFile.open(filename, std::ifstream::in);
	vector<int>* geneIds = geneExpression->genes;

	map<string, int>::iterator end = geneSymbolToId->end(); // for a gene not found in the network

	if (inFile.is_open()) {

		//read the first line (sample name)
		string samplesList;
		getline(inFile, samplesList);
		istringstream ss(samplesList);
		string s;
		getline(ss, s, delim);	//discard the first column (GENE header)
		while (ss) {	//for each column (sample)
			if (!getline(ss, s, delim))
				break;
			sampleIdToName->push_back(s);
		}

		while (inFile.good()) { //for each row (gene)
			string s;
			if (!getline(inFile, s))
				break;

			istringstream ss(s);
			vector<double> row;	//to save expression values of each gene

			int i = 0;	// for read gene symbol in the first column
			bool found = false;	// is a current gene found in the network

			while (ss) {	//for each column (sample)
				string s;

				if (!getline(ss, s, delim))
					break;
				if(i == 0){ //read gene symbols
					map<string, int>::iterator it = geneSymbolToId->find(s);
					if(it != end){	//found in the network
						found = true;
						geneIds->push_back(it->second);	//add the gene id to geneEx
					}else{
						found = false;
					}
				}else{		//read expression values
					row.push_back(atof(s.c_str()));
				}

				if(!found)	//not found in the network, so skip this row (gene)
					break;

				i++;	//go to the next column (sample)
			}

			if(found){
				geneExpression->matrix->push_back(row);
			}
		}
		inFile.close();
	} else {
		cerr << "Error opening file\n";
	}
}

void readMutations(const char* filename, Mutations* mutations,
		char delim, map<string, int>* geneSymbolToId){
	ifstream inFile;
	inFile.open(filename, std::ifstream::in);
	vector<int>* geneIds = mutations->genes;

	map<string, int>::iterator end = geneSymbolToId->end(); // gene not found in the network

	if (inFile.is_open()) {

		//read the first line (sample name)
		string samples;
		getline(inFile, samples);

		//read the output from the
		while (inFile.good()) {
			string s;
			if (!getline(inFile, s))
				break;

			istringstream ss(s);
			vector<int> row;	//to save mutation values of each gene

			int i = 0;	// for read gene in the first column;
			bool found = false;	// found in the network

			while (ss) {
				string s;
				if (!getline(ss, s, delim))
					break;
				if(i == 0){ //read gene symbols
					map<string, int>::iterator it = geneSymbolToId->find(s);
					if(it != end){	//found in the network
						found = true;
						geneIds->push_back(it->second);
					}else{
						found = false;
					}
				}else{		//read expression values
					row.push_back(atoi(s.c_str()));
				}
				if(!found)
					break;
				i++;
			}
			if(found){
				mutations->matrix->push_back(row);
			}
		}
		inFile.close();
	} else {
		cerr << "Error opening file\n";
	}
}

void readPointMutations(const char* filename, PointMutations* pointMutations,
		char delim, map<string, int>* geneSymbolToId){
	ifstream inFile;
	inFile.open(filename, std::ifstream::in);
	vector<int>* geneIds = pointMutations->genes;

	map<string, int>::iterator end = geneSymbolToId->end(); // gene not found in the network

	if (inFile.is_open()) {

		//read the first line (sample name)
		string samples;
		getline(inFile, samples);

		//read the output from the
		while (inFile.good()) {
			string s;
			if (!getline(inFile, s))
				break;

			istringstream ss(s);
			vector<int> row;

			int i = 0;	// for read gene in the first column;
			bool found = false;	// found in the network

			while (ss) {
				string s;
				if (!getline(ss, s, delim))
					break;
				if(i == 0){ //read gene symbols
					map<string, int>::iterator it = geneSymbolToId->find(s);
					if(it != end){	//found in the network
						found = true;
						geneIds->push_back(it->second);
					}else{
						found = false;
					}
				}else{		//read mutation values
					row.push_back(atoi(s.c_str()));
				}

				if(!found)
					break;

				i++;
			}

			if(found){
				pointMutations->matrix->push_back(row);
			}

		}
		inFile.close();
	} else {
		cerr << "Error opening file\n";
	}
}

void readCopyNumberVariation(const char* filename, CopyNumberVariation* copyNumberVariation,
		char delim, map<string, int>* geneSymbolToId){
	ifstream inFile;
	inFile.open(filename, std::ifstream::in);
	vector<int>* geneIds = copyNumberVariation->genes;

	map<string, int>::iterator end = geneSymbolToId->end(); // gene not found in the network

	if (inFile.is_open()) {

		//read the first line (sample name)
		string samples;
		getline(inFile, samples);

		//read the output from the
		while (inFile.good()) {
			string s;
			if (!getline(inFile, s))
				break;

			istringstream ss(s);
			vector<int> row;

			int i = 0;	// for read gene in the first column;
			bool found = false;	// found in the network

			while (ss) {
				string s;
				if (!getline(ss, s, delim))
					break;
				if(i == 0){ //read gene symbols
					map<string, int>::iterator it = geneSymbolToId->find(s);
					if(it != end){	//found in the network
						found = true;
						geneIds->push_back(it->second);
					}else{
						found = false;
					}
				}else{		//read mutation values
					row.push_back(atoi(s.c_str()));
				}

				if(!found)
					break;

				i++;
			}

			if(found){
				copyNumberVariation->matrix->push_back(row);
			}

		}
		inFile.close();
	} else {
		cerr << "Error opening file\n";
	}
}

int findIndex(vector<int>* geneIds, int currentGeneId){
	int index = -1;	//not found
	int size = geneIds->size();

	for (int i = 0; i < size; ++i) {
		if(geneIds->at(i) == currentGeneId){
			return i;
		}
	}
	return index;
}

void readGenesList(const char* filename, vector<int>* geneIds, map<string, int>* geneSymbolToId){
	ifstream inFile;
	inFile.open(filename, std::ifstream::in);

	if (inFile.is_open()) {

		//read the output from the
		while (inFile.good()) {
			string s;
			if (!getline(inFile, s))
				break;

			istringstream ss(s);
			map<string, int>::iterator it = geneSymbolToId->find(s);
			geneIds->push_back(it->second);

		}
		inFile.close();
	} else {
		cerr << "Error opening file\n";
	}
}


void readGenesListUpDown(const char* filename, vector<int>* geneIdsUpDown, map<string, int>* geneSymbolToId){
	ifstream inFile;
	inFile.open(filename, std::ifstream::in);

	int totalGenes = geneSymbolToId->size();

	if (inFile.is_open()) {

		//read the output from the
		while (inFile.good()) {
			string s;
			if (!getline(inFile, s))
				break;

			istringstream ss(s);
		  	size_t foundUp = s.find("UP");
		  	size_t foundDown = s.find("DOWN");

		  	trimStr(s, "_");
		  	string currentGeneName = s;

			map<string, int>::iterator it = geneSymbolToId->find(currentGeneName);

		  	if (foundUp!=std::string::npos){
		  		//the gene is upregulated
		  		geneIdsUpDown->push_back(it->second);
		  	}

		  	if (foundDown!=std::string::npos){
		  		//the gene is downregulated
		  		geneIdsUpDown->push_back(it->second + totalGenes);
		  	}

		}
		inFile.close();
	} else {
		cerr << "Error opening file\n";
	}
}

void readPhenotypeGenesFromFile(const char* filename, vector<int>* phenotypeGeneIdsUpDown, map<string, int>* geneSymbolToId){
	ifstream inFile;
	inFile.open(filename, ifstream::in);

	int totalGenes = geneSymbolToId->size();

	if (inFile.is_open()) {

		//read the output from the
		while (inFile.good()) {
			string line;
			if (!getline(inFile, line))
				break;
			istringstream lineStream(line);
			int i = 0;
			string geneSym;
			bool isPhenotype = false;

			while (lineStream) {	//for each column (parameter)
				string token;
				if (!getline(lineStream, token, '\t'))
					break;
				if(i == 0){	//gene symbol
					geneSym = token;
				}else if(i == 7){	//phenotype status
					if(token.compare("Y") == 0){
						isPhenotype = true;
					}
				}
				if(isPhenotype){
					size_t foundUp = geneSym.find("UP");
					size_t foundDown = geneSym.find("DOWN");
				  	trimStr(geneSym, "_");
					map<string, int>::iterator it = geneSymbolToId->find(geneSym);
					if (foundUp!=std::string::npos){
						//the gene is upregulated
						phenotypeGeneIdsUpDown->push_back(it->second);
					}
					if (foundDown!=std::string::npos){
						//the gene is downregulated
						phenotypeGeneIdsUpDown->push_back(it->second + totalGenes);
					}
				}
				i++;
			}
		}
		inFile.close();
	} else {
		cerr << "Error opening file of explained and phenotype gene list \n";
	}
}

void readDriverGenesFromFile(const char* filename, vector<DriverGeneFromFile>* driverGenesFromFile, map<string, int>* geneSymbolToId){
	ifstream inFile;
	inFile.open(filename, ifstream::in);

	int totalGenes = geneSymbolToId->size();
	int countDriver = 0;

	if (inFile.is_open()) {

		string header;
		if(inFile.good()){
			getline(inFile, header); //skip the header
		}

		//read the output from the
		while (inFile.good()) {
			string line;
			if (!getline(inFile, line))
				break;
			istringstream lineStream(line);
			int i = 0;
			string geneSym;
			bool isPhenotype = false;

			int currentDriverGeneId = -1;
			while (lineStream) {	//for each column (parameter)
				string token;
				if (!getline(lineStream, token, '\t'))
					break;

//				GENE	DRIVER_FREQUENCY	DRIVER_SNV_FREQUENCY	DRIVER_DELTION_FREQUENCY	DRIVER_AMPLIFICATION_FREQUENCY	CANCER_CENSUS	PAN_CANCER	IMPACT	MUTATION_FREQUENCY	SNV_FREQUENCY	DELTION_FREQUENCY	AMPLIFICATION_FREQUENCY
//				EGFR	0.393	0.262	0.000	0.241	Y	NA	120.1447552328	0.393	0.262	0.000	0.241
				if(i == 0){	//gene symbol
					map<string, int>::iterator it = geneSymbolToId->find(token);
					currentDriverGeneId = it->second;
					countDriver++;
				}else if(i == 1){
					driverGenesFromFile->at(currentDriverGeneId).driverFreq = atof(token.c_str());
				}else if(i == 2){
					driverGenesFromFile->at(currentDriverGeneId).driverSnpFreq = atof(token.c_str());
				}else if(i == 3){
					driverGenesFromFile->at(currentDriverGeneId).driverDelFreq = atof(token.c_str());
				}else if(i == 4){
					driverGenesFromFile->at(currentDriverGeneId).driverAmpFreq = atof(token.c_str());
				}else if(i == 5){
					if(token.compare("Y") == 0){
						driverGenesFromFile->at(currentDriverGeneId).isInCancerCensus = true;
					}else{
						driverGenesFromFile->at(currentDriverGeneId).isInCancerCensus = false;
					}
				}else if(i == 7){
					driverGenesFromFile->at(currentDriverGeneId).impactScore = atof(token.c_str());
				}else if(i == 8){
					driverGenesFromFile->at(currentDriverGeneId).mutFreq = atof(token.c_str());
				}else if(i == 9){
					driverGenesFromFile->at(currentDriverGeneId).snpFreq = atof(token.c_str());
				}else if(i == 10){
					driverGenesFromFile->at(currentDriverGeneId).delFreq = atof(token.c_str());
				}else if(i == 11){
					driverGenesFromFile->at(currentDriverGeneId).ampFreq = atof(token.c_str());
				}
				i++;
			}
		}
		inFile.close();
	} else {
		cerr << "Error opening file of explained and phenotype gene list \n";
	}

	cout << "Read " << countDriver << " from driver_list.txt file\n";
}



void readBenchmarkGeneList(string benchmarkGeneListFilename, vector<int>* cancerBenchmarkGenes, map<string, int>* geneSymbolToId){
	ifstream inFile;
	inFile.open(benchmarkGeneListFilename.c_str(), std::ifstream::in);

	map<string, int>::iterator end = geneSymbolToId->end();

	if (inFile.is_open()) {

		//for each row
		while (inFile.good()) {
			string s;
			if (!getline(inFile, s))
				break;

			istringstream rowStr(s);

			int i = 0;
			bool found = false;
			while (rowStr) {	//for each column
				string s;

				if (!getline(rowStr, s, '\t'))
					break;
				if(i == 0){ //read gene symbols at the first column

					//trim the whitespace at the end of the string
					s.erase(s.find_last_not_of(" \n\r\t")+1);

					map<string, int>::iterator it = geneSymbolToId->find(s);
					if(it != end){	//found in the network
						found = true;
						cancerBenchmarkGenes->push_back(it->second);	//add the gene id to geneEx
//						cout << s << endl;
					}else{
						found = false;
//						cout << s << " not found" << endl;
					}
				}
				else if(i == 7){	//read tumor type (somatic mutation)
				}

				if(!found)	//not found in the network, so skip this row (gene)
					break;

				i++;	//go to the next column (sample)
			}

		}
		inFile.close();
	} else {
		cerr << "Error opening file\n";
	}

}

