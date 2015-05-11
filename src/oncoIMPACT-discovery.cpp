//============================================================================
// Name        : oncoIMPACT-discovery.cpp
// Author      : Nok C Suphavilai
// Version     :
// Copyright   : 
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <stdio.h>
#include <stdlib.h>

#include "header/utilities.h"
#include "header/input.h"

#ifdef _WIN32
	#include <direct.h>
	#include <windows.h>
	#include <direct.h>
#else
	#include <sys/types.h>
	#include <sys/stat.h>
	#include <unistd.h>
#endif

using namespace std;

int main( int argc, char *argv[] ) {

	string configFilename;

	if ( argc != 2 ){ // argc should be 4 for correct execution
		cout<<"usage: "<< argv[0] <<" <config_file>\n";
		return 0;
	}else{
		configFilename = argv[1];
	}

	/*
	 * Read configuration from a file
	 */

	string outDir;

	string networkFilename;
	string expFilename;
	string snpFilename;
	string cnvFilename;
	string benchmarkGeneListFilename;

	string dbPath = "";

	string datatType = "";
	int numThreads;
	//0 = sensitive, 1 = stringent
	int mode = 0; //default mode = sensitive

	ifstream inConfigFile;
	inConfigFile.open(configFilename.c_str(), std::ifstream::in);

	int countArvInConfigFile = 0;
	if (inConfigFile.is_open()) {

		//read the output from the
		while (inConfigFile.good()) {
			string line;
			if (!getline(inConfigFile, line))
				break;

			//comment or empty lines
			if(line[0] == '#' || line.size() == 0){
				continue;
			}

			istringstream lineStrem(line);
			string configName;

			int i = 0;
			while (lineStrem) {	//for each column (sample)
				string token;

				if (!getline(lineStrem, token, '='))
					break;

				if(i == 0){ //read config name
					configName = token;
				}else if(i == 1){	//read config content
					if(configName.compare("outDir") == 0){
						outDir = token;
						countArvInConfigFile++;
						cout << configName << " = " << outDir << endl;
					}else if(configName.compare("network") == 0){
						countArvInConfigFile++;
						networkFilename = token;
						cout << configName << " = " << networkFilename << endl;
					}else if(configName.compare("numThreads") == 0){
						countArvInConfigFile++;
						numThreads = atoi(token.c_str());
						cout << configName << " = " << numThreads << endl;
					}else if(configName.compare("exp") == 0){
						countArvInConfigFile++;
						expFilename = token;
						cout << configName << " = " << expFilename << endl;
					}else if(configName.compare("snp") == 0){
						countArvInConfigFile++;
						snpFilename = token;
						cout << configName << " = " << snpFilename << endl;
					}else if(configName.compare("cnv") == 0){
						countArvInConfigFile++;
						cnvFilename = token;
						cout << configName << " = " << cnvFilename << endl;
					}else if(configName.compare("benchmarkGeneList") == 0){
						countArvInConfigFile++;
						benchmarkGeneListFilename = token;
						cout << configName << " = " << benchmarkGeneListFilename << endl;
					}else if(configName.compare("dbPath") == 0){
						countArvInConfigFile++;
						dbPath = token;
						cout << configName << " = " << dbPath << endl;
					}else if(configName.compare("dataType") == 0){
						if(token.compare("ARRAY") == 0){	//use sensitive mode
							mode = 0;
						}else if(token.compare("RNA_SEQ") == 0){
							mode = 1;
						}
						cout << configName << " = " << token << " use mode = " << mode << endl;
					}
				}

				i++;
			}

		}
		inConfigFile.close();
	} else {
		cerr << "Error opening configuration file\n";
	}

	if(countArvInConfigFile != 8){
		cout << "Incorrect configuration file\n";
	}else{
		cout << "Configuration file is read\n";
	}

	//create sub-directory to store output files
	#if defined(_WIN32)	//_WIN32 - Defined for applications for Win32 and Win64
		_mkdir((outDir + "/sensitive").c_str());
		_mkdir((outDir + "/stringent").c_str());
		_mkdir((outDir + "/sensitive/samples").c_str());
		_mkdir((outDir + "/stringent/samples").c_str());
	#else
		mkdir((outDir + "/sensitive").c_str(), 0777); // notice that 777 is different than 0777
		mkdir((outDir + "/stringent").c_str(), 0777); // notice that 777 is different than 0777
		mkdir((outDir + "/sensitive/samples").c_str(), 0777); // notice that 777 is different than 0777
		mkdir((outDir + "/stringent/samples").c_str(), 0777); // notice that 777 is different than 0777
	#endif

	/*
	 * Read all input file
	 */

	//Read Network File
	cout << "reading network file and create mapping <id, geneSymbol> ..."	<< endl;
	map<string, int> geneSymbolToId;
	vector<string> geneIdToSymbol;
	TIntAdjList network;
	readNetwork(networkFilename.c_str(), &network, &geneIdToSymbol, &geneSymbolToId, '\t');
	int totalGenes = network.size();
	int totalGenesUpDown = totalGenes * 2;


	//Read gene expression
	TDoubleMatrix originalGeneExpressionMatrix;
	cout << "reading gene expression matrix ..." << endl;
	//sample id mapping
	vector<string> sampleIdToName;
	vector<int> genesEx; // gene ids ; size = # DE genes (consider only genes that are in the network)
	vector<string> geneSymbols;
	GeneExpression geneExpression;	//all info of gene expression
	geneExpression.genes = &genesEx;
	geneExpression.matrix = &originalGeneExpressionMatrix;
	readGeneExpression(expFilename.c_str(), &geneExpression, '\t',
			&geneSymbolToId, &sampleIdToName);

	int totalSamples = originalGeneExpressionMatrix[0].size(); // = originalMutationMatrix.size()
	int numGenesEx = originalGeneExpressionMatrix.size();

	cout << "\ttotal genes in gene expression matrix is " << numGenesEx << endl;
	cout << "\ttotal samples in gene expression matrix is " << totalSamples
			<< endl;

	//Read SNP (row = gene, col = samples)
	cout << "reading point mutation matrix ..." << endl;
	vector<int> genesPointMut;	// gene ids ; size = # mutated genes (consider only genes that are in the network)
	TIntegerMatrix originalPointMutationsMatrix;
	PointMutations pointMutations;	//all info of point mutation
	pointMutations.genes = &genesPointMut;
	pointMutations.matrix = &originalPointMutationsMatrix;
	readPointMutations(snpFilename.c_str(), &pointMutations, '\t', &geneSymbolToId);

	totalSamples = originalPointMutationsMatrix[0].size(); // = originalMutationMatrix.size()
	int numGenesPointMut = originalPointMutationsMatrix.size();

	cout << "\ttotal genes in point mutation matrix is " << numGenesPointMut << endl;
	cout << "\ttotal samples in point mutation matrix is " << totalSamples
			<< endl;

	//Read CNV
	cout << "reading CNV matrix ..." << endl;
	vector<int> genesCNV;	// gene ids ; size = # mutated genes (consider only genes that are in the network)
	TIntegerMatrix originalCNVsMatrix;
	CopyNumberVariation CNVs;	//all info of CNV
	CNVs.genes = &genesCNV;
	CNVs.matrix = &originalCNVsMatrix;
	readCopyNumberVariation(cnvFilename.c_str(), &CNVs, '\t', &geneSymbolToId);

	totalSamples = originalCNVsMatrix[0].size(); // = originalMutationMatrix.size()
	int numGenesCNV = originalCNVsMatrix.size();

	cout << "\ttotal genes in CNV matrix is " << numGenesCNV << endl;
	cout << "\ttotal samples in CNV matrix is " << totalSamples
			<< endl;

	/*
	 * Combining point mutation and CNV matrix
	 */

	//create a list of mutated genes, which contain point mutation or CNV or both
	vector<int> genesMut; // gene ids ; size = # mutated genes
	vector<bool>* isMutated = new vector<bool>(totalGenes);
	//initialize isMutated, pointMutationCount, CNVCount
	for (int i = 0; i < totalGenes; ++i) {
		isMutated->at(i) = false;
	}
	//1. check point mutation
	for (int i = 0; i < numGenesPointMut; ++i) {
		isMutated->at(genesPointMut[i]) = true;
	}
	//2. check CNV
	for (int i = 0; i < numGenesCNV; ++i) {
		isMutated->at(genesCNV[i]) = true;
	}
	//create a list of all mutated genes
	for (int i = 0; i < totalGenes; ++i) {
		if(isMutated->at(i)){
			genesMut.push_back(i);
		}
	}
	delete isMutated;

	//combine point mutation matrix and CNV matrix (row = gene, column = sample)
	cout << "combining point mutation matrix and CNV matrix ..." << endl;

	int numGenesMut = genesMut.size();
	TIntegerMatrix originalMutationMatrix;
	//for each mutated gene i
	for (int i = 0; i < numGenesMut; ++i) {
		int currentGeneId = genesMut[i];

		//find an index of a current gene in point mutation matrix
		int pmRowIndex = findIndex(&genesPointMut, currentGeneId);

		//find an index of a current gene in CNV matrix
		int cnvRowIndex = findIndex(&genesCNV, currentGeneId);

		vector<int> hasAMutatedGene(totalSamples);
		//for each sample j
		for (int j = 0; j < totalSamples; ++j) {
			if(pmRowIndex >= 0 and originalPointMutationsMatrix[pmRowIndex][j] != 0){
				hasAMutatedGene[j] = 1;
			}
			if(cnvRowIndex >= 0 and originalCNVsMatrix[cnvRowIndex][j] != 0){
				hasAMutatedGene[j] = 1;
			}
		}
		originalMutationMatrix.push_back(hasAMutatedGene);
	}

	Mutations mutations;	//all info of mutations (point mutation and CNV)
	mutations.genes = &genesMut;
	mutations.matrix = &originalMutationMatrix;

	totalSamples = originalMutationMatrix[0].size();

	cout << "\ttotal genes in mutation matrix is " << numGenesMut << endl;
	cout << "\ttotal samples in mutation matrix is " << totalSamples
			<< endl;

	//TODO Read parameters from file
	int L = 0;
	int D = 0;
	double F = 0.0;

	ifstream inParameterFile;
	string parameterFilename = dbPath + "/parameters.dat";
	inParameterFile.open(parameterFilename.c_str(), ifstream::in);
	if (inParameterFile.is_open()) {
		if (inParameterFile.good()) {
			string line;
			if (!getline(inParameterFile, line)){
				cerr << "Error reading parameter file\n";
				return 0;
			}
			getline(inParameterFile, line);	//skip the header
			cout << line << endl;
			istringstream lineStream(line);

			int i = 0;
			while (lineStream) {	//for each column (parameter)
				string param;
				if (!getline(lineStream, param, '\t'))
					break;
				if(i == 0){	//L
					L = atoi(param.c_str());
				}else if(i == 1){	//D
					D = atoi(param.c_str());
				}else if(i == 2){	//F
					F = atof(param.c_str());
				}
				i++;
			}

		}
		inParameterFile.close();
	} else {
		cerr << "Error opening parameter file\n";
	}
	cout << "Parameters (L,D,F) are set to " << L << ", " << D << ", " << F << endl;

	//Read phenotype genes from file
	string phenotypeFileName = dbPath + "/exp_gene_freq.dat";
	vector<bool> isPhenotypeGenesUpDown(totalGenesUpDown);
	for (int i = 0; i < totalGenesUpDown; ++i) {
		isPhenotypeGenesUpDown[i] = false;
	}
	vector<int> phenotypeGeneIdsUpDown;	// phenotype gene ids
	vector<double> pValues(totalGenesUpDown);
	readPhenotypeGenesFromFile(phenotypeFileName.c_str(), &phenotypeGeneIdsUpDown, &geneSymbolToId);
	int numPhenotypeGene = phenotypeGeneIdsUpDown.size();
	cout << "There are " << numPhenotypeGene << " phenotype genes" << endl;
	for (int i = 0; i < numPhenotypeGene; ++i) {
		isPhenotypeGenesUpDown[phenotypeGeneIdsUpDown[i]] = true;
	}

	//Read driver list with the statistics
	vector<DriverGeneFromFile> driverGenesFromFile(totalGenes);
	string driverFilenameSensitive = dbPath + "/sensitive/driver_list.txt";
	readDriverGenesFromFile(driverFilenameSensitive.c_str(), &driverGenesFromFile, &geneSymbolToId);

	//TODO Read MODULE.dat


	//Read cancer benchmark gene list
	vector<int> cancerBenchmarkGenes;
	readBenchmarkGeneList(benchmarkGeneListFilename, &cancerBenchmarkGenes, &geneSymbolToId);
	vector<bool> isCancerBenchmarkGenes(totalGenes);
	int numBenchmarkGenes = cancerBenchmarkGenes.size();
	for (int i = 0; i < numBenchmarkGenes; ++i) {
		isCancerBenchmarkGenes[cancerBenchmarkGenes[i]] = true;
	}

	/*
	 * TODO Construct module for all input samples
	 */

	//TODO find explained genes for all input samples

	//TODO find mutated genes for all input samples

	/*
	 * Finding driver genes
	 */

	//TODO combined db modules with input sample's modules

	/*
	 * SENSITIVE
	 */

	//TODO create bipartite graph

	//TODO greedy set cover algorithm

	//TODO construct modules, merge, and trim

	//TODO calculate the IMPACT score for each input sample

	//TODO print out the result into outDir

	/*
	 * STRINGENT
	 */

	//TODO create bipartite graph

	//TODO greedy set cover algorithm

	//TODO construct modules, merge, and trim

	//TODO calculate the IMPACT score for each input sample

	//TODO print out the result into outDir

	return 0;
}
