#!/usr/bin/env python3
import os;
import re;
import sys;
import fileSplitter;
import createInchiKeys;
import CFMrunner;
import spectraNormalizer;
import tandemMS_Merger;
import spectraDataEditor;

def main():
	workplacePath = "/mnt/scratch/ozing003/CFM_workplace/";
	cfmPath = workplacePath + "cfm/bin/";
	cfmDataPath = workplacePath + "cfmData/";
	smilePath = cfmDataPath + "smileFile/";
	resultPath = cfmDataPath + "results/";

	#split files
	fileSplitter.main(smilePath + "Beta_lactams.csv");
	print("File splitting done!")
	fileName = fileSplitter.fileName;
	
	#create inchi keys and add then to a file containing a DB id and smiles string
	print("creating inchiKeys")
	#createInchiKeys.main(smilePath, fileName);
	print("inchiKeys created!")
	
	#run CFM_ID
	print("getting ready to run CFM-ID");
	#CFMrunner.main(cfmPath,cfmDataPath,fileName,resultPath);
	print("Done running CFM_ID!");
	
	#run spectraNormalizer
	print("normalizing spectra");
	spectraNormalizer.main(resultPath,fileName);
	print("done normalizing");
	
	#run spectraMerger
	print("merging the spectra")
	tandemMS_Merger.main(resultPath,fileName);
	print("merging done!")
	
	#run spectraDataEditor
	print("adding extra data to the spectra")
	spectraDataEditor.main(smilePath,resultPath,fileName);
	print("Complete the data is now ready for MS2LDA")

if __name__ == '__main__':
	main();
