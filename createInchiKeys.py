#!/usr/bin/env python3
import os;
import re;
import sys;


"""
creates 2 lists with the given file.
File should contain a ID and a smiles string.
"""
def createList(path):
	idList = [];
	smileList = [];
	for line in open(path):
		splittedLine = line.split();
		idList.append(splittedLine[0]);
		smileList.append(splittedLine[1]);
	return idList, smileList;

"""
creates a newFile from a list of smiles.
smileList = list of smiles strings
newPath = path for the new file to be placed
"""
def newFile(smileList,newPath):
	newFile = open(newPath,"w");
	for line in smileList:
		newFile.write(line + "\n");
	newFile.close();

"""
Runs the molconverter tool and takes the output file, reads it and puts the lines in a list and returns those.
The generated file gets deleted in the end.

"""
def createInchiKeys(molConvertPath,tempPath,newPath):
	newOut = outPath ;
	#run commandline program molconvert to turn a file of smiles strings in to a file of inchikeys.
	os.system("{0} -2:e inchikey {1} -o {2}".format(molConvertPath,newPath,newOut));
	inchiKeyList = [];
	for line in open(newOut,"r"):
		inchiKeyList.append(line);
	os.system("rm {}".format(newOut));
	return inchiKeyList;

"""

"""
def makeInchiSmileFile(idList,smileList,inchiKeyList, finalOutput):
	finalOut = open(finalOutput,"w");
	for i in range(len(idList)):
		newLine = idList[i] + " " + smileList[i] + " " + inchiKeyList[i];
		finalOut.write(newLine);
	finalOut.close();
	
def main(molConvertPath,filePath, fileName):
	newPath = filePath + "tempSmiles.txt";
	#for fileName (normaly something like "ethers_123_part_") in directory if file ends with 2 numbers followed by .txt save the file in a list.
	#match example   "smiles_1000_part_01.txt"
	foundFiles = [f for f in os.listdir(filePath) if re.search(re.escape(fileName) + r'[0-9]{2}.txt',f)];
	for aFile in foundFiles:
		splitName = aFile.split(".");
		path = filePath + aFile;
        
		#temporary files to store the smiles strings in and the inchi keyś
		newPath = filePath + "tempSmiles.txt";
        tempPath = filePath + "tempInchiKeys.txt";
        
        #new path for the output of the dataFile. Will contain the ID,SMILES,InchIKey.
		finalOutput = filePath + splitName[0] + "_dataFile." + splitName[1];
		idList, smileList = createList(path);
		newFile(smileList,newPath);
		inchiKeyList = createInchiKeys(molConvertPath,tempPath,newPath);
		makeInchiSmileFile(idList,smileList,inchiKeyList, finalOutput);
	os.system("rm {}".format(newPath));
	
if __name__ == '__main__':
	sys.exit(main())
