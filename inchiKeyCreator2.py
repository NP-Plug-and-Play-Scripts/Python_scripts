#!/usr/bin/env python3

import os
import sys
import re

"""
This script turns a file of NP-ID and SMILES strings combinations in to a file that contains NP-ID, SMILES and InchiKey. or NP-ID, neutral SMILES, SMILES, neutral InchiKey, InchiKey.
WARNING this requires molconvert a tool included in jchem. https://chemaxon.com/products/instant-jchem/download. 

CFM_pipeline step 3:
Second step in the pipeline. takes each of the splitted files and adds the inchi keys to it so they can later be added to the mgf files.
next script run is CFMrunner.py.

Made by: Rutger Ozinga 
Last edit: 22/11/2018
"""
import os;
import re;
import sys;


"""
creates 3 lists with the given file.
File should contain a ID, a neutral smiles string and an original smiles string.
"""
def createLists(path):
	idList = [];
	smileList = [];
	altSmileList = [];
	for line in open(path):
		line = line.strip();
		line =line.replace(" ", ",");
		splittedLine = line.split(",");
		idList.append(splittedLine[0]);
		if len(splittedLine) > 2:
			smileList.append(splittedLine[1]);
			altSmileList.append(splittedLine[2]);

		else:
			smileList.append(splittedLine[1]);
			altSmileList.append(splittedLine[1]);
	return idList, smileList, altSmileList;

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
creates inchiKeys with the path to molConvert, a tempPath for temporary files and a path to the smilestrings.
runs molconvert for each smile seperatly (takes ages).
"""
def createInchiKeys(molConvertPath,tempPath,smilePath):
	newOut = tempPath;
	print("Generating InchiKeys for " + smilePath)
		#run commandline program molconvert to turn a file of smiles strings in to a file of inchikeys.
	inchiKeyList = [];
	for smile in open(smilePath,"r"):
		os.system("{0} -2:e inchikey {1} -o {2}".format(molConvertPath + "molconvert" ,"'" + smile + "'",newOut));
		if os.path.exists(newOut) and os.path.getsize(newOut) > 0:
			for inchiKey in open(newOut,"r"):
				inchiKeyList.append(inchiKey);
			os.system("rm {}".format(newOut));
		else:
			inchiKeyList.append("None\n");
	print("InchiKeys for " + smilePath + " done!")
	return inchiKeyList;

"""
takes the id list (contains all the id's) smile list, alternative smile list, inchiKeyList, alternative inchikeyList and the path to 
the new file and writes the data in to the file
"""
def makeInchiSmileFile(idList,smileList,altSmileList,inchiKeyList, altInchiKeyList, finalOutput):
	finalOut = open(finalOutput,"w");
	print(len(idList), len(smileList),len(altSmileList), len(inchiKeyList), len(altInchiKeyList));
	for i in range(len(idList)):
		if smileList[i] != altSmileList[i]:
			#strip for the new inchi keys because molconvert realy likes to add line endings
			newLine = idList[i] + " " + smileList[i] + " " + altSmileList[i] + " " + inchiKeyList[i].strip() + " " + altInchiKeyList[i];
			finalOut.write(newLine);
		else:
			newLine = idList[i] + " " + smileList[i] + " " + inchiKeyList[i];
			finalOut.write(newLine);
	finalOut.close();



def main(molConvertPath,filePath, fileName, fileNumber):
	newAltPath = filePath + "tempAltSmiles.txt";
	print('working in ' + filePath + fileName);
	splitName = fileName.split(".");
	path = filePath + fileName;
	#temporary files to store the smiles strings in and the inchi keys
	newPath = filePath + "tempSmiles_" + fileNumber + ".txt";
	newAltPath = filePath + "tempAltSmiles_" + fileNumber + ".txt";
	tempPath = filePath + "tempInchiKeys_" + fileNumber + ".txt";
	tempAltPath = filePath + "tempAltInchiKeys_" + fileNumber + ".txt";
	#new path for the output of the dataFile. Will contain the ID,SMILES,InchIKey, Alternative Smile and alternative inchiKey.
	finalOutput = filePath + splitName[0] + "_dataFile_" + fileNumber + ".txt";
	idList, smileList, altSmileList = createLists(path);
	newFile(smileList,newPath);
	newFile(altSmileList, newAltPath);
	inchiKeyList = createInchiKeys(molConvertPath,tempPath,newPath);
	altInchiKeyList = createInchiKeys(molConvertPath,tempAltPath,newAltPath);
	makeInchiSmileFile(idList, smileList, altSmileList, inchiKeyList, altInchiKeyList, finalOutput);
	os.system("rm {}".format(newPath));
	
if __name__ == '__main__':
	main(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4]);
