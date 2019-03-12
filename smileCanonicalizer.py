#!/usr/bin/env python
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors
import sys;
import os;
import re;


def smileCanonMaker(smileFile):
	subStructureDict = {};
	comboList = [];
	for line in open(smileFile):
		splittedLine = line.split(",");
		#removes spaces and commas in the names and replaces them with _ and -.
		name = splittedLine[0].replace(" ","_").replace(",","-").replace(";","");
		smile = splittedLine[1];
		mol = Chem.MolFromSmiles(smile);
		cannonSmile = Chem.MolToSmiles(mol,isomericSmiles = False);
		comboList.append(name + "," + cannonSmile);
	return comboList;

def writeNewFile(newFile,comboList):
	newSmileFile = open(newFile,'w');
	for line in comboList:
		newSmileFile.write(line + "\n");
	newSmileFile.close();

def main(filePath, fileName):
	comboList = smileCanonMaker(filePath + fileName);
	writeNewFile(filePath + fileName.split(".")[0] + "_Cannonical.csv",comboList);

if __name__ == '__main__':
	main(sys.argv[1],sys.argv[2]);
