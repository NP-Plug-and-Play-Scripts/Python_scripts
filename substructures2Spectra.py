#!/usr/bin/env python3

import sys;
import os;
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors
 
def structureFileSplitter(structureFile):
	subStructureDict = {};
	for line in open(structureFile):
		comboList = [];
		splittedLine = line.split();
		smileID = splittedLine[1];
		substructureList = splittedLine[0].split(".");
		comboList.append(substructureList);
		weightList = []
		for smile in substructureList:
			mol = Chem.MolFromSmiles(smile);
			molWeight = Descriptors.ExactMolWt(mol);
			weightList.append(molWeight);
		comboList.append(weightList);
		subStructureDict[smileID] = comboList
	return subStructureDict;

"""
Creates a fake spectra with the list of lists containing substructures.
Each substructure is added to the spectra and the weight of a H+ atom is added to the substructure weight.
"""
def createSMGF(subStructureDict):
	spectraList = [];
	for key in subStructureDict.keys():
		spectra = [];
		spectra.append("BEGIN IONS");
		spectra.append("ID=" + key);
		spectra.append("TITLE=Fake Spectra of Molblocks generated Substructures with S. Stokman Rules");
		spectra.append("PEPMASS=" + str(subStructureDict[key][1][-1]));
		spectra.append("SMILES=" + subStructureDict[key][0][-1]);
		peakList = [];
		for weight in subStructureDict[key][1]:
			weight = float(weight);
			# adding the weight of 1 h+ atom 
			spectraWeight = str(round(weight + 1.008,8))
			# adding the weight of 2 h+ atoms
			spectraWeight2 = str(round(weight + 2.016,8))
			# removing the weight of 1 h+ atom 
			spectraWeight3 = str(round(weight - 1.008,8))
			# removing the weight of 2 h+ atoms
			spectraWeight4 = str(round(weight - 2.016,8))
			if weight != float(subStructureDict[key][1][-1]):
				peakList.append(spectraWeight4 + " 200");
				peakList.append(spectraWeight3 + " 300");
				peakList.append(spectraWeight + " 300");
				peakList.append(spectraWeight2 + " 200");
			else:
				peakList.append(spectraWeight + " 900");
		peakList.sort();
		for peak in peakList:
			spectra.append(peak);
		spectra.append("END IONS");
		spectraList.append(spectra);
	return spectraList;

"""
Creates a new file with a file path and a list of spectra
"""
def newFile(filePath,spectraList):
	newFile = open(filePath,"w");
	for spectra in spectraList:
		for line in spectra:
			newFile.write(line + "\n");
		newFile.write("\n");
	newFile.close();


def main(filePath):
	subStructureDict = structureFileSplitter(filePath);
	spectraList = createSMGF(subStructureDict);
	newFilePath = filePath.replace(".txt", "_spectra.mgf")
	newFile(newFilePath,spectraList);
	
if __name__ == '__main__':
	main(sys.argv[1]);
