#!/usr/bin/env python3
import os;
import sys;
import re;

"""
Turns the given spectra file in to a a list of lists. each list contains 1 spectra (from Begin IONS to End IONS)
first variable is th path to the spectra file (Mgf format).
"""
def makeSpectraList(spectraPath):
    spectraList = [];
    spectra = [];
    for line in open(spectraPath):
        #end of a spectra is indicated with END IONS so when it appears 
        #all data stored in the the list is added to the spectraList
        #and the list is emptied for the next spectra
        if line.startswith("END IONS"):
            if spectra != []:
                spectra.append(line);
                #sends the spectra (list of all the lines in the spectra) and the spectraZip to the method remakeSpectra.
                spectraList.append(spectra);
                spectra = [];
        else: 
            spectra.append(line.strip());
    return spectraList;

def addSpectraInfo(spectraList,spectraInfo):
    editedSpectraList = [];
    for spectra in spectraList:
        editedSpectra = [];
        cutUpTITLE = spectra[3].split(";");
        structure_id = cutUpTITLE[0].split("=")[1];
        ID = "ID=" + structure_id;
        newTitle = "TITLE=" + ";".join(cutUpTITLE[1:]);
        smile = "SMILES=" + spectraInfo[structure_id][0];
        inchiKey = spectraInfo[structure_id][1];
        for x in range(len(spectra)):
            if x < 3:
                if x == 1:
                    editedSpectra.append("IUPAC={}".format("Not_Added"));
                    editedSpectra.append(ID);
                    editedSpectra.append(newTitle);
                editedSpectra.append(spectra[x]);
            elif x == 3:
                editedSpectra.append(smile);
                editedSpectra.append(inchiKey);
            else:
                editedSpectra.append(spectra[x]);
        editedSpectraList.append(editedSpectra);
    return editedSpectraList;
            
"""
Turns the content of the given path in to a dictionary with the first entry of the line as key and the second and third as a list values
first vaiable contains the path to the file. file should contain lines with 3 tab separated values in this case ID, smiles string, inchikey
"""
def spectraInfoDict(infoPath):
    spectraInfo = {};
    for line in open(infoPath):
        splittedLine = line.split();
        spectraInfo[splittedLine[0]] = [splittedLine[1], splittedLine[2]];
    return spectraInfo;

"""
writes the info from the first variable in to the given file path.
first variable is a list of lists 
second variable is a path
"""
def writeFile(editedSpectra, newFile):
    editSpecFile = open(newFile,'w');
    for spectra in editedSpectra:
        for line in spectra:
            editSpecFile.write(str(line) + "\n");
    editSpecFile.close();

def main(smilePath,resultPath,fileName):
    infoPath = smilePath;
    spectraPath = resultPath;
    #for file in directory if file ends with 2 numbers followed by _output.mgf save the file in a list.
    #match example   "smiles_1000_part_01_output.mgf"
    foundInfoFiles = [f for f in os.listdir(infoPath) if re.search(re.escape(fileName) + r'[0-9]{2}_dataFile.txt',f)];
    foundSpectraFiles = [f for f in os.listdir(spectraPath) if re.search(re.escape(fileName) + r'[0-9]{2}_output_normalized_merged.mgf',f)];
    foundInfoFiles.sort();
    foundSpectraFiles.sort();
    for x in range(len(foundInfoFiles)):
        print(foundInfoFiles[x], foundSpectraFiles[x],x);
        spectraFilePath = spectraPath + foundSpectraFiles[x];
        infoFilePath = infoPath + foundInfoFiles[x];
        spectraList =makeSpectraList(spectraFilePath);
        spectraInfo = spectraInfoDict(infoFilePath);
        editedSpectra = addSpectraInfo(spectraList,spectraInfo);
        newFile = spectraFilePath.replace("normalized_merged","MS2LDA_ready");
        os.system("rm {}".format(infoPath + foundInfoFiles[x]));
        os.system("rm {}".format(spectraPath + foundSpectraFiles[x]));
        writeFile(editedSpectra,newFile);
    
if __name__ == '__main__':
    main();
