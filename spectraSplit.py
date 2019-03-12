#!/usr/bin/env python
import sys;
import os;
import re;

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

def writeNewFile(newFile,combiSpecCol):
    newSpecFile = open(newFile,'w');
    for spectra in combiSpecCol:
        for line in spectra:
            newSpecFile.write(line + "\n");
    newSpecFile.close();


def main(filePath, fileName):
    fullFilePath = filePath + fileName;
    spectraList = makeSpectraList(fullFilePath);
    spectraNum = len(spectraList);
    if spectraNum%2:
        middle = int((spectraNum/2) + 0.5);
        print(middle)
    else:    
        middle = int(spectraNum/2);
        print(middle)
        
    filePart1 = filePath + fileName.replace(".mgf","_part1.mgf");
    filePart2 = filePath + fileName.replace(".mgf","_part2.mgf");
    writeNewFile(filePart1,spectraList[:middle])
    writeNewFile(filePart2,spectraList[middle:-1])


if __name__ == '__main__':
    main(sys.argv[1],sys.argv[2]);
