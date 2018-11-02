#!/usr/bin/env python3

import os;
import re;
import sys;
import matplotlib.pyplot as plt;


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
        elif line == "\n":
            pass;
        else: 
            spectra.append(line.strip());
    return spectraList;


def makePeakOccuranceDict(spectraList):
    peakOccuranceDict = {}
    numOfIdentifiers = len([line for line in spectraList[0] if line[0].isalpha()]);
    for spectra in spectraList:
        totalPeaks = len(spectra)-numOfIdentifiers;
        if totalPeaks in peakOccuranceDict.keys():
            peakOccuranceDict[totalPeaks] += 1;
        else:
        # create a new array of intensities in the dict with the given weight as key
            peakOccuranceDict[totalPeaks] = 1;
    keyList = peakOccuranceDict.keys();
    keyList.sort();
    return keyList, peakOccuranceDict;

def makeDistributionGraph(keyList,peakOccuranceDict, path,editedName):
    valueList = [peakOccuranceDict[key] for key in keyList];
    plt.plot(keyList,valueList)
    plt.title(editedName + '_graph');
    plt.xlabel('Occurence', fontsize=18)
    plt.ylabel('Number of peaks', fontsize=16)
    plt.savefig(path + editedName + '_graph.png');

def getMean(peakOccuranceDict):
    totalValues = 0;
    totalPeaks = 0;
    for key in peakOccuranceDict.keys():
        totalPeaks += int(key) * int(peakOccuranceDict[key]);
        totalValues += int(peakOccuranceDict[key]);
        
    print(totalPeaks/totalValues);
def main(path, filename):
    #open new file
    editedName = filename.split(".")[0];
    newfile = open(path + editedName + "_PeakSummary.txt", "w");
    # for each file in the list run analyzeFile and then write the output to the new file
    spectraList = makeSpectraList(path + filename);
    keyList, peakOccuranceDict = makePeakOccuranceDict(spectraList);
    getMean(peakOccuranceDict);
    makeDistributionGraph(keyList,peakOccuranceDict,path,editedName);
    #close the file once done
    newfile.close();

if __name__ == '__main__':
    main(sys.argv[1],sys.argv[2]);
