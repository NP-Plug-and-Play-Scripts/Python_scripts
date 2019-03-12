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

"""
creates a dictionary out of a list of spectra.
the input lists consists of list of list each list contains 1 spectra.
It calculates the number of peaks in each spectra and adds the the weight to the dictionary allong with the occurance.
At the end it gets the list of keys in the dict and sorts it.
In = a list of lists containing spectra
out =  a sorted KeyList and a dict containing the weights as key and occurance a value.
"""
def makePeakOccuranceDict(spectraList):
    peakOccuranceDict = {}
    #gets the number of identifiers in the spectra. if the line of a spectra starts with a capitol letter the count goes up
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

"""
creates a graph based on the keyList the peak OccuranceDict and saves it at the given path with a given name.
the Graph is a bar graph with the number of peaks on the x axis and the occurance of these amounts on the y Axis.
"""
def makeDistributionGraph(keyList,peakOccuranceDict, path,editedName):
    valueList = [peakOccuranceDict[key] for key in keyList];
    plt.plot(keyList,valueList)
    plt.title(editedName + '_graph');
    plt.ylabel('Occurence', fontsize=18)
    plt.xlabel('Number of peaks', fontsize=16)
    plt.savefig(path + editedName + '_graph.png');

"""
Gets the mean value of a dictionary of weights and the amount of times these weights appear.
Requires a dict of weights as key and occurance as value.
"""
def getMean(peakOccuranceDict):
    totalValues = 0;
    totalPeaks = 0;
    for key in peakOccuranceDict.keys():
        totalPeaks += int(key) * int(peakOccuranceDict[key]);
        totalValues += int(peakOccuranceDict[key]);
    print(totalPeaks/totalValues);

"""
main method calls all the others
Requires a path and a file name.
path to the directory the file is in
name of the fiel you want analysed.
"""
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
