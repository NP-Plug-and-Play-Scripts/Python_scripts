#!/usr/bin/env python3
import math;
import re;
import os;

def makeSpectraList(filePath):
    spectraList = [];
    spectra = [];
    for line in open(filePath):
        if line.startswith("END IONS"):
            if spectra != []:
                spectra.append(line.strip());
                weightList, normList = normalizeSpectra(spectra);
                adv_normList = advancedNormalization(normList);
                spectraZip = zip(weightList,adv_normList);
                spectraList.append(remakeSpectra(spectra,spectraZip));
                spectra = [];
        else: 
            spectra.append(line.strip());
    return spectraList;
    
def normalizeSpectra(spectraList):
    total = 0;
    intensityList = [];
    weightList = []
    for line in spectraList[4:-1]:
        intensity = float(line.split()[1]);
        weightList.append(line.split()[0]);
        intensityList.append(intensity);
        total += intensity;
        n = len(spectraList[4:-1]);
    mean = total/n;
    meanSubbedIntensities = [x - mean for x in intensityList];
    std = math.sqrt(sum([x**2 for x in meanSubbedIntensities])/n);
    if(n == 1):
        normIntensities = [mean];
    else:
        normIntensities = [x/std for x in meanSubbedIntensities];
    return weightList, normIntensities;
"""
takes the normalized values and further normalizes it by putting it on a 0-900 scale. first it sets the values of 
everything above 2 and below -2 to 2 or -2 respectivly.
then it multiplies the values by 225. making the values go from -450 to 450 then by adding 450 the min value is now 0 and the highest values 900.
"""
def advancedNormalization(normIntensities):
    for x in range(len(normIntensities)):
        if normIntensities[x] > 2:
            normIntensities[x] = 2.0;
        elif normIntensities[x] < -2:
            normIntensities[x] = -2.0;
    magicInducedList = [x*225 for x in normIntensities];
    moreMagic = [x + 450 for x in magicInducedList];
    return moreMagic;
"""
takes the original spectra and the normalized one and creates a new spectra taking 
the headers from the original and the peaks from the normalized list of peaks.
"""
def remakeSpectra(originalSpectra, normalizedSpectra):
    newSpectra = [];
    for i in range(len(originalSpectra)):
        if(i < 3):
            newSpectra.append(originalSpectra[i]);
        elif(3 < i < len(originalSpectra)-1): 
            newSpectra.append(str(normalizedSpectra[i-4][0]) + " " + str(normalizedSpectra[i-4][1]));
        else:
            newSpectra.append(originalSpectra[i]);
    return newSpectra;
    
"""
takes a list of lists containing the normalized spectra and puts them in to a new files 
requires a list of lists and a file path.
"""
def writeNewFile(newFile,normSpecList):
    normSpecFile = open(newFile,'w');
    for spectra in normSpecList:
        for line in spectra:
            normSpecFile.write(str(line) + "\n");
    normSpecFile.close();


def main():
    filePath = "/mnt/scratch/ozing003/CFM_workplace/cfmData/results/";
    foundFiles = [f for f in os.listdir(filePath) if re.search(r'\w*[0-9]{2}_output.mgf',f)];
    
    for aFile in foundFiles:
        splitName = aFile.split("_");
        newName = "_".join(splitName[0:4]) + "_normalized_" + splitName[4];
        
        spectraList = makeSpectraList(filePath + "" + aFile);
        newFilePath = filePath + "" + newName;
        writeNewFile(newFilePath,spectraList);

main();
