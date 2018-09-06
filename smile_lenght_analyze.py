#!/usr/bin/env python

#from rdkit import Chem
#from rdkit.Chem import Descriptors
#from rdkit.Chem import rdMolDescriptors
import os;
import re;



def analyzeFile(filePath):
    '''
    returns a summary of the file in the form of the mean, median,
    largest, smallest , mean of the largest 25% and smallest 25% 
    of the smile string length..
    '''
    totalLines = 0;
    sortedList = [];
    sortList = [];
    mean = 0;
    median = 0;
    smallest = 0;
    #going through the file saving each line as a dict containing 2 enties for each line the length and the line.
    for line in open(filePath).readlines():
        line = line.strip();
        smileLength = len(line.split()[1]);
        lineDict = {"length":smileLength, "line":line};
        sortList.append(lineDict);
    #sorting the list from small to large
    sortedList = sorted(sortList, key=lambda k: k["length"]);
    
    #calculating the mean, median, etc
    mean, totalLines = calcMean(sortedList);
    median = calcMedian(sortedList, totalLines);
    smallMean,largeMean = calcQuarters(sortedList, totalLines);
    head, tail = os.path.split(filePath);
    output = "File name:" + tail + "\tMean:" + str(mean) + "\tMedian:" + str(median) + "\tMean of smallest 25%:" + str(smallMean) + "\tMean of the largest 25%:" + str(largeMean) + "\tSmallest smile:" + str(getSmallest(sortedList)) + "\tLargest smile:" + str(getLargest(sortedList)) + "\n";
    return output;
    
def calcMean(sortedList):
    total = 0;
    totalLines = 0;
    for smileDict in sortedList:
        totalLines += 1;
        total += int(smileDict["length"]);
        
    return total/totalLines, totalLines;
    
def calcMedian(sortedList, totalLines):
    sortedList;
    #if the modulo 2 of the total amount of lines is 0 
    if totalLines%2 == 0:
        #takes the middle 2 amounts and takes the average of it
        middle = totalLines /2
        sumOfMiddle = sortedList[middle]["length"] + sortedList[middle-1]["length"];
        return sumOfMiddle/2;
    else: 
        middle = totalLines/2 - 0.5
        middleLen = sortedList[middle]
        return middleLen
def calcQuarters(sortedList, totalLines):
    quarter = totalLines/4;
    smallQ = sortedList[0:quarter];
    largeQ = sortedList[-quarter:];
    smallTot=0
    largeTot=0
    for i in range(quarter):
        smallTot += smallQ[i]["length"];
        largeTot += largeQ[i]["length"];
    smallMean = smallTot/quarter;
    largeMean = largeTot/quarter;
    return smallMean, largeMean

def getLargest(sortedList):
    return sortedList[-1]["length"];

def getSmallest(sortedList):
    return sortedList[0]["length"];


def main():
    import glob, os
    path  = "/mnt/scratch/ozing003/CFM_workplace/cfmData/smileFile/";
    foundFiles = [f for f in os.listdir(path) if re.search(r'\w*[0-9]{2}.txt',f)]
    newfile = open(path + "summary.txt", "w");
    for file in foundFiles:
        newfile.write(analyzeFile(path + "/" + file));
    newfile.close();

if __name__ == "__main__":
    main();
