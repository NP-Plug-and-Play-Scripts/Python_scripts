#!/usr/bin/env python

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors
import os;
import re;
"""
This script analyzes all the files in a directory that end with 2 digits and .txt
its ment for the smile string files that are produced with the cfm-multithreaded bash script
for each of the file parts it looks at the average mol weight of the smile strings, the median, 
the smallest the largest the mean of the 25% longest smiles and the mean of the 25% smallest smiles

currently to change the path were you want to look you can change it in the main() function

Made by: Rutger Ozinga 
Last edit: 9/6/2018
"""


def analyzeFile(filePath):
    '''
    requires a file path.
    returns a summary of the file in the form of the mean, median,
    largest, smallest , mean of the largest 25% and smallest 25% 
    of the smile mol weight..
    '''
    totalLines = 0;
    sortedList = [];
    sortList = [];
    mean = 0;
    median = 0;
    smallest = 0;
    #going through the file saving each line as a dict containing 2 enties for each line, the mol weight and the line.
    for line in open(filePath).readlines():
        line = line.strip();
        mol = Chem.MolFromSmiles(line.split()[1]);
        molWeight = Descriptors.ExactMolWt(mol);
        lineDict = {"weight":molWeight, "line":line};
        sortList.append(lineDict);
    #sorting the list from small to large
    sortedList = sorted(sortList, key=lambda k: k["weight"]);
    
    #calculating the mean, median, etc
    mean, totalLines = calcMean(sortedList);
    median = calcMedian(sortedList, totalLines);
    lightMean,heavyMean = calcQuarters(sortedList, totalLines);
    #split the path in to the head and tail. head contains the entire path exept the last file and tail contains the last file.
    # example: /home/user/desktop/file.txt  will become   head = /home/user/desktop/   tail = file.txt
    head, tail = os.path.split(filePath);
    output = "File name:" + tail + "\tMean:" + str(mean) + "\tMedian:" + str(median) + "\tMean of lightest 25%:" + str(lightMean) + "\tMean of the heaviest 25%:" + str(heavyMean) + "\tLightest smile:" + str(getLightest(sortedList)) + "\tLargest smile:" + str(getHeaviest(sortedList)) + "\n";
    return output;
    
def calcMean(sortedList):
    """
    calculates the mean of the given list of dictionaries 
    also counts the number of lines in the file.
    returns the mean and total number of lines
    """
    total = 0;
    totalLines = 0;
    for smileDict in sortedList:
        totalLines += 1;
        total += int(smileDict["weight"]);
    #return mean weight and total lines
    return total/totalLines, totalLines;
    
def calcMedian(sortedList, totalLines):
    """
    calculates the median given a list of dicts containing smiles + weight and the total amount of lines
    returns the median
    """
    sortedList;
    #if the modulo 2 of the total amount of lines is 0 
    if totalLines%2 == 0:
        #takes the middle 2 amounts and takes the average of it
        middle = totalLines /2
        sumOfMiddle = sortedList[middle]["weight"] + sortedList[middle-1]["weight"];
        return sumOfMiddle/2;
    else: 
        #since it cant be divided by 2 it means 1 is left so the result of the division by 2 results in ??.5
        # by adding 0.5 we would get the middle number but since we work in python to get the middle number of a list its the number -1
        # so i remove 0.5 to get that number instantly.
        #example  5/2 = 2.5  the middle number is 3 so 2.5 + 0.5 but in python a list of 5 numbers is 0,1,2,3,4 so the middle one is 2 
        print(totalLines);
        middle = int(totalLines/2 - 0.5)
        middleLen = sortedList[middle]
        return middleLen
        
def calcQuarters(sortedList, totalLines):
    """
    calculates the mean of the heaviest 25% and lightest 25%
    returns those 2 means.
    """
    quarter = totalLines/4;
    #lightest 25%
    lightQ = sortedList[0:quarter];
    #heaviest 25%
    heavyQ = sortedList[-quarter:];
    #total length of all the light smiles
    lightTot=0
    #same for the heaviest smiles
    heavyTot=0
    for i in range(quarter):
        lightTot += lightQ[i]["weight"];
        heavyTot += heavyQ[i]["weight"];
    lightMean = lightTot/quarter;
    heavyMean = heavyTot/quarter;
    #returns the lightMean and heavyMean
    return lightMean, heavyMean

def getHeaviest(sortedList):
    #gets the heaviest smile
    return sortedList[-1]["weight"];

def getLightest(sortedList):
    #gets the lightest smile
    return sortedList[0]["weight"];


def main():
    path  = "/mnt/scratch/ozing003/CFM_workplace/cfmData/smileFile/";
    #for each file in the path check if it ends with 2 digits and .txt and save them in a list
    foundFiles = [f for f in os.listdir(path) if re.search(r'\w*[0-9]{2}.txt',f)]
    #open new file
    newfile = open(path + "Moll_summary.txt", "w");
    # for each file in the list run analyzeFile and then write the output to the new file
    for file in foundFiles:
        newfile.write(analyzeFile(path + "/" + file));
    #close the file once done
    newfile.close();

if __name__ == "__main__":
    main();
