#!/usr/bin/env python

#from rdkit import Chem
#from rdkit.Chem import Descriptors
#from rdkit.Chem import rdMolDescriptors
import os;
import re;
"""
This script analyzes all the files in a directory that end with 2 digits and .txt
its ment for the smile string files that are produced with the cfm-multithreaded bash script
for each of the file parts it looks at the average lenght of the smile strings, the median, 
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
    #split the path in to the head and tail. head contains the entire path exept the last file and tail contains the last file.
    # example: /home/user/desktop/file.txt  will become   head = /home/user/desktop/   tail = file.txt
    head, tail = os.path.split(filePath);
    output = "File name:" + tail + "\tMean:" + str(mean) + "\tMedian:" + str(median) + "\tMean of smallest 25%:" + str(smallMean) + "\tMean of the largest 25%:" + str(largeMean) + "\tSmallest smile:" + str(getSmallest(sortedList)) + "\tLargest smile:" + str(getLargest(sortedList)) + "\n";
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
        total += int(smileDict["length"]);
    #return mean and total lines
    return total/totalLines, totalLines;
    
def calcMedian(sortedList, totalLines):
    """
    calculates the median given a list of dicts containing smiles and the total amount of lines
    returns the median
    """
    sortedList;
    #if the modulo 2 of the total amount of lines is 0 
    if totalLines%2 == 0:
        #takes the middle 2 amounts and takes the average of it
        middle = totalLines /2
        sumOfMiddle = sortedList[middle]["length"] + sortedList[middle-1]["length"];
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
    calculates the mean of the largest 25% and smallest 25%
    returns those 2 means.
    """
    quarter = totalLines/4;
    #smallest 25%
    smallQ = sortedList[0:quarter];
    #largest 25%
    largeQ = sortedList[-quarter:];
    #total length of all the smallest smiles
    smallTot=0
    #same for the largest smiles
    largeTot=0
    for i in range(quarter):
        smallTot += smallQ[i]["length"];
        largeTot += largeQ[i]["length"];
    print(quarter);
    smallMean = smallTot/quarter;
    largeMean = largeTot/quarter;
    #returns the smallMean and largeMean
    return smallMean, largeMean

def getLargest(sortedList):
    #gets the longest smile
    return sortedList[-1]["length"];

def getSmallest(sortedList):
    #gets the smallest smile
    return sortedList[0]["length"];


def main():
    path  = "/mnt/scratch/ozing003/CFM_workplace/cfmData/smileFile/";
    #for each file in the path check if it ends with 2 digits and .txt and save them in a list
    foundFiles = [f for f in os.listdir(path) if re.search(r'smallSmile.txt',f)]
    #open new file
    newfile = open(path + "summary_small.txt", "w");
    # for each file in the list run analyzeFile and then write the output to the new file
    for file in foundFiles:
        newfile.write(analyzeFile(path + "/" + file));
    #close the file once done
    newfile.close();

if __name__ == "__main__":
    main();
