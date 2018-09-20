#!/usr/bin/env python3
"""
This script splits a given file in to smaller parts containing a equal amount of lines (or close to in case the total is uneven)
the cores is the amount of parts you want to split the file in to.

originaly made to split a file containing smile strings in to equal parts that would be ran through anothe program.
with the file being in parts, that way multiple can be running in the background on seperate cores.
"""
#the amount of cores you wish to use.
cores = 10

#name of the file that needs to be ran through cfm_id
smileFile = "/mnt/scratch/ozing003/CFM_workplace/cfmData/smileFile/Rutger_10000_SMILES.txt";
#path to the new File
newPath = "/mnt/scratch/ozing003/CFM_workplace/cfmData/smileFile/";
#length of the file. is calculated in lineDivider()
fileLength = 0;

"""
splits the given list in to a predivined number of parts and stores them in a dictionary
the way it divides the list is by taking the inner and outerr most enties 
of the length sorted List and adds them to the first entry.
It then moves up one spot at the start and down one spot a the end of the list and 
adds those to the next entry.
"""
def lineDivider(sortedList):
    global fileLength;
    fileLength = len(sortedList);
    fileNum = 1;
    partitionDict = {};
    for x in range(fileLength/2):
        # if the fileNum is already in the dict append the existing list belonging to the key
        if fileNum in partitionDict.keys():
            #add the Xth entry of a list to the dict
            partitionDict[fileNum].append(sortedList[x]);
            #add the Xth entry at the end of a list to the dict 
            #example: list = [1,2,3,4,5] index 0 would be 1 from the start and index -1 would be 5 so to get that its -x-1 (-0 -1 = -1) 
            partitionDict[fileNum].append(sortedList[-x-1]);
        # if the fileNum is not a key in the dict add at with the values at the given index as a list of values.
        else:
            partitionDict[fileNum] = [sortedList[x],sortedList[-x-1]];
            
        if fileNum == cores:
            fileNum = 1;
        else:
            fileNum += 1;
    #in case the list has a uneven number of lines add the remaining line to the next file.
    if fileLength%2 != 0:
        partitionDict[fileNum].append(sortedList[fileLength/2]);
    return partitionDict;
"""
takes  a file and adds all lines to a list. Then sorts the List based on length.
This method requires a fielPath.
"""
def fileLengthSorter(filePath):
    fileList = [];
    for line in open(filePath,"r"):
        fileList.append(line);
    fileList = sorted(fileList,key=len);
    return fileList;

"""
given a dictionary containing numbers as keys and a list of lines it puts the content of the library
in individual files. So the values belonging to Key 1 will be added a file and the values of Key 2 to another file.
"""
def createFiles(partDict):
    newFileName = "smiles_" + str(fileLength/cores) + "_part_";
    for part in partDict.keys():
        #formats the number in a way that instead of being 0,1,2,3 it will be 00,01,02,03
        partNum = "{0:02d}".format(part-1)
        filePart = open(newPath + newFileName + partNum + ".txt", "w"); 
        for line in partDict[part]:
            filePart.write(line);
"""Main Method calls all other methods"""
def main():
    sortedFileList = fileLengthSorter(smileFile);
    partDict = lineDivider(sortedFileList);
    createFiles(partDict);
    
    
if __name__ == "__main__":
    main();

