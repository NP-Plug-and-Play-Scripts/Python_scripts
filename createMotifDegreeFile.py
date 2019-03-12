#!/usr/bin/env python

import sys; 
def getExtraInfo(motifsPath):
	info = []
	count = 0
	for line in open(motifsPath):
		splitted = line.split(',');
		match = splitted[1];
		if "(None)" not in match:
			if "Best Match" not in match:
				print(count)
				print(match)
				count += 1
				info.append(match)


def getUniqueMotifList(motifsPath):
	uniqueMotifList = []
	for line in open(motifsPath):
		print(line)
		splittedLine = line.split(" ");
		motif = splittedLine[0].replace("\"","");
		if motif not in uniqueMotifList:
			if motif != "Mass2Motif":
				uniqueMotifList.append(motif);
	print(uniqueMotifList)
	print(len(uniqueMotifList))
	return uniqueMotifList;

def getDegreeDict(degreePath):
	degreeDict = {}
	for line in open(degreePath):
		splittedLine = line.split(",");
		degreeDict[splittedLine[0]] = splittedLine[1];
	return degreeDict
def makeMotifDegreeList(uniqueMotifList,degreeDict, newPath):
	newFile = open(newPath,'w');
	for motif in uniqueMotifList:
		newLine = motif + "," + degreeDict[motif] + "\n"
		newFile.write(newLine);
	newFile.close();

def main(motifPath,degreePath,newPath):
	#motifList = getUniqueMotifList(motifPath);
	#degreeDict = getDegreeDict(degreePath);
	#makeMotifDegreeList(motifList,degreeDict,newPath)
	getExtraInfo(motifPath)
if __name__ == "__main__":
	main(sys.argv[1],sys.argv[2],sys.argv[3]);
