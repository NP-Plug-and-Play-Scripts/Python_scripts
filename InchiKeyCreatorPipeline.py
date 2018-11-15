#!/usr/bin/env python3
import sys
import os
import re
import multiprocessing;
import fileSplitter;
import inchiKeyCreator3;
import InchiKeyFileMerger;


"""
Creates a fileList based on the given pattern and location looks for files that match the pattern.
dataPath = the path in which the wanted files are located
filePattern = the pattern of the file in this case name of the original unsplitted file so  "smiles" if the file was called "smiles.csv"
"""
def getFileList(dataPath,filePattern):
	pattern = re.escape(filePattern) + r'_[0-9]{1,8}_part_[0-9]{2}.txt';
	fileList = [f for f in os.listdir(dataPath) if re.search(pattern,f)];
	#sort it so the files go from 00 to 09;
	fileList.sort();
	return fileList;

def main(molconvertPath,filePath, csvSmileFile):
	#path to the csv file
	csvPath = filePath + csvSmileFile;
	#file name without .csv 
	filePattern = csvSmileFile.split(".")[0];
	#splits the file in to multiple parts
	fileSplitter.main(csvPath);
	#runs getFileList with the path and filePattern
	fileList = getFileList(filePath,filePattern);
	fileList.sort();
	jobs = [];
	#for each file in file list create a new job (multiprocess), add it to a list and start it.
	for fileNumber in range(len(fileList)):
		molconvertJob = multiprocessing.Process(target=inchiKeyCreator3.main, args=(molconvertPath,filePath,fileList[fileNumber],str(fileNumber)));
		jobs.append(molconvertJob);
		molconvertJob.start();
	#magic, goes through all the jobs and joins them. This causes the script to wait till all jobs are done before continuing
	for job in jobs:
		job.join();
	print('*** All jobs finished ***');
	InchiKeyFileMerger.main(filePath,filePattern);

if __name__ == '__main__':
	main(sys.argv[1],sys.argv[2],sys.argv[3])
