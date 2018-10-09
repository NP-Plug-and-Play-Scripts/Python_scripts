#!/usr/bin/env python3
import os;
import sys;
import re;

"""
Runs CFM ID with the given info
cfmPath = path to the cfm-id runable folder (../cfm/bin)
cfmData = path to the cfmData folder
inName = the name of the Input file
outPath = path to which the output is written
"""
def runCFM(cfmPath,cfmData,inName,outPath):
	cfmMode = cfmPath + "cfm-predict";
	cfmParam = cfmData + "params_metab_ce_cfm/param_output0.log";
	cfmConfig = cfmData + "param_config.txt";
	filePath = cfmData + "smileFile/" + inName;
	command = "{0} {1} 0.001 {2} {3} 0 {4}&".format(cfmMode,filePath,cfmParam,cfmConfig,outPath);
	os.system(command);
	#cfmPath/cfm-predict ${cfmData}${smileDir}${newFileName}${value}.txt 0.001 $cfmData/params_metab_ce_cfm/param_output0.log $cfmData/param_config.txt 0 $cfmData/results/${newFileName}${value}_output.mgf&

"""
get all the files in a given file path that corresponds to a given pattern.
cfmData = path to the cfmData file
filePattern = patern that matches to the first part of the file
"""
def getFileList(cfmData,filePattern):
	pattern = re.escape(filePattern) + r'[0-9]{2}.txt';
	fileList = [f for f in os.listdir(cfmData+"/smileFile/") if re.search(pattern,f)];
	#sort it so the files go from 00 to 09;
	fileList.sort();
	return fileList;

"""
Main method runs all the functions.
cfmPath = path to the cfm-d bin folder
cfmData = path to the cfmData folder
fileName = pattern that matches the first part of the file "ethers_131_part_" or something similair.
outputPath = path to the results folder.
"""
def main(cfmPath, cfmData,fileName,outputPath):
	fileList = getFileList(cfmData,fileName);
	for x in range(len(fileList)):
		newFileName = fileList[x].replace(".txt","_output.mgf");
		newFilePath = outputPath + newFileName;
		runCFM(cfmPath,cfmData,fileList[x],newFilePath)

if __name__ == '__main__':
    main(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4]);

