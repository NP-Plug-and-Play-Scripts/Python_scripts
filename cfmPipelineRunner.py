#!/usr/bin/env python3

"""
This script is the one you want to run if you have data from the NP-DB and want ran through CFM_predict so you get files ready for MS2LDA.
It calls all the other scripts and runs them after eachother.

IMPORTANT: its suggested to use the .sh instaler for CFM-ID found on this Github account. https://github.com/NP-Plug-and-Play-Scripts/Bash-scripts.git
This way the expected file structure is already in place. if not try to create a folder containing a folder with the following file structure.
(a tab means going in to a file)

-CFM_workplace (main folder for the workplace)
    -rdKit (library for cfm)
    -boost (library for cfm)
    -lpsolve (library for cfm)
    -cfm (is the cfm instalation and should contain a bin folder with runnable cfm tools).
        -bin (contains the runable cfm tools)
    -cfmData (is a folder that should contain 3 folders, and a cfm config file).
        -smileFile (location for the smiles data that you wish to run through this pipeline.)
        -results (location for the CFM_data results)
        -params_metab_ce_cfm  (contains the trained models for cfm_predict)
        -param_config.txt (contains the settings for cfm_predict)
    

MAKE IT EASY TO CHANGE PARAMETERS AND EXPLAIN THEM.
Made by: Rutger Ozinga 
Last edit: 10/10/2018
"""

import os;
import re;
import sys;
import fileSplitter;
import createInchiKeys;
import CFMrunner;
import spectraNormalizer;
import tandemMS_Merger;
import spectraDataEditor;
"""

"""
def main():
	workplacePath = "/mnt/scratch/ozing003/CFM_workplace/";
	cfmPath = workplacePath + "cfm/bin/";
	cfmDataPath = workplacePath + "cfmData/";
	smilePath = cfmDataPath + "smileFile/";
	resultPath = cfmDataPath + "results/";
	molConvertPath = "/mnt/scratch/ozing003/jchem/jchemsuite/bin/";
	#split files
    """STEP 1"""
	fileSplitter.main(smilePath + "AA_peptides_analogs.csv");
	print("File splitting done!")
	fileName = fileSplitter.fileName;
    """STEP 2"""
	#create inchi keys and add then to a file containing a DB id and smiles string
	print("creating inchiKeys")
	createInchiKeys.main(molConvertPath, smilePath, fileName);
	print("inchiKeys created!")
	"""STEP 3"""
	#run CFM_ID
	print("getting ready to run CFM-ID");
	CFMrunner.main(cfmPath,cfmDataPath,fileName,resultPath);
	print("Done running CFM_ID!");
	"""STEP 4"""
	#run spectraNormalizer
	print("normalizing spectra");
	spectraNormalizer.main(resultPath,fileName);
	print("done normalizing");
	"""STEP 5"""
	#run spectraMerger
	print("merging the spectra")
	tandemMS_Merger.main(resultPath,fileName);
	print("merging done!")
	"""STEP 6"""
	#run spectraDataEditor
	print("adding extra data to the spectra")
	spectraDataEditor.main(smilePath,resultPath,fileName);
	print("Complete the data is now ready for MS2LDA")

if __name__ == '__main__':
	main();
