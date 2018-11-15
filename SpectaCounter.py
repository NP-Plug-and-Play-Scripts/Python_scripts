#!/usr/bin/env python3
import sys;

def countSpectra(filePath):
	count = 0;
	for line in open(filePath):
		if line.startswith("BEGIN IONS"):
			count+=1;
	return count;

def main(filePath):
	specNum = countSpectra(filePath);
	print(specNum);

if __name__ == "__main__":
	main(sys.argv[1])
