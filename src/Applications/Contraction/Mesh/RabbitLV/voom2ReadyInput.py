#!/usr/bin/python
# This script converts node files from mm to cm 
#
# Written by Aditya Ponnaluri
# May 2016

import sys

InputFileStrings = ["" for x in range(len(sys.argv) - 1)];

for arg_iter in range(1,len(sys.argv)):
	InputFileStrings[arg_iter - 1] = sys.argv[arg_iter]


if len(sys.argv) - 1 == 0:
	print "No Input Arguments Provided"
	exit()

# First Read Node File:
nodeFile = open(InputFileStrings[0], 'r');
editedNodeFile = open(InputFileStrings[0].split('.')[0] + '_A.node', 'w')

FirstLine = nodeFile.readline().split()
NumNodes = FirstLine[0]

editedNodeFile.write(NumNodes + '\t3\n')

for nodeIter in range(int(NumNodes)):
	Line = nodeFile.readline().split()
	editedNodeFile.write(str(float(Line[0])/10) + '\t' + str(float(Line[1])/10) + '\t' + str(float(Line[2])/10) + '\n')

nodeFile.close()
editedNodeFile.close()
