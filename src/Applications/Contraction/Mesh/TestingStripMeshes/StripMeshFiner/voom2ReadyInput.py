#!/usr/bin/python
# This python script takes two inputs (Node, Element) files and makes them usable for voom2
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
springBCFile = open(InputFileStrings[0].split('.')[0] + '_A.SpringBCnode', 'w');

FirstLine = nodeFile.readline().split()
NumNodes = FirstLine[0]

editedNodeFile.write(NumNodes + '\t3\n')

for nodeIter in range(int(NumNodes)):
	Line = nodeFile.readline().split(',')
	editedNodeFile.write(Line[1] + '\t' + Line[2] + '\t' + Line[3] + '\n')

	if abs(float(Line[1])) <= 1.0E-9:
		springBCFile.write(str((nodeIter + 1)) + '\n');

nodeFile.close()
editedNodeFile.close()
springBCFile.close();


# Next Read Element File
elementFile = open(InputFileStrings[1], 'r');
editedElementFile = open(InputFileStrings[1].split('.')[0] + '_A.ele', 'w')

FirstLine = elementFile.readline().split()
NumElements = FirstLine[0]

editedElementFile.write(NumElements + '\t' + FirstLine[1] + '\n');

for eleIter in range(int(NumElements)):
	Line = elementFile.readline().split(',')
	editedElementFile.write(str(int(Line[1]) - 1) + '\t' + str(int(Line[2]) - 1) + '\t' + str(int(Line[3]) - 1) + '\t' + str(int(Line[4]) - 1) + '\n')
