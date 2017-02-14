#!/usr/bin/python
# Script finds the location of the top most node on the myocardium
# Argument 1: Inner Surface Nodes, Argument 2: Nodal Positions
# Written by Aditya Ponnaluri
# July 2016

# Usage: ./FindTopMostNode.py Small_B.innerSurfNode Small_B.node

import sys

InputFileStrings = ["" for x in range(len(sys.argv) - 1)];

for arg_iter in range(1,len(sys.argv)):
        InputFileStrings[arg_iter - 1] = sys.argv[arg_iter]


if len(sys.argv) - 1 == 0:
        print "No Input Arguments Provided"
        exit()

# First Read Node File:
innerSurfNodeFile = open(InputFileStrings[0], 'r');
nodeFile = open(InputFileStrings[1], 'r');

FirstLine = nodeFile.readline().split()
NumNodes = int(FirstLine[0])

nodePositions = [[0 for x in range(3)] for y in range(NumNodes)];

# Read all nodes into a table
for nodeIter in range(NumNodes):
	Line = nodeFile.readline().split()
	nodePositions[nodeIter][0] = float(Line[0]);
	nodePositions[nodeIter][1] = float(Line[1]);
	nodePositions[nodeIter][2] = float(Line[2]);

# Now read in inner surface node file and find min/max nodes
FirstLine = innerSurfNodeFile.readline().split()
NumNodes = int(FirstLine[0]);

minNode = 0;
maxNode = 0;
minZVal = float(1000.0);
maxZVal = -1000.0;

for nodeIter in range(NumNodes):
	Line = innerSurfNodeFile.readline().split();
	tempNode = int(Line[0]);
	if nodePositions[tempNode][2] < minZVal:
		minNode = tempNode;
		minZVal = nodePositions[tempNode][2];
	if nodePositions[tempNode][2] > maxZVal:
		maxNode = tempNode;
		maxZVal = nodePositions[tempNode][2];

print("Max Node: ", maxNode, " with Z-Position: ",maxZVal,"\n");
print("Min Node: ", minNode, " with Z-Position: ",minZVal,"\n");
