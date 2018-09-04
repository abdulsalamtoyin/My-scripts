#!/usr/bin/env python
 
import sys
 
in_file =sys.argv[1]
 
fh = open(in_file, 'r')
head = fh.readline()
lines = fh.readlines()
 
sigUp = 0
sigDown = 0
twoUp = 0
twoDown = 0
fourUp = 0
fourDown = 0
totalTested = 0
 
for line in lines:
 ln = line.strip().split()
 if ln[6] == 'OK':
 
  totalTested += 1
 
 ln[9] = float(ln[9])
 
 
 if ln[-1] == 'yes' and ln[9] > 1.5:
  sigUp += 1
 elif ln[-1] == 'yes' and ln[9] < -1.5:
  sigDown += 1
 
 if ln[9] > 1:
  twoUp += 1
 if ln[9] > 2:
  fourUp +=1
 elif ln[9] < -1:
  twoDown += 1
 if ln[9] < -2:
  fourDown += 1
 
fh.close()
 
print
print ("tested:"), totalTested
print ("sig increase:"), sigUp
print ("sig decrease:"), sigUp
print ("2x increase:"), twoUp
print ("2x decrease:"), twoDown
print ("4x increase:"), fourUp
print ("4x decrease:"), fourDown
print ()
