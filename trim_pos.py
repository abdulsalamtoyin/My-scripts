#!/usr/bin/python

############
#
#	test for input
#
#	Andreas Gisel December 2014
#
#	command: python script.py <filename> <score treshold>  (python script.py DRTA001A.fastq.stats 37)
#
############

import sys, re

cut5 = ''
cut3 = ''

#####Load arguments

inputfile = sys.argv[1]	
treshscore = sys.argv[2]

G_count = []
Q_mean = []

works = open(inputfile).readlines()

# column	count		min		max	sum			mean	Q1	med	Q3	IQR	lW	rW	A_Count	C_Count		G_Count	T_Count	N_Count	Max_count
# 1			35274720	2		34	1122616688	31.82	31	34	34	3	27	34	5876876	14288444	8319746	6728036	61618	35274720

for each in works:
 m = re.match("column", each)	# check header line
 if not m:
 				# print each
  data = each.split()
  G_count.append(int(data[-4]))	# load G content
  Q_mean.append(float(data[5]))	# load average score
 
G_count_red = G_count[20:]	#create partial list of last values

min = min(G_count_red)
min1 = min - int(min*0.1)	# accept 10% variation
max = max(G_count_red)
max1 = max + int(max*0.1)	# accept 10% variation

for pos, val in reversed(list(enumerate(G_count))):	# search value out of range - last in range is first nt of trimmed sequence
							# print pos, val, min1, max1
 if int(val) < min1 or int(val) > max1:
  cut5 = pos+2
  break

for pos, val in reversed(list(enumerate(Q_mean))):	# search value larger than threshold - first
							# print pos, val, treshscore
 if float(val) > float(treshscore):
  cut3 = pos+1
  break

print 'Cut 5prime is "', cut5
print 'Cut 3prime is "', cut3

   
