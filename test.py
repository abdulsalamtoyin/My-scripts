fh = open('/Users/Toyin/Desktop/diff4_out/gene_exp.diff')


fh.readline()
# # Get rid of the header line

fpkms1 = []
fpkms2 = []

for line in fh:
     data = line.split()
     # If the values are zero, then a loglog plot will fail.
     fpkm1 = float(data[7])
     fpkm2 = float(data[8])
     if fpkm1 == 0 or fpkm2 ==0:
         continue
     fpkms1.append(fpkm1)
     fpkms2.append(fpkm2)

from matplotlib import pyplot as mpl
from numpy import array

mpl.loglog(array(fpkms1), array(fpkms2), ',')

##Plot a diagonal line corresponding to equal expression
mpl.loglog([min(fpkms1), max(fpkms1)],
            [min(fpkms1), max(fpkms1)], 'r:')
mpl.show()

