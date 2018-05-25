import numpy as np
import matplotlib.pyplot as plt
import sys

fName = sys.argv[1]
if(len(sys.argv)==1):
	print "Please provide input filename"
	sys.exit()
fitFile = open(fName,"r")

lines = fitFile.readlines()
numpoints = int(lines[0])

fitFile.close()

fitArr = []

for line in lines:
	if(line=="\n"):
		continue
	fitArr.append([float(num) for num in line.split(",")[:-1]])

plt.figure()
imgplot=plt.imshow(np.array(fitArr))
imgplot.set_cmap('rainbow')
plt.colorbar()
plt.show()

