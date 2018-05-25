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

xdist = np.linspace(0.35,0.45,numpoints)
xinit = np.linspace(0.4,0.6,numpoints)

lines = lines[1:]
fitArr = []

for line in lines:
	if(line=="\n"):
		continue
	fitArr.append([float(num) for num in line.split(",")[:-1]])
fitArr=fitArr[:-1]
label_locs = np.linspace(0,numpoints-1,9)


plt.figure()
imgplot=plt.imshow(np.array(fitArr))
imgplot.set_cmap('rainbow')
plt.colorbar()
plt.xlabel("Initial Position of robot")
plt.xticks(label_locs,("0.4","0.425","0.45","0.475","0.5","0.525","0.55","0.575","0.6"))
plt.yticks(label_locs,("0.35","0.3625","0.375","0.3875","0.4","0.4125","0.425","0.4375","0.45"))
plt.ylabel("Position of distractor peak")
plt.show()

