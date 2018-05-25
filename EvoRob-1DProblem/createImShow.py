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

xdist = np.linspace(0,1,numpoints)
xinit = xdist

lines = lines[1:]
fitArr = []

for line in lines:
	if(line=="\n"):
		continue
	fitArr.append([float(num) for num in line.split(",")[:-1]])

print len(fitArr)

label_locs = np.linspace(0,numpoints-1,11)

plt.figure()
imgplot=plt.imshow(np.array(fitArr))
imgplot.set_cmap('viridis')
plt.colorbar()
plt.xlabel("Initial Position of robot")
plt.xticks(label_locs,("0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1.0"))
plt.yticks(label_locs,("0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1.0"))
plt.ylabel("Position of distractor peak")
plt.savefig("Selected/fitness.eps")
plt.show()