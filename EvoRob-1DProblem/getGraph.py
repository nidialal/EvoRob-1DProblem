import numpy as np
import matplotlib.pyplot as plt
import sys

fName = sys.argv[1]
if(len(sys.argv)==1):
	print "Please provide input filename"
	sys.exit()
fitFile = open(fName,"r")

lines = fitFile.readlines()

fitFile.close()
fitArr = []
fitArr=[float(num) for num in lines[0].split(",")[:-1]]

plt.figure()
plt.plot(fitArr)
plt.xlabel("No of Tournament-->")
plt.ylabel("Peak fitness of population-->")
plt.savefig("Selected/PeakFitness.eps")
plt.show()