import numpy as np
import matplotlib.pyplot as plt
import random

x = np.linspace(0,1,1000)
y = []

xdiff = 0.5

def testfn():
	test_var = (10.0/9) +xdiff
	test1 = False
	if(test_var>1):
		test1 = True

	test_var = xdiff - (1.0/9)
	test2=False
	if(test_var>0):
		test2 = True

	for x1 in x:
		val = 0
		val1 = 0
		val2 = 0
		val3 = 0

		xmin = 5.5/13
		xmax = 7.5/13
		if(x1>xmin and x1<xmax):
			val1 = np.sin(x1*np.pi*13)+1
			val1/=2

		xmin = 0
		xmax = (1.0/9) + xdiff
		if(test2==True):
			xmin = xdiff-(1.0/9)
		if(x1>xmin and x1<xmax):
			val2 = np.sin((x1+1-xdiff)*np.pi*9/2)
			val2 = 0 if val2 < 0 else val2
			if(val2>0.9998):
				print x1

		xmin = (8.0/9) +xdiff
		xmax = 1	
		if(test1==False):
			xmax = (10.0/9) + xdiff

		if(x1>xmin and x1<xmax):
			val3 = np.sin((x1-xdiff)*np.pi*9/2)
			val3 = 0 if val3 < 0 else val3
			if(val3>0.9998):
				print x1

		val = max(val1,val2,val3)

		y.append(val)

for xval in [0.0]:
	xdiff = xval
	print "Xdiff:"+str(xdiff)
	y =[]
	testfn()
	plt.plot(x,y)
	plt.xlabel("Position of robot")
	plt.ylabel("Sensor value")
	plt.show()
