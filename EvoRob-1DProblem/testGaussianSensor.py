import numpy as np
import matplotlib.pyplot as plt
import random
from math import *

x = np.linspace(0,1,1000)
y = []

def circular_distance(a,b):
	return min(abs(a-b),1-abs(a-b))

def fn_g(x, mu, sigma):
	return exp(-(x-mu)**2/(2*sigma**2))
	


def fn_s(distPos):

	mu_narrow = 0.5
	sig_narrow = 1.0/26
	mu_wide = distPos
	sig_wide = 1.0/18
	for x1 in x:

		g_narrow = fn_g(x1,mu_narrow,sig_narrow)
		g_wide = fn_g(circular_distance(x1,1.0),mu_wide,sig_wide)
		g = max(g_narrow,g_wide)

		y.append(g)

def testfn(xdiff):
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


		xmin = (8.0/9) +xdiff
		xmax = 1	
		if(test1==False):
			xmax = (10.0/9) + xdiff

		if(x1>xmin and x1<xmax):
			val3 = np.sin((x1-xdiff)*np.pi*9/2)
			val3 = 0 if val3 < 0 else val3

		val = max(val1,val2,val3)

		y.append(val)

for xval in [0.1]:
	distPos = xval
	y =[]
	fn_s(distPos)
	plt.plot(x,y)
	y = []
	testfn(distPos)
	plt.plot(x,y)
	plt.xlabel("Position of robot")
	plt.ylabel("Sensor value")

	plt.show()
