import numpy as np
import matplotlib.pyplot as plt
import random
from math import *

x = np.linspace(0,1,1000)
y = []

def circular_distance(a,b):
	return min(abs(a-b),1-abs(a-b))

def fn_gauss(x, mu, sigma):
	return exp(-(x-mu)**2/(2*sigma**2))
	
def fn_s(distPos):

	mu = distPos
	sig = 1.0/26
	for x1 in x:

		g_wide = max(fn_gauss(x1,mu,sig),fn_gauss(x1,mu+1,sig),fn_gauss(x1,mu-1,sig))
		#g_wide = fn_gauss(x1,mu,sig)
		y.append(g_wide)



for xval in [0.01,0.02,0.03,0.04,0.05,0.1,0.9,0.95,0.96,0.97,0.98,0.99]:
	distPos = xval
	y =[]
	fn_s(distPos)
	plt.plot(x,y,label=str(xval))
	plt.axvline(x=xval,color='r',linewidth=0.25)
	plt.xlabel("Position of robot")
	plt.ylabel("Sensor value")
	plt.legend(loc="upper right")
	plt.show()
