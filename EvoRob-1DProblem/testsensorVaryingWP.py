import numpy as np
import matplotlib.pyplot as plt
import random

x = np.linspace(0,1,1000)
y = []

xdiff = 0.5#random.uniform(-0.3,0.3)

st_pos = (8.0/9) +xdiff

for x1 in x:
	val = 0
	'''
	xmin = 5.5/13
	xmax = 7.5/13
	if(x1>xmin and x1<xmax):
		val = np.sin(x1*np.pi*13)+1
		val/=2
	'''
	xmin = 0
	xmax = (1.0/9) + xdiff
	if(x1>xmin and x1<xmax):
		val = np.sin((x1+1-xdiff)*np.pi*9/2)
		val = 0 if val < 0 else val

	xmin = (8.0/9) +xdiff
	xmax = 1
	if(x1>xmin and x1<xmax):
		val = np.sin((x1-xdiff)*np.pi*9/2)
		val = 0 if val < 0 else val

	y.append(val)

plt.plot(x,y)
plt.xlabel("Position of robot")
plt.ylabel("Sensor value")
plt.savefig('/Users/nidjac/Documents/Thesis/Poster/sensor.svg',dpi=100)
plt.show()