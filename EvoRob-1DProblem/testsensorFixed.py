import numpy as np
import matplotlib.pyplot as plt

x1 = np.linspace(0,1,1000)
y = []

for x in x1:
	val = 0

	val = 0
	xmin = 0
	xmax = (1.0/9) + xdiff
	if(x>xmin and x<xmax):
		val = np.sin((x+1-xdiff)*np.pi*9/2)
		val = 0 if val < 0 else val

	xmin = 5.5/13
	xmax = 7.5/13
	if(x>xmin and x<xmax):
		val = np.sin(x*np.pi*13)+1
		val = val/2

	xmin = (8.0/9) +xdiff
	xmax = 1
	if(x>xmin and x<xmax):
		val = np.sin((x-xdiff)*np.pi*9/2)
		val = 0 if val < 0 else val

	y.append(val)

plt.plot(x1,y)
plt.show()