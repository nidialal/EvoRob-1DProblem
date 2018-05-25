import numpy as np
import matplotlib.pyplot as plt
import random

x = np.linspace(0,1,1000)
y = []
dist =[]
xdiff = random.uniform(-0.2,0.2)

xmin = 0
xmax = (1.0/9) + xdiff
if(xmax < 0):
	st_pos = (8.0/9) +xdiff
	np_pos=random.uniform(0,st_pos-(1.0/13))
else:
	st_pos = (8.0/9) +xdiff
	np_pos=random.uniform(xmax,st_pos-(1.0/13))

x_target = np_pos +(1.0/26)

for x1 in x:
	value = 0
	xmin = 0
	xmax = (1.0/9) + xdiff

	if(x1>xmin and x1<xmax):
		value = np.sin((x1+1-xdiff)*np.pi*9/2)
		value = 0 if value < 0 else value

	xmin = np_pos
	xmax = xmin+(1.0/13)
	if(x1>xmin and x1<xmax):
		value = np.sin((x1-xmin)*np.pi*13)
		value = 0 if value < 0 else value

	xmin = (8.0/9) +xdiff
	xmax = 1
	if(x1>xmin and x1<xmax):
		value = np.sin((x1-xdiff)*np.pi*9/2)
		value = 0 if value < 0 else value

	y.append(value)
print "Target:" + str(x_target)
for pos in np.linspace(0,1,100):
	print pos
	print [abs(pos-x_target),abs(1-pos+x_target),abs(1-x_target+pos)]
	print "Distance = " + str(np.min([abs(pos-x_target),abs(1-pos+x_target),abs(1-x_target+pos)]))
plt.plot(x,y)
plt.show()