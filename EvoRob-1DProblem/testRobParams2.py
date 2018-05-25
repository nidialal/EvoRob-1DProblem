import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
from numpy import linalg as la
from math import log
import random
import os.path as op
import os


numNodes = 2#int(raw_input("Input number of nodes:"))
numEdges = 0
adjMat = [[0]*numNodes for x in range(numNodes) ]
edgeList=[]

t_rob = []
nodes = []

t1 = 500

init = []

A = 0.5
C = 2
D = 10
E = 4
F = 2

xdiff = 0
new = False

def set_params():
	global xdiff
	xdiff = random.uniform(-0.2,0.2)

def sum_square(list):
	return sum(map(lambda x: x*x, list))

def sum_Exp4(list):
	return sum(map(lambda x: x*x*x*x, list))

def g(y_k,l):
	return -y_k * (np.square(np.square(y_k)-1) + l)

def findNode(point): 
	for i in range(numNodes):
		arr = np.zeros(numNodes)
		arr.put(i,1)
		if(la.norm(point-arr) < 0.2):
			return i+1
	return 0

def fn_s(x,xdiff):

	val = 0
	xmin = 0
	xmax = (1.0/9) + xdiff
	if(x>xmin and x<xmax):
		val = np.sin((x+1-xdiff)*np.pi*9/2)
		val = 0 if val < 0 else val
		return val

	xmin = 5.5/13
	xmax = 7.5/13
	if(x>xmin and x<xmax):
		val = np.sin(x*np.pi*13)+1
		return val/2

	xmin = (8.0/9) +xdiff
	xmax = 1
	if(x>xmin and x<xmax):
		val = np.sin((x-xdiff)*np.pi*9/2)
		val = 0 if val < 0 else val
		return val

	return val


def model(t,init):
	
	p_vals = init[0:numNodes]
	y_vals = init[numNodes:-1]
	x_val = init[-1]
	p_squared = sum_square(p_vals)
	p_Exp4 = sum_Exp4(p_vals)
	y_squared = sum_square(y_vals)
	p_y_x_dot = []
	for i in range(numNodes):
		pdot = (float)(p_vals[i] * ( F * (1 - p_squared) + D * (( np.square(p_vals[i]) * p_squared ) - p_Exp4)))
		row = adjMat[i]
		col = [r[i] for r in adjMat]
		Eterm = 0
		for j,val in enumerate(row):
			if(val>0):
				Eterm -= p_vals[i] * p_vals[j] * np.square(y_vals[edgeList.index((i,j))])
		for j,val in enumerate(col):
			if(val>0):
				Eterm += np.square(p_vals[j]) * np.square(y_vals[edgeList.index((j,i))])
		pdot += E * Eterm
		p_y_x_dot.append(pdot)
	
	for i,edge in enumerate(edgeList):
		ydot = g(y_vals[i], 
			A - (beta[i] * np.square(p_vals[edge[0]])) + (C * (y_squared - np.square(y_vals[i])))) + (gamma[i]*fn_s(x_val%1,xdiff))
		p_y_x_dot.append(ydot)

	x_dot = 0
	for i in range(numNodes):
		x_dot += alpha[i] * p_vals[i]

	p_y_x_dot.append(x_dot)

	return p_y_x_dot

#__________________________________________________________________________

near_axis = random.choice(range(numNodes))
for i in range(numNodes):
	if(i==near_axis):
		init.append(random.uniform(0.89,0.99))
	else:
		init.append(random.uniform(0.01,0.11))

for i in range(numNodes):
	for j in range(numNodes):
		adjMat[i][j] = 0 if i==j else 1
		if(adjMat[i][j]==1):
			numEdges+=1
			edgeList.append((i,j))
			init.append(random.uniform(0.01,0.11))
init.append(random.random())

alpha = [0.0034601758140425966, -0.0033158891550608268]
beta = [1.4213619977314964, 1.571601421747769]
gamma = [0.39021643729082356, 0.3884054111263685]

set_params()

x_rob = []
x_cur = 0
r = ode(model).set_integrator('vode')
r.set_initial_value(init)
dt = 0.05
while r.successful() and r.t<t1:
	x_cur = r.y[-1]
	x_rob.append(x_cur%1)
	nodes.append(r.y[0:-1])
	t_rob.append(r.t)
	r.integrate(t1,step=True)

counter = 0
folderName = "Analysis_17"
'''
if(op.isfile("counter.txt")):
	f=open("counter.txt","r")
	counter = int(f.read())
	f.close()
else:
	counter = 1

folderName += str(counter)
if new:
	f=open("counter.txt","w")
	f.write(str(counter+1))
	f.close()

if not op.exists(folderName):
	os.makedirs(folderName)
'''
if(op.isfile(folderName+"/counter.txt")):
	f=open(folderName+"/counter.txt","r")
	counter = int(f.read())
	f.close()
else:
	counter = 1

f=open(folderName+"/counter.txt","w")
f.write(str(counter+1))
f.close()

plt.figure(figsize=(14.0, 7.0))
plt.subplot(411)
for i in range(numNodes):
	plt.plot(t_rob,[node[i] for node in nodes],label="p"+str(i+1))
plt.legend(loc="upper right")

plt.subplot(412)
for i in range(numEdges):
	plt.plot(t_rob,[node[numNodes+i] for node in nodes],label="t"+str(i+1))
plt.legend(loc="upper right")

usecolor=['w','r','c','m']
plt.subplot(413)
plt.scatter(t_rob,[x%1 for x in x_rob],c=[usecolor[findNode(node[:numNodes])] for node in nodes],s=0.1,label="Robot position")
x_sensor = np.linspace(0,1,200)
plt.plot([fn_s(x%1,xdiff)*50 for x in x_sensor],x_sensor)
plt.xlabel(str(alpha))
plt.legend(loc="upper right")

plt.subplot(414)
plt.scatter(t_rob,[fn_s(x%1,xdiff) for x in x_rob],s=0.1,label="Sensor value")
plt.legend(loc="upper right")

plt.savefig(folderName+"/Analysis_"+str(counter)+".png",dpi=500)


plt.show()

#_____________________________________________________________________________

