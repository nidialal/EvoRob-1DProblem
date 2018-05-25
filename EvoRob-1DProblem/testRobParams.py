import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
from numpy import linalg as la
from math import log
import random
import os.path as op

numNodes = 3#int(raw_input("Input number of nodes:"))
numEdges = 0
adjMat = [[0]*numNodes for x in range(numNodes) ]
edgeList=[]#(0,1),(1,2),(2,0)]

t_rob = []
nodes = []

t1 = 1000

init = []

A = 0.5
C = 2
D = 10
E = 4
F = 2


def findNode(point): 
	for i in range(numNodes):
		arr = np.zeros(numNodes)
		arr.put(i,1)
		if(la.norm(point-arr) < 0.2):
			return i+1
	return 0

def sum_square(list):
	return sum(map(lambda x: x*x, list))

def sum_Exp4(list):
	return sum(map(lambda x: x*x*x*x, list))

def g(y_k,l):
	return -y_k * (np.square(np.square(y_k)-1) + l)


def fn_s(x):
	val = 0
	
	xmin = 0
	xmax = 2.0/9
	if(x>xmin and x<xmax):
		val = np.sin((x+1)*np.pi*9/2)+1
	
	xmin = 5.5/13
	xmax = 7.5/13
	if(x>xmin and x<xmax):
		val = np.sin(x*np.pi*13)+1
	
	xmin = 7.0/9
	xmax = 1
	if(x>xmin and x<xmax):
		val = np.sin(x*np.pi*9/2)+1
	
	return val/2

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
			A - (beta[i] * np.square(p_vals[edge[0]])) + (C * (y_squared - np.square(y_vals[i])))) + (gamma[i]*fn_s(x_val%1))
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

#print("Enter the adjacency matrix")
for i in range(numNodes):
	for j in range(numNodes):
		adjMat[i][j] = 0 if i==j else 1 #int(raw_input("Edge present from "+str(i+1) + " to " + str(j+1) + "(0/1)?"))
		if(adjMat[i][j]==1):
			numEdges+=1
			edgeList.append((i,j))
			init.append(random.uniform(0.01,0.11))
init.append(random.random())
'''
alpha=[0.7356660000436958, 0.1, -0.008458952250947949]
beta=[1.4, 1.4497456923125227, 1.4, 1.5773978430922075, 1.4, 1.4]
gamma=[0.59466123110656666, 0.69324424654526795, 0.62227499426877142, 0.68009445513350686, 0.12354814306528494, 0.0017866266678713382]
'''
'''
alpha=[-0.0012502653330501243, -0.15967661027763147, 0.0049301736501063606]
beta=[1.4, 1.4, 1.4, 1.4340929920998298, 1.5893428711495408, 1.4133171146514312]
gamma=[0.004771831905898464, 0, 0.20717444666277246, 0, 0, 0.14540689721666192]
'''

alpha=[-0.006277055900152817,-0.00035155286902427332,0.3772863170106737]
beta=[1.5530286259797004,1.4,1.5577080978345177,1.4791469776983825,1.5198986146670015,1.4]
gamma=[0.013473934985301018,0.006742951732529554,0,0.00927179515277148,0.083542553578403969,0]
'''
alpha=[0.06833431970523246, -0.0076494309320067405, -0.0019873691694322184]
beta=[1.4140015939664383, 1.5340264942314312, 1.572868739745367, 1.4, 1.4, 1.46334275197439]
gamma=[0.28803369126798928, 0.5114578598105987, 0.6837488749320128, 0.804521164215374, 0.1091332730920701, 0.20453770854105999]
'''
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
if(op.isfile("Analysis_1/counter.txt")):
	f=open("Analysis_1/counter.txt","r")
	counter = int(f.read())
	f.close()
else:
	counter = 1

f=open("Analysis_1/counter.txt","w")
f.write(str(counter+1))
f.close()

plt.figure("Analysis",dpi=300,figsize=(14.0, 7.0))
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
plt.xlabel(str(alpha))
plt.legend(loc="upper right")

plt.subplot(414)
plt.scatter(t_rob,[fn_s(x%1) for x in x_rob],s=0.1,label="Sensor value")
plt.legend(loc="upper right")
plt.savefig("Analysis_1/Analysis"+str(counter)+".png",dpi=500)

#_____________________________________________________________________________

