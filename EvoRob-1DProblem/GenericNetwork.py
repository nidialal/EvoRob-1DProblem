import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import odeint
from numpy import linalg as la
from math import log

def sum_square(list):
	return sum(map(lambda x: x*x, list))

def sum_Exp4(list):
	return sum(map(lambda x: x*x*x*x, list))

def g(y_k,l):
	return -y_k * (np.square(np.square(y_k)-1) + l)

def model(init,t,A,B,C,D,E,F,numNodes,numEdges,adjMat,egdeList):
	
	p_vals = init[0:numNodes]
	y_vals = init[numNodes:]
	p_squared = sum_square(p_vals)
	p_Exp4 = sum_Exp4(p_vals)
	y_squared = sum_square(y_vals)
	p_y_dot = []
	print adjMat
	for i in range(numNodes):
		pdot = (float)(p_vals[i] * ( F * (1 - p_squared) + D * (( np.square(p_vals[i]) * p_squared ) - p_Exp4)))
		row = adjMat[i]
		col = [r[i] for r in adjMat]
		Eterm = 0
		print "i:" + str(i)
		print adjMat[i]
		print "Row: " + str(row)
		print "Col: " + str(col)
		for j,val in enumerate(row):
			if(val>0):
				Eterm -= p_vals[i] * p_vals[j] * np.square(y_vals[edgeList.index((i,j))])
		for j,val in enumerate(col):
			if(val>0):
				Eterm += np.square(p_vals[j]) * np.square(y_vals[edgeList.index((j,i))])
		pdot += E * Eterm
		p_y_dot.append(pdot)
	
	for i,edge in enumerate(edgeList):
		ydot = g(y_vals[i], 
			A - (B * np.square(p_vals[edge[0]])) + (C * (y_squared - np.square(y_vals[i]))))
		p_y_dot.append(ydot)

	return tuple(p_y_dot)

def findNode(point):
	'''
	******* Checks if the solution is in a node
	'''

	if(la.norm(point-np.array((1,0,0,0))) < 0.2):
		return 1
	elif(la.norm(point-np.array((0,1,0,0))) < 0.2):
		return 2
	elif(la.norm(point-np.array((0,0,1,0))) < 0.2):
		return 3
	elif(la.norm(point-np.array((0,0,0,1))) < 0.2):
		return 4
	else:
		return 0

t = np.linspace(0,500,10000)
A = 0.5
B = 2.5
C = 2
D = 10
E = 4
F = 2

init = []

numNodes = int(raw_input("Enter number of nodes:"))
numEdges = int(raw_input("Enter number of edges:"))
edgeList = []
adjMat = [[0]*numNodes for x in range(numNodes) ]

for i in range(numNodes):
	init.append(float(raw_input("Input initial p"+str(i+1)+":")))

print("Enter the nodes connecting edges(0 to n-1)")
for i in range(numEdges):
	start = int(raw_input("Input start node:"))
	end = int(raw_input("Input end node:"))
	edgeList.append((start,end))
	adjMat[start][end] = 1
	init.append(float(raw_input("Input initial y"+str(i+1)+":")))

print adjMat

print edgeList

nodes = odeint(model,init,t,args=(A,B,C,D,E,F,numNodes,numEdges,adjMat,edgeList,),atol = 0.0000000000000001,rtol=0.00000000001)

file = open('nodes','w')

nodelist = []
for node in nodes:
	x = findNode(node[:numNodes])
	if(x>0):
		if(len(nodelist)==0):
			nodelist.append(x)
		else:
			if(nodelist[len(nodelist)-1]!=x):
				nodelist.append(x)


print "NodeList: " + str(nodelist)

plt.figure("P series")
for i in range(numNodes):
	plt.plot(t,[node[i] for node in nodes],label="p"+str(i+1))
plt.legend(loc="upper right")

plt.figure("Time Series")

for i in range(numEdges):
	plt.plot(t,[node[numNodes+i] for node in nodes],label="t"+str(i+1))
plt.legend(loc="upper right")
plt.show()

