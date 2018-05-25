import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
from numpy import linalg as la
from math import log
import random
from datetime import datetime
import os.path as op
import os


numNodes = 3
numEdges = 0
adjMat = [[0]*numNodes for x in range(numNodes) ]
edgeList=[]

t_rob = []
nodes = []

t1 = 200

init = []

A = 0.5
C = 2
D = 10
E = 4
F = 2

xdiff = 0

#Analysis14
alpha= [0.066854506661658208, -0.016054975161765565, -0.0099148250331080774]
beta= [1.5820192906249351, 1.40011004603279, 1.4947056086663628, 1.4, 1.539039341360442, 1.4372235419325299]
gamma= [1, 1, 0.3066188749894846, 0.23354678297422582, 0.3899771328992091, 0.17733083804406102]
'''
alpha = [0.015653494680099689, -0.017189886807194774, -0.002266398348500176]
beta = [1.5544651638807674, 1.5417802524752053, 1.5013509311970452, 1.5604450164122607, 1.5312561696848859, 1.402247596625195]
gamma = [0.6092215257899631, 1, 0.3371989927354151, 1, 0.5819841347180218, 0.5942477151132515]
'''
def set_params(x):
	global xdiff
	xdiff = x

def findNode(point): 
	for i in range(numNodes):
		arr = np.zeros(numNodes)
		arr.put(i,1)
		if(la.norm(point-arr) < 0.2):
			return i
	return numNodes+1

def findY(yvals):
	for i,y in enumerate(yvals):
		if(y>0.99):
			return i+1
	return 0

def sum_square(list):
	return sum(map(lambda x: x*x, list))

def sum_Exp4(list):
	return sum(map(lambda x: x*x*x*x, list))

def g(y_k,l):
	return -y_k * (np.square(np.square(y_k)-1) + l)


def fn_s(x,xdiff,wid_peak=True):
	nar_peak=True
	val = 0
	if(nar_peak):
		xmin = 5.5/13
		xmax = 7.5/13
		if(x>xmin and x<xmax and nar_peak):
			val = np.sin(x*np.pi*13)+1
			return val/2
	if(wid_peak):
		xmin = 0
		xmax = (1.0/9) + xdiff
		if(x>xmin and x<xmax):
			val = np.sin((x+1-xdiff)*np.pi*9/2)
			val = 0 if val < 0 else val
			return val
		
		
		xmin = (8.0/9) +xdiff
		xmax = 1
		if(x>xmin and x<xmax):
			val = np.sin((x-xdiff)*np.pi*9/2)
			val = 0 if val < 0 else val
			return val
	
	return val

def model(t,init,wid_peak):
	
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
			A - (beta[i] * np.square(p_vals[edge[0]])) + (C * (y_squared - np.square(y_vals[i])))) + (gamma[i]*fn_s(x_val%1,xdiff,wid_peak))
		p_y_x_dot.append(ydot)

	x_dot = 0
	for i in range(numNodes):
		x_dot += alpha[i] * p_vals[i]

	p_y_x_dot.append(x_dot)

	return p_y_x_dot

numpoints = 401
fit_arr =[]
xdist = np.linspace(0,0.5,int(numpoints/2),endpoint = False)
xdist = np.append(xdist,np.linspace(-0.5,0,int(numpoints/2)+1))

xinit = np.linspace(0,1,numpoints)

for i in range(numNodes):
	for j in range(numNodes):
		adjMat[i][j] = 0 if i==j else 1
		if(adjMat[i][j]==1):
			numEdges+=1
			edgeList.append((i,j))

for i in xdist:
	wid_peak = True
	scores = []

	if(i>(2.75/9) or i<(-2.75/9)):
		wid_peak = False
	
	for j in xinit:
		init = []
		near_axis = random.choice(range(numNodes))
		for y in range(numNodes):
			if(y==near_axis):
				init.append(random.uniform(0.89,0.99))
			else:
				init.append(random.uniform(0.01,0.11))
		for y in range(numNodes):
			for z in range(numNodes):
				if(adjMat[y][z]==1):
					init.append(random.uniform(0.01,0.11))
		init.append(0)
		init[-1] = j
		set_params(i)
		x_rob = []
		t_rob = []
		nodes = []
		node_index = []
		y_index = []
		x_cur = 0

		r = ode(model).set_integrator('vode',atol=1e-6)
		r.set_initial_value(init).set_f_params(wid_peak)
		dt = 0.05
		while r.successful() and r.t<t1:
			x_cur = r.y[-1]
			x_rob.append(x_cur%1)
			r.integrate(t1,step=True)

		rob_dist = [abs(pos-0.5) for pos in x_rob[int(0.8*len(x_rob)):]]
		score = np.max(rob_dist)
		scores.append(score)

	fit_arr.append(scores)

np.savetxt('fit_arr.out',fit_arr)

'''
plt.figure()
imgplot=plt.imshow(np.array(fit_arr))
imgplot.set_cmap('YlGnBu')
plt.colorbar()
plt.xlabel("Initial Position of robot")
plt.ylabel("Position of distractor peak")
fig.tight_layout()
plt.savefig("Graph_1.svg")
'''
