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

t1 = 400

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


def fn_s(x,xdiff):
	test_var = (10.0/9) +xdiff
	test1 = False
	if(test_var>1):
		test1 = True

	test_var = xdiff - (1.0/9)
	test2=False
	if(test_var>0):
		test2 = True

	val = 0
	val1 = 0
	val2 = 0
	val3 = 0

	#Narrow Peak
	xmin = 5.5/13
	xmax = 7.5/13
	if(x>xmin and x<xmax):
		val1 = np.sin(x*np.pi*13)+1
		val1/=2

	#WidePeak
	xmin = 0
	xmax = (1.0/9) + xdiff
	if(test2==True):
		xmin = xdiff-(1.0/9)
	if(x>xmin and x<xmax):
		val2 = np.sin((x+1-xdiff)*np.pi*9/2)
		val2 = 0 if val2 < 0 else val2

	#WidePeak
	xmin = (8.0/9) +xdiff
	xmax = 1	
	if(test1==False):
		xmax = (10.0/9) + xdiff
	if(x>xmin and x<xmax):
		val3 = np.sin((x-xdiff)*np.pi*9/2)
		val3 = 0 if val3 < 0 else val3

	val = max(val1,val2,val3)
	
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

numpoints = int(raw_input("Enter the resolution"))
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
fName = "fitness_"+datetime.now().strftime("%Y-%m-%d--%H:%M:%S")+".txt"
fit_file = open(fName,"a+")
fit_file.write(str(numpoints)+"\n")
fit_file.close()
for i in xdist:
	i=round(i,5)
	wid_peak = True
	scores = []

	if(i>(2.75/9) or i<(-2.75/9)):
		wid_peak = False
	
	for j in xinit:
		j = round(j,5)
		print (i,j)

		init[-1] = j
		set_params(i)
		x_rob = []
		t_rob = []
		nodes = []
		node_index = []
		y_index = []
		x_cur = 0

		r = ode(model).set_integrator('vode',atol=1e-6)
		r.set_initial_value(init)
		dt = 0.05
		while r.successful() and r.t<t1:
			x_cur = r.y[-1]
			x_rob.append(x_cur%1)
			r.integrate(t1,step=True)

		rob_dist = [abs(pos-0.5) for pos in x_rob[int(0.8*len(x_rob)):]]
		score = np.max(rob_dist)
		scores.append(score)
		fit_file = open(fName,"a+")
		fit_file.write(str(score)+",")
		fit_file.close()
	fit_file = open(fName,"a+")
	fit_file.write("\n\n")
	fit_file.close()

