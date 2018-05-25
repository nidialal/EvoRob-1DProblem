import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from numpy import linalg as la
from math import log
import random
from datetime import datetime
import os.path as op
import os
import datetime


numNodes = 3
numEdges = 0
adjMat = [[0]*numNodes for x in range(numNodes) ]
edgeList=[]

t_rob = []
nodes = []

t1 = 2000

init = []

A = 0.5
C = 2
D = 10
E = 4
F = 2

xdiff = 0

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

def model(init,t):
	
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


folderName = "../../Reports/figures/1D/"

for i in range(numNodes):
	for j in range(numNodes):
		adjMat[i][j] = 0 if i==j else 1
		if(adjMat[i][j]==1):
			numEdges+=1
			edgeList.append((i,j))

alpha= [0.066854506661658208, -0.016054975161765565, -0.0099148250331080774]
beta= [1.5820192906249351, 1.40011004603279, 1.4947056086663628, 1.4, 1.539039341360442, 1.4372235419325299]
gamma= [1, 1, 0.3066188749894846, 0.23354678297422582, 0.3899771328992091, 0.17733083804406102]
init = [0.9307040982988228, 0.07834107744994431, 0.016283740224093738, 0.040660717110571114, 0.09700957946018898, 0.09701901081633599, 0.07674411748194522, 0.010574052587478358, 0.021614691301056338, 0]
near_axis = 0
t=np.linspace(0,200,2000)

conditions = [[0.1,0.3,"areaA"],[0.42,0.2,"AreaB"],[0.42,0.5,"AreaC"],[0.55,0.2,"AreaD"],[0.56,0.6,"AreaE"],[0.8,0.2,"AreaF"]] #[0.9,0.427,"areaH"],[0.2,0.427,"areaG"]]#,

color_list = plt.cm.jet(np.linspace(0, 1, 7))
color_list=np.delete(color_list,1,0)

for c in conditions:

	xdist = c[0]
	init[-1] = c[1]
	desc = c[2] #("at"+str(xdist)+"_starting"+str(init[-1])).replace(".","_")

	set_params(xdist)

	x_rob = []
	t_rob = t
	nodes = []
	node_index = []
	y_index = []
	x_cur = 0

	sol = odeint(model,init,t)

	nodes = sol#[0]

	x_rob = [node[-1] for node in nodes]

	fig=plt.figure("Analysis_"+desc,figsize=(12.0, 8.0))

	plt.subplot(411)
	for i in range(numNodes):
		plt.plot(t_rob,[node[i] for node in nodes],label="p"+str(i+1),color=color_list[i])
	plt.legend(loc="upper right")
	plt.xlabel('Time units')
	plt.ylabel('Value of p cells')
	plt.xlim(0, 200)

	plt.subplot(412)
	for i in range(numEdges):
		plt.plot(t_rob,[node[numNodes+i] for node in nodes],label="y"+str(i+1),color=color_list[i])
	plt.legend(loc="upper right")
	plt.xlabel('Time units')
	plt.ylabel('Value of y cells')
	plt.xlim(0, 200)

	plt.subplot(413)
	plt.scatter(t_rob,[x%1 for x in x_rob],c=[color_list[findNode(node[:numNodes])] for node in nodes],s=0.1)
	x_sensor = np.linspace(0,1,200)
	plt.plot([fn_s(x,xdiff)*10 for x in x_sensor],x_sensor,color='#ca3587')
	#plt.axvspan(t_rob[int(0.8*len(t_rob))],t_rob[-1], color='yellow', alpha=0.5)
	plt.xlabel('Time units')
	plt.ylabel('Position of robot')
	plt.xlim(0, 200)

	plt.subplot(414)
	plt.scatter(t_rob,[fn_s(x%1,xdiff) for x in x_rob],s=0.1,color='#ca3587')
	plt.xlabel('Time units')
	plt.ylabel('Sensor value')
	plt.xlim(0, 200)

	fig.savefig(folderName+"/Analysis_"+desc+".pdf",bbox_inches='tight',pad_inches=0.0,dpi=75)
	'''
	plt.figure("RobMotion"+desc,figsize=(10,2))
	plt.scatter(t_rob,[x%1 for x in x_rob],s=0.1,color='#3fb7e9')
	x_sensor = np.linspace(0,1,200)
	plt.plot([fn_s(x,xdiff)*40 for x in x_sensor],x_sensor,color='#ca3587')
	plt.xlabel("Integration time/Sensor Value")
	plt.ylabel("Robot Position")
	plt.subplots_adjust(bottom=.25, left=.25)
	#plt.savefig(folderName+"/RobMotion_"+desc+".eps")
	'''
	#plt.show()


#_____________________________________________________________________________