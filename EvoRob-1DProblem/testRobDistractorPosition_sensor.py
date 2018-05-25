import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
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

t1 = 1000

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


def fn_s(x,xdiff,nar_peak=True,wid_peak=True):

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
	if(nar_peak):
		xmin = 5.5/13
		xmax = 7.5/13
		if(x>xmin and x<xmax):
			val1 = np.sin(x*np.pi*13)+1
			val1/=2

	#WidePeak
	if(wid_peak):
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

def model(t,init,nar_p,wide_p):
	
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
			A - (beta[i] * np.square(p_vals[edge[0]])) + (C * (y_squared - np.square(y_vals[i])))) + (gamma[i]*fn_s(x_val%1,xdiff,nar_p,wide_p))
		p_y_x_dot.append(ydot)

	x_dot = 0
	for i in range(numNodes):
		x_dot += alpha[i] * p_vals[i]

	p_y_x_dot.append(x_dot)

	return p_y_x_dot

#__________________________________________________________________________
'''
counter = 0
folderName = "Analysis_"
if(op.isfile("counter.txt")):
	f=open("counter.txt","r")
	counter = int(f.read())
	f.close()
else:
	counter = 2
counter = 15
folderName += str(counter-1)

if not op.exists(folderName):
	os.makedirs(folderName)
'''
'''
folderName = "Selected"

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

alpha= [0.066854506661658208, -0.016054975161765565, -0.0099148250331080774]
beta= [1.5820192906249351, 1.40011004603279, 1.4947056086663628, 1.4, 1.539039341360442, 1.4372235419325299]
gamma= [1, 1, 0.3066188749894846, 0.23354678297422582, 0.3899771328992091, 0.17733083804406102]
'''


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
conditions = [[0.0,True,True,"at_0"]]#,[-1.5/9,True,True,"leftDist"],[-0.35,True,True,"leftDist_Merged"],[0.5,True,True,"mid_overlap"],[1.5/9,True,True,"rightDist"],[0.35,True,True,"rightDist_Merged"],[0.0,True,False,"noDist"],[0.0,False,True,"noTarget"]]

for c in conditions:

	xdist = c[0]
	nar_p = c[1]
	wide_p = c[2]
	desc = c[3]

	set_params(xdist)

	x_rob = []
	t_rob = []
	nodes = []
	node_index = []
	y_index = []
	x_cur = 0
	t_1=datetime.datetime.now()

	r = ode(model).set_integrator('vode',atol=1e-6)
	r.set_initial_value(init).set_f_params(nar_p,wide_p)
	dt = 0.05

	while r.successful() and r.t<t1:
		x_cur = r.y[-1]
		x_rob.append(x_cur%1)
		nodes.append(r.y[0:-1])
		node = findNode(r.y[:numNodes])
		if(node!=0):
			if(len(node_index)==0):
				node_index.append(node)
			else:
				if(node_index[len(node_index)-1]!=node):
					node_index.append(node)
		yval = findY(r.y[numNodes:-1])
		if(yval!=0):
			if(len(y_index)==0):
				y_index.append(yval)
			else:
				if(y_index[len(y_index)-1]!=yval):
					y_index.append(yval)
		t_rob.append(r.t)
		r.integrate(t1,step=True)

	t2=	datetime.datetime.now()
	print (t2-t_1).total_seconds()

	color_list = plt.cm.cool(np.linspace(0, 1, 7))
	color_list=np.delete(color_list,1,0)

	plt.figure("Analysis_"+desc,figsize=(12.0, 6.0))

	plt.subplot(411)
	for i in range(numNodes):
		plt.plot(t_rob,[node[i] for node in nodes],label="p"+str(i+1),color=color_list[i])
	plt.legend(loc="upper right")

	plt.subplot(412)
	for i in range(numEdges):
		plt.plot(t_rob,[node[numNodes+i] for node in nodes],label="y"+str(i+1),color=color_list[i])
	plt.legend(loc="upper right")

	plt.subplot(413)
	plt.scatter(t_rob,[x%1 for x in x_rob],c=[color_list[findNode(node[:numNodes])] for node in nodes],s=0.1,label="Robot position")
	x_sensor = np.linspace(0,1,200)
	plt.plot([fn_s(x,xdiff,nar_p,wide_p)*40 for x in x_sensor],x_sensor,color='#ca3587')
	plt.xlabel("Robot Position")
	plt.legend(loc="upper right")

	plt.subplot(414)
	plt.scatter(t_rob,[fn_s(x%1,xdiff,nar_p,wide_p) for x in x_rob],s=0.1,label="Sensor value",color='#ca3587')
	plt.legend(loc="upper right")

	plt.savefig(folderName+"/Analysis_"+desc+".eps")

	plt.figure("RobMotion"+desc,figsize=(10,2))
	plt.scatter(t_rob,[x%1 for x in x_rob],s=0.1,color='#3fb7e9')
	x_sensor = np.linspace(0,1,200)
	plt.plot([fn_s(x,xdiff,nar_p,wide_p)*40 for x in x_sensor],x_sensor,color='#ca3587')
	plt.xlabel("Integration time/Sensor Value")
	plt.ylabel("Robot Position")
	plt.subplots_adjust(bottom=.25, left=.25)
	plt.savefig(folderName+"/RobMotion_"+desc+".eps")



#_____________________________________________________________________________

