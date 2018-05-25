import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
from pprint import pprint

t1=10
nodes = []
t_nodes = []
initial = 0

def model(t,init):
	'''
	******* Differential equations to solve
	'''
	x = init
	xdot = 1.3 - fn_s(x)
	return xdot

def fn_s(x):
	val = np.sin(1*x) + np.cos(3*x) -0.5
	return val if val > 0 else 0

def solout(t,y):
	global nodes
	print y
	nodes.append(y[0])
	t_nodes.append(t)
	return 0


r = ode(model).set_integrator('dopri5',nsteps=10000)
r.set_initial_value(initial)
r.set_solout(solout)
r.integrate(t1)

print nodes
plt.figure("Robot motion")
plt.plot(t_nodes,nodes)

plt.figure("Sensor output")
plt.plot(np.linspace(0,2,100),[fn_s(x) for x in np.linspace(0,2,100)])
plt.show()
