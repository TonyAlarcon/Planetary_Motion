###################################################################################
#	Program determines solar planetary trajectory from Newton's law of 
#	universal gravitation, F = G m1 * m2 / r^2 
# 
#	**This is a work in process program, recommended updates suggested on
#	comment**
#
###################################################################################

import numpy as np
import math as math 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from decimal import Decimal


#Gravitational Constant SI Units
G = 6.67408E-11

###################################################################################
#	Class describes an astronomical body by defining position=(x,y), velocity = (vx,vy)
# 	radius (spherical approximation is assumed), mass, and orbit
#
#	Orbit defines the time (in Solar days) to complete an orbit around the Sun.This 
#	parameter is optional, default is zero. 
#
#	*Update should include switching position and velocity to arrays to reduce parameters.
###################################################################################


class Body:
	def __init__(self, x, y, vx, vy, mass, radius, orbit = 0): 
		self.x = x
		self.y = y
		self.vx = vx
		self.vy = vy
		self.mass = mass
		self.radius = radius
		self.orbit = orbit
		
	
		
		
#Initializing Sun and planets (No pluto)
	
sun = Body(0.0 , 0.0, 0.0 , 0.0, 1.989E30, 695700000)
earth = Body(-147095000000 , 0.0, 0.0 , -30300.0 , 5.972E24, 695700000, 365)
mercury = Body(-46000000000, 0.0, 0.0 , -58980, 0.33011E24, 2439700, 88)
venus = Body(-107480000000, 0.0, 0.0, -35260, 4.8675E24, 6051800, 225)
mars = Body(-206620000000, 0.0, 0.0, -26500, 6.4171E23, 3389500, 687)
jupiter = Body(-740520000000, 0.0, 0.0, -13720, 1898.19E24, 71492000, 4332)
saturn = Body(-1352550000000, 0.0, 0.0, -10180, 568.34E24, 54364000, 10760)
uranus = Body(-2741300000000, 0.0, 0.0, -7110, 86.813E24, 24973000, 30700)
neptune = Body(-4444450000000, 0.0, 0.0, -5500, 102.413E24, 24341000, 60200)


###################################################################################
#	Functions uses Newton's law of universal gravitation, F = G m1 * m2 / r^2 
#	to determine acceleration in x and y components. 
#
#	Functions will be used in RK4 to determine trajectory.	
#
#	*Update might include merging into one function/class
#	*Update should take into account foces due to planets, which is ignored
#	at the moment. 
###################################################################################

def ax(x,y):
	a = - G * sun.mass * x/ (x**2 + y**2)**(3/2.0)
	return a
	
def ay(x,y):
	a = - G * sun.mass * y/ (x**2 + y**2)**(3/2.0)
	return a
	


#Time Step is one Earth Day in seconds
h = 86400


###################################################################################
#	RK4 returns 2-element position array (x,y) 
#
#	Functions will be used in RK4 to determine trajectory.	
#
#	*Update should include control flow statement to return velocity if desired
#	an example can be seen in Cyclotron repo 
#
###################################################################################

def rk4(p0x,p0y,v0x,v0y, orbit):
	h = 86400
	r = []
	vx = []
	vy = []	
	
	for i in range(0,orbit):
		p_ix = p0x
		v_ix = v0x
		p_iy = p0y
		v_iy = v0y

		k1vx = h * ax( p_ix, p_iy )
		k1vy = h * ay( p_ix, p_iy )
		l1x = h * v_ix
		l1y = h * v_iy

		k2vx = h * ax( p0x + l1x * 0.5, p0y + l1y * 0.5 )
		k2vy = h * ay( p0x + l1x * 0.5, p0y + l1y * 0.5 )
		l2x = h * (v0x + k1vx * 0.5)
		l2y = h * (v0y + k1vy * 0.5)
	
		k3vx = h * ax( p0x + l2x * 0.5, p0y + l2y * 0.5 )
		k3vy = h * ay( p0x + l2x * 0.5, p0y + l2y * 0.5 )
		l3x = h * (v0x + k2vx * 0.5)
		l3y = h * (v0y + k2vy * 0.5)

		k4vx = h * ax( p0x + l3x , p0y + l3y )
		k4vy = h * ay( p0x + l3x , p0y + l3y )
		l4x = h * (v0x + k3vx )
		l4y = h * (v0y + k3vy )


		p0x = p0x + (l1x + 2.0 * (l2x + l3x) + l4x) / 6.0
		p0y = p0y + (l1y + 2.0 * (l2y + l3y) + l4y) / 6.0

		v0x = v0x + (k1vx + 2.0 * (k2vx + k3vx) + k4vx) / 6.0
		v0y = v0y + (k1vy + 2.0 * (k2vy + k3vy) + k4vy) / 6.0

		r.append(p0x)
		r.append(p0y)
		vx.append(v0x)
		vy.append(v0y)
		
	return r

# returns the x coordinate from RK4 array. 
def get_x(body):
	x = []

	results = rk4(body.x, body.y, body.vx, body.vy, body.orbit)
	
	even = np.arange(0,len(results),2)

	for i in even:
		x.append(results[i])

	return x
	
# returns the y coordinate from RK4 array. 
def get_y(body):

	y = []
	
	results = rk4(body.x, body.y, body.vx, body.vy, body.orbit)
	
	odd = np.arange(1,len(results),2)

	for i in odd:
		y.append(results[i])

	return y
	

plt.plot(get_x(earth), get_y(earth) )
plt.plot(get_x(mercury), get_y(mercury) )
plt.plot(get_x(venus), get_y(venus) )
plt.plot(get_x(mars), get_y(mars) )
# plt.plot(get_x(jupiter), get_y(jupiter) )
# plt.plot(get_x(saturn), get_y(saturn) )
# plt.plot(get_x(uranus), get_y(uranus) )
# plt.plot(get_x(neptune), get_y(neptune) )
plt.show()


	
	
	
	
	
