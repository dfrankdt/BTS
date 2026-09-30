#!/usr/bin/env python3
"""
Agent Based Chemotaxis: Run and tumble where probability of switching to 
tumbling (koff) may be a function of an external gradient.

Produces
 - Figure 1: Position of Np particles after Nt timesteps in case of zero (blue)
   and active (red) chemotaxis, see Figure 11.1(a)
 - Figure 2: Variance as a function of time, illustrating expected effective
   diffusion (blue) and saturated variance due to chemotaxis (red), see 
   Figure 11.1(b)
"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as manimation
rng = np.random.default_rng()

# =============================================================================
# External Chemical Concentration
# =============================================================================
def c(x, y):
	Lx = 5
	Ly = 3
	z = -( (x - Lx)**2 + (y - Ly)**2 )
	return z

# =============================================================================
# Statistics
# =============================================================================
def doMSD(x, y):
	"""
	Do the mean-square-distance calculation on each timestep
	"""
	Nt = np.size(x, 0) - 1
	Np = np.size(x, 1)
	msd = np.zeros(Nt+1)
	for kt in range(Nt+1):
		xm = sum(x[kt,:])/Np
		ym = sum(y[kt,:])/Np
		msd[kt] = sum( (x[kt,:] - xm)**2 + (y[kt,:] - ym)**2 )/Np
	return msd

# =============================================================================
# Run and Tumble
# =============================================================================
def xyTrajectory(v, kon, k0, chi, dt, Nt, Np):
	"""
	Pass Np particles through a run and tumble process.  All particles start
	moving (state 1).
	"""
	x = np.zeros( (Nt+1, Np) )
	y = np.zeros( (Nt+1, Np) )
	s = np.ones(Np)
	theta = rng.uniform(0, 2*np.pi, Np)
	for kt in range(Nt):
		# --- Move
		dx = v*np.cos(theta)*s*dt
		dy = v*np.sin(theta)*s*dt
		
		# --- Concentration Gradient
		dc = c(x[kt,:] + dx, y[kt,:] + dy) - c(x[kt,:], y[kt,:])
		koff = k0 - chi * dc/dt
		
		# --- Determine whether particle switches direction
		R = rng.uniform(0, 1, Np)
		pswitch = (s==0)*(kon*dt > R) + (s==1)*(koff*dt > R)
		ds = ((s==0)*(1) + (s==1)*(-1))*pswitch
		dtheta = (s==0)*(pswitch) * rng.uniform(0, 2*np.pi, Np)
		
		# --- Update everything		
		theta = np.mod(theta + dtheta, 2*np.pi)
		s = s + ds
		x[kt+1, :] = x[kt, :] + dx
		y[kt+1, :] = y[kt, :] + dy
		
	return x, y

# =============================================================================
# Main Simulation Function
# =============================================================================
def chemotax_agents():
	# --- Parameters
	Np = 1000
	kon = 10
	k0 = 1
	v = 1
	Tmax, Nt = 50, 1000
	dt = Tmax/Nt
	
	# --- Simulations
	chi = 0
	x, y = xyTrajectory(v, kon, k0, chi, dt, Nt, Np)
	
	chi = 0.1
	xc, yc = xyTrajectory(v, kon, k0, chi, dt, Nt, Np)
	
	# --- Plotting
	fig, ax = plt.subplots()
	ax.plot(x[-1,:], y[-1,:], 'b.')
	ax.plot(xc[-1, :], yc[-1, :], 'r.')
	ax.set(title = f'Final Position of {Np:4d} particles after {Nt:4d} timesteps')
	ax.set(xlabel = 'x', ylabel = 'y')
	
	# --- Statistics
	msd = doMSD(x, y)
	msdc = doMSD(xc, yc)
	
	fig2, ax2 = plt.subplots()
	ax2.plot(np.linspace(0, Tmax, Nt+1), msd, 'b')
	ax2.plot(np.linspace(0, Tmax, Nt+1), msdc, 'r')
	ax2.set(title = f'Drift from mean')
	ax2.set(xlabel = 't', ylabel = 'Variance')
	plt.show()
	


# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    chemotax_agents()
