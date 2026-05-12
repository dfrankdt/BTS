#!/usr/bin/env python3
"""
Plots: We create the plots shown in Figures 1.1 - 1.5 of the text.

"""

# --- To Do: Figure out how to display reasonably sized arrows in quiver

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# =============================================================================
# Figure 1.1
# =============================================================================
def figure1_1():
	"""
	We create Figure 1.1, illustrating du/dt vs f(u) for f(u) = au (u - 1)(alpha - u)
	for a = 10 and alpha = 0.25.
	"""
	# --- Set parameters, identify function
	a, alpha = 10, 0.25
	u = np.linspace(-0.2, 1.2, 2**9 + 1)
	f = a*u * (u - 1) * (alpha - u)
	f0 = np.zeros(u.size)

	# --- Initialize the plot	
	fig1, ax = plt.subplots()
	p_ep = ax.plot([0, alpha, 1], [0, 0, 0], 'ok')
	p_ax = ax.plot(u, f0, 'r--')
	p_fn = ax.plot(u, f, 'b')

	# --- Add arrows
	ax.annotate("", xytext=(0.25, 0.5), xy = (0.05, 0.5),
					arrowprops = dict(arrowstyle="->"))
	ax.annotate("", xytext=(0.65, 0.5), xy = (0.85, 0.5),
					arrowprops = dict(arrowstyle="->"))
	
	# --- Personalize the labels, axis ranges
	ax.set(ylim = (-1, 1))
	ax.set(xlabel = 'u', ylabel = 'du/dt')
	plt.show()

# =============================================================================
# Figure 1.2
# =============================================================================
def figure1_2():
	"""
	We create Figure 1.2, illustrating the 'trick' of plotting F(u), an integral
	of f(u) as a function of u, then switching the axis to create the desired graph.
	"""
	
	# --- Set parameters, identify function
	a, al = 10, 0.25
	u = np.linspace(-0.2, 1.2, 2**9 + 1)
	F = (al * np.log(np.abs(1-u)) 
		- np.log(np.abs(u-al)) 
		+ (1-al)*np.log(np.abs(u))
		)/(al*(al-1)*a)
	F0 = np.zeros(u.size)
	Fal = al*np.ones(u.size)
	F1 = np.ones(u.size)

	# --- Initialize the plot shown in subfigure (a)
	fig2a, ax2a = plt.subplots()
	ax2a.plot(u, F)
	ax2a.plot(u, F0, '--r')

	ax2a.set(xlim = (-0.2, 1.2), ylim = (-3, 2))
	ax2a.set(xlabel = 'u', ylabel = 'F(u)')
	
	# --- Initialize the plot shown in subfigure (b)
	fig2b, ax2b = plt.subplots()
	ax2b.plot(F, u)
	ax2b.plot([-3, 2], [0, 0], '--')
	ax2b.plot([-3, 2], [al, al], '--')
	ax2b.plot([-3, 2], [1, 1], '--')

	ax2b.set(xlim = (-3, 2), ylim = (-0.2, 1.2))
	ax2b.set(xlabel = 't', ylabel = 'u(t)')
	
	plt.show()

# =============================================================================
# Figure 1.3
# =============================================================================
def figure1_3():
	"""
	We create Figure 1.3, illustrating the phase portrait for the DE
	
	  u'' + f(u) = 0
	  
	with f(u) as given in the above figures.
	"""
	
	# --- Set parameters, identify function
	a, al = 10, 0.25
	u = np.linspace(-0.2, 1.2, 2**9 + 1)
	f = a*u*(1 - u)*(u - al)

	# --- Initialize the plot shown in subfigure (a)	
	fig3a, ax3a = plt.subplots()

	ax3a.set(xlim = (-0.2, 1.2), ylim = (-1, 1))
	ax3a.set(xlabel = 'u', ylabel = 'v')
	ax3a.plot([0, al, 1], [0, 0, 0], 'og')
	ax3a.plot([0, 0], [-1, 1], '--b')
	ax3a.plot([al, al], [-1, 1], '--b')
	ax3a.plot([1, 1], [-1, 1], '--b')
	ax3a.plot([-0.2, 1.2], [0, 0], '--r')
	
	# --- Add arrows, leveraging the vector field
	u = [0, al, 1]
	v = [-0.5, 0.5]
	U, V = np.meshgrid(u, v)
	dU = V
	dV = -a*U*(1 - U)*(U - al)
	ax3a.quiver(U, V, dU, dV, 
				pivot='mid',
				scale_units='inches', scale=4, width=0.004,
				headlength=2, headaxislength=2
				)

	u = [-0.1, al/2, (al + 1)/2, 1.1]
	v = [0]
	U, V = np.meshgrid(u, v)
	dU = V
	dV = -a*U*(1 - U)*(U - al)
	ax3a.quiver(U, V, dU, dV, 
				pivot='mid',
				scale_units='inches', scale=4, width=0.004,
				headlength=2, headaxislength=2
				)

	u = [al/2, (al + 1)/2]
	v = [-0.25, 0.25]
	U, V = np.meshgrid(u, v)
	dU = V
	dV = -a*U*(1 - U)*(U - al)
	ax3a.quiver(U, V, dU, dV, 
				pivot='mid',
				scale_units='inches', scale=4, width=0.004,
				headlength=2, headaxislength=2
				)

	# --- Initialize the plot shown in subfigure (b)
	fig3b, ax3b = plt.subplots()

	ax3b.set(xlim = (-0.2, 1.2), ylim = (-1, 1))
	ax3b.set(xlabel = 'u', ylabel = 'v')
	ax3b.plot([0, al, 1], [0, 0, 0], 'og')
	ax3b.plot([0, 0], [-1, 1], '--b')
	ax3b.plot([al, al], [-1, 1], '--b')
	ax3b.plot([1, 1], [-1, 1], '--b')
	ax3b.plot([-0.2, 1.2], [0, 0], '--r')
	
	# --- Construct the Hamiltonian
	def H(u, v):
		F = -a * (u**4/4 - (1+al)*u**3/3 + al*u**2/2)
		z = 1/2 * v**2 + F
		return z
	
	# --- Identify the level curves
	u_vals = np.array([-0.1, 0, al/1, (al + 1)/2, 1])
	levels = H(u_vals, 0)
	levels = np.sort(levels)
	
	# --- Do the contour plotting
	u = np.linspace(-0.2, 1.2, 2**8+1)
	v = np.linspace(-1, 1, 2**8+1)
	U, V = np.meshgrid(u, v)
	Z = H(U, V)
	ax3b.contour(U, V, Z, levels, colors=['black'], linestyles=['solid'])

	# --- Add a couple arrows
	u = [0.6]
	v = [-0.5, 0.5]
	U, V = np.meshgrid(u, v)
	dU = V
	dV = -a*U*(1 - U)*(U - al)
	ax3b.quiver(U, V, dU, dV, 
				pivot='mid',
				scale_units='inches', scale=4, width=0.004,
				headlength=2, headaxislength=2
				)
	plt.show()

# =============================================================================
# Figure 1.4
# =============================================================================
def figure1_4():
	"""
	We create Figure 1.4, illustrating the phase portrait for a typical saddle
	i.e. for
	
	  u'' + u = 0
	  
	"""
	# --- Initialize the plot
	fig4, ax4 = plt.subplots()

	ax4.set(xlim = (-1, 1), ylim = (-1, 1))
	ax4.set(xlabel = 'u', ylabel = 'v')
	ax4.plot(0, 0, 'ok')
	
	# -- Add the nullclines
	ax4.plot([0, 0], [-1, 1], '--')
	ax4.plot([-1, 1], [0, 0], '--')

	# --- Add arrows, leveraging the vector field
	u = [-0.5, 0.5]
	v = [0]
	U, V = np.meshgrid(u, v)
	dU = V
	dV = U
	ax4.quiver(U, V, dU, dV, 
				pivot='mid',
				scale_units='inches', scale=4, width=0.004,
				headlength=2, headaxislength=2
				)
	u = [0]
	v = [-0.5, 0.5]
	U, V = np.meshgrid(u, v)
	dU = V
	dV = -U
	ax4.quiver(U, V, dU, dV, 
				pivot='mid',
				scale_units='inches', scale=4, width=0.004,
				headlength=2, headaxislength=2
				)

	u = [-0.25, 0.25]
	v = [-0.25, 0.25]
	U, V = np.meshgrid(u, v)
	dU = V
	dV = U
	ax4.quiver(U, V, dU, dV, 
				pivot='mid',
				scale_units='inches', scale=4, width=0.004,
				headlength=2, headaxislength=2
				)

		# --- Construct the Hamiltonian
	def H(u, v):
		z = 1/2 * v**2 - 1/2 * u**2
		return z
	
	# --- Identify the level curves
	levels = np.array([H(0.25, 0), H(0, 0.25)])
	levels = np.sort(levels)
	print(levels)
	
	# --- Do the contour plotting
	u = np.linspace(-1, 1, 2**8+1)
	v = np.linspace(-1, 1, 2**8+1)
	U, V = np.meshgrid(u, v)
	Z = H(U, V)
	ax4.contour(U, V, Z, levels, colors=['black'], linestyles=['solid'])

	plt.show()

# =============================================================================
# Figure 1.5
# =============================================================================
def figure1_5():
	"""
	We create Figure 1.5, illustrating the phase portrait for a typical (a)
	node and (b) spiral i.e. for
	
	  u'' - 1.3 u' + 0.24 u = 0     (a)
	  u'' - 0.2 u' + 1.01 u = 0 	(b)
	  
	"""
	# --- Initialize the plot shown in subfigure (a)
	fig5a, ax5a = plt.subplots()

	ax5a.set(xlim = (-1, 1), ylim = (-1, 1))
	ax5a.set(xlabel = 'u', ylabel = 'v')
	ax5a.plot(0, 0, 'ok')
	
	# --- Add the nullclines
	ax5a.plot([-1, 1], [-5, 5], '--')
	ax5a.plot([-1, 1], [-1, 1], '--')

	# --- Add arrows, leveraging the vector field
	u = [-0.15]
	v = [-.75]
	U, V = np.meshgrid(u, v)
	dU = -U + 0.2*V
	dV = 0.3*U - 0.3*V
	ax5a.quiver(U, V, dU, dV, 
				pivot='mid',
				scale_units='inches', scale=1/2, width=0.004,
				headlength=2, headaxislength=2
				)

	# --- Identify approximate solution trajectories
	def de_rhs(t, x):
		A = np.array([ [-1, 0.2], [0.3, -0.3] ])
		dx = A@x
		return dx
	t_final = 20
	x0_list = np.array([[1, 0], [-1, 0], [0, 1], [0, -1]])
	k_tests = np.size(x0_list, 0)
	for k in np.arange(0, k_tests):
		x0 = x0_list[k,:]
		soln = solve_ivp(de_rhs, [0, t_final], x0, dense_output=True)
		tt = np.linspace(0, t_final, 2**9+1)
		xx = soln.sol(tt).T
		ax5a.plot(xx[:,0], xx[:,1],'k')
	

	# --- Initialize the plot shown in subfigure (b)
	fig5b, ax5b = plt.subplots()

	ax5b.set(xlim = (-1, 1), ylim = (-1, 1))
	ax5b.set(xlabel = 'u', ylabel = 'v')
	ax5b.plot(0, 0, 'ok')
	
	# --- Add the nullclines
	ax5b.plot([-1, 1], [1/10, -1/10], '--')
	ax5b.plot([-1, 1], [-10, 10], '--')

	# --- Identify approximate solution trajectories
	# --- Note: here we integrate backward in time to capture the unstable EP
	def de_rhs(t, x):
		A = np.array([ [0.1, 1], [-1, 0.1] ])
		dx = A@x
		return dx
	t_final = -40
	x0_list = np.array([ [-1, 1/10] ])
	k_tests = np.size(x0_list, 0)
	for k in np.arange(0, k_tests):
		x0 = x0_list[k,:]
		soln = solve_ivp(de_rhs, [0, t_final], x0, dense_output=True)
		tt = np.linspace(0, t_final, 2**9+1)
		xx = soln.sol(tt).T
		ax5b.plot(xx[:,0], xx[:,1],'k')

	plt.show()
	
# =============================================================================
# Main Simulation Function
# =============================================================================
def de_chapt_plots():
	#figure1_1()
	#figure1_2()
	figure1_3()
	#figure1_4()
	#figure1_5()



# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    de_chapt_plots()
