# Notes
# - fix the initial profile
# - the CN should probably have true zero Dirichlet conditions

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as manimation
# -- Nonlinearity
def F(u):
	y = u * (1-u)
	return y

def doCN(u0, x, t):
	dx = x[1]-x[0]
	dt = t[1]-t[0]
	Nx = np.size(x) - 1
	Nt = np.size(t) - 1

	U = np.zeros( (Nx+1, Nt+1) )
	U[:,0] = u0

	# -- Matrices for performing CN
	gam = dt/(2*dx**2)
	D2 = -2*np.eye(Nx+1)
	D2 = D2 + np.eye(Nx+1, k=1) + np.eye(Nx+1, k=-1)

	Acn = np.eye(Nx+1) - gam*D2
	Bcn = np.eye(Nx+1) + gam*D2
	
	Acn[0,:] = np.zeros(Nx+1)
	Acn[0, 0] = 1
	Acn[-1,:] = np.zeros(Nx+1)
	Acn[-1, -1] = 1

	# -- Initialization of CN method
	uk = u0
	for kt in range(Nt):
		y = Bcn@uk + dt*F(uk)
		y[0], y[-1] = 0, 0
		ukp1 = np.linalg.solve(Acn, y)
		uk = ukp1
		U[:,kt+1] = ukp1
	return U

# -- Parameter
Y = 5

# -- Spatial and Temporal Scales
Nt, Nx = 2**6, 2**8
dt, dx = 0.2,  Y/Nx
x = np.linspace(0, Y, Nx+1)
t = np.linspace(0, dt*Nt, Nt+1)

# -- Initial profile and array for state variable
u0 = 0.2*np.exp(- (x-Y/2)**2/(Y/4))*np.sin(np.pi*x/Y)
U = doCN(u0, x, t)

# Initialize movie
fig, ax = plt.subplots()
p_update = ax.plot([], [], 'b', label='Time Evolution')[0]
p_init = ax.plot(x, u0, '--r', label='Initial Profile')
ax.set(ylim=(0,1))
ax.set(xlabel='x', ylabel='u(x, t)')
ax.legend(loc='upper left')

def update(frame):
    tk = t[frame]
    u = U[:, frame]
    p_update.set_xdata(x)
    p_update.set_ydata(u)
    ax.set(title=f'Time t = {tk:.2f} s')
    ax.legend(loc='upper left')
    return(p_update)
        
ani = manimation.FuncAnimation(fig=fig, func=update, frames=range(0, Nt+1), interval=100)
plt.show()



