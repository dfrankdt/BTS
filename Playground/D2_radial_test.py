#!/usr/bin/env python3
"""
Second difference operator leveraging radial symmetry

I'd like to be sure that the second difference operator in polar coordinates
with radial symmetry is correct.  So I'll construct such an operator and test
it on the function

 u(r, t) = exp(-r^2)
 
for which the second derivative in space satisfies

  nabla^2 u = 1/r d/dr ( r du/dr ) = 4(r^2 - 1) * exp(-r^2)

"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as manimation


# =============================================================================
# Second Difference Operator
# =============================================================================
def D2_polar(r, u):
	dr = r[1] - r[0]
	Nr = len(r) - 1
	
	# --- Interior discretizations
	rp = (r[2:Nr+1] + r[1:Nr])/2
	rm = (r[1:Nr] + r[0:Nr-1])/2

	# --- Second derivative operator
	D2 = -2*np.eye(Nr-1)
	D2p = rp[0:Nr-2]/r[1:Nr-1]
	D2m = rm[1:Nr-1]/r[2:Nr]
	D2 = D2 + np.diag( D2p, k=1) + np.diag(D2m, k=-1)

	# --- Adjust at boundaries
#	D2[0, 0] = D2[0, 0] + rm[0]/r[1]
	D2u = 1/dr**2 * (D2@u[1:Nr])
	D2u[-1] = D2u[-1] + 1/dr**2 * (rp[-1]/r[-2]*u[-1])
	D2u[0] = D2u[0] + 1/dr**2 * (rm[0]/r[1]*u[0])
	print(np.linalg.cond(D2))
	return D2u

# =============================================================================
# Main Simulation Function
# =============================================================================
def D2_radial_test():
	# --- Parameters
	R = 2
	Nr = 2**9
	r = np.linspace(0, R, Nr+1)
	
	u = np.exp(-r**2)
	upp = D2_polar(r, u)

	fig, (ax1, ax2) = plt.subplots(1, 2, figsize= (12.8, 4.8) )
	ax1.plot(r[1:Nr], upp)
	ax1.plot(r[1:Nr], 4*(r[1:Nr]**2 - 1)*u[1:Nr], '--r')
	ax2.plot(r[1:Nr], upp - 4*(r[1:Nr]**2 - 1)*u[1:Nr])
	plt.show()

# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    D2_radial_test()
