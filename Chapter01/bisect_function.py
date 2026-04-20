#!/usr/bin/env python3
"""
Bisection: Illustration of the method of bisection.

This code is based on bisect_function.m

"""
# =============================================================================
# Target Function
# =============================================================================
def f(u):
	"""
	This is the function we would like to zero. Feel free to replace this function.
	"""
	y = 1 - 2 * u
	return y

# =============================================================================
# Bisection Method
# =============================================================================
def bisect_method(f, a, b, atol):
	# Inputs: 
	#  f - function which we wish to find a root
	#  a, b - interval on which we seek a root
	#  atol - absolute tolerance
	#
	# Outputs: 
	#  u0 - estimate of the root with |u - u0| < atol
	"""
	DO THE WORK HERE
	"""
	# Do the rest of the work
	check = 1
	while check > atol:
		fa = f(a)
		fb = f(b)
		
		u0 = (a + b)/2
		f0 = f(u0)
		
		check = abs(f0)

		ftest = ((f0*fa) > 0)

		a = ftest*u0 + (1 - ftest)*a
		fa = ftest*f0 + (1 - ftest)*fa
		
		b = (1 - ftest)*u0 + ftest*b
		fb = (1 - ftest)*f0 + ftest*fb
		
	return u0

# =============================================================================
# Main Simulation Function
# =============================================================================
def bisect_function():
	a, b = 0, 10
	fa, fb = f(a), f(b)
	# --- Check once to see whether a root is guaranteed on the interval
	if fa*fb>0:
		print('Error in bisection method')
		print('The IVT does not guarantee a root on [a, b]')
		print('Exiting')
		exit()

	atol = 1e-12

	u0 = bisect_method(f, a, b, atol)
	print(u0)


# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    bisect_function()
