# Functions

import math

def octant(angle):
	"""Determine what octant an angle is in

	Parameter:
	----------
	angle:      float; angle in radians

	Returns:
	--------
	quadrant:   int; the octant {1, 2, 3,...8}
				which the angle inhabits
	"""
	angle %= (2*math.pi)
	if angle == 0:
		return 1
	else:
		quad = math.pi/4.0
		return int(math.ceil(angle/quad))

def rad(degrees):
	return math.pi/180*degrees

def deg(radians):
	return 180.0/math.pi*radians

def l2norm_2d(new_psi, old_psi):
	"""Compare the L2 engineering norms of two 2-dimentional arrays
	
	Parameters:
	-----------
	new_psi:        array of the latest values
	old_psi:        array of the reference values.
					Must be the same shape as new_psi.
	
	Returns:
	--------
	fluxdiff:       float; the L2 norm of the arrays
	"""
	assert new_psi.shape == old_psi.shape, \
		"Matrix dimensions do not match!"
	fluxdiff = 0
	ny, nx = new_psi.shape
	for j in range(ny):
		for i in range(nx):
			fluxdiff += (new_psi[j, i] - old_psi[j, i])**2
	fluxdiff = math.sqrt(fluxdiff/nx*ny)
	return fluxdiff


def is_used(val, used_set):
	"""Check whether a value is used in a set, within tolerance
	
	Parameters:
	-----------
	val:        float; value to check for
	used_set:   Iterable of floats; values to check in.
	
	Returns:
	--------
	Boolean; whether val is in used_set, within tolerance
	"""
	for u in used_set:
		if math.isclose(val, u):
			return True
	return False
