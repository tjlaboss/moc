# Quadrature
#
# Quadrature classes: Tabuchi-Yamamoto, Uniform Distributed, and Leonard.

import numpy


class Quadrature(object):
	"""Base class for a generic quadrature

	Parameters:
	-----------
	phis:       array of azimuthal angles (from TrackGenerator)
	np:         int; number of polar angles to use
	name:       str; name of this quadrature model

	Attributes:
	-----------
	na:         number of azimuthal angles and weights
	wa:         array of azimuthal weights
	thetas:     array of polar angles
	wp:         array of polar weights
	"""
	def __init__(self, phis, np, name = "Generic"):
		self.phis = phis
		self.np = np
		self.name = name
		# Azimuthal
		self.na = len(phis)
		self.wa = numpy.zeros(self.na)
		# Polar
		self.wp = numpy.zeros(np)
		self.sinthetas = numpy.zeros(np)
	
	def __str__(self):
		return self.name + " Quadrature"


class AreaQuadrature(Quadrature):
	"""A very simple quadrature used for testing area calculations"""
	def __init__(self, phis, np, name="Area"):
		super().__init__(phis, np, name)
		self.wa = numpy.ones(self.na)
		self.wa /= (4*self.na)
		self.wp = numpy.ones(self.np)/self.np
		self.sinthetas = numpy.ones(self.np)


class AzimuthalQuadrature(Quadrature):
	"""Base class providing a good azimuthal angle quadrature,
	decoupled from the polar angle quadrature.

	Adapted from the OpenMOC azimuthal quadrature; see
	https://mit-crpg.github.io/OpenMOC/methods/track_generation.html#azimuthal-angle-quadrature

	Parameters:
	-----------
	phis:       array of azimuthal angles (from TrackGenerator)
	np:         int; number of polar angles to use
	name:       str; name of this quadrature model

	Attributes:
	-----------
	na:         number of azimuthal angles and weights
	wa:         array of azimuthal weights
	thetas:     array of polar angles
	wp:         array of polar weights
	"""
	def __init__(self, phis, np, name = "Generic Azimuthal"):
		super().__init__(phis, np, name)
		# Set the azimuthal angles and weights
		if self.na == 1:
			self.wa = numpy.ones(1)/8
		else:
			c = 1 / (2 * numpy.pi)
			self.wa[0] = c*(phis[1] + phis[0]) / 2
			self.wa[-1] = c*(numpy.pi - phis[-1] - phis[-2]) / 2
			for m in range(1, self.na - 1):
				self.wa[m] = c*(phis[m+1] - phis[m-1]) / 2


class UniformDistributedQuadrature(AzimuthalQuadrature):
	"""An equal angle, equal weight polar quadrature

	Parameters:
	-----------
	phis:       array of azimuthal angles (from TrackGenerator)
	np:         int; number of polar angles to use
	name:       str; name of this quadrature model

	Attributes:
	-----------
	na:         number of azimuthal angles and weights
	wa:         array of azimuthal weights
	thetas:     array of polar angles
	wp:         array of polar weights
	"""
	def __init__(self, phis, np, name = "Equal Angle"):
		super().__init__(phis, np, name)
		# Set the polar angles and weights
		dtheta = numpy.pi/(2*np)
		theta_bar_0 = 0
		for i in range(np):
			theta_bar_1 = theta_bar_0 + dtheta
			u0 = numpy.cos(theta_bar_0)
			u1 = numpy.cos(theta_bar_1)
			sintheta = numpy.sin(numpy.arccos((u0 + u1)/2))
			self.sinthetas[i] = sintheta
			self.wp[i] = (u0 - u1)*sintheta
			theta_bar_0 = theta_bar_1


class LeonardQudadrature(AzimuthalQuadrature):
	"""Leonard's optimum quadrature set"""
	def __init__(self, phis, np, name = "Leonard"):
		assert 0 < np <= 3, \
			"LO Quadrature is only available for 1-3 polar angles."
		super().__init__(phis, np, name)
		# Actual leonard values
		if np == 1:
			self.wp = numpy.array([1.0])
			self.sinthetas = numpy.array([0.752244])
		elif np == 2:
			self.wp = numpy.array([0.139473, 0.860527])
			self.sinthetas = numpy.array([0.273658, 0.865714])
		else:
			self.wp = numpy.array([0.020530, 0.219161, 0.760309])
			self.sinthetas = numpy.array([0.103840, 0.430723, 0.905435])
		self.wp *= self.sinthetas


class YamamotoQuadrature(AzimuthalQuadrature):
	"""A Tabuchi-Yamamoto quadrature set"""
	def __init__(self, phis, np, name = "Tabuchi-Yamamoto"):
		assert 0 < np <= 3, \
			"TY Quadrature is only available for 1-3 polar angles."
		super().__init__(phis, np, name)
		# Actual tabuchi-yamamoto values
		if np == 1:
			self.wp = numpy.array([1.0])
			self.sinthetas = numpy.array([0.798184])
		elif np == 2:
			self.wp = numpy.array([0.212854, 0.787146])
			self.sinthetas = numpy.array([0.363900, 0.899900])
		else:
			self.wp = numpy.array([0.046233, 0.283619, 0.670148])
			self.sinthetas = numpy.array([0.166648, 0.537707, 0.932954])
		self.wp *= self.sinthetas



if __name__ == "__main__":
	# test
	from cell import test_cell
	from trackgenerator import TrackGenerator
	
	t = TrackGenerator(test_cell, 8*2, dtarget = .01)
	qua1 = YamamotoQuadrature(t.phis[0, :], 3)
	print("\nThetas:")
	print([numpy.arcsin(theta) for theta in qua1.sinthetas])
	print("Polar weights:")
	print(qua1.wp, qua1.wp.sum())
	print("Azimuthal weights:")
	print(qua1.wa, qua1.wa.sum())