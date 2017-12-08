# Quadrature
#
# Quadrature classes: Tabuchi-Yamamoto, and TODO: more to come

import numpy


class Quadrature(object):
	"""Generic quadrature

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


class YamamotoQuadrature(Quadrature):
	"""A Tabuchi-Yamamoto quadrature set"""
	
	def __init__(self, phis, np, name = "Tabuchi-Yamamoto"):
		super().__init__(phis, np, name)
		# Set the azimuthal angles and weights
		c = 1/(2*numpy.pi)
		self.wa[0] = c*(phis[1] + phis[0])/2
		self.wa[-1] = 1 - c*(phis[-2] + phis[-1])/2
		for m in range(1, self.na - 1):
			self.wa[m] = c*(phis[m - 1] + phis[m + 1])/2
		
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
		# TODO: precalculate sin(theta)*wt_theta


if __name__ == "__main__":
	# test
	from cell import test_cell
	from trackgenerator import TrackGenerator
	
	t = TrackGenerator(test_cell, 4, dtarget = .5)
	qua1 = YamamotoQuadrature(t.phis[0, :], 3)
	print("\nThetas:")
	print([numpy.arcsin(theta) for theta in qua1.sinthetas])
	print("Polar weights:")
	print(qua1.wp, qua1.wp.sum())
	print("Azimuthal weights:")
	print(qua1.wa, qua1.wa.sum())