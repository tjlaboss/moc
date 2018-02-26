# FSR
#
# Module for Flat Source Regions

from math import pi

class FlatSourceRegion(object):
	"""FSR for MOC
	One group, two dimensions

	Parameters:
	-----------
	area:		float, cm^2; surface area of the FSR
	xs:			float, cm^-1; total macroscopic xs in the FSR (one-group)
	q:			float, optional  [Default: 0.0]
	name:		str, optional; useful name of FSR

	Attributes:
	-----------
	source:		`q` divided by sigma_t and 4pi
	"""
	def __init__(self, area, xs, q=0.0, name=""):
		self.area = area
		self.xs = xs
		self.name = name
		self.q = q
		if q == 0:
			self.source = 0
		elif q > 0 and xs > 0:
			self.source = q / (4*pi*xs)
		else:
			raise ValueError("Cannot have nonzero source with 0 cross section")

	def __str__(self):
		return self.name

	def report(self):
		"""Get a human-readable report

		Returns:
		--------
		string of the name, area, and source
		"""
		rep = """\
FSR {}
------
Area:	   {:.4f} cm^2
Q/Sigma_t: {:.4f}""".format(self.name, self.area, self.source)
		return rep

