# Ray
#
# Class for the Ray (Track)

import warnings
from functions import octant
from pylab import *
from decimal import Decimal

D = 4   # Number of decimal places to round to

def dround(number, numd = D):
	"""Advanced decimal rounding
	
	Parameters:
	-----------
	number:         float; the number you want to round
	numd:           int; the number of decimal places to round to
	"""
	return float(Decimal(str(number)).quantize(Decimal(str(10**-numd))))


class Ray(object):
	"""A particle ray to track

	Parameters:
	-----------
	x0:         float; where the ray starts on the x-axis (cm)
	y0:         float;    "   on the y-axis (cm)
	phi:        float; azimuthal angle (radians)
	pin:        Cell;  fuel pin cell to trace the ray over
	
	Attributes:
	-----------
	x1:         float; where the ray stops on the x-axis (cm)
	y1:         float;    "   on the y-axis (cm)
	length:     float; total length of the ray, from start to finish
	index0;     int;   index in the angular flux vector at (x0, y0),
	                   assigned by the generator
	index1:     int;   index in the angular flux vector at (x1, y1)
	phindex:    int; index of the azimuthal angle phi
	
	deprecated:
		key0:       tuple; key of ray's starting coordinates (rounded)
		key1:       tuple; key of ray's stopping coordinates (rounded)
		octant:     int; octant in which Ray.phi falls
	"""
	
	def __init__(self, x0, y0, phi, phindex, pin):
		assert pin.xmin <= x0 <= pin.xmax, \
			"x0={:.2} is not in the pin pitch!".format(x0)
		assert pin.ymin <= y0 <= pin.ymax, \
			"y0={:.2} is not in the pin pitch!".format(y0)
		self.x0 = x0
		self.y0 = y0
		self.phi = phi
		self.octant = octant(phi)
		self.phindex = phindex
		self.pin = pin
		self.index0 = None
		self.index1 = None
		self.comp = None
		self.x1, self.y1 = self.get_edge_coordinates()
		if self.x1 is not None and self.y1 is not None:
			# print("x1: {:.3}, y1: {:.3}".format(self.x1, self.y1))
			self.length = sqrt((self.x1 - self.x0)**2 +
			                   (self.y1 - self.y0)**2)
		else:
			warnings.warn("Track's end point is the same as its start point; "
			              "length will be 0")
			self.length = 0
		self.key0 = (dround(self.x0, D), dround(self.y0, D))
		self.key1 = (dround(self.x1, D), dround(self.y1, D))
	
	def __str__(self):
		rep = "Ray:"
		rep += "\n{} -> {}".format(self.index0, self.index1)
		rep += "\nkey0: " + str(self.key0)
		rep += "\nkey1: " + str(self.key1)
		return rep
	
	def swap_direction(self):
		"""Reverse the order of the track.
		Will swap index0, index1, key0, key1
		"""
		self.index0, self.index1 = self.index1, self.index0
		self.key0, self.key1     = self.key1, self.key0
		
	
	def get_indices(self, fwd = True):
		"""Get the starting and stopping indices of the ray
		according to the direction it is traveling.
		
		Paramters:
		----------
		fwd:        Boolean; whether the track is going in the
					default direction.
					[Default: True]
		
		Returns:
		--------
		start:      int; index of the flux at track start
		stop:       int; index of the flux at track start
		"""
		assert (self.index0 is not None) and (self.index1 is not None), \
			"Track indices have not been set yet."
		if fwd:
			return self.index0, self.index1
		else:
			return self.index1, self.index0
		
	
	def get_edge_coordinates(self):
		"""Get the coordinates at which the ray hits the edge

		TODO: I replaced some clever-sounding trig with a bunch of
		clunky if logic. Now it works, but it's repetitive and inelegant.
		Now that the tracer works as intended, rewrite this.

		Returns:
		--------
		x1, y1:     x- and y-coordinates of the end of the ray
		"""
		w = abs(tan(self.phi))
		if self.x0 == self.pin.xmin:
			# going to the right
			dx = self.pin.xmax - self.pin.xmin
			if self.octant <= 4:
				# going up
				dy = self.pin.ymax - self.y0
				# check if it hits the x or y edge first
				if dy/w > dx:
					# then it hits the right edge first
					#print("left edge -> right edge")
					x1 = self.pin.xmax
					y1 = self.y0 + dx*w
					return x1, y1
				else:
					# then it hits the top first
					#print("left edge -> top edge")
					x1 = self.x0 + dy/w
					y1 = self.pin.ymax
					return x1, y1
			else:
				# going down
				dy = self.y0 - self.pin.ymin
				if dy/w > dx:
					# then it hits the right edge first
					#print("left edge -> right edge")
					x1 = self.pin.xmax
					y1 = self.y0 - dx*w
					return x1, y1
				else:
					# then it hits the bottom edge first
					#print("left edge -> bottom edge")
					x1 = self.x0 + dy/w
					y1 = self.pin.ymin
					return x1, y1
		
		elif self.x0 == self.pin.xmax:
			# going to the left
			dx = self.pin.xmax - self.pin.xmin
			if self.octant <= 4:
				# going up
				dy = self.pin.ymax - self.y0
				# check if it hits the x or y edge first
				if dy/w > dx:
					# then it hits the right edge first
					#print("right edge -> left edge")
					x1 = self.pin.xmin
					y1 = self.y0 + dx*w
					return x1, y1
				else:
					# then it hits the top first
					#print("right edge -> top edge")
					x1 = self.x0 - dy/w
					y1 = self.pin.ymax
					return x1, y1
			else:
				# going down
				dy = self.y0 - self.pin.ymin
				if dy/w > dx:
					# then it hits the left edge first
					#print("right edge -> left edge")
					x1 = self.pin.xmin
					y1 = self.y0 - dx*w
					return x1, y1
				else:
					# then it hits the bottom edge first
					#print("right edge -> bottom edge")
					x1 = self.x0 - dy/w
					y1 = self.pin.ymin
					return x1, y1
		
		elif self.y0 == self.pin.ymax:
			# going down
			dy = self.pin.ymax - self.pin.ymin
			if 6 >= self.octant > 2:
				# going left
				dx = self.x0 - self.pin.xmin
				# check if it hits the x or y edge first
				if dx*w > dy:
					# then it hits the bottom first
					#print("top edge -> bottom edge")
					x1 = self.x0 - dy/w
					y1 = self.pin.ymin
					return x1, y1
				else:
					# then it hits the left first
					#print("top edge -> left edge")
					x1 = self.pin.xmin
					y1 = self.y0 - dx*w
					return x1, y1
			else:
				# going right
				dx = self.pin.xmax - self.x0
				# check if it hits the x or y edge first
				if dx*w > dy:
					# then it hits the bottom first
					#print("top edge -> bottom edge")
					x1 = self.x0 + dy/w
					y1 = self.pin.ymin
					return x1, y1
				else:
					# then it hits the right first
					#print("top edge -> right edge")
					x1 = self.pin.xmax
					y1 = self.y0 - dx*w
					return x1, y1
		
		elif self.y0 == self.pin.ymin:
			# going up
			dy = self.pin.ymax - self.pin.ymin
			if 6 >= self.octant > 2:
				# going left
				dx = self.x0 - self.pin.xmin
				# check if it hits the x or y edge first
				if dx*w > dy:
					# then it hits the top first
					#print("bottom edge -> top edge")
					x1 = self.x0 - dy/w
					y1 = self.pin.ymax
					return x1, y1
				else:
					# then it hits the left first
					#print("top edge -> left edge")
					x1 = self.pin.xmin
					y1 = self.y0 + dx*w
					return x1, y1
			else:
				# going right
				dx = self.pin.xmax - self.x0
				# check if it hits the x or y edge first
				if dx*w > dy:
					# then it hits the top first
					#print("bottom edge -> top edge")
					x1 = self.x0 + dy/w
					y1 = self.pin.ymax
					return x1, y1
				else:
					# then it hits the right first
					#print("top edge -> right edge")
					x1 = self.pin.xmax
					y1 = self.y0 + dx*w
					return x1, y1
		
		else:
			raise NotImplementedError("Ray is not starting on the cell edge!")
	
	def get_dist_to_collision(self, forward):
		u = cos(self.phi)
		v = sin(self.phi)
		R = self.pin.radius
		
		if forward:
			x = self.x0
			y = self.y0
			sign_ = +1
		else:
			x = self.x1
			y = self.y1
			sign_ = -1
		
		# Test whether the ray goes through the fuel cell at all
		try:
			if R > abs((self.x1 - self.x0)*self.y0 - (self.y1 - self.y0)*self.x0)/ \
					sqrt((self.x1 - self.x0)**2 + (self.y1 - self.y0)**2):
				# Then an intersection expected
				radicand = R**2 + u**2*(x**2 - y**2) + \
				           2.0*u*v*x*y - x**2
				s = abs(u*x + v*y + sign_*sqrt(abs(radicand)))
				return s
			else:
				return forward*self.length
		except ZeroDivisionError as err:
			errstr = "{}; track length may be zero.".format(err)
			warnings.warn(errstr)
			return forward*self.length
	
	def trace(self, forward = True):
		# s1: distance before fuel pin
		s1 = self.get_dist_to_collision(forward)
		# s2: distance after fuel pin
		s2 = self.get_dist_to_collision(not forward)
		# sf: distance in fuel pin
		sf = self.length - s2 - s1
		# print("S1={:.3}; S2={:.3}; Sf={:.3};\tlength={:.3}".format(s1, s2, sf, self.length))
		return s1, sf, s2

