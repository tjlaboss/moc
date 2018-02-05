# Track Generator
#
# Class for the track laydown

import cell
import ray
from functions import *
import pylab
from warnings import warn

DEBUG = False

class TrackGenerator(object):
	"""Thing that lays down the tracks over some cell geometry
	
	Parameters:
	-----------
	cell:           cell.Cell; asdf
	nazim:          int; number of azimuthal angle pairs
					(equivalent to _num_azim_2)
	dtarget:        float; target distance between tracks (cm)
					Will be modified in order to make tracks cyclic
	
	Attributes:
	-----------
	many, to be described
	"""
	def __init__(self, cell_, nazim, dtarget):
		self.cell = cell_
		self.nazim = nazim
		
		n = nazim//2
		
		self.target_angles = pylab.zeros(n)
		self.dazim = pylab.zeros(n)
		
		self.phis = pylab.zeros((2, n))
		
		self.nxs = pylab.zeros(n, dtype = int)
		self.nys = pylab.zeros(n, dtype = int)
		self.integral_nxs = pylab.zeros(n, dtype = int)
		self.integral_nys = pylab.zeros(n, dtype = int)
		self.dxs = pylab.zeros(n)
		self.dys = pylab.zeros(n)
		
		
		delta_x = cell_.xmax - cell_.xmin
		delta_y = cell_.ymax - cell_.ymin
		for i in range(n):
			# New method, from OpenMOC docs
			ang = math.pi*(i + 0.5)/nazim
			self.target_angles[i] = ang
			self.nxs[i] = math.fabs(delta_x*math.sin(ang)/dtarget) + 1
			self.nys[i] = math.fabs(delta_y*math.cos(ang)/dtarget) + 1
			phi_eff = math.atan2(delta_y*self.nxs[i], delta_x*self.nys[i])
			self.phis[0, i] = phi_eff
			self.phis[1, i] = math.pi - phi_eff
			self.dxs[i] = delta_x/self.nxs[i]
			self.dys[i] = delta_y/self.nys[i]
			self.dazim[i] = self.dxs[i]*math.sin(phi_eff)
		
		#print("Angles:", [round(deg(ang)) for ang in self.phis])
		self._tracks = []  # TODO: Initialize as list of the right length
		# Integrate nxs, nys
		for i in range(1, n):
			self.integral_nxs[i] = self.nxs[:i].sum()
			self.integral_nys[i] = self.nys[:i].sum()
		# Dictionaries of generated tracks using coordinates as keys:
		self.track_starts = dict()  # key = (x0, y0)
		self.track_stops = dict()   # key = (x1, y1)
		self.track_phis = dict()    # key = phi
		#
		self.__flux_index = -1
	
	def _increment(self):
		self.__flux_index += 1
		return self.__flux_index
	
	
	def _ij(self, i, j):
		"""Given the x index and the y index,
		return a single index of the two in the same vector.

		Parameters:
		-----------
		i:      int; x index
		j:      int; y index

		Returns:
		--------
		k:      int; index of x and y
		"""
		k = i + j*2*self._nx
		return k
	
	
	def _assign_index(self, x, y, a, i = None, j = None, dx = None, dy = None):
		"""Given the coordinates of a track endpoint, assign the index
		in the angular flux matrix
		
		If you are on the x-edge, you must specify the index j,
		or it can be calculated from the spacing dy.
		
		If you are on the y-edge, you must specify the index i,
		or it can be calculated from the spacing dx.
		 
		Required Parameters:
		--------------------
		x:          float, cm; x-coordinate of endpoint
		y:          float; cm; y-coordinate of endpoint
		
		Optional Parameters:
		--------------------
		i:          int; x-index on the y-edge
		j:          int; y-index on the x-edge
		dx:         float, cm; x-spacing on the y-edge,
					if the index i is not given
		dy:         float, cm; y-xpacing on the x-edge,
					if the index j is not given
		
		Returns:
		--------
		k:          int; index in the psi vector
		"""
		stepx = self.integral_nxs[a]
		stepy = self.integral_nys[a]
		if j is None and dy is not None:
			j = self._jndex(y, dy)
		if i is None and dx is not None:
			i = self._index(x, dx)
		if x == self.cell.xmin:
			return 0*self.nazim + stepy + j
		elif x == self.cell.xmax:
			return 1*self.nazim + stepy + j
		elif y == self.cell.ymin:
			return 2*self.nazim + stepx + i
		elif y == self.cell.ymax:
			return 3*self.nazim + stepx + i
		else:
			raise ValueError("Invalid coordinates: ({}, {})".format(x, y))
	
	def _jndex(self, y, dy):
		j = (self.cell.ymax - y)/dy  # - 0.5
		return math.floor(j)
	
	def _index(self, x, dx):
		i = (self.cell.xmax + x)/dx  # - 0.5
		return math.floor(i)
	
	def _track_by_quadrant(self, x, y, phi, phindex):
		"""Temporary function used for debugging ray spacing
		
		This used to give different parameters to rays, but now
		it has been simplified. At this point, it's just a wrapper
		for Ray() creation and plotting.
		"""
		track = ray.Ray(x, y, phi, phindex, self.cell)
		
		# Plotting, dbug
		if track is not None:
			#phideg = int(180*phi/math.pi)
			#print("({}, {}) to ({:.2}, {:.2}) @ {} deg".format(track.x0, track.y0, track.x1, track.y1, phideg))
			
			# TODO: Turn this off. Plots should be done as part of track.trace()
			dists = track.trace()
			
			self.cell.plot_track(track, dists)
			
			if DEBUG:
				pylab.show()
				self.cell.figure, self.cell.axis = self.cell._set_plot()
		
		return track
			
	
	def _add_track_to_dict(self, track):
		"""Use a track's starting and stopping coordinates as keys
		in this TrackGenerator's dictionaries for easy lookup.
		
		Parameter:
		----------
		track:      Ray; track to add to the dictionaries
		"""
		self._tracks.append(track)
		if track.index0 not in self.track_starts:
			self.track_starts[track.index0] = {}
			self.track_starts[track.index1] = {}
		self.track_starts[track.index0][track.index1] = track
		self.track_starts[track.index1][track.index0] = track
		if track.index1 not in self.track_stops:
			self.track_stops[track.index1] = {}
		self.track_stops[track.index1][track.index0] = track
		
		if track.phindex not in self.track_phis:
			self.track_phis[track.phindex] = []
		self.track_phis[track.phindex].append(track)
	
	
	def __link_tracks(self):
		"""Form cyclic tracks"""
		for track in self._tracks:
			starts = self.track_starts[track.index1]
			for i0 in starts:
				if i0 != track.index0:
					other_track = starts[i0]
			if track.index1 == other_track.index1:
				other_track.swap_direction()
			#print(track, other_track, "\n\n")
			track.next_track = other_track
			other_track.last_track = track
			if track not in self.track_phis[track.phindex]:
				self.track_phis[track.phindex].append(track)
			
	def _cycle_track(self, track0, a, b=0):
		"""Lay cyclic tracks from some starting point
		
		Parameters:
		-----------
		track0:         ray.Ray; track to start from
		a:              int; index of azimuthal angle
						TODO: Figure out if this is always `track0.phindex`
		b:              int; index (0 or 1) of angular direction.
						If 0, it is an angle in the first quadrant.
						If 1, it is the complementary angle in quadrant 2.
						[Default: 0]
		
		Returns:
		--------
		count:          int; the number of tracks laid, including track0
		yleft:          list of y-values hit on the left edge
		"""
		old_track = track0
		d = 1
		count = 1
		yleft = set()
		yright = set()
		while True:
			x = old_track.x1
			y = old_track.y1
			if pylab.isclose(x, self.cell.xmin):
				yleft.add(round(y, 5))
			elif pylab.isclose(x, self.cell.xmax):
				yright.add(round(y, 5))
			if pylab.isclose(x, track0.x0) and pylab.isclose(y, track0.y0):
				# then we're back where we started
				break
			# Reflection logic
			if y == self.cell.ymax:
				d = -1
			elif y == self.cell.ymin:
				d = +1
			else:
				b = 1*(not b)
			phi = self.phis[b, a]*d
			#print("x0: {:.2f}\ty0: {:.2f}\tphi: {}".format(x, y, round(deg(phi))))
			track = self._track_by_quadrant(x, y, phi, a)
			track.last_track = old_track
			old_track.next_track = track
			track.index0 = self._increment()
			old_track.index1 = track.index0
			self._add_track_to_dict(track)
			old_track = track
			count += 1
		# link the first and last tracks
		track.index1 = track0.index0
		track.next_track = track0
		track0.last_track = track
		return count, yleft, yright
	
	def generate(self):
		"""Actually live up to its name and make the tracks"""
		print("Generating tracks...")
		for a in range(self.nazim//2):
			dy = self.dys[a]
			ny = self.nys[a]
			
			count0 = 0
			yused_xmin = set()
			yused_xmax = set()
			nxy = 2*(self.nxs[a] + self.nys[a])
			
			for n in range(ny):
				for x00, yused in zip((self.cell.xmin, self.cell.xmax), (yused_xmin, yused_xmax)):
					y00 = self.cell.ymax - (0.5 + n)*dy
					if not is_used(y00, yused) and count0 < nxy and x00==self.cell.xmin:
						b = 0
						phi = self.phis[b, a]
						if True:
							print("X00: {:.2f}\tY00: {:.2f}\tphi: {}".format(x00, y00, round(deg(phi))))
						track00 = self._track_by_quadrant(x00, y00, phi, a)
						track00.index0 = self._increment()
						self._add_track_to_dict(track00)
						count1, ynewmin, ynewmax = self._cycle_track(track00, a, b)
						count0 += count1
						yused_xmin.update(ynewmin)
						yused_xmax.update(ynewmax)
				
			# Sometimes I get exactly 2 two many tracks on one angle.
			# I don't know why.
			if count0 > nxy:
				warn("Possibly too many tracks ({} out of {})".format(count0, nxy))
			elif count0 < nxy:
				warn("Possibly too few tracks ({} out of {})".format(count0, nxy))
			
		print("...done.\n")

					
if __name__ == "__main__":
	# test
	t = TrackGenerator(cell.test_cell, 120, dtarget = .1)
	t.generate()
	if not DEBUG:
		pylab.show()
	for trac in t._tracks:
		print(trac.index0, "\t", trac.index1, " \t ", trac.key0, trac.key1)
