# Track Generator
#
# Class for the track laydown

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
	ntotal:            int; total number of nodes on x and y edges
	"""
	def __init__(self, cell_, nazim, dtarget):
		self.cell = cell_
		self.nazim = nazim
		if self.cell.axis:
			titext = "{} Azimuthal angles\t$\delta_z = {}$ cm".format(
				2*nazim, dtarget)
			self.cell.axis.set_title(titext, fontsize="14")
		
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
		self.ntotal = 2 * (self.nxs.sum() + self.nys.sum())
		
		#print("Angles:", [round(deg(ang)) for ang in self.phis[0,:]])
		self._tracks = []
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
		self.generated = False

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
	
	def _new_track(self, x, y, phi, phindex):
		"""Wrapper for Ray() creation and plotting.

		Parameters:
		-----------
		x:			float, cm; x0 of the new ray
		y:			float, cm; y0 of the new ray
		phi:		float, radians; azimuthal angle
		phindex:	int; azimuthal index
		"""
		track = ray.Ray(x, y, phi, phindex, self.cell)
		dists = track.trace()
		self.cell.plot_track(track, dists)
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
		while True:
			x = old_track.x1
			y = old_track.y1
			if pylab.isclose(x, track0.x0) and pylab.isclose(y, track0.y0):
				# then we're back where we started
				break
			elif pylab.isclose(x, self.cell.xmin):
				# then this hits a new y-ordinate
				yleft.add(round(y, 5))
			# Reflection logic
			if y == self.cell.ymax:
				d = -1
			elif y == self.cell.ymin:
				d = +1
			else:
				b = 1*(not b)
			phi = self.phis[b, a]*d
			#print("x0: {:.2f}\ty0: {:.2f}\tphi: {}".format(x, y, round(deg(phi))))
			track = self._new_track(x, y, phi, a)
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
		return count, yleft

	def generate(self, plot_each_cycle=False):
		"""Actually live up to its name and make the tracks

		Algorithm:
		----------
		We know the number of angles per edge and the required spacing
		from the correction calculations in __init__(). Therefore, we know
		where all the endpoints of the tracks need to be. By following
		reflection logic, we can start a track at one of the required
		locations and follow it until it returns to its starting coordinates.
		For some special angles	(e.g., 90 degrees), the cycle will terminate
		before all required sites have been hit. Then, start a new cycle
		at the next unused site.

		For most problems, this function is responsible for most of the runtime.

		Parameter:
		----------
		plot_each_cycle:		Boolean, optional; whether to plot the complete
								cycle of tracks for each azimuthal angle.
								Useful for debugging.
								[Default: False]
		"""
		print("Generating tracks...")
		for a in range(self.nazim // 2):
			dy = self.dys[a]
			ny = self.nys[a]

			count0 = 0
			yused_xmin = set()
			nxy = 2 * (self.nxs[a] + self.nys[a])
			x00 = self.cell.xmin

			for n in range(ny):
				y00 = self.cell.ymax - (0.5 + n) * dy
				if not is_used(y00, yused_xmin) and count0 < nxy:
					b = 0
					phi = self.phis[b, a]
					track00 = self._new_track(x00, y00, phi, a)
					track00.index0 = self._increment()
					self._add_track_to_dict(track00)
					count1, ynew = self._cycle_track(track00, a, b)
					count0 += count1
					yused_xmin.update(ynew)
			if count0 != nxy:
				warn("Wrong number of tracks ({} out of {})".format(count0, nxy))
			if self.cell.figure and plot_each_cycle:
				# Plot the full cycle for this azimuthal angle
				pylab.show()
				self.cell.figure, self.cell.axis = self.cell._set_plot()

		print("...done.\n")
		self.generated = True

					
if __name__ == "__main__":
	# test
	import cell as c
	test_cell = c.Cell(c.PITCH, c.RADIUS, c.SIGMA_NF, c.SIGMA_A)
	t = TrackGenerator(test_cell, 16, dtarget = .25)
	t.generate()
	if not DEBUG:
		pylab.show()
	for trac in t._tracks:
		print(trac.index0, "\t", trac.index1, " \t ", trac.key0, trac.key1)
