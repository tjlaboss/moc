# Cell
#
# Class and methods for a fuel pin cell in MOC

from functions import *
import fsr
import pylab
import random

# Parameters for the Cell in the assignment
PITCH = 1.26                        # cm; pin pitch
RADIUS = 0.4                        # cm; fuel pin radius
AREA_FUEL = pylab.pi*RADIUS**2		# cm^2; fuel pin area
AREA_MOD = PITCH**2 - AREA_FUEL		# cm^2; moderator area
AREA_RATIO = AREA_MOD/AREA_FUEL
N238 = 2.2E-2                       # 10^24 atoms/cm^3; number density of U238
SIGMA_P238 = 11.4                   # b; Potential scatter xs of U238
SIGMA_PO = 4.0                      # b; Potential scatter xs of O16
SIGMA_P = SIGMA_P238 + 2*SIGMA_PO   # b; Potential scatter xs of fuel
SIGMA_NF = N238*SIGMA_P             # cm^-1; fuel potential scatter xs
SIGMA_AS = [1E-6, .25, 1.0, 5.0, 1E6]  # cm^-1; moderator absorption xs
MODS = [fsr.FlatSourceRegion(AREA_FUEL, a, 0.0, "Mod") for a in SIGMA_AS]
FUEL = fsr.FlatSourceRegion(AREA_FUEL, SIGMA_NF, 1.0, "Fuel")
MOD = MODS[1]
BOUNDARY_CONDITIONS = {"reflective", "periodic", "vacuum"}

class Cell(object):
	"""A resonant pin cell surrounded by moderator
	
	Parameters:
	-----------
	pitch:          float; lattice pitch (cm)
	radius:         float; pin radius (cm)
	fuel:			fsr.FlatSourceRegion; FSR for the fuel
	sigma_y_mod:    fsr.FlatSourceRegion; FSR for the moderator
	boundary:       str; one of {"reflective", "periodic", "vacuum"}
	plot:           Boolean; whether to produce a plot of the model
					[Default: True]
	
	"""
	def __init__(self, pitch, radius, fuel, mod,
				 boundary = "reflective", plot = True):
		self.pitch = pitch
		self.radius = radius
		boundary = boundary.lower()
		assert boundary in BOUNDARY_CONDITIONS, \
			"Boundary must be reflective, periodic, or vacuum"
		self.boundary = boundary
		# For a square universe...
		self.xmin = -pitch/2.0
		self.xmax = +pitch/2.0
		self.ymin = -pitch/2.0
		self.ymax = +pitch/2.0

		if fuel is not None:
			assert isinstance(fuel, fsr.FlatSourceRegion), \
				"You probably forgot to change SIGMA_NF to fuel_fsr."
			self.fuel = fuel
		if mod is not None:
			assert isinstance(mod,  fsr.FlatSourceRegion), \
				"You probably forgot to change SIGMA_A to mod_fsr."
			self.mod = mod
		
		if plot:
			self.figure, self.axis = self._set_plot()
		else:
			self.figure, self.axis = None, None
		
	def _set_plot(self):
		"""Set up a base plot"""
		fig, ax = pylab.subplots()
		ax.add_artist(pylab.Circle((0, 0), self.radius, color = "y"))
		ax.grid()
		ax.set_xlim([self.xmin, self.xmax])
		ax.set_ylim([self.ymin, self.ymax])
		ax.set_aspect("equal")
		return fig, ax
	
	def plot_track(self, track, s, npoints = 15, col = None):
		"""After tracing a ray, add its track to the plot
		
		
		Parameters:
		-----------
		track:      Ray() to plot
		s: tuple of the following:
			s1:         float; distance from the edge to the cell
			sf:         float; distance through the cell
			s2:         float; distance from the cell to the edge
		col:        str or tuple(len=3); color to use for the plot
					[Default: None; will choose randomly]
		
		Returns:
		--------
		plotted:    Boolean; False if there is not plot, True if there is
		"""
		if self.axis is None:
			return False
		
		try:
			s1, sf, s2 = s
		except TypeError:
			if isinstance(s, float):
				print("Warning: expected tuple (s1, sf, s2); got s1")
				s1 = s
				sf = 0
				s2 = 0
			else:
				raise
		
		if col is None:
			col = (random.random(), random.random(), random.random())
		u = pylab.cos(track.phi)
		v = pylab.sin(track.phi)
		
		# Instantiate these; they'll be updated as we go
		xstop = track.x0
		ystop = track.y0
		
		if s1:
			x1start = track.x0
			xstop = x1start + s1*u
			x1vals = pylab.linspace(x1start, xstop, npoints)
			y1start = track.y0
			ystop = y1start + s1*v
			y1vals = pylab.linspace(y1start, ystop, npoints)
			self.axis.plot(x1vals, y1vals, '-', color = col)
		
		if sf:
			xfstart = xstop
			xfstop = xfstart + sf*u
			xfvals = pylab.linspace(xfstart, xfstop, npoints)
			yfstart = ystop
			yfstop = yfstart + sf*v
			yfvals = pylab.linspace(yfstart, yfstop, npoints)
			self.axis.plot(xfvals, yfvals, ':', color = col)
		
		if s2:
			x2start = track.x1
			x2stop = x2start - s2*u
			x2vals = pylab.linspace(x2start, x2stop, npoints)
			y2start = track.y1
			y2stop = y2start - s2*v
			y2vals = pylab.linspace(y2start, y2stop, npoints)
			self.axis.plot(x2vals, y2vals, '-', color = col)
		
		return True
			

if __name__ == "__main__":
	test_cell = Cell(PITCH, RADIUS, FUEL, MOD)
	import ray
	track = ray.Ray(-.25, -PITCH/2, rad(60), None, test_cell)
	segments = track.trace()
	if test_cell.plot_track(track, segments):
		pylab.show()
	else:
		raise SystemError("Test Failed!")
