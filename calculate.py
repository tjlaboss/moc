# Calculate
#
# A script to do flux calculations

import pylab
import cell
import quadrature
from functions import l2norm_2d
from trackgenerator import TrackGenerator

# Global constants
PLOT = False
MAXITER = 1000
EPS = 1E-5       # Convergence criterion
N_AZIM_2 = 16    # number of azimuthal angle pairs
D_AZIM = 0.1   # target spacing between parallel tracks (cm)
QFSR = 1.0       # flat fixed source magnitude


class Calculator(object):
	"""Docstring
	
	Parameters:
	-----------
	model:          Cell; model of the cell geometry
	generator:      TrackGenerator; object to lay down the MOC tracks
	quad:           Quadrature; which quadrature model to use
	eps:            float; convergence criterion for flux solution
					[Default: TBD]
	plot:           Boolean; whether to produce a plot during calculation
					[Default: True]
	
	Attributes:
	-----------
	psi:            numpy.array; angular flux, size = (np, ntotal)
						psi[p,:]: flux vector for polar angle p
						psi[p,i]: angular flux at node i for angle p
	ntotal:         int; the number of nodes on the x and y boundaries.
						Serves as the length of psi[p,:]
	modr_flux:      float; scalar flux in the moderator
	fuel_flux:      float; scalar flux in the fuel
	
	deprecated
	-----------
	flux_dict_0:    nested dictionaries:
					{phi : {(x0, y0) : {(x1, y1) : flux } } }
	flux_dict_1:    nested dictionaries:
					{phi : {(x1, y1) : {(x0, y0) : flux } } }
	"""
	def __init__(self, model, generator, quad, source = QFSR, eps = EPS, plot = True):
		self.model = model
		self.generator = generator
		self.quad = quad
		#self.source = source / (4*pylab.pi)
		self.eps = eps
		self.plot = plot
		self.psi = pylab.zeros((self.quad.np, self.generator.ntotal, 2))
		self.modr_flux = 1.0
		self.fuel_flux = 1.0
		# Area calculation
		self.modr_area_tally = None
		self.fuel_area_tally = None

	def transport_sweep(self, calculate_area=False):
		"""Sweep over all angles, tracks, and segments

		Parameters:
		-----------
		q:					float; source in the fuel
		calculate_area:		Boolean, optional; whether to sample the surface
							areas of the FSRs.
							[Default: False]

		Returns:
		--------
		fuel_flux:		float; scalar flux in the fuel
		modr_flux:		float; scalar flux in the moderator
		"""
		assert self.generator.generated, \
			"You must generate tracks before sweeping!"
		# Shorthand
		sig_fuel = self.model.fuel.xs
		sig_mod = self.model.mod.xs
		# Initialize the fluxes in each region
		modr_flux = 0.0
		fuel_flux = self.model.fuel.q / (4 * pylab.pi)
		if calculate_area:
			self.modr_area_tally = 0
			self.fuel_area_tally = 0
		for a in range(self.quad.na):
			wa = self.quad.wa[a]
			wk = self.generator.dazim[a]  # effective track spacing
			track_list = self.generator.track_phis[a]
			for track in track_list:
				for fwd in (True, False):
					index0, index1 = track.get_indices(fwd)
					dists = track.trace(fwd)
					s1, sf, s2 = dists
					if calculate_area:
						self.fuel_area_tally += sf*wk*wa
						self.modr_area_tally += (s1 + s2)*wk*wa
					for p in range(self.quad.np):
						if self.model.boundary == "reflective":
							psi = self.psi[p, index0, fwd*1]	# angular flux --> pull from array
						elif self.model.boundary == "vacuum":
							psi = 0.0
						else:
							errstr = "Unknown Boundary Condition: {}"
							raise NotImplementedError(errstr.format(self.model.boundary))
						sintheta = self.quad.sinthetas[p]
						wp = self.quad.wp[p]  		# includes sintheta
						weight = wp * wa * wk
						# Moderator before cell (absorption, no source)
						if s1:
							attenuation = pylab.exp(-sig_mod * s1 / sintheta)
							delta_psi1 = (psi - self.model.mod.source) * (1 - attenuation)
							modr_flux += weight * delta_psi1 / self.model.mod.area
							psi -= delta_psi1
						# Fuel inside cell (scatter, with source)
						if sf:
							attenuation = pylab.exp(-sig_fuel * sf / sintheta)
							delta_psif = (psi - self.model.fuel.source) * (1 - attenuation)
							fuel_flux += weight * delta_psif / self.model.fuel.area
							psi -= delta_psif
						# Moderator after cell (absorption, no source)
						if s2:
							attenuation = pylab.exp(-sig_mod * s2 / sintheta)
							delta_psi2 = (psi - self.model.mod.source) * (1 - attenuation)
							modr_flux += weight * delta_psi2 / self.model.mod.area
							psi -= delta_psi2
						try:
							self.psi[p, index1, fwd*1] = psi
						except IndexError:
							print(self.psi.shape)
							raise

		fuel_flux *= 4*pylab.pi/sig_fuel
		if modr_flux:
			modr_flux *= 4*pylab.pi/sig_mod
		return fuel_flux, modr_flux


	def calculate(self):
		"""Do the flux calculation
		
		Tells the TrackGenerator to generate the tracks, which is an intensive
		process. Then, proceeds to iterate until the boundary flux and/or the
		
		Return:
		-------
		converged:      Boolean; whether the calculation converged
		fluxdiff:           float; the l2norm of the psi vector from the most recent iteration
		"""
		# Lay down the tracks and initialize the fluxes
		if not self.generator.generated:
			self.generator.generate()
		# And then iterate.
		count = 0
		fluxdiff, fdiff, mdiff = [1 + EPS]*3
		print("Sweeping...")
		while fluxdiff > EPS or fdiff > EPS or mdiff > EPS:
			psi_old = pylab.array(self.psi)
			fuel_flux, modr_flux = self.transport_sweep(calculate_area=not count)
			print("\nIteration", count)
			# Check for convergence
			if count > 0:
				fluxdiff = l2norm_2d(self.psi.sum(axis=2), psi_old.sum(axis=2))
				fdiff = abs(fuel_flux - self.fuel_flux)/self.fuel_flux
				mdiff = abs(modr_flux - self.modr_flux)/self.modr_flux
				print("norm:", fluxdiff)
				print("fuel diff: {:.4%} \t mod diff: {:.4%}".format(fdiff, mdiff))
				print("fuel flux: {:.5f} \t mod flux: {:.5f}".format(self.fuel_flux, self.modr_flux))
			self.fuel_flux = fuel_flux
			self.modr_flux = modr_flux
			# Don't allow infinite loops
			count += 1
			if count >= MAXITER:
				print("FAILED TO CONVERGE AFTER", count, "ITERATIONS")
				return False, fluxdiff

		print("...done.")
		print("Converged after", count, "iterations.")
		#print(self.psi.sum(axis=2).sum(axis=0))
		return True, fluxdiff
				
		
		
		
if __name__ == "__main__":
	import fsr
	fuel_fsr = fsr.FlatSourceRegion(cell.AREA_FUEL, cell.SIGMA_NF, QFSR, "Fuel")
	for sigma_a in cell.SIGMA_AS[1:2]:
		mod_fsr = fsr.FlatSourceRegion(cell.AREA_MOD, sigma_a, 0.0, "Mod")
		cell0 = cell.Cell(cell.PITCH, cell.RADIUS, fuel_fsr, mod_fsr,
						  plot = PLOT, boundary="reflective")
		trackgen = TrackGenerator(cell0, N_AZIM_2, D_AZIM)
		tabuchi = quadrature.YamamotoQuadrature(trackgen.phis[0,:], 3)
		calc = Calculator(cell0, trackgen, tabuchi, plot = PLOT)
		calc.calculate()
		flux_ratio = calc.fuel_flux / calc.modr_flux
		print("Flux in fuel: {:.5f}\tFlux in mod: {:.5f}\tRatio: {:.5f}".format(
			calc.fuel_flux, calc.modr_flux, flux_ratio))

		def errea(a1, a2):
			return abs(a1 - a2)/a1
		print()
		print("-" * 60)
		print("Fuel area")
		print("\tActual: {:.6f}\tEffective: {:.6f}".format(cell.AREA_FUEL, calc.fuel_area_tally))
		print("\tError: {:.4%}".format(errea(cell.AREA_FUEL, calc.fuel_area_tally)))
		print()
		print("Moderator area")
		print("\tActual: {:.6f}\tEffective: {:.6f}".format(cell.AREA_MOD, calc.modr_area_tally))
		print("\tError: {:.4%}".format(errea(cell.AREA_MOD, calc.modr_area_tally)))

		# Wigner-Bell math!
		'''
		print("Fuel XS: {:.4f} cm^-1; \t Mod XS: {:.4f} cm^-1".format(cell.SIGMA_NF, sigma_a))
		mkl = 2*cell.RADIUS		# mean Kord length
		pff = mkl*cell.SIGMA_NF / (1 + mkl*cell.SIGMA_NF)
		print("pff: {:.5f} \t 1 - pff: {:.5f}".format(pff, 1-pff))
		ratio = sigma_a / cell.SIGMA_NF / (1/pff - 1)*cell.AREA_RATIO
		print("Wigner prediction of flux ratio: {:.5f}".format(ratio))
		'''
		if PLOT:
			pylab.show()


