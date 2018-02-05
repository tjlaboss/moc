# Calculate
#
# A script to do flux calculations

import math
import numpy
import cell
import quadrature
from functions import deg, l2norm_2d
from trackgenerator import TrackGenerator

# Global constants
PLOT = True
MAXITER = 100
EPS = 1E-3      # Convergence criterion  # TODO: Change back to 1E-5
N_AZIM_2 = 12    # number of azimuthal angle pairs
D_AZIM = 0.25    # target spacing between parallel tracks (cm)
QFSR = 1.0      # flat fixed source magnitude


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
	psi:            numpy.array; angular flux, size = (np, nxy)
						psi[:,p]: flux vector for polar angle p
						psi[i,p]: angular flux at node i for angle p
	nxy:            int; the number of nodes on the x and y boundaries.
						Serves as the length of psi[:,p]
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
		self.source = source
		self.eps = eps
		self.plot = plot
		self.psi = numpy.zeros((self.quad.np, 4*self.generator.nazim))
		self.modr_flux = -1.0
		self.fuel_flux = -1.0
		
		
	def calculate(self):
		"""Do the flux calculation
		
		DESCRIBE WHAT I'VE DONE HERE
		
		Return:
		-------
		converged:      Boolean; whether the calculation converged
		diff:           float; the l2norm of the psi vector from the most recent iteration
		"""
		sig_mod = self.model.sigma_y_mod
		sig_fuel = self.model.sigma_n_fuel
		# Precalculations
		q = QFSR/sig_fuel/(4*math.pi)
		# Lay down the tracks and initialize the fluxes
		self.generator.generate()
		# And then iterate.
		diff = 1 + EPS
		count = 0
		while diff > EPS:
		#while count < 10:    # restore the 'while' loop
			# Reset the fluxes before sweeping
			psi_old = numpy.empty(self.psi.shape)
			psi_old[:,:] = self.psi
			fuel_flux = QFSR/sig_fuel
			modr_flux = 0.0  # No source in fuel
			for a in range(self.quad.na):
				wa = self.quad.wa[a]
				dazim = self.generator.dazim[a] # effective track spacing
				track_list = self.generator.track_phis[a]
				#track0 = track_list[0]
				#track = track_list[1]
				#while track is not track0:
				# TODO: Determine if it matters where we start??
				for track in track_list:
					for fwd in (True, False):
						index0, index1 = track.get_indices(fwd)
						if count == -1:#0: # debug, delete this
							print("\n", int(deg(track.phi)), "deg")
							if fwd:
								print(index0, "->", index1, " \t ", track.key0, track.key1)
							else:
								print(index0, "->", index1, " \t ", track.key1, track.key0)
						# Precalculate some stuff before the polar quadrature loop
						dists = track.trace(fwd)
						if self.plot:
							self.model.plot_track(track, dists)
						s1, sf, s2 = dists
						for p in range(self.quad.np):
							psi = self.psi[p, index0]  # angular flux --> pull from array
							sintheta = self.quad.sinthetas[p]
							wp = self.quad.wp[p]    # includes sintheta
							weight = 2*math.pi*wp*wa*dazim
							# Moderator before cell (absorption, no source)
							if s1:
								attenuation = math.exp(-sig_mod*s1/sintheta)
								delta_psi1 = psi*(1 - attenuation)
								modr_flux += weight*delta_psi1
								psi -= delta_psi1
							# Fuel inside cell (scatter, with source)
							if sf:
								attenuation = math.exp(-sig_fuel*sf/sintheta)
								delta_psif = (psi - q/sig_fuel)*(1 - attenuation)
								fuel_flux += weight*delta_psif
								psi -= delta_psif
							# Moderator after cell (absorption, no source)
							if s2:
								attenuation = math.exp(-sig_mod*s2/sintheta)
								delta_psi2 = psi*(1 - attenuation)
								modr_flux += weight*delta_psi2
								psi -= delta_psi2
							try:
								self.psi[p, index1] = psi
							except IndexError:
								print(self.psi.shape)
								raise

					#print(track, track.next_track)
					#track = track.next_track
			
			print(count)
			# Check for convergence
			if count > 0:
				diff = l2norm_2d(self.psi, psi_old)
				fdiff = abs(fuel_flux - self.fuel_flux)/self.fuel_flux
				mdiff = abs(modr_flux - self.modr_flux)/self.modr_flux
				print("norm:", diff)
				print("fuel diff: {:.4%} \t mod diff: {:.4%}".format(fdiff, mdiff))
				#print(self.psi)
			self.fuel_flux = fuel_flux
			self.modr_flux = modr_flux
			# Don't allow infinite loops
			count += 1
			if count >= MAXITER:
				print("FAILED TO CONVERGE AFTER", count, "ITERATIONS")
				return False, diff
		
		print("Converged after", count, "iterations.")
		#print(self.psi)
		print(self.psi[0,:] + self.psi[1,:] + self.psi[2,:])
		print("Flux in fuel: {}\tFlux in mod: {}".format(self.fuel_flux, self.modr_flux))
		return True, diff
				
		
		
		

for sigma_a in cell.SIGMA_AS[0:2]:
	cell0 = cell.Cell(cell.PITCH, cell.RADIUS,
	                  cell.SIGMA_NF, sigma_a, plot = PLOT)
	trackgen = TrackGenerator(cell0, N_AZIM_2, D_AZIM)
	tabuchi = quadrature.YamamotoQuadrature(trackgen.phis[0,:], 3)
	calc = Calculator(cell0, trackgen, tabuchi, plot = PLOT)
	calc.calculate()
