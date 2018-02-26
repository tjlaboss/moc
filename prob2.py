# Problem 2
#
# 22.212/PSet03/Problem 2
#   by Travis J. Labossiere-Hickman (email: tjlaboss@mit.edu)

import math
from cell import *
from trackgenerator import TrackGenerator
import quadrature
import calculate
import fsr

AZIM_PAIRS = 16		# Half the number of azimuthal angles
D_AZIM = 0.1		# Azimuthal ray spacing target
PLOT = False		# Whether to plot the tracks
EPS = 1E-5			# Numerical tolerance
SIGMA_F = 0.1/EPS	# Fuel XS for Dancoff factors

# All parts of the problem will be done on the same pincell geometry.
# Initialize a single track generator and quadrature based on that geometry.
cell0 = Cell(PITCH, RADIUS, None, None, plot=PLOT)
trackgen = TrackGenerator(cell0, AZIM_PAIRS, D_AZIM)
trackgen.generate()
ty3 = quadrature.YamamotoQuadrature(trackgen.phis[0, :], 3)
report = """
Problem 2:
----------"""
fuel_fsr = fsr.FlatSourceRegion(AREA_FUEL, SIGMA_F, 1.0, "Fuel (Dancoff)")
for sigma_a in SIGMA_AS[:]:
	header = "Mod XS = {:.3f} cm^-1".format(sigma_a)
	report += '\n\t' + header
	mod_fsr = fsr.FlatSourceRegion(AREA_MOD, sigma_a, 0.0, "Moderator")
	pincell_r = Cell(PITCH, RADIUS, fuel_fsr, mod_fsr, plot=PLOT,
					 boundary="reflective")
	pincell_v = Cell(PITCH, RADIUS, fuel_fsr, mod_fsr, plot=PLOT,
					 boundary="vacuum")
	sigma_f = fuel_fsr.xs
	calc_r = calculate.Calculator(pincell_r, trackgen, ty3, plot=PLOT, eps=EPS)
	calc_v = calculate.Calculator(pincell_v, trackgen, ty3, plot=PLOT, eps=EPS)
	print(header)
	# Get the isolated (vacuum) flux
	flux0 = calc_v.transport_sweep()[0]
	# Get the infinite lattice (reflective) flux
	calc_r.calculate()
	flux1 = calc_r.fuel_flux
	flux_ratio_r = calc_r.fuel_flux / calc_r.modr_flux
	# Wigner Rational Approximation
	l = 2 * pincell_r.radius	# mean chord length
	pff_wigner = l * sigma_f / (1 + l * sigma_f)
	# Dancoff correction
	x = fuel_fsr.q/sigma_f
	c = 1 - (x - calc_r.fuel_flux) / (x - flux0)
	xf = mod_fsr.xs * mod_fsr.area / (flux_ratio_r * fuel_fsr.xs * fuel_fsr.area)
	pff_actual = 1 / (1 + xf)
	report += "\n\t\tFuel/Mod flux ratio: {:.5f}".format(flux_ratio_r)
	report += "\n\t\tP_F->F (wigner):     {:.5f}".format(pff_wigner)
	report += "\n\t\tP_F->F (actual):     {:.5f}".format(pff_actual)
	report += "\n\t\tDancoff Correction:  {:.5f}".format(c)
	report += "\n"

print(report)
