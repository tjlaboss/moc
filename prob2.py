# Problem 2
#
# 22.212/PSet03/Problem 2
#   by Travis J. Labossiere-Hickman (email: tjlaboss@mit.edu)

import math
from cell import *
from trackgenerator import TrackGenerator
import quadrature
import calculate

AZIM_PAIRS = 16		# Half the number of azimuthal angles
D_AZIM = 0.1		# Azimuthal ray spacing target
PLOT = False		# Whether to plot the tracks
EPS = 1E-5			# Numerical tolerance
SIGMA_F = 1/EPS		# Fuel XS for Dancoff factors

report = """
Problem 2:
----------"""
for sigma_a in SIGMA_AS[1:3]:
	header = "Mod XS = {:.2f} cm^-1".format(sigma_a)
	report += '\n\t' + header
	pincell = Cell(PITCH, RADIUS, SIGMA_F, sigma_a, plot=PLOT)
	sigma_f = pincell.sigma_n_fuel
	trackgen = TrackGenerator(pincell, AZIM_PAIRS, D_AZIM)
	trackgen.generate()
	ty3 = quadrature.YamamotoQuadrature(trackgen.phis[0, :], 3)
	calc = calculate.Calculator(pincell, trackgen, ty3, plot=PLOT, eps=EPS)
	print(header)
	# Get the vacuum flux from the first iteration
	calc.transport_sweep(calc.source)
	flux0 = calc.fuel_flux
	# Then converge
	calc.calculate()
	flux_ratio = calc.fuel_flux / calc.modr_flux
	report += "\n\t\tFuel/Mod flux ratio: {:.5f}".format(flux_ratio)
	l = 2*pincell.radius	# mean chord length
	#pff = l*sigma_f / (1 + l*sigma_f)
	pff = sigma_f / (sigma_f + 1/(l*N238) )
	report += "\n\t\tP_F->F:              {:.5f}".format(pff)
	#d = 0 # TODO: implement
	x = 4*math.pi*calc.source/sigma_f
	c = 1 - (x - calc.fuel_flux) / (x - flux0)
	report += "\n\t\tDancoff Correction:  {:.5f}".format(c)
	report += "\n"

print(report)
