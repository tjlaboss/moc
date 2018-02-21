# Problem 1
#
# 22.212/PSet03/Problem 1
#   by Travis J. Labossiere-Hickman (email: tjlaboss@mit.edu)

from cell import *
from trackgenerator import TrackGenerator
import quadrature
import calculate

AZIM_PAIRS = 16		# Half the number of azimuthal angles
D_AZIM = 0.1		# Azimuthal ray spacing target
PLOT = False		# Whether to plot the tracks
EPS = 1E-5			# Numerical tolerance

report = """
Problem 1:
----------"""
for sigma_a in SIGMA_AS[1:2]:
	report += "\n\tMod XS = {:.2f} cm^-1".format(sigma_a)
	pincell = Cell(PITCH, RADIUS, SIGMA_NF, sigma_a, plot=PLOT)
	trackgen = TrackGenerator(pincell, AZIM_PAIRS, D_AZIM)
	ty3 = quadrature.YamamotoQuadrature(trackgen.phis[0, :], 3)
	calc = calculate.Calculator(pincell, trackgen, ty3, plot=PLOT)
	calc.calculate()
	flux_ratio = calc.fuel_flux / calc.modr_flux

print(report)
