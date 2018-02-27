# Problem 1
#
# 22.212/PSet03/Problem 1
#   by Travis J. Labossiere-Hickman (email: tjlaboss@mit.edu)

from cell import *
from trackgenerator import TrackGenerator
import quadrature
import calculate

AZIM_PAIRS = 128	# Half the number of azimuthal angles
D_AZIM = 0.005		# Azimuthal ray spacing target
PLOT = False		# Whether to plot the tracks
EPS = 1E-5			# Numerical tolerance

report = """
Problem 1:
----------"""
fuel_fsr = fsr.FlatSourceRegion(AREA_FUEL, SIGMA_NF, 1.0, "Fuel")
for sigma_a in SIGMA_AS:
	report += "\n\tMod XS = {:.2f} cm^-1".format(sigma_a)
	mod_fsr = fsr.FlatSourceRegion(AREA_MOD, sigma_a, 0.0, "Moderator")
	pincell = Cell(PITCH, RADIUS, fuel_fsr, mod_fsr, plot=PLOT)
	trackgen = TrackGenerator(pincell, AZIM_PAIRS, D_AZIM)
	ty3 = quadrature.YamamotoQuadrature(trackgen.phis[0, :], 3)
	calc = calculate.Calculator(pincell, trackgen, ty3, plot=PLOT)
	calc.calculate()
	flux_ratio = calc.fuel_flux / calc.modr_flux
	report += "\n\t\tFuel-to-mod flux ratio: {:.5f}".format(flux_ratio)
	report += "\n"

print(report)
print("DAZIM: {}, NAZIM2: {}".format(D_AZIM, AZIM_PAIRS))
