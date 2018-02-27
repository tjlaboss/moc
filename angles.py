# Angles
#
# Converge in azimuthal angle

from cell import *
from trackgenerator import TrackGenerator
import quadrature
import calculate

D_AZIM = 0.001		# Azimuthal ray spacing target
PLOT = False		# Whether to plot the tracks

report = """
Problem 1:
----------"""
fuel_fsr = fsr.FlatSourceRegion(AREA_FUEL, SIGMA_NF, 1.0, "Fuel")
last = -1
for naz2 in [2, 4, 8, 16, 32, 64, 128, 256, 512]:
	sigma_a = SIGMA_AS[1]
	report += "\n\tN_Azim = {}".format(naz2*2)
	mod_fsr = fsr.FlatSourceRegion(AREA_MOD, sigma_a, 0.0, "Moderator")
	pincell = Cell(PITCH, RADIUS, fuel_fsr, mod_fsr, plot=PLOT)
	trackgen = TrackGenerator(pincell, naz2, D_AZIM)
	ty3 = quadrature.YamamotoQuadrature(trackgen.phis[0, :], 3)
	calc = calculate.Calculator(pincell, trackgen, ty3, plot=PLOT)
	calc.calculate()
	flux_ratio = calc.fuel_flux / calc.modr_flux
	diff = abs(flux_ratio - last) / last
	last = flux_ratio
	report += "\n\t\tFuel-to-mod flux ratio: {:.5f}".format(flux_ratio)
	if naz2 != 2:
		report += "  (diff: {:.3%})".format(diff)
	report += "\n"

print(report)
