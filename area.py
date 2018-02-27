# Area calculator
#
# Example area calculation

from calculate import *
import fsr


def errea(a1, a2):
	return abs(a1 - a2) / a1


fuel = fsr.FlatSourceRegion(cell.AREA_FUEL, cell.SIGMA_NF, 1.0, "fuel")
modr = fsr.FlatSourceRegion(cell.AREA_MOD,            0.2, 0.0, "moderator")
cell0 = cell.Cell(cell.PITCH, cell.RADIUS, fuel, modr, plot = PLOT)

report = ""
for dazim in (0.5, 0.1, 0.05, 0.01, 0.005, 0.001):
	trackgen = TrackGenerator(cell0, nazim=2, dtarget=dazim)
	#area_quad = quadrature.AreaQuadrature(trackgen.phis[0, :], 1)
	area_quad = quadrature.YamamotoQuadrature(trackgen.phis[0, :], 1)
	calc = Calculator(cell0, trackgen, area_quad, plot=PLOT, eps=1E-5)
	calc.calculate()

	report += "\n\n" + "=" * 60
	report += "\nD_AZIM = {} cm".format(dazim)
	report += "\n" + "-" * 60
	report += "\nFuel area"
	report += "\n\tActual: {:.6f}\tEffective: {:.6f}".format(cell.AREA_FUEL, calc.fuel_area_tally)
	report += "\n\tError: {:.4%}".format(errea(cell.AREA_FUEL, calc.fuel_area_tally))
	report += "\n\nModerator area"
	report += "\n\tActual: {:.6f}\tEffective: {:.6f}".format(cell.AREA_MOD, calc.modr_area_tally)
	report += "\n\tError: {:.4%}".format(errea(cell.AREA_MOD, calc.modr_area_tally))

print(report)

