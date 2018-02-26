# Ratios
#
# Test flux ratio calculations

import cell as c
import fsr
import quadrature
import calculate
from trackgenerator import TrackGenerator

cell0 = c.Cell(c.PITCH, c.RADIUS, None, None)
gen = TrackGenerator(cell0, 8, 0.1)
gen.generate()
ty1 = quadrature.YamamotoQuadrature(gen.phis[:,0], 1)
report = ""
# Effectively homogenous
report += "\n\nCell 1: Fuel XS=1.0, Mod XS=1.0, Qfuel=1.0, Qmodr=1.0"
fuel1 = fsr.FlatSourceRegion(c.AREA_FUEL, xs=1.0, q=1.0)
mod1 = fsr.FlatSourceRegion(c.AREA_FUEL, xs=1.0, q=1.0)
cell1 = c.Cell(c.PITCH, c.RADIUS, fuel1, mod1, plot=False)
calc1 = calculate.Calculator(cell1, gen, ty1, plot=False)
calc1.calculate()
flux1 = calc1.fuel_flux / calc1.modr_flux
report += "\n\tFuel flux: {:.4f} \t Mod flux: {:.4f}".format(calc1.fuel_flux, calc1.modr_flux)
report += "\n\tFuel-to-mod flux ratio: {:.3f}".format(flux1)

# Double everything in fuel
report += "\n\nCell 2: Fuel XS=2.0, Mod XS=1.0, Qfuel=2.0, Qmodr=1.0"
fuel2 = fsr.FlatSourceRegion(c.AREA_FUEL, xs=2.0, q=2.0)
cell2 = c.Cell(c.PITCH, c.RADIUS, fuel2, mod1, plot=False)
calc2 = calculate.Calculator(cell2, gen, ty1, plot=False)
calc2.calculate()
flux2 = calc2.fuel_flux / calc2.modr_flux
report += "\n\tFuel flux: {:.4f} \t Mod flux: {:.4f}".format(calc2.fuel_flux, calc2.modr_flux)
report += "\n\tFuel-to-mod flux ratio: {:.3f}".format(flux2)

# Double source in fuel
report += "\n\nCell 3: Fuel XS=1.0, Mod XS=1.0, Qfuel=2.0, Qmodr=1.0"
fuel3 = fsr.FlatSourceRegion(c.AREA_FUEL, xs=1.0, q=2.0)
cell3 = c.Cell(c.PITCH, c.RADIUS, fuel3, mod1, plot=False)
calc3 = calculate.Calculator(cell3, gen, ty1, plot=False)
calc3.calculate()
flux3 = calc3.fuel_flux / calc3.modr_flux
report += "\n\tFuel flux: {:.4f} \t Mod flux: {:.4f}".format(calc3.fuel_flux, calc3.modr_flux)
report += "\n\tFuel-to-mod flux ratio: {:.3f}".format(flux3)

# Double source everywhere
report += "\n\nCell 4: Fuel XS=1.0, Mod XS=1.0, XS=1.0, Qfuel=2.0, Qmodr=2.0"
fuel4 = fsr.FlatSourceRegion(c.AREA_FUEL, xs=1.0, q=2.0)
mod4 = fsr.FlatSourceRegion(c.AREA_MOD, xs=1.0, q=2.0)
cell4 = c.Cell(c.PITCH, c.RADIUS, fuel4, mod4, plot=False)
calc4 = calculate.Calculator(cell4, gen, ty1, plot=False)
calc4.calculate()
flux4 = calc4.fuel_flux / calc4.modr_flux
report += "\n\tFuel flux: {:.4f} \t Mod flux: {:.4f}".format(calc4.fuel_flux, calc4.modr_flux)
report += "\n\tFuel-to-mod flux ratio: {:.3f}".format(flux4)

# Double source in fuel, no source in mod
report += "\n\nCell 5: Fuel XS=1.0, Mod XS=1.0, Qfuel=2.0, Qmodr=0.0"
mod5 = fsr.FlatSourceRegion(c.AREA_MOD, xs=1.0, q=0.0)
cell5 = c.Cell(c.PITCH, c.RADIUS, fuel4, mod5, plot=False)
calc5 = calculate.Calculator(cell5, gen, ty1, plot=False)
calc5.calculate()
flux5 = calc5.fuel_flux / calc5.modr_flux
report += "\n\tFuel flux: {:.5f} \t Mod flux: {:.5f}".format(calc5.fuel_flux, calc5.modr_flux)
report += "\n\tFuel-to-mod flux ratio: {:.3f}".format(flux5)

# Source in fuel, no source in mod
report += "\n\nCell 6: Fuel XS=1.0, Mod XS=1.0, Qfuel=1.0, Qmodr=0.0"
cell6 = c.Cell(c.PITCH, c.RADIUS, fuel1, mod5, plot=False)
calc6 = calculate.Calculator(cell6, gen, ty1, plot=False)
calc6.calculate()
flux6 = calc6.fuel_flux / calc6.modr_flux
report += "\n\tFuel flux: {:.5f} \t Mod flux: {:.5f}".format(calc6.fuel_flux, calc6.modr_flux)
report += "\n\tFuel-to-mod flux ratio: {:.3f}".format(flux6)


print(report)



