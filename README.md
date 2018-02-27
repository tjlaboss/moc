# moc

Method-of-Characteristics solver for 22.212

Containing:

### cell.py

Module with problem parameters and plotter. Run this to view the geometry and a single ray.

### fsr.py

Module for Flat Source Regions.

### ray.py

Module containing the ray tracer.

### trackgenerator.py

Module with a generator for cyclic tracks. Run this to view them.

### quadrature.py

Module with quadrature sets. 
* Azimuthal: Equal angle
* Polar: Uniform Distributed, Leonard's Optimum, Tabuchi-Yamamoto

### calculate.py

Module containing the transport solver. Run this to get fuel and moderator fluxes.

### functions.py

Module with a couple of generic functions

### area.py

Test out quadratures and check for spatial convergence by calculating the effective areas of FSRs.

### angles.py

Check for angular convergence by running an increasingly high number of azimuthal angles. SLOW.

### ratios.py

Test flux ratios as a function of source strength and region cross section.

### prob1.py

22.212/PSet03/Problem 1: Flux Ratios

### prob2.py

22.212/PSet03/Problem 2: Dancoff Factors

### Report/

Project report and LaTeX source

