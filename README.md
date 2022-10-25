# Seasonality-in-Sea-Surface-Salinity-Balance-Equation
Data and code to reproduce results in "Evaluation of Seasonality in Sea Surface Salinity Balance Equation via Function Registration"

## Set up
Before running the code, please download [the fdasrvf package](https://github.com/jdtuck/fdasrvf_MATLAB) and [the functional data analysis package](https://www.psych.mcgill.ca/misc/fda/downloads/FDAfuns/R/inst/Matlab/fdaM/).

## Code description:
- [ArgoSIOForcingAdvectionVerticalMixing.nc](ArgoSIOForcingAdvectionVerticalMixing.nc): Salinity of the sea surface and forcing minus advection and vertical entrainme data from the following sources
    - Roemmich & Gilson SSS: https://sio-argo.ucsd.edu/RG_Climatology.html
    - OAFlux evaporation: https://climatedataguide.ucar.edu/climate-data/oaflux-objectively-analyzed-air-sea-fluxes-global-oceans
    - MIMOC mixed layer depth: https://www.pmel.noaa.gov/mimoc/
    - GPM IMERG precipitation: https://gpm.nasa.gov/data/imerg
    - OSCAR surface velocity: https://podaac.jpl.nasa.gov/dataset/OSCAR_L4_OC_third-deg
    - NOAA Coastwatch Ekman Upwelling: https://coastwatch.pfeg.noaa.gov/erddap/griddap/erdQAstressmday.html

- [This code](SSS_FMAV_registration.m) reproduces results presented in the paper.
