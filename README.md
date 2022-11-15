# WD_Gaia_Quick_Radius_Estimate
Estimate white dwarf radius from Gaia DR3 astrometry.
    
Queries for Gaia DR3 data and distances from Bailer-Jones et al. 
(2021, AJ, 161, 147).
Scales Bergeron et al. DA or DB cooling models 
(http://www.astro.umontreal.ca/~bergeron/CoolingModels/)
until the radius reproduces the Gaia magnitude at the parallactic distance.
Corrects for extinction if coefficient A_G is provided in magnitudes, 
following the prescription of Gentile-Fusillo et al. (2021, MNRAS, 508, 3877).
Uses https://github.com/keatonb/BergeronToPandas to read in Bergeron models.
If error on Teff given, includes contribution roughly in quadrature.
Returns measurement and lower, upper bounds defined by Bailer-Jones et al. 
confidence interval.

ra, dec should be in decimal degrees

teff,teff_err in Kelvin

searchradius is in arcseconds

modelmass in solar units must match a file from the Bergeron models
            that is new enough to include synthetic Gaia DR3 magnitudes.
