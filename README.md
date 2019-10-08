# WD_Gaia_Quick_Radius_Estimate
Estimate WD Radius from Gaia astrometry.
  
Queries for Gaia data and distances from Bailer-Jones et al. (2018, ApJ, 156, 2).
Scales Bergeron et al. DA (no DBs currently) cooling models (http://www.astro.umontreal.ca/~bergeron/CoolingModels/)
until the radius reproduces the Gaia magnitude at the parallactic distance.
Does not account for extinction or reddening.  Uses https://github.com/keatonb/BergeronToPandas to read in Bergeron models. 
To estimate uncertainty from Teff error, perturb input Teff by errorbars.
    Returns measurement and lower, upper bounds defined by Bailer-Jones et al. 
    confidence interval.

ra, dec should be in decimal degrees

teff in Kelvin from spectroscopy/colors

searchradius is in arcseconds

modelmass in solar units must match a file from the Bergeron models that is new enough to include synthetic Gaia magnitudes.
