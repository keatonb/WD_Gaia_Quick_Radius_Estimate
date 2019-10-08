#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 13:39:24 2019

https://github.com/keatonb/WD_Gaia_Quick_Radius_Estimate

@author: keatonb
"""

from __future__ import print_function, division

import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia
from astroquery.vizier import Vizier
import numpy as np
from scipy.interpolate import interp1d
import os

from BergeronToPandas import BergeronToPandas

def WDGaiaRadius(ra,dec,teff,teff_err=None,searchradius=5.0,modelmass=0.6,modelpath='AllModels/'):
    """Estimate WD Radius from Gaia astrometry.
    
    Queries for Gaia data and distances from Bailer-Jones et al. (2018, ApJ, 156, 2).
    Scales Bergeron et al. DA (no DBs currently) cooling models 
    (http://www.astro.umontreal.ca/~bergeron/CoolingModels/)
    until the radius reproduces the Gaia magnitude at the parallactic distance.
    Does not account for extinction or reddening.  
    Uses https://github.com/keatonb/BergeronToPandas to read in Bergeron models.
    If error on Teff given, includes contribution roughly in quadrature.
    Returns measurement and lower, upper bounds defined by Bailer-Jones et al. 
    confidence interval.
    
    ra, dec should be in decimal degrees
    teff,teff_err in Kelvin from spectroscopy
    searchradius is in arcseconds
    modelmass in solar units must match a file from the Bergeron models
                that is new enough to include synthetic Gaia magnitudes.
    """
    #First query Gaia sources
    coord = SkyCoord(ra=ra, dec=dec, unit=(u.degree, u.degree), frame='icrs')
    searchradius = u.Quantity(searchradius, u.arcsec)
    j = Gaia.cone_search_async(coord, searchradius)
    r = j.get_results()
    print("Gaia sources within {}: {}".format(searchradius,len(r)))
    if len(r) > 0:
        r.sort("dist")
        #Accept nearest
        sourceid = r["source_id"][0]
        Gmag = r["phot_g_mean_mag"][0]
        print("nearest has Gaia G mag of {}".format(Gmag))
        #Then get Bailer-Jones distance
        v = Vizier(columns=['*'],
                column_filters={"Source":str(sourceid)}).get_catalogs('I/347')
        dist = np.array([v[0]["rest"][0],v[0]["b_rest"][0],v[0]["B_rest"][0]])

        #Referencing Bergeron cooling models
        modelfile = os.path.join(modelpath, "Table_Mass_{:.1f}".format(modelmass))
        print("Referencing Bergeron model file at "+modelfile)
        
        models = BergeronToPandas(modelfile).values()[0]
        
        #Constants
        G = 6.67259e-8 #cm3 g-1 s-2
        Msun = 1.99e33 #g
        Rsun = 6.96e10 #cm
        
        modelteff = models['Teff'].values
        modellogg = models['logg'].values
        modelradius = np.sqrt(G*modelmass*Msun/(10.**modellogg))/Rsun #solar radii
        
        modelGmag = models['G'].values
        
        #Define interpolation functions
        teff2radius = interp1d(modelteff,modelradius)
        teff2Gmag = interp1d(modelteff,modelGmag)
        
        distancemodulus = 5.*np.log10(dist)-5.
        
        #Start with values assuming along loaded model track
        trackradius = teff2radius(teff)
        trackGmag = teff2Gmag(teff)
        estGmag = trackGmag + distancemodulus 
        
        #Scale to measured magnitude
        remaining = estGmag - Gmag
        fluxratio = 2.512**remaining
        radiusratio = np.sqrt(fluxratio)
        gaiaradius = trackradius * radiusratio
        
        #Include effect of Teff_err if given
        if teff_err is not None:
            testteff = np.array([teff,teff+teff_err,teff-teff_err])
            trackradius = teff2radius(testteff)
            trackGmag = teff2Gmag(testteff)
            estGmag = trackGmag + distancemodulus[0]
            remaining = estGmag - Gmag
            fluxratio = 2.512**remaining
            radiusratio = np.sqrt(fluxratio)
            gr2 = trackradius * radiusratio
            #combine roughly in quadrature
            gaiaradius = gaiaradius[0] + np.sqrt((gaiaradius - gaiaradius[0])**2. + 
                                                 (gr2 - gr2[0])**2.)*np.array([1,-1,1])
        return gaiaradius
    else: #No Gaia results?
        print("Returning 'None's")
        return None,None,None
    