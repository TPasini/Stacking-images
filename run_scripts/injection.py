import astropy.io.fits as fits
import os
import numpy as np
from lib.lib_fits import flatten
from astropy.wcs import WCS
from astropy.cosmology import FlatLambdaCDM

def Lorentz2D(flux, re_as, center_ra, center_dec, x, y, pixscale, beam):
    #pixscale = float(1.5)
    re_px = re_as / pixscale
    flux_freq = flux #* (freq / ref_freq) ** spidx
    I_0 = flux_freq / (2. * np.pi * re_as ** 2)  # Jy/arcsec^2
    r = np.sqrt((x - center_ra) ** 2 + (y - center_dec) ** 2)
    lorentz = I_0 * np.exp(- r / re_px) * pixscale ** 2  # Jy/pixel
    lorentz = I_0 * np.exp(- r / re_px) * 1.133 * beam ** 2
    #print((flux_freq / (2. * np.pi * re_px ** 2)), np.max(lorentz))

    return lorentz

def inject_images(file, ra, dec, z, flux):

    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

    root_name = str(file[0:22])

    if os.path.exists(root_name+'_iminj.fits'):
        os.remove(root_name+'_iminj.fits')

    hdul = fits.open(file)

    pixspacing = ((hdul[0].header['CDELT2']) * 60) #in arcmin
    pixscale = pixspacing*60
    beam = ((hdul[0].header['BMAJ']) * 3600)

    head, data = flatten(file)
    w = WCS(head)

    dataorig = hdul[0].data
    header_orig = hdul[0].header

    rapix, decpix = w.wcs_world2pix(ra, dec, 0)

    min_x = int(rapix) - 150
    max_x = int(rapix) + 150
    min_y = int(decpix) - 150
    max_y = int(decpix) + 150
    y, x = np.mgrid[min_y:max_y, min_x:max_x]
    conv = cosmo.kpc_proper_per_arcmin(z)
    re=167 #kpc
    re_as = (re / (conv.value/60))
    fakesource = Lorentz2D(flux, re_as, rapix, decpix, x, y, pixscale, beam)
    dataorig[0, 0, min_y:max_y, min_x:max_x] += fakesource
    fits.writeto(root_name+'_iminj.fits', dataorig, header=header_orig)
