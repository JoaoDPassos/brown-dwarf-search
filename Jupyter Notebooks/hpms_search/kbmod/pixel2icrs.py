from pathlib import Path
from astropy.io import fits
from astropy.wcs import WCS

mas_p_pixel =  249.85968
img_path = Path("./kbmod_fits_examples/g_band/rings.v3.skycell.2327.023.wrp.g.55334_37545.combined.fits")

with fits.open(img_path.open('rb')) as hdul:
    hdr = hdul[1].header
    wcs = WCS(hdr)

def pix2mas(pix):
    return pix * mas_p_pixel

def mas2pix(mas):
    return mas / mas_p_pixel

def pix2as(pix):
    return pix2mas(pix)/1000

def as2pix(ars):
    return mas2pix(ars)*1000

def ppd2masyr(ppd):
    # Takes a veloctiy in pixels per day (ppd) and returns that velocity in mas/yr
    return ppd * 365 * mas_p_pixel

def masyr2ppd(masyr):
    # Takes a veloctiy in pixels mas/yr and returns that velocity in pixels per day (ppd)
    return masyr/(365*mas_p_pixel)
