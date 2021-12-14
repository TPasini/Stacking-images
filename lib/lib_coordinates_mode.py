#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2019 - Francesco de Gasperin
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

# LOFAR Transients Key Project
#
# General purpose astronomical coordinate handling routines.
#
# Modified (simplified) by Casey to work with pyDAL scripts

import math, numpy, scipy, logging, datetime, pytz

core_lat = 52.9088
core_lon = -6.8689 # Note that wcstools takes a +ve longitude as WEST.

wcstools_path= '/data/users/swinbank/app/lib/libwcstools.so.1'

def julian_date(time=None, modified=False):
    """Return the Julian Date: the number of days (including fractions) which
    have elapsed between noon, UT on 1 January 4713 BC and the specified time.
    If modified is True, return the Modified Julian Date: the number of days
    (including fractions) which have elapsed between the start of 17 November
    1858 AD and the specified time. Takes a datetime.datetime object as
    input."""
    if not time:
        time = datetime.datetime.now(pytz.utc)
    mjdstart = datetime.datetime(1858, 11, 17, tzinfo=pytz.utc)
    mjd = time - mjdstart
    mjd_daynumber = mjd.days + mjd.seconds / (24.0*60**2) + mjd.microseconds / (24.0*60**2*1000**2)
    if modified:
        return mjd_daynumber
    return 2400000.5 + mjd_daynumber

def jd2lst(jd, lon=core_lon):
    """Return the Local Sidereal Time in seconds.  Starts with Julian
    Date.  Default logitude is that of LOFAR core: note that wcstools takes a
    positive longitude as being WEST."""
    import ctypes
    wcstools = ctypes.cdll.LoadLibrary(wcstools_path)
    wcstools.setlongitude(ctypes.c_double(lon))
    wcstools.jd2lst.restype = ctypes.c_double
    c_jd = ctypes.c_double(jd)
    return wcstools.jd2lst(c_jd)

def mjd2lst(mjd, lon=core_lon):
    """Modified Julian Date to Local Sidereal Time; uses jd2lst"""
    return jd2lst(mjd + 2400000.5, lon=lon)

def sec2deg(seconds):
    """Seconds of time to degrees of arc"""
    return 15.0 * seconds / 3600.0

def sec2days(seconds):
    """Seconds to number of days"""
    return seconds / (24.0 * 3600)

def sec2hms(seconds):
    """Seconds to hours, minutes, seconds"""
    hours, seconds = divmod(seconds, 60**2)
    minutes, seconds = divmod(seconds, 60)
    return (int(hours), int(minutes), seconds)

def mjds2lst(mjds, lon=core_lon):
    """Return the local sidereal time in degrees.  Starts with mean julian
    date in seconds, as is the time column from a measurement set.  Default
    logitude is that of LOFAR core."""
    # That is, this should be a drop-in replacement for old_mjds2lst(), below.
    return sec2deg(mjd2lst(sec2days(mjds), lon=lon))

def altaz(mjds, ra, dec, lat=core_lat):
    """Calculates the azimuth and elevation of source from time and position
    on sky.  Takes MJD in seconds and ra, dec in degrees.  Returns (alt,az) in
    degrees."""

    #compute hour angle in degrees
    ha = mjds2lst( mjds ) - ra
    if (ha < 0): ha = ha + 360

    #convert degrees to radians
    ha, dec, lat = list(map(math.radians, (ha, dec, lat)))

    #compute altitude in radians
    sin_alt = math.sin(dec)*math.sin(lat) + math.cos(dec)*math.cos(lat)*math.cos(ha)
    alt = math.asin(sin_alt)
    
    #compute azimuth in radians
    #divide by zero error at poles or if alt = 90 deg
    cos_az = (math.sin(dec) - math.sin(alt)*math.sin(lat))/(math.cos(alt)*math.cos(lat))
    az = math.acos(cos_az)

    #convert radians to degrees
    hrz_altitude, hrz_azimuth = list(map(math.degrees, (alt, az)))

    #choose hemisphere
    if (math.sin(ha) > 0): hrz_azimuth = 360 - hrz_azimuth;

    return hrz_altitude, hrz_azimuth

def ratohms(radegs):
    """Convert RA in decimal degrees format to hours, minutes,
    seconds format.

    Keyword arguments:
    radegs -- RA in degrees format

    Return value:
    ra -- tuple of 3 values, [hours,minutes,seconds]

    """

    radegs %= 360
    raseconds = radegs * 3600 / 15.0
    return sec2hms(raseconds)

def dectodms(decdegs):
    """Convert Declination in decimal degrees format to deg, minutes,
    seconds format.

    Keyword arguments:
    decdegs -- Dec. in degrees format

    Return value:
    dec -- list of 3 values, [degrees,minutes,seconds]

    """
    
    if abs(decdegs) > 90:
        raise ValueError
    decd = int(decdegs)
    decm = int((decdegs-decd)*60)
    decs = (((decdegs-decd)*60)-decm)*60
    if decd < 0:
        decm = -1*decm
        decs = -1*decs
    dec = (decd,decm,decs)
    return dec

def hmstora(rah,ram,ras):
    """Convert RA in hours, minutes, seconds format to decimal
    degrees format.

    Keyword arguments:
    rah,ram,ras -- RA values (h,m,s)

    Return value:
    radegs -- RA in decimal degrees

    """
    hrs = (float(rah)+(float(ram)/60)+(float(ras)/3600.0)) % 24

    return 15*hrs
# Converts an hms format RA to decimal degrees

def dmstodec(decd,decm,decs):
    """Convert Dec in degrees, minutes, seconds format to decimal
    degrees format.

    Keyword arguments:
    decd,decm,decs -- list of Dec values (d,m,s)

    Return value:
    decdegs -- Dec in decimal degrees

    """
    if float(decd) < 0:
        decm = -1*float(decm)
        decs = -1*float(decs)
        
    decdegs = float(decd)+(float(decm)/60.)+(float(decs)/3600.0)
    
    if abs(decdegs) > 90:
        raise ValueError

    return decdegs
# Converts a dms format Dec to decimal degrees  

def getCoord(ra, dec):
    """ Convert 2 strings to a RA and DEC in degrees.
	The strings can be in ##:##:## ##:##:## format
	or in ##d##m##s ##h##m##s format
	or in ###.### ###.###
	"""
    try:
        ra = float(ra)
        dec = float(dec)
    except:
        import re
        # split ra and dec using : or [h|d]sm format
        if ':' in ra and ':' in dec:
            rah, ram, ras = ra.split(':')
            ra = hmstora(float(rah), float(ram), float(ras))
            decd, decm, decs = dec.split(':')
            dec = dmstodec(float(decd), float(decm), float(decs))
        else:
            rah, ram, ras = re.split("h|m|s", ra)[0:3]
            ra = hmstora(float(rah), float(ram), float(ras))
            decd, decm, decs = re.split("d|m|s", dec)[0:3]
            dec = dmstodec(float(decd), float(decm), float(decs))
        if re.search('-0*$',decd): dec *= -1 # -0 bug fix

    return ra, dec

def angsep(ra1,dec1,ra2,dec2):
    """Find the angular separation of two sources, in arcseconds,
    using the proper spherical trig formula

    Keyword arguments:
    ra1,dec1 - RA and Dec of the first source, in decimal degrees
    ra2,dec2 - RA and Dec of the second source, in decimal degrees

    Return value:
    angsep - Angular separation, in arcsec

    """
    if ra1 == ra2 and dec1 == dec2: return 0

    b = (numpy.pi/2)-numpy.radians(dec1)
    c = (numpy.pi/2)-numpy.radians(dec2)

    return 3600*numpy.degrees(numpy.arccos((numpy.cos(b)*numpy.cos(c))+(numpy.sin(b)*numpy.sin(c)*numpy.cos(numpy.radians(ra1-ra2)))))

def angsep2(ra1deg, dec1deg, ra2deg, dec2deg):
    """Returns angular separation between two coordinates (all in degrees)"""
    import math

    if ra1deg == ra2deg and dec1deg == dec2deg: return 0

    ra1rad=ra1deg*math.pi/180.0
    dec1rad=dec1deg*math.pi/180.0
    ra2rad=ra2deg*math.pi/180.0
    dec2rad=dec2deg*math.pi/180.0

    # calculate scalar product for determination
    # of angular separation
    x=math.cos(ra1rad)*math.cos(dec1rad)*math.cos(ra2rad)*math.cos(dec2rad)
    y=math.sin(ra1rad)*math.cos(dec1rad)*math.sin(ra2rad)*math.cos(dec2rad)
    z=math.sin(dec1rad)*math.sin(dec2rad)

    if x+y+z >= 1: rad = 0
    else: rad=math.acos(x+y+z)

    # Angular separation
    deg=rad*180/math.pi
    return deg

# Find angular separation of 2 positions, in arcseconds

def alphasep(ra1,ra2,dec1,dec2):
    """Find the angular separation of two sources in RA, in arcseconds

    Keyword arguments:
    ra1,dec1 - RA and Dec of the first source, in decimal degrees
    ra2,dec2 - RA and Dec of the second source, in decimal degrees

    Return value:
    angsep - Angular separation, in arcseconds

    """

    Dra = abs(ra1-ra2)
    if Dra > 180: Dra = 360-Dra
    return 3600*Dra*math.cos(math.radians((dec1+dec2)/2.0))

# Find angular separation in RA of 2 positions, in arcseconds

def deltasep(dec1,dec2):
    """Find the angular separation of two sources in Dec, in arcseconds

    Keyword arguments:
    dec1 - Dec of the first source, in decimal degrees
    dec2 - Dec of the second source, in decimal degrees

    Return value:
    angsep - Angular separation, in arcseconds

    """

    return 3600*abs(dec1-dec2)

# Find angular separation in Dec of 2 positions, in arcseconds

def alpha(l,m,alpha0,delta0):
    """Convert a coordinate in l,m into an coordinate in RA

    Keyword arguments:
    l,m -- direction cosines, given by (offset in cells) x cellsi (radians)
    alpha_0, delta_0 -- centre of the field

    Return value:
    alpha -- RA in decimal degrees

    """
    return alpha0+(math.degrees(math.atan2(l,((math.sqrt(1-(l*l)-(m*m))*math.cos(math.radians(delta0)))-(m*math.sin(math.radians(delta0)))))))

# Find the RA of a point in a radio image, given l,m and field centre

def delta(l,m,alpha0,delta0):
    """Convert a coordinate in l,m into an coordinate in Dec

    Keyword arguments:
    l,m -- direction cosines, given by (offset in cells) x cellsi (radians)
    alpha_0, delta_0 -- centre of the field

    Return value:
    delta -- Dec in decimal degrees

    """
    return math.degrees(math.asin(m*math.cos(math.radians(delta0))+(math.sqrt(1-(l*l)-(m*m))*math.sin(math.radians(delta0)))))

# Find the declination of a point in a radio image, given l,m and field centre

def l(ra,dec,cra,cdec,incr):
    """Convert a coordinate in RA,Dec into a direction cosine m

    Keyword arguments:
    ra,dec -- Source location
    cra,cdec -- centre of the field
    incr -- number of degrees per pixel (negative in the case of RA)

    Return value:
    l -- Direction cosine

    """
    return (math.cos(math.radians(dec))*math.sin(math.radians(ra-cra)))/(math.radians(incr))

# Find the l direction cosine in a radio image, given an RA and Dec and the
# field centre

def m(ra,dec,cra,cdec,incr):
    """Convert a coordinate in RA,Dec into a direction cosine m

    Keyword arguments:
    ra,dec -- Source location
    cra,cdec -- centre of the field
    incr -- number of degrees per pixel

    Return value:
    m -- direction cosine

    """
    return ((math.sin(math.radians(dec))*math.cos(math.radians(cdec)))-(math.cos(math.radians(dec))*math.sin(math.radians(cdec))*math.cos(math.radians(ra-cra))))/math.radians(incr)

# Find the l direction cosine in a radio image, given an RA and Dec and the
# field centre

def lm_to_radec(ra0,dec0,l,m):
    print('This function should be the inverse of radec_to_lmn, but it is not.  There is likely an error here.')
    sind0=math.sin(dec0)
    cosd0=math.cos(dec0)
    dl=l
    dm=m
    d0=dm*dm*sind0*sind0+dl*dl-2*dm*cosd0*sind0
    sind=math.sqrt(abs(sind0*sind0-d0))
    cosd=math.sqrt(abs(cosd0*cosd0+d0))
    if (sind0>0):
     sind=abs(sind)
    else:
     sind=-abs(sind)

    dec=math.atan2(sind,cosd)

    if l!=0:
     ra=math.atan2(-dl,(cosd0-dm*sind0))+ra0
    else:
     ra=math.atan2((1e-10),(cosd0-dm*sind0))+ra0

# Calculate RA,Dec from l,m and phase center.
# Note: As done in Meqtrees, which seems to differ from l, m functions above.
# Meqtrees equation may have problems, judging from my difficulty fitting a fringe to L4086 data.  Pandey's equation is now used in radec_to_lmn

    return (ra,dec)

def radec_to_lmn(ra0,dec0,ra,dec):
    l=math.cos(dec)*math.sin(ra-ra0)
#    l=-math.cos(dec)*math.sin(ra-ra0)  # old
    sind0=math.sin(dec0)
    if sind0 != 0:
#     m=-(math.cos(ra-ra0)*math.cos(dec)-math.cos(dec0))/math.sin(dec0)   # from sarod -- wrong for L4086 simulations!
     m=math.sin(dec) * math.cos(dec0) - math.cos(dec) * math.sin(dec0) * math.cos(ra - ra0)  # from pandey;  gives same results for casa and cyga
    else:
     m=0
    n = math.sqrt(1-l**2-m**2)
    return (l,m,n)

# Calculate l,m,n from RA,Dec and phase center.
# Note: As done in Meqtrees, which seems to differ slightly from l,m functions above.

def eq_to_gal(ra,dec):
    """Find the Galactic co-ordinates of a source given the equatorial
    co-ordinates

    Keyword arguments:
    (alpha,delta) -- RA, Dec in decimal degrees

    Return value:
    (l,b) -- Galactic longitude and latitude, in decimal degrees

    """

    R = [[-0.054875539726,-0.873437108010,-0.483834985808],[0.494109453312,-0.444829589425,+0.746982251810],[-0.867666135858,-0.198076386122,+0.455983795705]]

    s = [math.cos(math.radians(ra))*math.cos(math.radians(dec)),math.sin(math.radians(ra))*math.cos(math.radians(dec)),math.sin(math.radians(dec))]
    
    sg = []
    
    sg.append(s[0]*R[0][0]+s[1]*R[0][1]+s[2]*R[0][2])
    sg.append(s[0]*R[1][0]+s[1]*R[1][1]+s[2]*R[1][2])
    sg.append(s[0]*R[2][0]+s[1]*R[2][1]+s[2]*R[2][2])

    b = math.degrees(math.asin(sg[2]))
    l = math.degrees(math.atan2(sg[1],sg[0]))

    if l<0:
        l = l+360

    return (l,b)

# Return the Galactic co-ordinates of a point given in equatorial co-ordinates

def gal_to_eq(l,b):
    """Find the Galactic co-ordinates of a source given the equatorial
    co-ordinates

    Keyword arguments:
    (l,b) -- Galactic longitude and latitude, in decimal degrees

    Return value:
    (alpha,delta) -- RA, Dec in decimal degrees

    """

    Rinv = [[-0.054875539692115144, 0.49410945328828509, -0.86766613584223429], [-0.87343710799750596, -0.44482958942460415, -0.19807638609701342], [-0.4838349858324969, 0.74698225182667777, 0.45598379574707293]]

    sg = [math.cos(math.radians(l))*math.cos(math.radians(b)),math.sin(math.radians(l))*math.cos(math.radians(b)),math.sin(math.radians(b))]
    
    s = []
    
    s.append(sg[0]*Rinv[0][0]+sg[1]*Rinv[0][1]+sg[2]*Rinv[0][2])
    s.append(sg[0]*Rinv[1][0]+sg[1]*Rinv[1][1]+sg[2]*Rinv[1][2])
    s.append(sg[0]*Rinv[2][0]+sg[1]*Rinv[2][1]+sg[2]*Rinv[2][2])

    dec = math.degrees(math.asin(s[2]))
    ra = math.degrees(math.atan2(s[1],s[0]))

    if ra<0:
        ra = ra+360

    return (ra,dec)

class CoordSystem(object):
    """A container for constant strings representing different coordinate
    systems."""
    FK4 = "B1950 (FK4)"
    FK5 = "J2000 (FK5)"

def coordsystem(name):
    """Given a string, return a constant from class CoordSystem."""
    mappings = {'j2000': CoordSystem.FK5,
        'fk5': CoordSystem.FK5,
        CoordSystem.FK5.lower(): CoordSystem.FK5,
        'b1950': CoordSystem.FK4,
        'fk4': CoordSystem.FK4,
        CoordSystem.FK4.lower(): CoordSystem.FK4}

    return mappings[name.lower()]

def convert_coordsystem(ra, dec, insys, outsys):
    """Convert RA & dec (given in decimal degrees) between equinoxes"""

    #Uses the wcstools C library: you'll need that in
    #/usr/local/lib/libwcstools.so.1, or, if it's elsewhere, to adjust the path
    #below. (Note that this builds as libwcs.a by default; you'll need to make
    #it shared & change the name to avoid conflict with wcslib. Build all the
    #individual object files with -fPIC, then run:
    #
    #  gcc -shared -Wl,-soname,libwcstools.so.1 -o libwcstools.so.1 *o -lc
    #
    #to make a shared library.)

    import ctypes
    wcstools = ctypes.cdll.LoadLibrary(wcstools_path)
    c_ra, c_dec = ctypes.c_double(ra), ctypes.c_double(dec)
    p_ra, p_dec = ctypes.pointer(c_ra), ctypes.pointer(c_dec)

    if insys == CoordSystem.FK4:
        if outsys == CoordSystem.FK5:
            wcstools.fk425(p_ra, p_dec)
    elif insys == CoordSystem.FK5:
        if outsys == CoordSystem.FK4:
            wcstools.fk524(p_ra, p_dec)

    return c_ra.value, c_dec.value
