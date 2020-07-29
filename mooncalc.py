""" Library to compute moon position, azimuth, elevation, rise/set times 
    and phases.


    Heavily based on the suncalc.js javascript library by Vladimir
    Agafonkin at https://github.com/mourner/suncalc, but missing
    sun-only calculations.  Minimal dependencies so should be able to
    run on micropython as well.
"""

import time



try: # make things work both in micropython and regular python:
    from collections import namedtuple
except:
    from ucollections import namedtuple

from math import pi, sin, cos, tan, asin, atan2, acos, sqrt
RAD = pi/180

# sun calculations are based on https://www.aa.quae.nl/en/reken/zonpositie.html formulas

# date/time constants and conversions
DAYS_SEC = 60 * 60 * 24

# Julian dates: see
# https://docs.kde.org/trunk5/en/extragear-edu/kstars/ai-julianday.html
# the formulas below are designed to work in the neighborhood of the
# year 2000 a.D., so we'll convert to an epoch starting on 1/1/2000
J1970 = 2440588
J2000 = 2451545

def to_julian(date):
    return date / DAYS_SEC - 0.5 + J1970

def from_julian(j):
    return time.localtime((j + 0.5 - J1970) * DAYS_SEC)

# convert date to (days since 2000). Including fractional days.
def to_days(date):
    return to_julian(date) - J2000

# general calculations for position
E = RAD * 23.4397 # obliquity of the Earth

def right_ascension(l, b):
    return atan2(sin(l) * cos(E) - tan(b) * sin(E), cos(l))
 

def declination(l, b):
    return asin(sin(b) * cos(E) + cos(b) * sin(E) * sin(l))


def azimuth(H, phi, dec):
    return atan2(sin(H), cos(H) * sin(phi) - tan(dec) * cos(phi))


def altitude(H, phi, dec):
    return asin(sin(phi) * sin(dec) + cos(phi) * cos(dec) * cos(H))


def sidereal_time(d, lw):
    return RAD * (280.16 + 360.9856235 * d) - lw


def astro_refraction(h):
    if h < 0: # the following formula works for positive altitudes only.
        h = 0 # if h = -0.08901179 a div/0 would occur.

    # formula 16.4 of "Astronomical Algorithms" 2nd edition
    # by Jean Meeus (Willmann-Bell, Richmond) 1998.
    #
    # 1.02 / tan(h + 10.26 / (h + 5.10)) h in degrees,
    # result in arc minutes -> converted to RAD:
    return 0.0002967 / tan(h + 0.00312536 / (h + 0.08901179))


#
# general sun calculations
# I left out a bunch of functions (sunrise etc) as they are not pertinent to
# the moon specifically. Ping me and I can add them
#

def solar_mean_anomaly(d):
    """
    d: fractional days since 1/1/2000
    """
    return RAD * (357.5291 + 0.98560028 * d)

def ecliptic_longitude(M):
    C = RAD * (1.9148 * sin(M) + 0.02 * sin(2 * M)
               + 0.0003 * sin(3 * M)) # equation of center
    P = RAD * 102.9372                # perihelion of the Earth
    return M + C + P + pi

SunCoordsStruct = namedtuple('SunCoordsStruct', 'dec ra')

def sun_coords(d):
    """
    d: fractional days since 1/1/2000
    """
    M = solar_mean_anomaly(d)
    L = ecliptic_longitude(M)
    
    return SunCoordsStruct(
        dec = declination(L, 0),
        ra = right_ascension(L, 0))

# adds a custom time to the times config
def add_time (angle, rise_name, set_name):
    TIMES.append([angle, rise_name, set_name])

    
# calculations for sun times
J0 = 0.0009

def julian_cycle(d, lw): 
    return round(d - J0 - lw / (2 * pi))

def approx_transit(Ht, lw, n):
    return J0 + (Ht + lw) / (2 * pi) + n

def solar_transit(ds, M, L):
    return J2000 + ds + 0.0053 * sin(M) - 0.0069 * sin(2 * L)

def hour_angle(h, phi, d):
    return acos((sin(h) - sin(phi) * sin(d)) / (cos(phi) * cos(d)))

#
# moon calculations, based on
# http://aa.quae.nl/en/reken/hemelpositie.html formulas
#

MoonCoordsStruct = namedtuple('MoonCoordsStruct', 'ra dec dist')


def moon_coords(d): # geocentric ecliptic coordinates of the moon
    """
    d: fractional days since 1/1/2000
    """
    L = RAD * (218.316 + 13.176396 * d) # ecliptic longitude
    M = RAD * (134.963 + 13.064993 * d) # mean anomaly
    F = RAD * (93.272 + 13.229350 * d)  # mean distance
    
    l  = L + RAD * 6.289 * sin(M)       # longitude
    b  = RAD * 5.128 * sin(F)           # latitude
    dt = 385001 - 20905 * cos(M)        # distance to the moon in km
    
    return MoonCoordsStruct(ra=right_ascension(l, b),
                            dec=declination(l, b),
                            dist=dt)

MoonPositionStruct=namedtuple('MoonPositionStruct',
                              'azimuth altitude distance parallactic_angle')

def get_moon_position(date, lat, lng):
    """
    date in seconds since 1/1/1970
    lat from -90 to +90
    lng from -180 to +180
    """
    lw  = RAD * -lng
    phi = RAD * lat
    d   = to_days(date)

    c = moon_coords(d)
    H = sidereal_time(d, lw) - c.ra
    h = altitude(H, phi, c.dec)
    # formula 14.1 of "Astronomical Algorithms" 2nd edition by Jean Meeus
    # (Willmann-Bell, Richmond) 1998.
    pa = atan2(sin(H), tan(phi) * cos(c.dec) - sin(c.dec) * cos(H))

    h = h + astro_refraction(h)  # altitude correction for refraction

    return MoonPositionStruct(
        azimuth=azimuth(H, phi, c.dec),
        altitude=h,
        distance=c.dist,
        parallactic_angle=pa)


# calculations for illumination parameters of the moon,
# based on http://idlastro.gsfc.nasa.gov/ftp/pro/astro/mphase.pro
# formulas and Chapter 48 of "Astronomical Algorithms" 2nd edition
# by Jean Meeus (Willmann-Bell, Richmond) 1998.

MoonIlluminationStruct = namedtuple('MoonIlluminationStruct',
                                    'fraction phase angle')


# TODO: either make adjustment in d+=TIME_EPOCH_DIFF or pass that in from caller
def get_moon_illumination(date):
    d=to_days(date)

    s = sun_coords(d)
    m = moon_coords(d)

    sdist = 149598000 # distance from Earth to Sun in km
    
    phi = acos(sin(s.dec) * sin(m.dec)
               + cos(s.dec) * cos(m.dec) * cos(s.ra - m.ra))
    inc = atan2(sdist * sin(phi), m.dist - sdist * cos(phi))
    angle = atan2(cos(s.dec) * sin(s.ra - m.ra),
                  sin(s.dec) * cos(m.dec) - cos(s.dec)
                  * sin(m.dec) * cos(s.ra - m.ra))
  
    return MoonIlluminationStruct(
        fraction=(1 + cos(inc)) / 2,
        phase=0.5 + 0.5 * inc * (-1 if angle < 0 else 1) / pi,
        angle=angle)


def hours_later(date, h):
    return date + 60 * 60 * h


#
# calculations for moon rise/set times are based on
# http://www.stargazing.net/kepler/moonrise.html article
#

MoonTimesStruct=namedtuple('MoonTimesStruct',
                           'rise set, alwaysup alwaysdown')

def get_moon_times(date, lat, lng, inUTC=True):
    """ date is seconds since epoch
    lat from -90 to +90,
    lng from -180 to 180
    """
    # We look at the whole day, in 2h intervalsm so start at hour 0:
    t=list(time.localtime(date + (time.timezone if inUTC else 0)))
    t[3:6]=[0,0,0]

    t = time.mktime(tuple(t)) - (time.timezone if inUTC else 0)

    hc = 0.133 * RAD
    h0 = get_moon_position(t, lat, lng).altitude - hc
    h1 = h2 = rise = set = a = b = xe = ye = d = roots = x1 = x2 = dx = 0

    # go in 2-hour chunks, each time seeing if a 3-point quadratic curve
    # crosses zero (which means rise or set)
    for i in [2*i + 1 for i in range(12)]:
        h1 = get_moon_position(hours_later(t, i), lat, lng).altitude - hc
        h2 = get_moon_position(hours_later(t, i + 1), lat, lng).altitude - hc

        a = (h0 + h2) / 2 - h1
        b = (h2 - h0) / 2
        xe = -b / (2 * a)
        ye = (a * xe + b) * xe + h1
        d = b * b - 4 * a * h1
        roots = 0
        
        if (d >= 0):
            dx = sqrt(d) / (abs(a) * 2)
            x1 = xe - dx
            x2 = xe + dx
            if abs(x1) <= 1: roots+=1
            if abs(x2) <= 1: roots+=1
            if x1 < -1: x1 = x2

        if roots == 1:
            if (h0 < 0): rise = i + x1
            else: set = i + x1

        elif roots == 2:
            rise = i + (x2 if ye < 0 else x1)
            set = i + (x1 if ye < 0 else x2)

        if rise and set: break

        h0 = h2

    return MoonTimesStruct(hours_later(t, rise) if rise else -1,
                           hours_later(t, set) if set else -1,
                           not rise and not set & ye > 0,
                           not rise and not set and ye <= 0)
