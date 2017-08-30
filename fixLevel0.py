Usage="Usage: fixLevel0(<input Level0 FITS file name>, <output FITS name>,\
         checkLevel0=True)"
Version="0.6.3"                 # fixLevel0 version number

# Reads the specified FITS file and corrects the FITS header/table according to
# the scan number and the dataLog file entries for it and the bracketting
# entries.  The corrections are:

# Object name should be the name in the dataLog except that REF observations use
# the corresponding SRC or OTF name.  An updated object name will also get
# updated coordinates and an updated doppler correction.

# The CTYPEn keywords will be corrected based on the COORDCD and VELFRAME
# For GALACTIC coordinates, CTYPE2 = 'GLON--GLS' and CTYPE3 = 'GLAT--GLS'
# All other coordinates have CTYPE2 = 'RA--GLS' and CTYPE3 = 'DEC--GLS'
# CTYPE1 = 'VELO-<frm>' where <frm> is the value of VELFRAME
# (VELFRAME = 'LSR' for objects in the Galaxy and 'HEL' for extragalactic)
# (Class adopts this convention.)

# The coordinates of the catalog object are already stored in CRVAL2 and CRVAL3
# so the keywords RA & DEC are redundant (and misnamed when COORDCD = GALACTIC).
# So RA and DEC keywords are deleted

# The [C II] synthesizer was fixed at 17.584 at the beginning of the flight (we
# think).  When the SCIENCE2 PDU cycle was done, it would have reverted to
# 17.583333 (prior to scan 7868).

# The [C II] rest frequency will be corrected to 1900.5369 GHz (was 1900.5469)

# The [N II] synthesizer was probably set to 13.515 prior to scan 1007, then it
# was set to 13.5151.  For scan 7868, it was reset to 13.5173 and at scan 9498
# it was set to 13.5155.

# The synthesizer settings in the corrected FITS files are as above.
# Using the revised values, change the CRPIX1 column. (This leaves the
# CRVAL1 value in the header to be the catalog velocity of the Object so that
# it will be used as the object velocity in the Class file header, for example.)

# Change the Axial Focus position to reflect the change made prior to scan
# 3082 (Position 39992 changed to -14007).

# Put in the RA & Dec of Venus in its spiral maps.

# Line describing the changes made will be returned.  It will also be placed in
# in a HISTORY entry of the output FITS header.  A new keyword, LEVEL0, contains
# the current version number of fixLevel0 and the date it fixed the FITS file.

# An "x" preface in the FITS file name is removed from the Level 0.6 FITS file

# Receiver header data are corrected.
# Level 0.5 columns for mixers #1 and 2 are swapped for both C+ and N+

# Fix spelling of APPARENT

import os
import sys
from math import *
import numpy as np
import time
import commands
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.time import Time 
import getCat

STALE = 300    # Max time difference (sec) for using bracketting scan's name
STODIR = os.getenv('STODIR')
STODATA = os.getenv('STODATA')
if not STODIR or not STODATA:
    print "Needs environment variables STODIR & STODATA"
    raise LookupError

dFrq = 1.1e9/1024               # Channel spacing in Hz
previous = [0,0,0,0]

# Get the STO2 object catalog
catFile = STODIR + '/STO-arizona/DAC/STO2.cat'
objectDict = getCat.getCat(catFile)

# Read in the ephemeris table
ephFile = STODATA + '/../scripts/ephem.txt'
eph = file(ephFile)
ephText = eph.readlines()
eph.close()
ephObj = []                     # Planet names
ephTab = []                     # Ephemeris tables
ephSiz = []                     # Number of entries in each table
ephSort = []                    # Non-zero means needs sorting
initSize = 5
ind = None

for line in ephText:
    item = line.split()
    if len(item) == 1:
        # Planet name
        try:
            ind = ephObj.index(item[0]) # Index of the that planet
            # Appending to a previous planet's table
        except:
            # Defining a new planet
            ephObj.append(item[0]) # Add the planet name
            ephSiz.append(0)       # Indicate an empty table
            ephSort.append(0)      # Non-zero if tstamps not chronological
            # Table contains  tstamp, RA, Dec, Vel, hPar for initSize entries
            ephTab.append(np.zeros((5, initSize), dtype=np.float))
            ind = len(ephObj) - 1  # Index of the new object
        continue
    if item[0] == 'Time':
        # Column header
        continue
    if len(item) != 5:
        print "Ignored " + line
        continue
    if ind is None:
        print "Ignoring entry with no planet name: " + line
        continue
    if ephSiz[ind] >= ephTab[ind].shape[1]:
        # Resize the table to contain initSize more entries
        tmp = np.zeros((5, ephTab[ind].shape[1] + initSize), dtype=np.float)
        for i in range(5):
            tmp[i][:ephTab[ind].shape[1]] = ephTab[ind][i]
        ephTab[ind] = tmp
        del tmp
    # Append the new entry
    if ephSiz[ind] and float(item[0]) < ephTab[ind][0][ephSiz[ind]-1]:
        ephSort[ind] = 1        # Indicate tstamp not chronological
    for i in range(5):
        ephTab[ind][i][ephSiz[ind]] = float(item[i])
    ephSiz[ind] += 1            # New entry count

for i in range(len(ephSort)):
    if ephSort[i]:
        # Put ephTab[i] into chronological order
        indexArr = ephTab[i][0].argsort()
        tmp = np.zeros((ephTab[i].shape[1], ephTab[i].shape[0]), dtype=np.float)
        unsorted = ephTab[i].transpose() # Now 1st index gets [tstamp, RA, ...] entry
        for j in range(len(indexArr)):
            tmp[j] = unsorted[indexArr[j]]
        ephTab[i] = tmp.transpose()
        del (tmp, unsorted, indexArr)
del ephText

# Find the planet/moon named in the ephemeris stored in ephObj & ephTab lists
# and return (RA, Dec, Velocity, horParallax) interpolated for the timestamp
# specified (sec).  RA & Dec & horParallax are in degrees, Velocity in km/s.
def interpCoords(name, tstamp):
    try:
        ind = ephObj.index(name)
    except:
        print "%s not in ephemeris list" % name
        return (-999., -999., 0., 0.)
    num = ephSiz[ind]
    if tstamp < ephTab[ind][0][0] - STALE or tstamp + STALE > ephTab[ind][0][num-1]:
        print "Time stamp %.0f for %s is outside the ephemeris loaded" % (tstamp, name)
        return (-999., -999., 0., 0.)
    lon = np.interp(tstamp, ephTab[ind][0][:num], ephTab[ind][1][:num])
    lat = np.interp(tstamp, ephTab[ind][0][:num], ephTab[ind][2][:num])
    vel = np.interp(tstamp, ephTab[ind][0][:num], ephTab[ind][3][:num])
    hpar = np.interp(tstamp, ephTab[ind][0][:num], ephTab[ind][4][:num])
    return (lon, lat, vel, hpar)

# Return  RA & Dec corrected for parallax given:
#   obslat: Observatory Latitude (deg)
#   lst: Local Mean Sidereal Time (hours)
#   hpar: Horizontal parallax (deg)
#   alt: Altitude of the observatory (km)
#   ra: Right Ascension (deg)
#   dec: Declination (deg)
def horPar(obslat, lst, hpar, alt, ra, dec):

    a = 6378.1600                       # Equatorial radius of earth (km)
    r = 6378.137 / sin(radians(hpar))   # Geocentric distance to the object

    lstd = 15.0 * lst                   # Sidereal time in degrees
    had = lstd - ra                     # Hour angle in degrees
    har = radians(had)                  # Hour angle in radians
    lstr = radians(lstd)
    obslatr = radians(obslat)
    cos2l = cos(2.0 * obslatr)
    cos4l = cos(4.0 * obslatr)
    cos6l = cos(6.0 * obslatr)
    sin2l = sin(2.0 * obslatr)
    sin4l = sin(4.0 * obslatr)
    
    # Radius Vector from center of the earth to the observatory
    rho = a * (0.998327073 + 0.001676438 * cos2l - 0.000003519 * cos4l \
               + 0.000000008 * cos6l) + alt

    latgc = obslat - 0.192429 * sin2l + 0.000323 * sin4l # Geocentric latitude
    latgc = radians(latgc)
    rar = radians(ra)
    decr = radians(dec)
    cosdec = cos(decr)
    sindec = sin(decr)
    cosra  = cos(rar)
    sinra  = sin(rar)
    coslgc = cos(latgc)
    sinlgc = sin(latgc)
    coslst = cos(lstr)
    sinlst = sin(lstr)
    sinha  = sin(har)
    cosha  = cos(har)

# Use the relations from Smart, Section 120, p. 205.

# Make the RA calculation:

    dr = rho / r

    t1 = sinha * coslgc
    t2 = cosdec - dr * cosha * coslgc
    t3 = dr * t1 / t2
    delrad = degrees(atan(-t3))
    radps  = ra + delrad

    if radps <  0.0:
        radps = radps + 360.
    if radps >= 360.:
        radps = radps - 360.
        
# Make the DEC calculation: 

    haprmd = had - delrad
    coshap = cos(radians(haprmd))

    t4 = coshap * (sindec - dr * sinlgc)
    t5 = cosdec * cosha - dr * coslgc
    t6 = t4 / t5;
    decdps = degrees(atan(t6))

    return (radps, decdps)


# Return the key name in altDict (default: objectDict) corresponding
# to the case insensitive <name>.  If not found, return None
def noCaseKey(name, altDict=objectDict):
    name = name.lower()
    for i in altDict.keys():
        if i and i.lower() == name:
            return i
    return None


# Read the specified dataLog file specified:

# Returns a tuple of two lists: The first is a list of scan numbers in the file
# and the second is a list of the split lines w/o the scan number element.

# Each line begins with mmddyy hh:mm:ss scan#. The timestamp in the first two
# columns is converted into Unix time (sec) which is the first element of the
# second list, the fourth column is the second element, ...

def getDataLog(dataLogName=None):

    if dataLogName == None:
        print "Must specify the dataLog file name"
        return ([], [])

    try:
        dataLog = file(dataLogName, 'r')
    except Exception as e:
        print "Could not open %s\n%s" % (dataLogName, e)
        return ([], [])

    log = dataLog.readlines()
    if len(log) == 0:
        print "File %s is empty" % dataLogName
        return ([], [])

    scanNo = []
    rest = []
    for line in log:
        items = line.split()
        scanNo.append(int(items[2]))
        tstamp = items[0] + " " + items[1]
        items[2] = time.mktime(time.strptime(tstamp, '%m%d%y %H:%M:%S'))
        rest.append(items[2:])
    return(scanNo, rest)

dataLogFile = STODATA + '/logs/dataLog.log'
dataLog = getDataLog(dataLogFile)

# Do the actual work (see documentation above)
def fixLevel0(fitsFile=None, fitsFileOut=None, checkLevel0=True):
    # Unless checkLevel0 is False, will only process Level 0 FITS files (& only once)
    # If fixLevel0 aborts, its return message begins with "Error"
    
    global dataLog
    global objectDict
    global previous
    
    if not fitsFile or not fitsFileOut:
        print Usage
        return 'Error: Must specify both an input and output FITS name'

    if fitsFile == fitsFileOut and not checkLevel0:
        inMode = 'update'
    else:
        inMode = 'readonly'        
        if os.access(fitsFileOut, os.F_OK):
            return 'Error: Output file %s already exists' % fitsFileOut
    
    try:
        if not dataLog:
            return "Error: No dataLog"
    except Exception as (e):
        print e
        return "Error: No dataLog"
        
    if not objectDict:
        return "Error: No objectDict"
        
    try:
        hduList = fits.open(fitsFile, mode=inMode)
    except Exception as (e):
        print "Error opening FITS file %s\n%s" % (fitsFile, e)
        return 'Error: Could not open %s' % fitsFile

    if checkLevel0:
        try:
            unit = hduList[1].columns['DATA'].unit
        except Exception as (e):
            print 
            msg = "Error FITS file %s has no DATA column" % fitsFile
            print msg
            print e
            return msg
        if unit != 'counts':
            msg = 'Error: %s is not a Level 0 FITS file' % fitsFile
            print msg
            return msg
        else:
            try:
                msg = hduList[1].header['level0']
                msg = "%s is a Level 0 file already fixed by %s" % (fitsFile, msg)
                print msg
                return msg
            except:
                pass

    # FITS file needs fixing
    try:
        tbdata = hduList[1].data
    except Exception as e:
        msg = "Error: FITS file %s does not have a binary table" % fitsFile
        print msg
        hduList.close()
        return msg

    # Get scanID number
    try:
        scan = hduList[1].header['scan']
    except Exception as e:
        msg = "Error: FITS file %s does not have a keyword SCAN" % fitsFile
        print msg
        hduList.close()
        return msg

    # Get dataLog index
    try:
        iLog = dataLog[0].index(scan)
    except:
        msg = "Error: Scan number %d (%s) not found in dataLog" % (scan, fitsFile)
        hduList.close()
        return msg

    msg = fitsFile + ':'
    pref = ' fixed'
    tstamp = dataLog[1][iLog][0]
    mode = dataLog[1][iLog][1]
    kind = dataLog[1][iLog][2]
    name = dataLog[1][iLog][-1]
    otRow = dataLog[1][iLog][-2]

    if name.find('_REF') >= 0 or mode == 'REF':
        # REF objects aren't used for doppler correction.
        # Use the name of a bracketting SRC or OTF
        # Get entries for prior scan and next scan
        pDt = nDt = STALE + 1
        for i in range(1, 5):
            pLog = dataLog[0].index(scan+i)
            if dataLog[1][pLog][1] != 'REF':
                pDt = dataLog[1][pLog][0] - tstamp
                break
        for i in range(1, 5):
            nLog = dataLog[0].index(scan-i)
            if dataLog[1][nLog][1] != 'REF':
                nDt = tstamp - dataLog[1][nLog][0]
                break
        if min(pDt, nDt) > STALE:
            msg += ' No alternate name'
            pref = ','
        elif pDt < nDt:
            name = dataLog[1][pLog][-1]
            kind = dataLog[1][pLog][2]
        else:
            name = dataLog[1][nLog][-1]
            kind = dataLog[1][nLog][2]
    lName = name.lower()

    # Get Object name and coordinates as stored in the Level 0.5 
    fitsObj =  hduList[1].header['object']
    fitsLong =  hduList[1].header['crval2']
    fitsLat =  hduList[1].header['crval3'] 
    fitsEqui = hduList[1].header['equinox']
    fitsVel = hduList[1].header['crval1']   # aka VELOCITY
    fitsCoord = hduList[1].header['coordcd']
    fitsFrame = hduList[1].header['velframe']
    fitsDop = hduList[1].header['gon_vel']

    unixTime = hduList[1].header['unixtime']
    gonLong =  hduList[1].header['gon_long']
    gonLat =  hduList[1].header['gon_lat']
    gonAlt = hduList[1].header['altitude']
    parallactic = hduList[1].header['parlactc']

    if kind == 'SPIRAL' or kind == 'BASIC':
        # Try to guess the object name from the SPIRAL/BASIC name
        for i in '3576 etaCar2 etaCar5 etaCar venus moon sky-bssrc'.split():
            if lName.find(i.lower()) != -1:
                if i == '3576':
                    name = 'NGC3576'
                elif i == 'etaCar':
                    name = 'etaCar2' # One spiral map omitted the "2"
                elif i == 'sky-bssrc':
                    name = 'moon'    # Blank sky near the moon
                    lName = 'moon'
                else:
                    name = i
                    break

    planet = False
    # Check for a planet (Only two ephemeris objects observed: Venus & Moon)
    for i in 'venus moon'.split():
        if lName.find(i) != -1:
            name = i            # Got one
            planet = True
            # Get coordinates from ephemeris table
            (objLong, objLat, horParallax, velocity) = interpCoords(name, tstamp)
            if objLong < 0:
                humTime = time.strftime("%d-%b-%Y %T", time.gmtime(tstamp))
                msg += ' Error: No ephemeris for %s @ %s UTC' % (name, humTime)
                return msg
            coordcd = 'APPARENT'
            frame = 'GEO'
            catName = name
            # Get local mean sidereal time (hours)
            timeObj = Time(tstamp, format='unix', \
                           location=(gonLong, gonLat, gonAlt))
            lmst = timeObj.sidereal_time('mean').value
            (pRa, pDec) = horPar(gonLat, lmst, horParallax, \
                                 gonAlt/1000., objLong, objLat)
            dRa =  (pRa - objLong)*3600.
            dDec = (pDec - objLat)*3600.
            break

    if not planet:
        # Case insensitive search for catalog name
        catName = noCaseKey(name, objectDict)
        if catName is None:
            # Not in catalog
            msg += ' Error: %s not in catalog' % name
            return msg
        else:
            (objLong, objLat, coordcd, velocity, frame) = objectDict[catName]

###    print 'catName %s long %.3f lat %.3f coordcd %s velocity %.2f frame %s' \
###                           % (catName, objLong, objLat, coordcd, velocity, frame)
    velocity  *= 1000.          # Convert to m/s

    
    cmd = STODIR + "/STO-arizona/dop/getDop %.0f %.4f %.4f %.4f %.4f %.4f %s %s"\
          % (unixTime, gonLong, gonLat, gonAlt, objLong, objLat, coordcd, frame)
    if coordcd == 'B1950':
        objEqui = 1950.
    elif coordcd == 'J2000':
        objEqui = 2000.
    else:
        objEqui = 0.

    try:
        gonVel = 1000.*float(commands.getoutput(cmd)) # Doppler velocity (m/sec)
    except Exception as e:
        print "Error executing %s\n%s" % (cmd, e)
        msg += ' Error in running getDop'
        return msg

    nII1row = list(tbdata.field('telescop')).index('NII_1')
    nII2row = list(tbdata.field('telescop')).index('NII_2')
    cII1row = list(tbdata.field('telescop')).index('CII_1')
    cII2row = list(tbdata.field('telescop')).index('CII_2')

    # Beam offsets wrt center of the beam diamond pattern on the sky (arcsec)
    offAz = np.array([99.3, 17.5, -99.3, -17.5], dtype=np.float32)
    offEl = np.array([17.5, -99.3, -17.5, 99.3], dtype=np.float32)
    rows = (nII1row, nII2row, None, cII2row) # Corresponding FITS row numbers

    # Make the offsets relative to the fiducial beam
    fidBeam = 3          # Zero-based index for the fiducial ([C II] #2 beam)
    offAz -= offAz[fidBeam]
    offEl -= offEl[fidBeam]
    beamOffID = 'Nominal offsets with C+ #2 being the fiducial'
    
    # Convert to RA & Dec offsets in degrees
    sinp = sin(radians(parallactic))
    cosp = cos(radians(parallactic))
    longOff = [0, 0, 0, 0, 0, 0]
    latOff = [0, 0, 0, 0, 0, 0]
    ra1=[]
    dec1=[]
    
    for i in range(2):
        longOff[i] = (-offAz[i]*cosp - offEl[i]*sinp)/3600. # RA
        latOff[i] = (-offAz[i]*sinp + offEl[i]*cosp)/3600.  # Dec
#        print "i = %d RA-DEC offsets = %.0f %.0f arcsec" % (i, 3600.*longOff[i], 3600.*latOff[i]) 
        if coordcd == 'GALACTIC':
            # Convert offsets to galactic
            cosd = cos(radians(hduList[1].header['udp_dec']))
            cosLat = cos(radians(objLat)) 
            if objLong != previous[0] or objLat != previous[1]:
                # new object.  Need fresh computation from catalog
                cGal = SkyCoord(frame="galactic", l=objLong, b=objLat, unit='deg')
                ra1.append(cGal.fk5.ra.value + longOff[i]/cosd)
                dec1.append(cGal.fk5.dec.value + latOff[i])
                previous = [objLong, objLat, cGal.fk5.ra.value, cGal.fk5.dec.value]
                pref = ' fixed uncached'
            else:
                # catalog entry has not changed.  Use cached value
                ra1.append(previous[2] + longOff[i]/cosd)
                dec1.append(previous[3] + latOff[i])
    if coordcd == 'GALACTIC':
        cOff = SkyCoord(frame='fk5', ra=np.array(ra1), dec=np.array(dec1), unit='deg')
        glongOff = (cOff.galactic.l.value - objLong)*cosLat
        glatOff = cOff.galactic.b.value - objLat
        longOff=glongOff.tolist() + [0,0,0,0]
        latOff=glatOff.tolist() + [0,0,0,0]

    # Now bulk-apply the offsets
    tbdata.field('cdelt2')[:] = longOff
    tbdata.field('cdelt3')[:] = latOff
    
    hduList[1].header['bmoff_id'] = beamOffID
    msg += pref + ' Beam offsets'
    pref = ','
    # Fix [C II] rest frequency
    tbdata.field('restFreq')[cII2row] = 1900.5369e9 # (Hz)

    # Fix [N II] and [C II] synthesizers
    if scan < 1007:
        hduList[1].header['syn_n_ii'] = (-13.515, 'N+ synthesizer freq (<0: set manually)')
    elif scan < 7868:
        hduList[1].header['syn_n_ii'] = (-13.5151, 'N+ synthesizer freq (<0: set manually)')

    elif scan < 9498:
        hduList[1].header['syn_n_ii'] = (-13.5173, 'N+ synthesizer freq (<0: set manually)')
    else:
        hduList[1].header['syn_n_ii'] = (-13.5155, 'N+ synthesizer freq (<0: set manually)')

    msg += pref + ' [N II] syn freq'
                    
    if scan < 7867:
        hduList[1].header['syn_c_ii'] = (-17.584, 'C+ synthesizer freq (<0: set manually)')
    else:
        hduList[1].header['syn_c_ii'] = (-17.583333, 'C+ synthesizer freq (<0: set manually)')
    msg += pref + ' [C II] syn & rest freq'
    
    if scan >= 3082:
        hduList[1].header['axfocpos'] = (-14007, 'Axial focus (encoder units)')
    else:
        hduList[1].header.comments['axfocpos'] = 'Axial focus (encoder units)'

    if fitsObj != catName:
        hduList[1].header['object'] = (catName, 'Catalog object name')
        msg += pref + ' name'
        
    if fitsLong != objLong or fitsLat != objLat or fitsCoord != coordcd or fitsEqui != objEqui:
        hduList[1].header['crval2'] = (objLong, 'Long. of Object in catalog (deg)')
        hduList[1].header['crval3'] = (objLat, 'Lat. of Object in catalog (deg)')
        hduList[1].header['equinox'] = objEqui
        hduList[1].header['coordcd'] = coordcd
        msg += pref + ' coords'
    else:
        # Just fix the comments
        hduList[1].header.comments['crval2'] = 'Long. of Object in catalog (deg)'
        hduList[1].header.comments['crval3'] = 'Lat. of Object in catalog (deg)'

    # Fix CTYPE2 & CTYPE3 according to the coordinate system
    if coordcd == 'GALACTIC':
        hduList[1].header['ctype2'] = 'GLON--GLS'
        hduList[1].header['ctype3'] = 'GLAT--GLS'
    else:
        hduList[1].header['ctype2'] = 'RA---GLS'
        hduList[1].header['ctype3'] = 'DEC--GLS'
    msg += pref + ' ctype2&3'
            
    # Get rid of duplicates/misnamed
    del hduList[1].header['ra']
    del hduList[1].header['dec']
    if fitsVel != velocity or fitsFrame != frame:
        hduList[1].header['crval1'] = hduList[1].header['velocity'] = velocity
        hduList[1].header.comments['crval1'] = "Object %s velocity (m/s)" % frame
        hduList[1].header['velframe'] = frame
        hduList[1].header['ctype1'] = 'VELO-%s' % frame[:3]
        msg += pref + ' velocity'
    if abs(fitsDop - gonVel) > 100:
        hduList[1].header['gon_vel'] = gonVel
        msg += pref + ' doppler'
    hduList[1].header['level0'] = ('Version %s' % Version, \
                                   'Applied %s' % time.strftime('%d %b %Y', time.gmtime()))

    # Fix up the per-pixel header (pixel 1 and 2 receiver data swapped)
    tbdata.field('BiasVolt')[nII1row], tbdata.field('BiasVolt')[nII2row] = tbdata.field('BiasVolt')[nII2row], tbdata.field('BiasVolt')[nII1row]
    tbdata.field('BiasVolt')[cII1row], tbdata.field('BiasVolt')[cII2row] = tbdata.field('BiasVolt')[cII2row], tbdata.field('BiasVolt')[cII1row]

    tbdata.field('BiasCurr')[nII1row], tbdata.field('BiasCurr')[nII2row] = tbdata.field('BiasCurr')[nII2row], tbdata.field('BiasCurr')[nII1row]
    tbdata.field('BiasCurr')[cII1row], tbdata.field('BiasCurr')[cII2row] = tbdata.field('BiasCurr')[cII2row], tbdata.field('BiasCurr')[cII1row]

    tbdata.field('IFPower')[nII1row], tbdata.field('IFPower')[nII2row] = tbdata.field('IFPower')[nII2row], tbdata.field('IFPower')[nII1row]
    tbdata.field('IFPower')[cII1row], tbdata.field('IFPower')[cII2row] = tbdata.field('IFPower')[cII2row], tbdata.field('IFPower')[cII1row]

    tbdata.field('BiasPot')[nII1row], tbdata.field('BiasPot')[nII2row] = tbdata.field('BiasPot')[nII2row], tbdata.field('BiasPot')[nII1row]
    tbdata.field('BiasPot')[cII1row], tbdata.field('BiasPot')[cII2row] = tbdata.field('BiasPot')[cII2row], tbdata.field('BiasPot')[cII1row]

    msg += pref + ' rx hdr'
    
    # Calculate CRPIX1
    niiSyn = abs(hduList[1].header['syn_n_ii'])
    ciiSyn = abs(hduList[1].header['syn_c_ii'])

    nii_crpix1 = 1 + (tbdata.field('restFreq')[nII1row] - 1.e9*(108*niiSyn+1.1))/dFrq + \
            (hduList[1].header['crval1'] + gonVel)/tbdata.field('cdelt1')[nII1row]
    cii_crpix1 = 1 + (tbdata.field('restFreq')[cII2row] - 1.e9*(108*ciiSyn+1.1))/dFrq + \
            (hduList[1].header['crval1'] + gonVel)/tbdata.field('cdelt1')[cII2row]
    tbdata.field('crpix1')[nII1row] = tbdata.field('crpix1')[nII2row] = nii_crpix1
    tbdata.field('crpix1')[cII2row] = cii_crpix1

    i = msg.find(': ')          # Skip past file name
    hisMsg = Version + msg[i+1:]
    hduList[1].header['history'] = hisMsg
                                   
    # Choose the [N II] rows and the CII_2 mixer row
    mask = (tbdata['line'] == '[N II]') + (tbdata['telescop'] == 'CII_2')
    newData = tbdata[mask]      # Selected rows
    hdu = fits.BinTableHDU(data=newData, header=hduList[1].header)
    newList = fits.HDUList([hduList[0], hdu])
    if fitsFileOut == fitsFile:
        newList.flush()
    else:
        newList.writeto(fitsFileOut)

    # Clean up
    newList.close()
    hduList.close()
    # Make sure no mmap's to the FITS files remain open (probably not needed):
    del (newList, hdu, newData, tbdata, hduList)
    return msg


# Call fixLevel0 for scans <begScan> through <endScan>.  The input scans will be
# in subdirectories of <inPath> named nnnnn where nnnnn is the 5-digit scanID
# number (with leading zeros if needed).  The output will be in similarly named
# subdirectories of <outPath>.  A <logFile> will be appended with the message
# returned by fixLevel0 for each FITS file.  Unless checkLevel0 is False,
# fixLevel0 will ignore any FITS files found that are not Level 0 or have
# already been processed by fixLevel0
def fixLoop(inPath, outPath, begScan, endScan, logFile='fixLevel0.log', checkLevel0=True):
    Usage = 'Usage: fixLoop(inPath, outPath, begScan, endScan, logFile="fixLevel0.log", checkLevel0=True)'

    if not os.path.isdir(inPath) or not os.access(inPath, os.R_OK):
        print "%s is not a directory or has no read access" % inPath
        print Usage
        return

    if not os.path.isdir(outPath) or not os.access(outPath, os.W_OK):
        print "%s is not a directory or has no write access" % outPath 
        print Usage
        return

    if inPath[-1] != '/':
        inPath += '/'
    if outPath[-1] != '/':
        outPath += '/'

    try:
        begScan = int(begScan)
        endScan = int(endScan)
    except:
        print "begScan and/or endScan are not numbers: %s & %s" % (begScan, endScan)
        print Usage
        return
    if begScan > endScan or begScan < 0:
        print "Illegal scan numbers for begScan & endScan: %d, %d" % (begScan, endScan)
        print Usage
        return

    try:
        log = file(logFile, 'a')
    except Exception as e:
        print "Could not write to log file: %s" % logFile
        print Usage
        return
    log.write("Version %s, outPath: %s\n" % (Version, outPath))

    numFixed = 0
    for scan in range(begScan, endScan+1):
        print '\rProcessing scan %d...' % (scan)
        zscan = "%05d" % scan
        # Get the file names in sequence number order from the ls command
        (status, fitsNames) = commands.getstatusoutput( \
                              'ls %s%s/*.fits | sort -t _ -k %d -n' % \
                              (inPath, zscan, inPath.count('_') + 2))
        fitsNames = fitsNames.split()
        if status != 0 or fitsNames[0] == 'ls:':
            msg = "No files found for scan #%d\n" % scan
            log.write(msg)
            print '\n' + msg
        else:
            outSubdir = "%s/%s" % (outPath, zscan)
            if not os.access(outSubdir, os.F_OK):
                os.mkdir(outSubdir) # Create the output subdirectory
            for fits in fitsNames:
                n = fits.rfind('/') + 1 # FITS file name in its subdirectory
                fname = "%s/%s" % (zscan, fits[n:])
                if fits[n] == 'x':
                    oFname = "%s/%s" % (zscan, fits[n+1:]) # Discard leading "x"
                else:
                    oFname = fname
                msg = fixLevel0(inPath + fname, outPath + oFname, checkLevel0)
                log.write(msg + '\n')
                n = msg.find(': ') + 2
                if msg[:5].lower() != 'error' and msg[n:n+5].lower() != 'error':
                    numFixed += 1
                else:
                    print '\n' + msg
                sys.stdout.flush()

    log.close()
    print "\nFixed %d files" % numFixed
