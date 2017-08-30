#LATEST update 23 June 2017
# %load 'STO2_v2.py'
import os
import numpy as np
import glob
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from scipy import signal
from scipy.signal import butter, filtfilt
from pylab import *
import matplotlib.mlab as mlab
import warnings
import math



warnings.filterwarnings("ignore")
 
# convenience routine to load data from file
print("readSTO2")
def readSTO2(ifile, lin, retdd=False, retcl=False, verbose=False, trim=None, rclip=None, tnorm=True, badpix=None):
    """
    Convenience routine for class rsto2 to read STO-2 Level 0.5 data of a single line
    Input:
        ifile:      data file name (full path might be required)
        lin:        desired Level 0.5 line: 0: NII_1; 1: NII_2; 2: CII_1; 3: CII_2; 4: OI_FFT; 5: OI_Freq_LK 
                           Level 0.6 line: 0: NII_1; 1: NII_2; 2: CII_2 (other, not observed lines were removed!)
        retdd:      if True, raw data table dd1 is returned (default: False)
        retcl:      if True, class object is returned (default: False)
        verbose:    print some info to screen
        trim:       trim both ends of the velocity and intensity arrays by the given number of 
                    pixels (default: None) Note: better use rclip!
        rclip:      replaces the given number of pixels at both ends of the velocity and intensity arrays by the value 
                    of the next adjoining pixels (default: None)
        tnorm:      flag to return spectrum divided by integration time: spectrum/tint
        badpix:     pixels with bad values to be replaced by interpolated value from neighboring pixels
    Output:
        separate arrays for velocity, and intensity, pos (astropy SkyCoord object), FITS header, and (optional) the raw data table
    """
    rs = rsto2(ifile, verbose=False, trim=None, rclip=rclip, tnorm=tnorm, badpix=None)
    vv, spec = rs.getData(lin)
    
    if retcl:
        if retdd: return vv, spec, rs.getPos(lin), rs.getHeader(), rs.getRawData(), rs
        else: return vv, spec, rs.getPos(lin), rs.getHeader(), rs
    else:
        if retdd: return vv, spec, rs.getPos(lin), rs.getHeader(), rs.getRawData()
        else: return vv, spec, rs.getPos(lin), rs.getHeader()



class rsto2:
    
    def __init__(self, ifile, trim=None, rclip=None, verbose=False, tnorm=True, badpix=None):
        '''
        reading STO-2 Level 0 and 1 data
        the difference is that some entries in the data table dd (see retdd=True) are not
        available in the Level 0 data, e.g. Tsys
        
        ifile:       data file name (full path might be required)
        lin:         desired line: 0: NII_1; 1: NII_2; 2: CII_1; 3: CII_2; 4: OI_FFT; 5: OI_Freq_LK
                     available are only 0, 1, and 3 
        verbose:     printing some info  (default: False)
        trim:        trim both ends of the velocity and intensity arrays by the given number of 
                     pixels (default: None) Note: better use rclip!
        rclip:       replaces the given number of pixels at both ends of the velocity and intensity arrays by the value 
                     of the next adjoining pixels (default: None)
        tnorm:       flag to return spectrum divided by integration time: spectrum/tint
        badpix:      pixels with bad values to be replaced by interpolated value from neighboring pixels
        '''
        if rclip==0: rclip=None
        if trim==0: trim = None
        
        
        hl = fits.open(ifile)
        self.hd = hl[0].header
        self.hd1 = hl[1].header
        self.dd1 = hl[1].data
        hl.close()
        off2 = np.array(self.dd1['CDELT2'], dtype=np.float)
        off3 = np.array(self.dd1['CDELT3'], dtype=np.float)
        RA =np.float(self.hd1['UDP_RA'])*u.deg
        DEC=np.float(self.hd1['UDP_DEC'])*u.deg
        self.pos = [SkyCoord(RA + off2[0]*u.deg/cos(radians(DEC)) , DEC + off3[0]*u.deg, frame='icrs'),
                    SkyCoord(RA + off2[1]*u.deg/cos(radians(DEC)) , DEC + off3[1]*u.deg, frame='icrs'),
                    SkyCoord(RA + off2[2]*u.deg/cos(radians(DEC)) , DEC + off3[2]*u.deg, frame='icrs')]
        self.spec = self.dd1['DATA']
        n_pix = self.dd1['MAXIS1'][0]
        n_row = self.hd1['NAXIS2']
        self.vv = np.ndarray((n_row,n_pix), dtype=np.float)
        for i in range(n_row):
            n_pixl = self.dd1['MAXIS1'][i]
            self.vv[i,:n_pixl] = (np.float(self.hd1['CRVAL1']) + (1 + np.arange(n_pixl) - self.dd1['CRPIX1'][i]) * self.dd1['CDELT1'][i])
            
        self.vv/= 1000.0
        self.hd1['TUNIT3'] = 'km/s' 

    
        if tnorm:
            tint = np.float(self.hd1['OBSTIME'])
            self.spec = self.spec/tint
            self.hd1['TUNIT24'] = 'counts/sec'

            
#         if verbose: 
#             #     # 3076  4014   REF    33    CII_2   2.226  322.99963    1.90015  -50.0   36.8  17.585350 -0.00370 -57.09700  2760264
#             print(' scan obsid  Type ot_row   line  obstime     l          b        vLSR     v0    syn_CII biasvolt biascurr totpower')
#             print('%5s %5s %5s %5s %8s %7.3f %10.5f %10.5f %6.1f %6.1f %10.6f %8.5f %8.3f %8.0f'%(
#                 hd1['SCAN'], hd1['OBSID'], hd1['TYPE'], hd1['OT_ROW'], dd1['TELESCOP'], 
#                 np.float(hd1['OBSTIME']), pos.galactic.l.deg, pos.galactic.b.deg, 
#                 np.float(hd1['VELOCITY'])/1000., vv[0], np.float(hd1['SYN_C_II']),
#                 np.float(dd1['BIASVOLT'][lin]),np.float(dd1['BIASCURR'][lin]),np.float(dd1['TOTPOWER'][lin])))
    
        if trim!=None:
            # cut away the outmost pixels - but better to use rclip below!
            vv = vv[trim:-trim]
            spec = spec[trim:-trim]
    
        if rclip!=None:
            # replace the values of the outermost pixels [0,rclip-1] and [-rclip,last_pixel] with the value 
            # of the first inside pixel [rclip] and [-rclip-1] respectively
            # to remove edge effects, but keep the number of pixels
            if verbose: print('STO2_v1: clipping edge values....', rclip)
            for rcl in range(rclip):
                self.spec[:,rcl]   = self.spec[:,rclip]
                self.spec[:,-rcl-1] = self.spec[:,-rclip]
        
        if badpix!=None:
            self.badpix = badpix
            self.spec = repPix(self.spec, badpix)
            
        
        
    def getData(self, lin):
        """
            Function to retrieve STO2 spectra for a single line:
            Input:
                lin: line number (0: [NII] pixel 1, 1: [NII] pixel 2, and 2: [CII])
            Output:
                velocity array, intensity array
        """
        return self.vv[lin,:], self.spec[lin,:]
    
    def getPos(self, lin):
        """
            return position as astropy SkyCoord object
        """
        return self.pos[lin]
    
    def getHeader(self):
        """
            return full FITS header of observation
        """
        return self.hd1
    
    def getRawData(self):
        """
            return raw data array for all lines
            Keywords: MAXIS1, CRPIX1, CDELT1, CDELT2, CDELT3, BiasVolt, BiasCurr, TotPower, IFPower, 
                      BiasPot, PIDstate, LNABias, LNACurr, BeamFac, TSYS, Trx, OutRange, IntTime,
                      NumOff, Telescop, Line, RESTFREQ, IMAGFREQ, DATA
        """
        return self.dd1
    
    def getTint(self):
        """
            return the integration time for this observation
        """
        return np.float(self.hd1['OBSTIME'])
    
    def getSpecUnit(self, apflag=True):
        """
            return unit of spectrum as string
            apflag: convert string to astropy units (default: True)            
        """
        spunit = self.hd1['TUNIT24'].replace('counts','count').replace('sec','s')
        print(spunit)
        if apflag: spunit = u.Unit(spunit) 
        return spunit
    
    def getVeloUnit(self, apflag=True):
        """
            return units of velocity axis as string
            apflag: convert string to astropy units (default: True)
        """
        vunit = self.hd1['TUNIT3']
        if apflag: vunit = u.Unit(vunit) 
        return vunit


##################################################################################
# misc. functions to manipulate spectra
##################################################################################

def writeSTO2(ofile, dat, header):
    """
    Writing STO-2 style Level 1 FITS files that should be CLASS readable
    
    ofile:  name of output file; should end in .fits
    dat:    data structure for output data; updated version of the Level 0 data structure
    header: header for file; should be updated version of Level 0 header
    
    with respect to Lev ):
    updated header: OFFRA, OFFDEC, TCAL
    added header:   CALSCAN1, OFFSCAN1, CALSCAN2, OFFSCAN2, SEQRASTER, HOTSEQ1, HOTSEQ2
    updated bin table: CDELT2, CDELT3, TSYS,  INTTIME (->0 for OI), NUMOFF, DATA (now calibrated) 
    """
    fits.writeto(ofile, dat, header)



def cleanSpec(dat, thresh = 500., verbose=False, getReport=False):
    """
    Remove spikes from a spectrum
    Better use cleanSpecPeak() below!
    """
    
    ndat = dat.copy()
    fdat = butter_lowpass_filtfilt(ndat, 100, 10000)

    # try to clip the bad data    
    Ddat = ndat - fdat
    tag = mlab.find(np.abs(Ddat)<thresh)
    if verbose: print('median ndat/fdat/Ddat: ', np.median(ndat), np.median(fdat), np.median(Ddat))
    if verbose: print('tag: ', tag.tolist())
    
    report = []
    if tag.shape[0]>0:
        J = find(diff(tag)>1)
        nn = ndat.shape[0]
        for K in split(tag, J+1):
            # this is the array of data to be changed!
            # we are assuming isolated spikes and not spike city
            # the ends can be trimmed!
            if verbose: print('affected pixels: ', K)
            if ((K.min()>5)&(K.max()<nn-5)&(np.mean(Ddat[K])>thresh*0.7)):
                lw = np.arange(K.min()-4,K.min(),1)
                up = np.arange(K.max()+1,K.max()+5,1)
                updw = np.hstack([lw,up])
                # interpolate
                it = np.interp(K,updw,dat[updw])
                ndat[K] = it
                for kk in K:
                    rstr = 'cleaned pixel: %i %12.4f %14.2f'%(kk, ndat[kk], dat[kk])
                    report.append(rstr)
                if verbose: print('updw:  ', updw, dat[updw])
                if verbose: print('split: ', K, ndat[K], Ddat[K], it)

    if getReport: return ndat, report
    else: return ndat
    

def cleanSpecPeak(dat, thresh = 1E9, verbose=False, getReport=False, bflag=True, boff=1, vv=None):
    """
    Remove spikes from a spectrum
    
    getReport: return a short message about what was cleaned.
    bflag: also clean the next pixel around the spike
    vv: provide velocity information (velocity of affected pixels returned in report)
    """
    
    ndat = dat.copy()

    # try to clip the bad data    
    tag = mlab.find(np.abs(ndat)>thresh)
    if verbose: print('tag: ', tag.tolist())
    
    # correcting neighboring pixels
    if bflag: 
        brange = np.arange(boff) + 1
    
    report = np.array([(0., 0., '', -1, 0.0, 0.0, 0.0)], dtype=[('scan', 'i4'),('obsid', 'i4'),('type', 'S10'),('pixel', 'i4'),('nval', 'f8'), ('val', 'f8'), ('velo', 'f4')])
    nrep = 0
    
    if tag.shape[0]>0:
        J = find(diff(tag)>1)
        nn = ndat.shape[0]
        for K in split(tag, J+1):
            # this is the array of data to be changed!
            # we are assuming isolated spikes and not spike city
            # the ends can be trimmed!
            if verbose: print('affected pixels: ', K)
            if ((K.min()>7)&(K.max()<nn-7)&(np.mean(ndat[K])>thresh*0.9)):
                if bflag:
                    K = np.hstack((K.min()-brange,K))
                    K = np.hstack((K,K.max()+brange))
                lw = np.arange(K.min()-4,K.min(),1)
                up = np.arange(K.max()+1,K.max()+5,1)
                updw = np.hstack([lw,up])
                # interpolate
                it = np.interp(K,updw,dat[updw])
                if verbose: print('updw:  ', updw, dat[updw])
                if verbose: print('split: ', K, ndat[K], it)
                ndat[K] = it
                for kk in K:
                    if vv==None:
                        #rstr = 'cleaned pixel: %i %12.4f %14.2f %7.2f'%(kk, ndat[kk], dat[kk], 9999.0)
                        rsa = np.array([(0., 0., '', kk, ndat[kk], dat[kk], 9999.0)], dtype=[('scan', 'i4'),('obsid', 'i4'),('type', 'S10'),('pixel', 'i4'),('nval', 'f8'), ('val', 'f8'), ('velo', 'f4')])
                        nrep += 1
                    else:
                        #rstr = 'cleaned pixel: %i %12.4f %14.2f %7.2f'%(kk, ndat[kk], dat[kk], vv[kk])
                        rsa = np.array([(0., 0., '', kk, ndat[kk], dat[kk], vv[kk])], dtype=[('scan', 'i4'),('obsid', 'i4'),('type', 'S10'),('pixel', 'i4'),('nval', 'f8'), ('val', 'f8'), ('velo', 'f4')])
                        nrep += 1
                    report = np.vstack((report,rsa))

    if nrep > 0: report = report[1:]
    else: report=None
    
    if getReport: return ndat, report
    else: return ndat



def removeSpike(idat, thresh):
    """
    simple spike removal - experimental!
    thresh: threshold to determine pixels values to be replaced with median

    """
    dat = idat.copy()
    med = np.median(dat)
    sel = np.where(dat>med+thresh)
    if len(sel[0])>0: dat[sel[0]] = med
    sel = np.where(dat<med-thresh)        
    if len(sel[0])>0: dat[sel[0]] = med
    
    return dat




def repPix(dat, rpix):
    """
    replacing pixel values with interpolated values from neighboring pixels
    rpix: array of pixels to be "repaired"
    """
    nn = dat.shape[0]
    ndat = dat.copy()
    
    mask = np.ones(nn, dtype=bool)
    mask[rpix] = False

    # interpolate
    it = np.interp(rpix,np.arange(nn)[mask],dat[mask])
    ndat[rpix] = it
    
    return ndat




def resample(vv, spec, step, verbose=False):
    """
    resampling the velocity axis (vv) and the intensity (spec) in increments of channels (step) by 
    averaging the velocity channels (average is center of new channel) and
    averaging the intensity channels (average to preserve the integrated intensity
    
    Input:
      vv: 1-D array of velocities v (with constant delta v)
      spec: 1-D array of intensities spec (unit does not matter)
      step: number of channels averaged for resampled spectrum
      
    Output: nvv, nspec
      nvv: new 1-D velocity array
      nspec: new 1-D spectrum
    """
    step = int(step)
    nvv = np.zeros(int(spec.shape[0]//step))
    nspec = np.zeros(int(spec.shape[0]//step))
    
    if verbose: print('', step, spec.shape[0]-step)
    for i in np.arange(step,spec.shape[0]-step, step):
        if verbose: print(i//step, i-step, i)
        nvv[i//step] = np.mean(vv[i-step:i])
        nspec[i//step] = np.mean(spec[i-step:i])
        
    return nvv[1:-1], nspec[1:-1]


def decontamSpec(repspec, cleanspec, xx, reprange, repborderwidth, gap=0.):
    """
       Function to remove emission from spectrum by replacing from rescaled, clean spectrum
       input:
            repspec: 1-d spectrum to be repaired
            cleanspec: 1-d clean spectrum
            xx: x-axis valid for both spectra (can be anything, velocity, frequency, or wavelength
            reprange: xx-range to be repaired
            gap:      unused xx-range width between reprange and borderwidth
            repborderwidth: xx-range outside of reprange used for scaling
        output:
            1-d emission-removed spectrum
    """
    ruef_sm = smoothData(repspec)
    huf_sm = smoothData(cleanspec)
    # range for replacement is between -40 and 5 km/s
    # ruef needs to be repaired
    # use huf for repair
    if reprange==None: reprange = [-40.,0.]
    if repborderwidth==None: repborderwidth = 30.
    vsel = np.where(((xx>reprange[0]-gap-repborderwidth)&(xx<reprange[0]-gap))|((xx>reprange[1]+gap)&(xx<reprange[1]+gap+repborderwidth)))[0]
    druef = ruef_sm[vsel]
    dhuf = huf_sm[vsel]
    fxx = xx[vsel]
    isel = np.where((xx>reprange[0])&(xx<reprange[1]))[0]
    nxx = xx[isel]
    rat = druef/dhuf
    # linear interpolation
    #f = np.interp(nxx,fxx,rat)
    # non-linear interpolation
    from scipy.interpolate import interp1d
    f2 = interp1d(fxx, rat, kind='slinear')
    f = f2(nxx)
    newspec = repspec.copy()
    newspec[isel] = cleanspec[isel] * f
    
    return newspec



def smoothData(dat, cutoff=100, fs=10000):
    """
        Butterworth filter data smoothing
        
        for more info see notebook: STO2_cleaning_data
        or scipy documentation
    """
    return butter_lowpass_filtfilt(dat, cutoff, fs)


def butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a

def butter_lowpass_filtfilt(data, cutoff, fs, order=5):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = filtfilt(b, a, data)
    return y



def getScanList():
    converters = {6: lambda s: float(s or 0)}
    ifile = '/rdat/Projects/STO-2/Data/data_new_all_v3.txt'
    dat = np.genfromtxt(ifile, dtype=None, names=['scan','obsid','target','type','ot_row','obstime','l','b',
                                                  'l_off','b_off','date','synCII','biasv','biasc','tpower'],
        delimiter=[6,7,16,8,7,12,12,11,9,10,25,11,9,11,9], autostrip=True, skip_header=2)

    return dat



def mkRep(bp, vel, spe, hd, typ):
    """
    bp: bad pixels for report
    vel: velocity array
    spe: spectrum
    hd:  header
    typ: spectrum type
    """
    rep0 = np.array([(0., 0., '', -1, 0.0, 0.0, 0.0)], dtype=[('scan', 'i4'),('obsid', 'i4'),('type', 'S10'),('pixel', 'i4'),('nval', 'f8'), ('val', 'f8'), ('velo', 'f4')])

    # create first report row
    if len(bp)>0:
        rep = [np.array([(hd['scan'], hd['obsid'], typ, bp[0], spe[bp[0]], 0.0, vel[bp[0]])], 
              dtype=[('scan', 'i4'),('obsid', 'i4'),('type', 'S10'),('pixel', 'i4'),('nval', 'f8'), ('val', 'f8'), ('velo', 'f4')])]

    # if there are more report rows, append them to the first row
    if len(bp)>1:
        for i in range(1,len(bp),1):
            rep.append(np.array([(hd['scan'], hd['obsid'], typ, bp[i], spe[bp[i]], 0.0, vel[bp[i]])], 
              dtype=[('scan', 'i4'),('obsid', 'i4'),('type', 'S10'),('pixel', 'i4'),('nval', 'f8'), ('val', 'f8'), ('velo', 'f4')]))
            
    return rep




def top3perc(dat, lim=0.09, verbose=False):
    """
    sorting and returning the args for the highest 3 percent (or what set with lim)
    of values in dat. Does not work on 2-d arrays!
    dat: should be a 1d-array, e.g. a single spectrum
    """
    fdat = dat.flatten()
    nn = fdat.shape[0]
    
    # sorting in increasing order, highest values at end
    ss = np.argsort(fdat, axis=None)
    en = int(-lim*nn)
    if verbose: print('higest entries in data array (indices/values):       ', en, dat[ss[en:]])
    
    # return indices of highest values
    return ss[int(-lim*nn):]




                

