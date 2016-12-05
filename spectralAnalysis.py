import numpy as np
from astropy.io import ascii, fits
import pydis
import glob


#disDir = "/Users/cbrasseur/Documents/ObsAstro/APO/Q4JH03/UT161120/DIS/"


# apparently we didn't take any bias frames??? wtf????
#bias = np.zeros((1028, 2048)) 

# opens 2 matplotlib graphs
#flat, fmask_out = pydis.flatcombine(disDir + 'rflat.lis', bias,
#                                   trim=True, mode='spline', response=True,
#                                   display=True, output=disDir+'FLAT.fits')


#HeRed = disDir + 'He.0007r.fits'
#ArRed = disDir + 'Ar.0008r.fits'
#NeRed = disDir + 'Ne.0009r.fits'

# using the neon frame for calibration
#wfit = pydis.HeNeAr_fit(NeRed, trim=True, fmask=fmask_out,
#                          interac=True, mode='poly', fit_order=5)



def makeSpectrum(dataFile,stdDataFile,
                 fitFile,flatFile,fmaskFile,
                 traceFile = None,stdTraceFile = None,
                 interactive=False,spectrumFile = None):

    # NOTE: we did not take bias frames, so bias is assumed 0
    bias = np.zeros((1028, 2048))
    flat = fits.getdata(flatFile)
    fmask_out = np.loadtxt(fmaskFile).astype(int)

    wfit = np.loadtxt(fitFile)

    
    ## reduce the science image ##

    # open the image, get the data and a couple important numbers
    img = pydis.OpenImg(dataFile, trim=True)
    raw = img.data
    exptime = img.exptime
    airmass = img.airmass

    # remove bias and flat, divide by exptime
    data = ((raw - bias) / flat) / exptime


    ## Read in and reduce the flux standard star ##
    img = pydis.OpenImg(stdDataFile, trim=True)
    stdraw = img.data
    stdexptime = img.exptime
    stdairmass = img.airmass
    stddata= ((stdraw - bias) / flat) / stdexptime


    if traceFile:
        trace = np.loadtxt(traceFile)
    else:
        trace = pydis.ap_trace(data, fmask=fmask_out, nsteps=7, interac=False, display=True) 

    if stdTraceFile:
        stdtrace = np.loadtxt(stdTraceFile)
    else:
        stdtrace = pydis.ap_trace(stddata, fmask=fmask_out, nsteps=7, interac=False, display=True)
    
    ext_spec, sky, fluxerr = pydis.ap_extract(data, trace, apwidth=5,skysep=1,
                                              skywidth=7, skydeg=0)
    ext_std, stdsky, stderr = pydis.ap_extract(stddata, stdtrace, apwidth=5,
                                               skysep=1, skywidth=7, skydeg=0)

    # subtract the sky from the 1-d spectrum
    flux_red = (ext_spec - sky) # the reduced object
    flux_std = (ext_std - stdsky) # the reduced flux standard
    
    # map the wavelength using the HeNeAr fit
    wfinal = pydis.mapwavelength(trace, wfit, mode='poly')
    wfinalstd = pydis.mapwavelength(stdtrace, wfit, mode='poly')

    # correct the object and flux std for airmass extinction
    flux_red_x = pydis.AirmassCor(wfinal, flux_red, airmass,
                                  airmass_file='apoextinct.dat')
    flux_std_x = pydis.AirmassCor(wfinalstd, flux_std, stdairmass,
                                  airmass_file='apoextinct.dat')

    sensfunc = pydis.DefFluxCal(wfinalstd, flux_std_x, mode='spline',
                                stdstar='spec50cal/feige34.dat',display=False)

    # final step in reduction, apply sensfunc
    ffinal,efinal = pydis.ApplyFluxCal(wfinal, flux_red_x, fluxerr, 
                                       wfinalstd, sensfunc)

    # Saving the calebrated spectrum

    if not spectrumFile:
        spectrumFile = dataFile[:-5]+"_CALSPEC.fits"
    wavelengthCol = fits.Column(name='Wavelength', format='D', array=wfinal)
    fluxCol = fits.Column(name='Flux', format='D', array=ffinal)
    errorCol = fits.Column(name='Flux Error', format='D', array=efinal)

    cols = fits.ColDefs([wavelengthCol,fluxCol,errorCol])
    
    tbhdu = fits.BinTableHDU.from_columns(cols)

    tbhdu.writeto(spectrumFile)
    


# NOTE: don't call this from jupyter notebook, need interactive functionality
def makingTraceFile(dataFile,flatFile,fmaskFile,nsteps=7,traceFile=None):

    # NOTE: we did not take bias frames, so bias is assumed 0
    bias = np.zeros((1028, 2048))
    flat = fits.getdata(flatFile)
    fmask_out = np.loadtxt(fmaskFile).astype(int)

    # open the image, get the data and a couple important numbers
    img = pydis.OpenImg(dataFile, trim=True)
    raw = img.data
    exptime = img.exptime
    airmass = img.airmass

    # remove bias and flat, divide by exptime
    data = ((raw - bias) / flat) / exptime

    trace = pydis.ap_trace(data, fmask=fmask_out, nsteps=nsteps, interac=True, display=True)

    if not traceFile:
        traceFile = dataFile[:-5]+"_TRACE.txt"

    np.savetxt(traceFile,trace)

    return traceFile

    
# NOTE: don't call this from jupyter notebook, need interactive functionality
# band should be r or b (red or blue)
def makingFit(disDir,band,calFle):

    # making a list of the flat fields
    with open (band+'flat.lis','w') as FLATFLE:
        FLATFLE.write('\n'.join(glob.glob(disDir + "*flat.0*"+band+".fits")))


    bias = np.zeros((1028, 2048)) # apparently we didn't take any bias frames??? wtf????

    flat, fmask_out = pydis.flatcombine(band+'flat.lis', bias, trim=True, mode='spline', 
                                        response=True, display=True, output=band+'_FLAT.fits')

    calFle = disDir + calFle

    wfit = pydis.HeNeAr_fit(calFle, trim=True, fmask=fmask_out,
                            interac=True, mode='poly', fit_order=5)

    np.savetxt(band+"_WFIT.txt",wfit)
    np.savetxt(band+"_FMASK.txt",fmask_out)

    return band+'_FLAT.fits',band+"_FMASK.txt",band+"_WFIT.txt"
