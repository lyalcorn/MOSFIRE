#!kpython3
description = '''
acq_long2pos -- Acquire mask nod image sets at two different longslit mask
                positions.

This script is used in conjunction with the long2pos CSU mask to obtain spectra
of a target in two longslits which are offset from the field center in the
wavelength direction in order to sample the full range of wavelengths in the
passband.

                    |x|
                    | | Upper left position, covers long wavelengths
                    |x|
                                 |  | Alignment position
                                              |x|
                                              | | Lower right, short wavelengths
                                              |x|


The normal behavior, which applies when the CSU is configured with the standard
"long2pos" or "long2pos (align)" masks is to acquire a nod pair in each slit. 
If the CSU is configured with the "long2pos_specphot" or
"long2pos_specphot (align)", we acquire spectra in the center of each slit as
well as at the ends.  This is typically used when the observer wants both a
spectrum from both the "wide" portion of the slit and the narrow portion.  The
spectrum from the "wide" part can be used as a spectrophotometric observation
(no slit loss from the extra-wide slit).

Modification History:
    2021-04-23 [JoshW] Initial adaptation from the c-shell script.  Incorporates
                       changes to improve positioning accuracy.
'''

## Import General Tools
from pathlib import Path
import argparse
import logging
from datetime import datetime
from time import sleep
import numpy as np
import subprocess

import ktl

import mosfire


##-------------------------------------------------------------------------
## Parse Command Line Arguments
##-------------------------------------------------------------------------
## create a parser object for understanding command-line arguments
p = argparse.ArgumentParser(description=description)
## add flags
p.add_argument("-v", "--verbose", dest="verbose",
    default=False, action="store_true",
    help="Be verbose! (default = False)")
p.add_argument("-i", "--imaging", dest="imaging", type=bool,
    default=False, action="store_true",
    help="Take data in imaging mode rather than spectroscopy? (for testing)")
## add options
p.add_argument("-f", "--filter", dest="filter", type=str, default='',
    choices=['Y', 'J', 'H', 'K', 'Ks', 'J2', 'J3', 'nb1061', ''],
    help="The filter to take data in.")
p.add_argument("-e", "--exptime", dest="exptime", type=float, default=0,
    help="The exposure time (defaults to 0 which means use existing setting).")
args = p.parse_args()


##-------------------------------------------------------------------------
## Create logger object
##-------------------------------------------------------------------------
log = logging.getLogger('long2pos')
log.setLevel(logging.DEBUG)
LogFormat = logging.Formatter('%(asctime)s %(levelname)8s: %(message)s')
## Set up console output
LogConsoleHandler = logging.StreamHandler()
LogConsoleHandler.setLevel(logging.DEBUG)
LogConsoleHandler.setFormatter(LogFormat)
log.addHandler(LogConsoleHandler)
## Set up file output
now_str= datetime.utcnow().strftime('%Y%m%d_%H%M%SUT')
LogFileName = Path('~/acq_long2pos_{now_str}.log').expanduser()
LogFileHandler = logging.FileHandler(LogFileName)
LogFileHandler.setLevel(logging.DEBUG)
LogFileHandler.setFormatter(LogFormat)
log.addHandler(LogFileHandler)


##-------------------------------------------------------------------------
## Cache keyword services
##-------------------------------------------------------------------------
mosfireGS = ktl.cache(service='mosfire')
mds = ktl.cache(service='mds')
magiq = ktl.cache(service='magiq')
dcs = ktl.cache(service='dcs')


##-------------------------------------------------------------------------
## acq_long2pos
##-------------------------------------------------------------------------
def acq_long2pos(wide=False, waittime=10):
    log.info('Starting long2pos acquisition')
    log.debug('Setting SCRIPTRUN=1')
    mosfire.start_scriptrun()
    log.info('Marking base')
    mosfire.markbase()
    log.debug('Setting PATTERN=long2pos')
    mosfireGS['PATTERN'].write('long2pos')

    pixel_scale = mosfireGS['pscale'].read(binary=True)
    log.info('Offsetting to upper left position')
    mosfire.mxy(250.25*pixel_scale, -89.0*pixel_scale)
    log.info(f'Waiting {waittime} seconds for guider')
    sleep(waittime)
    log.debug('Setting XOFFSET=-45 YOFFSET=0')
    mosfireGS['XOFFSET'].write(-45.0)
    mosfireGS['YOFFSET'].write(0.0)

    if wide is True:
        log.debug('Setting FRAMEID=A')
        mosfireGS['FRAMEID'].write('A')
        log.info('Taking wide slit data')
        mosfire.take_exposure()

    log.info('acquiring slit nod data in upper left position')
    mosfire.sltmov(-7)
    log.debug('Setting YOFFSET=-7 FRAMEID=B')
    mosfireGS['YOFFSET'].write(-7.0)
    mosfireGS['FRAMEID'].write('B')
    mosfire.take_exposure()
    mosfire.sltmov(14)
    log.debug('Setting YOFFSET=+7 FRAMEID=A')
    mosfireGS['YOFFSET'].write(7.0)
    mosfireGS['FRAMEID'].write('A')
    mosfire.take_exposure()

    log.info(f'Returning to base before offsetting to lower right position')
    mosfire.gotobase()
    log.info(f'Waiting {waittime} seconds for guider')
    sleep(waittime)

    log.info('Offsetting to lower right position')
    mosfire.mxy(-250.25*pixel_scale, 89.0*pixel_scale)
    log.info(f'Waiting {waittime} seconds for guider')
    sleep(waittime)
    log.debug('Setting XOFFSET=+45 YOFFSET=0')
    mosfireGS['XOFFSET'].write(+45.0)
    mosfireGS['YOFFSET'].write(0.0)

    if wide is True:
        log.debug('Setting FRAMEID=A')
        mosfireGS['FRAMEID'].write('A')
        log.info('Taking wide slit data')
        mosfire.take_exposure()

    log.info('acquiring slit nod data in lower right position')
    mosfire.sltmov(-7)
    log.debug('Setting YOFFSET=-7 FRAMEID=B')
    mosfireGS['YOFFSET'].write(-7.0)
    mosfireGS['FRAMEID'].write('B')
    mosfire.take_exposure()
    mosfire.sltmov(14)
    log.debug('Setting YOFFSET=+7 FRAMEID=A')
    mosfireGS['YOFFSET'].write(7.0)
    mosfireGS['FRAMEID'].write('A')
    mosfire.take_exposure()

    log.info(f'Returning to base')
    mosfire.gotobase()
    log.debug('Setting XOFFSET=0 YOFFSET=0')
    mosfireGS['XOFFSET'].write(0.0)
    mosfireGS['YOFFSET'].write(0.0)
    log.debug('Setting SCRIPTRUN=0')
    mosfire.stop_scriptrun()


def acq_long2pos_specphot():
    acq_long2pos(wide=True)


def run_acq_long_2pos(min_guide_cycles=5):
    '''Adopted from acq_long2pos shell script and modified based on advice from
    Shui and based on the experience of the KCWI nod and shuffle mode.
    
    - Check exposure time of instrument and guider and decide whether guiding
      is advised.
    - Go to base between the two positions (posA and posC)
    - After offsetting, wait for guider centroid to converge before starting
      exposure.
    '''
    # Check that MOSFIRE is the selected instrument
    instrument_is_MOSFIRE() # raises exception if not

    log.info(f"Science Exposure Parameters:")
    # Read exposure time
    ITIME = mosfire.exptime()
    log.info(f"  Exposure Time = {ITIME:.0f} s")
    # Read coadds
    COADDS = mosfire.coadds()
    log.info(f"  coadds = {COADDS}")
    # Read sampling mode & number of reads
    sampmode = mosfire.sampmode()
    log.info(f"  Sampling Mode = {sampmode}")
    # Check the mask
    maskname = mosfireGS['maskname'].read()
    if maskname[:9] != 'long2pos':
        log.error(f"Mask {maskname} is not a long2pos mask")
        raise Exception(f"Mask {maskname} is not a long2pos mask")

    # Check guider and science exposure times
    log.info(f"Guider Exposure Parameters:")
    guider_exptime = magiq['TTIME'].read(binary=True)
    log.info(f"  TTIME    = {guider_exptime:.0f} s")
    guider_xypos = magiq['XYPOS'].read()
    log.info(f"  XYPOS    = {guider_xypos}")
    guider_centroid = magiq['CENTROID'].read()
    log.info(f"  CENTROID = {guider_centroid}")
    guider_fwhm = magiq['FWHM'].read(binary=True)
    log.info(f"  FWHM     = {guider_fwhm:.0f}")
    guider_snr = magiq['SNR'].read(binary=True)
    log.info(f"  SNR      = {guider_snr:.0f}")
#     guider_counts = magiq['STARCNT'].read(binary=True)
#     log.info(f"  STARCT   = {guider_counts:.0f} s")

    # Should we guide?
    if ITIME > min_guide_cycles*guider_exptime:
        we_should_guide = True
        


    # Set the specified observing parameters
    if args.exptime > 0:
        log.info(f'Setting exposure time to {exptime}')
        mosfire.set_exptime(exptime)
        if args.exptime < 30:
            log.info(f'Setting sampling mode to CDS')
            mosfire.set_sampmode('CDS')

    if args.filter != '':
        obsmode = f'{args.filter}-spectroscopy' if args.imaging is False\
                  else f'{args.filter}-imaging'
        log.info(f'Setting OBSMODE = {obsmode}')
        mosfire.set_obsmode(obsmode)

    specphot = (maskname=='long2pos_specphot')
    acq_long2pos(wide=specphot, wait=3*guider_exptime)


if __name__ == '__main__':
    run_acq_long_2pos()
