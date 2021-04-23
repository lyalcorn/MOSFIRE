#!kpython3

## Import General Tools
import inspect
from datetime import datetime
from pathlib import Path
import argparse
import logging

from mosfire.mask import Mask


description ='''This script takes an xml file containing a mask design
(generated by MAGMA) and evaluates when the rotator will be in one of MOSFIRE's
"bad angles" (e.g. -180, 0, +180) on the rotator.
'''


##-------------------------------------------------------------------------
## Parse Command Line Arguments
##-------------------------------------------------------------------------
## create a parser object for understanding command-line arguments
p = argparse.ArgumentParser(description=description)
## add flags
p.add_argument("-v", "--verbose", dest="verbose",
    default=False, action="store_true",
    help="Be verbose! (default = False)")
p.add_argument("-p", "--plot", dest="plot",
    default=False, action="store_true",
    help="Generate plot")
## add options
p.add_argument('maskfile', type=str,
               help="The XML file containing your mask")
p.add_argument("--night", dest="night", type=str,
    help="The UT night to check (in YYYY-MM-DD format).  Defaults to today.")
args = p.parse_args()


##-------------------------------------------------------------------------
## Create logger object
##-------------------------------------------------------------------------
log = logging.getLogger('check_mask')
log.setLevel(logging.DEBUG)
## Set up console output
LogConsoleHandler = logging.StreamHandler()
if args.verbose:
    LogConsoleHandler.setLevel(logging.DEBUG)
else:
    LogConsoleHandler.setLevel(logging.INFO)
LogFormat = logging.Formatter('%(asctime)s %(levelname)8s: %(message)s',
                              datefmt='%Y-%m-%d %H:%M:%S')
LogConsoleHandler.setFormatter(LogFormat)
log.addHandler(LogConsoleHandler)
## Set up file output
# LogFileName = None
# LogFileHandler = logging.FileHandler(LogFileName)
# LogFileHandler.setLevel(logging.DEBUG)
# LogFileHandler.setFormatter(LogFormat)
# log.addHandler(LogFileHandler)


##-------------------------------------------------------------------------
## Check Mask Angles
##-------------------------------------------------------------------------
def check_mask_angles(maskfile, night=None, plot=False,
                      skipprecond=False, skippostcond=True):
    this_script_name = inspect.currentframe().f_code.co_name
    log.debug(f"Executing: {this_script_name}")

    ##-------------------------------------------------------------------------
    ## Pre-Condition Checks
    if skipprecond is True:
        log.debug('Skipping pre condition checks')
    else:
        if night in [None, '']:
            now = datetime.utcnow()
            night = now.strftime('%Y-%m-%d')
    
    ##-------------------------------------------------------------------------
    ## Script Contents

    mask = Mask(maskfile)
    mask.find_bad_angles(night=night, plot=plot)
    
    ##-------------------------------------------------------------------------
    ## Post-Condition Checks
    if skippostcond is True:
        log.debug('Skipping post condition checks')
    else:
        pass


if __name__ == '__main__':
    check_mask_angles(args.maskfile, night=args.night, plot=args.plot)
