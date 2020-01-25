#!/usr/env/python

## Import General Tools
import sys
from pathlib import Path
import logging
import yaml
import random
import re

from datetime import datetime as dt
from datetime import timedelta as tdelta
from time import sleep
import numpy as np
import subprocess

from astropy.table import Table, Column
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.modeling import models, fitting
from scipy import ndimage

from instruments import connect_to_ktl, create_log

try:
    from ktl import Exceptions as ktlExceptions
except:
    pass


services = {}


##-------------------------------------------------------------------------
## Define Common Functions
##-------------------------------------------------------------------------
def connect(service):
    if type(service) is str and service.lower() == 'all':
        services.update(connect_to_ktl(name, serviceNames))
    else:
        services.update(connect_to_ktl(name, [service]))


def get(keyword, service='mosfire', mode=str):
    """Generic function to get a keyword value.  Converts it to the specified
    mode and does some simple parsing of true and false strings.
    """
    if service not in services.keys():
        log.error(f'Not connected to service "{service}"')
        return None

    log.debug(f'Querying {service} for {keyword}')
    if services == {}:
        return None
    assert mode in [str, float, int, bool]
    kwresult = services[service][keyword].read()
    log.debug(f'  Got result: "{kwresult}"')

    # Handle string versions of true and false
    if mode is bool:
        if kwresult.strip().lower() == 'false':
            result = False
        elif kwresult.strip().lower() == 'true':
            result = True
        else:
            try:
                result = bool(int(kwresult))
            except:
                result = None
        if result is not None:
            log.debug(f'  Parsed to boolean: {result}')
        else:
            log.error(f'  Failed to parse "{kwresult}"')
        return result
    # Convert result to requested type
    try:
        result = mode(kwresult)
        log.debug(f'  Parsed to {mode}: {result}')
        return result
    except ValueError:
        log.warning(f'Failed to parse {kwresult} as {mode}, returning string')
        return kwresult


def set(keyword, value, service='mosfire', wait=True):
    """Generic function to set a keyword value.
    """
    if service not in services.keys():
        log.error(f'Not connected to service "{service}"')
        return None

    log.debug(f'Setting {service}.{keyword} to "{value}" (wait={wait})')
    if services == {}:
        return None
    services[service][keyword].write(value, wait=wait)
    log.debug(f'  Done.')


##-------------------------------------------------------------------------
## Rotator Control
##-------------------------------------------------------------------------
def _set_rotpposn(rotpposn):
    '''Set the rotator position in stationary mode.
    
    This only tries to set the position once, use `set_rotpposn` in practice
    as it makes multiple attempts which seems to be more reliable.
    '''
    log.info(f'Setting ROTPPOSN to {rotpposn:.1f}')
    set('rotdest', rotpposn, service='dcs')
    sleep(1)
    set('rotmode', 'stationary', service='dcs')
    sleep(1)
    done = get('rotstat', service='dcs') == 'in position'
    while done is False:
        sleep(1)
        done = get('rotstat', service='dcs') == 'in position'


def set_rotpposn(rotpposn):
    '''Set the rotator position in stationary mode.
    '''
    try:
        _set_rotpposn(rotpposn)
    except ktlExceptions.ktlError as e:
        log.warning(f"Failed to set rotator")
        log.warning(e)
        sleep(2)
        log.info('Trying again ...')
        set_rotpposn(rotpposn)


##-------------------------------------------------------------------------
## Software
## * check for camera connection (ICE?)
##-------------------------------------------------------------------------
def watchrot_ok():
    return get('WATCHROTOK', mode=bool)

def watchslew_ok():
    return get('WATCHSLEWOK', mode=bool)

def watchfcs_ok():
    return get('WATCHFCSOK', mode=bool)

def watchdarcorr_ok():
    return get('DARENABL', mode=bool)

def autodisplay_ok():
    return get('AUTODISPLAYOK', mode=bool)

def watchprocesses_ok():
    ok = watchrot_ok()
    ok &= watchslew_ok()
    ok &= watchfcs_ok()
    ok &= watchdarcorr_ok()
    ok &= autodisplay_ok()
    return ok
