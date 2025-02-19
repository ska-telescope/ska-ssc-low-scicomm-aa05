#!/usr/bin/env python

"""
Basic visibility calibration for AA0.5

The basic approach is:
- Basic flagging
- Delay calibration
- Averaging in time and frequency
- Visibility (de-)rotation
- Bandpass estimate
- Gain phase solution
- Bandpass refinement
- Gain solution

Currently this uses CASA as the main toolkit but we
may change to DP3 later depending on outcomes.
"""

from casatasks import gaincal, bandpass, applycal, split
from argparse import ArgumentParser
from subprocess import call
from astropy.time import Time
from time import strftime
from astropy.coordinates import AltAz, HADec, EarthLocation, SkyCoord
import astropy.units as u
import numpy as np
from casacore.tables import table

def main(args):

    # Step 0: initial items
    # Set the location of SKA-Low
    # This is the location of C1
    obsloc = EarthLocation(lat=-26.826161995090644*u.deg,
                           lon=116.7656960036497*u.deg,
                           height=347.0376321554915*u.m)
    t = table(args.msname, ack=False, readonly=True)
    # Get the coordinates of the target field
    tfield = table(t.getkeyword('FIELD'), ack=False, readonly=True)
    phasedir = np.squeeze(tfield.getcol('PHASE_DIR'))*180./np.pi
    target = SkyCoord(ra=phasedir[0], dec=phasedir[1], unit="deg")
    print('Target coords:', target.to_string(style='hmsdms'))

    # Step 1: flagging
    # Skip for now until we have a good aoflagger strategy
    #cmd = f'aoflagger {args.msname}'
    #call(cmd.split())

    # Step 2: delay calibration
    gaincal(vis=args.msname,
            caltable='delay.tb',
            uvrange='>500m',
            antenna='*&',
            solint='inf',
            refant='s10-3',
            refantmode='strict',
            minblperant=2,
            gaintype='K',
           )
    applycal(vis=args.msname,
             gaintable=['delay.tb'],
            )

    # Step 3: averaging
    split(vis=args.msname,
          outputvis='avg.ms',
          width=6,
          timebin='5s',
         )

    # Step 4: derotation
    # Mitch provided branch...
    cmd = f'../vis_rotate.py -c DATA -m eep --eepcorr avg.ms'
    call(cmd.split())

    # Step 5: Bandpass estimate
    # For this we want to select 10 minutes at the highest elevations
    print('Computing HA, elevation for target track')
    t = table('avg.ms', ack=False, readonly=True) # Smaller one is faster
    times = np.unique(t.getcol('TIME'))
    timesobj = Time(times/(24.*3600.), format='mjd')
    target_azel = target.transform_to(AltAz(obstime=timesobj, location=obsloc))
    target_hadec = target.transform_to(HADec(obstime=timesobj, location=obsloc))
    tippytop = np.max(target_azel.alt)
    ttloc = np.where(target_azel.alt==tippytop)[0]
    ttha = target_hadec.ha[ttloc].value
    print(f'Target reaches maximum elevation ({tippytop}) at HA={ttha}')
    bpstarttime = timesobj[ttloc]-5.*u.min
    bpendtime = timesobj[ttloc]+5.*u.min
    bptimerange = f"{bpstarttime.strftime('%H:%M:%S')[0]}~{bpendtime.strftime('%H:%M:%S')[0]}"
    print(f'Time range for bandpass solution: {bptimerange}')
    bandpass(vis='avg.ms',
             caltable='bandpass1.tb',
             uvrange='>500m',
             antenna='*&',
             timerange=bptimerange,
             solint='inf',
             refant='s10-3',
             minblperant=2,
             bandtype='B',
             smodel=[56.,0.,0.,0.],
            )
    
    # Step 6: Gain phase estimate
    # This is to make sure the next bandpass estimate is high quality
    gaincal(vis='avg.ms',
            caltable='gain_bp.tb',
            uvrange='>500m',
            antenna='*&',
            timerange=bptimerange,
            solint='2min',
            refant='s10-3',
            refantmode='strict',
            minblperant=2,
            gaintype='G',
            calmode='p',
            smodel=[56.,0.,0.,0.],
            gaintable=['bandpass1.tb'],
           )
    
    # Step 7: Refine bandpass
    bandpass(vis='avg.ms',
             caltable='bandpass2.tb',
             uvrange='>500m',
             antenna='*&',
             timerange=bptimerange,
             solint='inf',
             refant='s10-3',
             minblperant=2,
             bandtype='B',
             smodel=[56.,0.,0.,0.],
             gaintable=['gain_bp.tb'],
            )

    # Step 8: Time dependent gains
    gaincal(vis='avg.ms',
            caltable='gain.tb',
            uvrange='>500m',
            antenna='*&',
            solint='2min',
            refant='s10-3',
            refantmode='strict',
            minblperant=2,
            gaintype='G',
            calmode='ap',
            smodel=[56.,0.,0.,0.],
            gaintable=['bandpass2.tb'],
           )
    
    # Finally: apply the bandpass and gains
    applycal(vis='avg.ms',
             gaintable=['bandpass2.tb', 'gain.tb'],
            )

    print('Done.')

if __name__ == '__main__':
    ap = ArgumentParser()
    ap.add_argument('msname', help='Input MS')
    args = ap.parse_args()
    main(args)
