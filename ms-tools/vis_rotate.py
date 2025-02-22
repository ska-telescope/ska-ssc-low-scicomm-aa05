#!/usr/bin/env python

from casacore.tables import table
import numpy as np
from numpy.linalg import inv
import matplotlib.pyplot as plt
from argparse import ArgumentParser
from tqdm import tqdm
from astropy.time import Time
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
import astropy.units as u
from glob import glob
from ska_low_mccs_calibration import eep

"""
Define the station rotation angles for AA0.5
"""
angles = {}
angles['s8-1'] = 251.3 * np.pi/180.
angles['s8-6'] = 193.6 * np.pi/180.
angles['s9-2'] = 56.1 * np.pi/180.
angles['s10-3'] = 44.2 * np.pi/180.

# Define some colors and legend strings for plotting
pc = ['k', 'g', 'b', 'r']
corr = ['XX', 'XY', 'YX', 'YY']

# Set the location of SKA-Low
# This is the location of C1
obsloc = EarthLocation(lat=-26.826161995090644*u.deg,
                       lon=116.7656960036497*u.deg,
                       height=347.0376321554915*u.m)

"""
This is the rotation function used if mode = 'ground'
"""
def rotmat(ang, do_inv):
    mat = np.array([[np.cos(-ang), -np.sin(-ang)],
                    [np.sin(-ang), np.cos(-ang)]], dtype=np.complex64)
    if do_inv:
        return inv(mat), inv(mat.conj().T)
    else:
        return mat, mat.conj().T

"""
Rotation function used if mode = 'analytic'
Based on Randall's getJonesAnalyticSimple
Updated 21 Feb 2025 with extra cos(theta) factor
"""
def rot_analytic(time, loc, radec, ang, do_inv):
    altaz = radec.transform_to(AltAz(obstime=time, location=loc))
    tr = altaz.zen.rad
    azr = altaz.az.rad
    pr = (np.pi/2.) - azr
    ##print(tr,pr,ang)
    ct = np.cos(tr)
    jpp = -np.sin(pr-ang)*ct
    jpt = np.cos(pr-ang)*ct*ct
    jqp = np.cos(pr-ang)*ct
    jqt = np.sin(pr-ang)*ct*ct
    # include 0.5 factor to match scaling of EEPs approximately
    J=np.array([[jpp,jpt],[jqp,jqt]],dtype=np.complex64)
    if do_inv:
        return inv(J), inv(J.conj().T)
    else:
        return J, J.conj().T


"""
Rotation function used if mode = 'eep'
"""
def rot_eep(time, loc, ra, dec, ang, feep, do_inv, do_pacorr):
    # We will get a Jones matrix for each time step
    J = eep.station_beam_matrix(feep,
                                ang, # -ve seems necessary from testing
                                loc,
                                time,
                                ra,
                                dec,
                                pa_correction=do_pacorr)
    # Swap sign of X-axis (test)
    for i in range(J.shape[0]):
        J[i][0][1] = -J[i][0][1]
        J[i][1][0] = -J[i][1][0]
    # Prepare the arrays to return
    Jarr = np.zeros_like(J)
    JarrT = np.zeros_like(J)
    if do_inv: # inv is the other way around for this Jones matrix...
        for tstep in range(len(J)):
            Jarr[tstep] = J[tstep]
            JarrT[tstep] = J[tstep].conj().T
    else:
        for tstep in range(len(J)):
            Jarr[tstep] = inv(J[tstep])
            JarrT[tstep] = inv(J[tstep].conj().T)
    return Jarr, JarrT


def rotate_ms(msname,
              column='CORRECTED_DATA',
              write_to_column='__NONE__',
              mode='analytic',
              invert=False,
              plots=False,
              dryrun=False,
              eepdir='/shared/eep-data/Perturbed_Vogel_HARP/Average_EEPs/',
              eepbase='HARP_SKALA41_randvogel_avg_',
              eepsuff='.npz',
              eepcorr=False):
    # First open up the MS
    if dryrun:
        print('DRY RUN ... opening MS readonly')
        t = table(msname, readonly=True)
    else:
        t = table(msname, readonly=False)

    # Check that the requested column exists
    if column not in t.colnames():
        print(f'Error: MS column {column} does not exist in {msname}')
        print(t.colnames())
        t.close()
        exit(1)

    # Grab the times (relative to first timestamp)
    timecol = np.sort(np.unique(t.getcol('TIME')),axis=None)
    times = Time(timecol/(24.*3600.), format='mjd')
    # for plotting purposes:
    timecol -= timecol[0] # seconds since start of obs

    # Get the coordinates of the target field
    tfield = table(t.getkeyword('FIELD'), readonly=True)
    phasedir = np.squeeze(tfield.getcol('PHASE_DIR'))*180./np.pi
    target = SkyCoord(ra=phasedir[0], dec=phasedir[1], unit="deg")
    print('Target coords:', target.to_string(style='hmsdms'))
    # For use with EEPs we only need RA,Dec in Degrees
    ra = phasedir[0]
    dec = phasedir[1]

    # Open the antenna table and obtain antenna IDs
    tant = table(t.getkeyword('ANTENNA'), readonly=True)
    antnames = tant.getcol('NAME')
    print('Stations:', antnames)

    # Open the spectral window table and get the channel frequencies
    tspw = table(t.getkeyword('SPECTRAL_WINDOW'), readonly=True)
    chan_freq = np.squeeze(tspw.getcol('CHAN_FREQ'))
    # If we are using EEP mode, get the EEPs
    if mode == 'eep':
        print('Getting EEPs ...')
        chan_freq_round_mhz = np.round(chan_freq/1.e6)
        chan_freq_to_get = np.unique(chan_freq_round_mhz)
        feeps = {} # For storing the EEP ndarrays
        for fmhz in chan_freq_to_get:
            feeps[fmhz] = eep.load_eeps(fmhz,
                                        eepdir,
                                        filebase=eepbase,
                                        suffix=eepsuff)
        print(f'Collected {len(feeps)} EEPs for {len(chan_freq)} channels')

    # Iterate over the unique (cross correlation baselines)
    for tsub in t.iter(('ANTENNA1','ANTENNA2')):
        # Which baseline are we dealing with?
        ant1 = np.unique(tsub.getcol('ANTENNA1'))[0]
        ant2 = np.unique(tsub.getcol('ANTENNA2'))[0]
        if ant1 == ant2:
            print(f'Skipping autocorrelation {ant1} ({antnames[ant1]})')
            continue
        print(f'Processing baseline {ant1}&{ant2} ({antnames[ant1]}&{antnames[ant2]}), using {column}')
        # If we're just rotating the ground coordinates we can already
        # create the rotation matrix for each station
        # Take into account whether we are rotating or unrotating
        if mode=='ground':
            rot_ant1, rot_ant1_T = rotmat(angles[antnames[ant1]], invert)
            rot_ant2, rot_ant2_T = rotmat(angles[antnames[ant2]], invert)

        # Get the relevant data from the MS and make an array to store results
        unrot_data = tsub.getcol(column)
        rot_data = np.ones(unrot_data.shape, dtype=unrot_data.dtype)

        # Rotate the data for this baseline
        prev_fmhz = 0.
        for j in tqdm(range(unrot_data.shape[1])): # frequency
            if mode == 'eep':
                fmhz = round(chan_freq[j]/1.e6)
                if fmhz != prev_fmhz: # skip making new EEP if same freq as last
                    eepJ1, eepJ1T = rot_eep(times, obsloc,
                                            ra, dec,
                                            angles[antnames[ant1]]*180./np.pi,
                                            feeps[fmhz], invert,
                                            eepcorr)
                    eepJ2, eepJ2T = rot_eep(times, obsloc,
                                            ra, dec,
                                            angles[antnames[ant2]]*180./np.pi,
                                            feeps[fmhz], invert,
                                            eepcorr)
                prev_fmhz = fmhz*1.
            for i in range(unrot_data.shape[0]): # time
                if mode == 'analytic':
                    rot_ant1, rot_ant1_T = rot_analytic(times[i], obsloc, target,
                                                        angles[antnames[ant1]],
                                                        invert)
                    rot_ant2, rot_ant2_T = rot_analytic(times[i], obsloc, target,
                                                        angles[antnames[ant2]],
                                                        invert)
                elif mode == 'eep':
                    rot_ant1 = eepJ1[i]
                    rot_ant1_T = eepJ1T[i]
                    rot_ant2 = eepJ2[i]
                    rot_ant2_T = eepJ2T[i]
                # At this point the operation is the same no matter
                # where the Jones matrices came from
                rot_data[i,j,:] = np.matmul(rot_ant1,
                                            np.matmul(unrot_data[i,j,:].reshape((2,2)),
                                                      rot_ant2_T)).reshape((4))
                # Alternative formulation:
                #vrecmat = unrot_data[i,j,:].reshape((2,2))
                #rvrec = rot_ant1_T@vrecmat@rot_ant2
                #rot_data[i,j,:] = rvrec.reshape((4))

        # Place the rotated data back where it came from, or in DATA,
        # depending on what the user asked for
        if not dryrun:
            if write_to_column != '__NONE__':
                print(f'Writing rotated data to {write_to_column} column')
                tsub.putcol(write_to_column, rot_data)
            else:
                print('Writing rotated data to',column,'column')
                tsub.putcol(column, rot_data)

        # If requested, make a plot that shows the before and after data
        # TO DO : Take into account the FLAG column
        if plots:
            out_file = f'visdata_{antnames[ant1]}_{antnames[ant2]}.png'
            print('Plotting to', out_file)
            # Make the plot, 2x2 with amp/phase for the baseline, before and after
            fig, ax = plt.subplots(2,2,figsize=(12,10))
            for k in range(4):
                ax[0][0].plot(timecol, np.mean(np.abs(unrot_data[:,:,k]), axis=1), f'{pc[k]}.', label=corr[k])
                ax[1][0].plot(timecol, np.mean(np.angle(unrot_data[:,:,k])*180./np.pi, axis=1), f'{pc[k]}.', label=corr[k])
                ax[0][1].plot(timecol, np.mean(np.abs(rot_data[:,:,k]), axis=1), f'{pc[k]}.', label=corr[k])
                ax[1][1].plot(timecol, np.mean(np.angle(rot_data[:,:,k])*180./np.pi, axis=1), f'{pc[k]}.', label=corr[k])
                ax[0][0].set_xlabel('Time (sec)')
                ax[1][0].set_xlabel('Time (sec)')
                ax[0][1].set_xlabel('Time (sec)')
                ax[1][1].set_xlabel('Time (sec)')
                ax[0][0].set_ylabel('Amp unrotated')
                ax[1][0].set_ylabel('Phase unrotated (deg)')
                ax[0][1].set_ylabel('Amp rotated')
                ax[1][1].set_ylabel('Phase rotated (deg)')
                ax[1][0].set_ylim((-200,200))
                ax[1][1].set_ylim((-200,200))
                ax[0][0].legend(loc='best')
                ax[0][1].legend(loc='best')
                ax[1][0].legend(loc='best')
                ax[1][1].legend(loc='best')
                ax[0][1].yaxis.set_label_position("right")
                ax[0][1].yaxis.tick_right()
                ax[1][1].yaxis.set_label_position("right")
                ax[1][1].yaxis.tick_right()
                fig.suptitle(f'{msname} - {column}')
                plt.tight_layout()
                plt.savefig(out_file, dpi=200, bbox_inches='tight')

    t.flush()
    t.close()
    tant.close()


if __name__ == '__main__':
    ap = ArgumentParser()
    ap.add_argument('msname',
                    help='Name of MS. This can be a wildcard string. If more than one MS name matches, the process will run for each of them one at a time. [no default]')
    ap.add_argument('-c','--column',
                    help='MS column name to read. This is also the column that will be written to if not --dryrun, and unless --write_to_column is used. [default %(default)s]', default='CORRECTED_DATA')
    ap.add_argument('-w','--write_to_column',
                    help='Write result to another column? (USE WITH CAUTION) [default value of COLUMN argument]',
                    default='__NONE__')
    ap.add_argument('-m','--mode',
                    help='Conversion mode. Choose from "ground" rotation (not recommended), "analytic", or "eep". [default %(default)s]',
                    default='analytic', choices=['ground', 'analytic', 'eep'])
    ap.add_argument('-i','--invert',
                    help='Invert sense of rotation? Use this option if you are forward modeling (applying to MODEL_DATA) or undoing a previously-applied rotation. [default %(default)s]',
                    default=False, action='store_true')
    ap.add_argument('-p', '--plots',
                    help='Generate plots to show before-and-after data? These will be named like "visdata_STATION1_STATION2.png" [default %(default)s]',
                    default=False, action='store_true')
    ap.add_argument('-d', '--dryrun',
                    help='Only compute the rotations but do not write back into the MS? Useful with -p. [default %(default)s]',
                    default=False, action='store_true')
    ap.add_argument('--eepdir',
                    help='EEP base directory. Ignored unless MODE is EEP. [default %(default)s]',
                    default='/shared/eep-data/Perturbed_Vogel_HARP/Average_EEPs/')
    ap.add_argument('--eepbase',help='EEP filename base. Ignored unless MODE is EEP. [default %(default)s]',
                    default='HARP_SKALA41_randvogel_avg_')
    ap.add_argument('--eepsuff',
                    help='EEP filename suffix. Ignored unless MODE is EEP. [default %(default)s]',
                    default='.npz')
    ap.add_argument('--eepcorr',
                    help='Apply the EEP PA correction? Ignored unless MODE is EEP. [default %(default)s]',
                    default=False, action='store_true')
    args = ap.parse_args()

    ms_name_list = glob(args.msname)
    if len(ms_name_list)==0:
        print(f'Error: no MS by that name ({args.msname}) here')
    del args.msname
    vargs = vars(args)
    if len(ms_name_list)==1:
        rotate_ms(ms_name_list[0], **vargs)
    else:
        print(f'Running sequentially for {len(ms_name_list)} MSs')
        for this_ms in ms_name_list:
            rotate_ms(this_ms, **vargs)

