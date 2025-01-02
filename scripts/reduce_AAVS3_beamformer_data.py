#!/usr/bin/env python3


import argparse
import glob
import subprocess
import time

import numpy as np
from astropy.time import Time


def find_nearest(array, value, return_index=True):
    """Find value of an array element closest to provided value.

    Parameters
    ----------
    array : input array
    value : value being searched for
    return_index : also return index of found value

    Returns
    -------
    Value in the array that is closest to the searched value
    or if return_index is set to True, return tuple: index and
    searched value.
    """
    array = np.asarray(array)
    index = (np.abs(array - value)).argmin()
    if return_index:
        return index, array[index]
    else:
        return array[index]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        usage="reduce_AAVS3_beamformer_data.py --dir <input_dir> [options]",
        description="Fold AAVS3 raw beamformer data. To allow greater\n"
        "flexibility each DSPSR option has to be specified in full, i.e. flag\n"
        'and a value. E.g "--bin 1024" for 1024 profile phase bins.',
        formatter_class=lambda prog: argparse.HelpFormatter(
            prog, max_help_position=100, width=250
        ),
        epilog="Copyright (C) 2024 by SKA Observatory, Maciej Serylak (maciej.serylak@skao.int)",
    )
    parser.add_argument(
        "--dir",
        dest="directory",
        metavar="<input_dir>",
        required=True,
        help="specify directory with PSRDADA files",
    )
    parser.add_argument(
        "--outdir",
        dest="output_directory",
        metavar="<output_dir> ",
        default="./",
        help="specify output directory (default: ./)",
    )
    parser.add_argument(
        "--eph",
        dest="ephem_file",
        metavar="<ephem_file>",
        required=True,
        help="specify ephemeris option",
    )
    parser.add_argument(
        "--source",
        dest="source",
        metavar="<source>",
        help="specify source name",
    )
    parser.add_argument(
        "--pfb",
        dest="second_PFB",
        metavar="<second_PFB>",
        default="-F 256:D",
        help='specify secondary PFB option (default: "-F 256:D")',
    )
    parser.add_argument(
        "--bin",
        dest="phase_bins",
        metavar="<phase_bins>",
        default="-b 1024",
        help='specify number or profile phase bins option (default: "-b 1024")',
    )
    parser.add_argument(
        "--t_int",
        dest="integration",
        metavar="<integration>",
        default="-L 8",
        help='specify integration length option (default: "-L 8")',
    )
    parser.add_argument(
        "--seek",
        dest="seek",
        metavar="<seek>",
        default="-S 0",
        help='specify processing start at t=seek seconds option (default: "-S 0")',
    )
    parser.add_argument(
        "--total",
        dest="total",
        metavar="<total>",
        default=" ",
        help="specify processing length of t=total seconds option",
    )
    parser.add_argument(
        "--telescope",
        dest="telescope",
        metavar="<telescope>",
        default="-k MWA",
        help='specify telescope name (default: "-k MWA")',
    )
    parser.add_argument(
        "--set_key",
        dest="set_key",
        metavar="<set_key>",
        default=" ",
        help='set observation attributes option (example: "-set coord=09:53:09.3097+07:55:35.75")',
    )
    parser.add_argument(
        "--psrsh",
        dest="psrshell",
        metavar="<psrshell>",
        default='-j "delete chan 236-255,0-19"',
        help="specify psrsh command option to be run by DSPSR before output (default: '-j \"delete chan 236-255,0-19\"')",
    )
    parser.add_argument(
        "--threads",
        dest="threads",
        metavar="<threads>",
        default=" ",
        help="specify number of CPU processor threads option",
    )
    parser.add_argument(
        "--cuda",
        dest="cuda",
        # metavar="<cuda>",
        action="store_true",
        help="use the GPU when running DSPSR, disables --threads option (default: False)",
    )
    parser.add_argument(
        "--correct",
        dest="correct",
        action="store_true",
        help="only correct PSRDADA header",
    )
    parser.add_argument(
        "-v", "--verbose", dest="verbose", action="store_true", help="verbose mode"
    )
    args = parser.parse_args()

    # Start timing the script.
    script_start_time = time.perf_counter()

    # Prepare AAVS frequency table.
    ADC_clock = 800  # Station ADC clock frequency in MHz.
    raw_nchannels = 512  # Number of raw channels. Station uses 384 channels (64 - 447).
    frequency_table = [i * ADC_clock / (raw_nchannels * 2) for i in range(raw_nchannels)]

    # Find all PSRDADA files within directory.
    dada_file_list = glob.glob(args.directory + "*.dada", recursive=False)

    # Loop through files and correct the PSRDADA header.
    for dada_file in dada_file_list:
        print("Processing {0:}".format(dada_file))
        completed_process_1 = subprocess.run(
            "dada_edit " + dada_file, shell=True, check=True, capture_output=True
        )
        dada_header = completed_process_1.stdout.decode().split(
            "\n"
        )  # Read PSRDADA header.
        dada_header.remove("")  # Remove empty lines in the header.
        # Remove multiple white spaces from the header.
        for i, dada_header_field in enumerate(dada_header):
            dada_header[i] = " ".join(dada_header_field.split())
            # Turn the header into dictionary.
        dada_header_dictionary = {}
        for dada_header_field in dada_header:
            key, value = dada_header_field.split(" ")
            dada_header_dictionary[key] = value
        # Correct header values.
        # Correct total bandwidth in MHz. Here it is a single channel.
        completed_process_2 = subprocess.run(
            "dada_edit -c BW=0.9259259375 " + dada_file,
            shell=True,
            check=True,
            capture_output=True,
        )
        # Correct individual channel bandwidth in Hz.
        completed_process_3 = subprocess.run(
            "dada_edit -c BANDWIDTH=925925.937500 " + dada_file,
            shell=True,
            check=True,
            capture_output=True,
        )
        # Correct/add source name if specified.
        if args.source:
            completed_process_4 = subprocess.run(
                "dada_edit -c SOURCE=" + args.source + " " + dada_file,
                shell=True,
                check=True,
                capture_output=True,
            )
        # Correct channel frequency value.
        channel_index, channel_frequency = find_nearest(
            frequency_table, float(dada_header_dictionary["FREQ"])
        )
        completed_process_5 = subprocess.run(
            "dada_edit -c FREQ=" + str(channel_frequency) + " " + dada_file,
            shell=True,
            check=True,
            capture_output=True,
        )
        # Correct resolution of data (2048 * 2 * 2).
        try:
            dada_resolution = int(dada_header_dictionary["RESOLUTION"])
            if dada_resolution != 8192:
                completed_process_6 = subprocess.run(
                    "dada_edit -c RESOLUTION=8192 " + dada_file,
                    shell=True,
                    check=True,
                    capture_output=True,
                )
        except KeyError:
            completed_process_6 = subprocess.run(
                "dada_edit -c RESOLUTION=8192 " + dada_file,
                shell=True,
                check=True,
                capture_output=True,
            )
        # Correct timestamp offset.
        sample_time = float(dada_header_dictionary["TSAMP"]) * 1e6
        unix_timestamp_offset = dada_header_dictionary["UNIXTIME_MSEC"]
        picoseconds_offset = int(float(unix_timestamp_offset) * 1e9)
        try:
            if int(dada_header_dictionary["PICOSECONDS"]) == 0:
                completed_process_7 = subprocess.run(
                    "dada_edit -c PICOSECONDS="
                    + str(picoseconds_offset)
                    + " "
                    + dada_file,
                    shell=True,
                    check=True,
                    capture_output=True,
                )
        except KeyError:
            completed_process_7 = subprocess.run(
                "dada_edit -c PICOSECONDS=" + str(picoseconds_offset) + " " + dada_file,
                shell=True,
                check=True,
                capture_output=True,
            )

    # Loop through files again and fold the data with DSPSR.
    if args.correct:
        pass
    else:
        # Loop through files and run DSPSR.
        for dada_file in dada_file_list:
            print("Folding {0:}".format(dada_file))
            completed_process_8 = subprocess.run(
                "dada_edit " + dada_file, shell=True, check=True, capture_output=True
            )
            dada_header = completed_process_8.stdout.decode().split(
                "\n"
            )  # Read PSRDADA header.
            dada_header.remove("")  # Remove empty lines in the header.
            # Remove multiple white spaces from the header.
            for i, dada_header_field in enumerate(dada_header):
                dada_header[i] = " ".join(dada_header_field.split())
            # Turn the header into dictionary.
            dada_header_dictionary = {}
            for dada_header_field in dada_header:
                key, value = dada_header_field.split(" ")
                dada_header_dictionary[key] = value
            # Create output PSRCHIVE filename.
            # Get observing time.
            unix_time = dada_header_dictionary["UNIXTIME"]
            unix_time = Time(unix_time, format="unix")
            filename_time = unix_time.strftime("D%Y%m%dT%H%M%S")
            # Get_channel_number.
            channel_index, channel_frequency = find_nearest(
                frequency_table, float(dada_header_dictionary["FREQ"])
            )
            # Get source name.
            source = dada_header_dictionary["SOURCE"]
            # Create PSRCHIVE output filename.
            archive_filename = "{0:}_{1:}_channel_{2:03d}".format(
                source, filename_time, channel_index
            )
            # Fold the data with DSPSR.
            if args.cuda:
                args.threads = " --cuda 0 "
            completed_process_9 = subprocess.run(
                "dspsr -a PSRFITS "
                + args.set_key
                + " "
                + args.phase_bins
                + " -A "
                + args.integration
                + " "
                + args.telescope
                + " "
                + args.second_PFB
                + " -E "
                + args.ephem_file
                + " -O "
                + archive_filename
                + " "
                + args.seek
                + " "
                + args.total
                + " "
                + args.psrshell
                + " "
                + args.threads
                + " "
                + dada_file,
                shell=True,
                check=True,
                capture_output=True,
            )

    # End timing the script and output running time.
    script_end_time = time.perf_counter()
    print(
        "\nScript running time: {:.1f} s.\n".format(script_end_time - script_start_time)
    )
