#!/usr/bin/env python3


import argparse
import ctypes
import io
import socket
import struct
import sys
import time

import astropy.units as u
import jinja2
import numpy as np
import pcapyplus
from astropy.time import Time


def format_mac_address(mac_string):
    return ":".join("%02x" % b for b in bytearray(mac_string))


def read_ethII_header(stream):
    header = {}
    header["dest_mac_addr"] = format_mac_address(stream.read(6))
    header["src_mac_addr"] = format_mac_address(stream.read(6))
    header["frame_length"] = struct.unpack("!H", stream.read(2))[0]
    return header


def read_ipv4_header(stream):
    header = {}
    stream.read(2)
    header["udp_length"] = struct.unpack("!H", stream.read(2))[0]
    stream.read(8)
    header["src_addr"] = socket.inet_ntoa(stream.read(4))
    header["dest_addr"] = socket.inet_ntoa(stream.read(4))
    return header


def read_udp_header(stream):
    header = {}
    header["src_port"] = struct.unpack("!H", stream.read(2))[0]
    header["dest_port"] = struct.unpack("!H", stream.read(2))[0]
    header["length"] = struct.unpack("!H", stream.read(2))
    header["checksum"] = struct.unpack("!H", stream.read(2))
    return header


def read_spead_header(stream):
    header = {}
    header["magic_number"] = struct.unpack("b", stream.read(1))[0]
    header["version"] = struct.unpack("b", stream.read(1))[0]
    header["item_pointer_width"] = struct.unpack("b", stream.read(1))[0]
    header["heap_addr_width"] = struct.unpack("b", stream.read(1))[0]
    header["reserved"] = struct.unpack("!H", stream.read(2))[0]
    header["num_items"] = struct.unpack("!H", stream.read(2))[0]
    return header


class Descriptor(ctypes.BigEndianStructure):
    _fields_ = [
        ("is_value", ctypes.c_uint8, 1),
        ("id", ctypes.c_uint16, 15),
        ("value", ctypes.c_uint64, 48),
    ]


class DescriptorUnion(ctypes.Union):
    _fields_ = [("struct", Descriptor), ("uint64", ctypes.c_uint64)]


# Descriptor of the LFAA SPEAD packet. This descriptor is based on
# SKA1 LOW LFAA to CSP Interface Control Document (Number 100-000000-004
# Rev. 3 and later) and deviates from the SPEAD standard significantly.
# The details of the changes to the packet are described in
# https://confluence.skatelescope.org/display/SE/SP-3800+CBF+SPEAD+format+according+to+ICD+version
descriptor_map = {
    1: "heap_counter",
    4: "pkt_len",
    12304: "scan_id",
    45056: "csp_channel_info",
    45057: "csp_antenna_info",
    13056: "samples",
}


def parse_descriptor_id(desc_id):
    return descriptor_map.get(desc_id, "unknown")


def read_spead_descriptor(stream):
    u = DescriptorUnion()
    u.uint64 = struct.unpack("Q", stream.read(8))[0]
    descriptor_id = u.struct.id
    value_or_address = u.struct.value
    return {parse_descriptor_id(descriptor_id): value_or_address}


NSAMPS_PER_PACKET = 2048
NPOL = 2


def read_spead_data(stream, n):
    bad_spead_data = False
    try:
        spead_data = np.array(struct.unpack(n * "b", stream.read(n))).reshape(
            NSAMPS_PER_PACKET, NPOL, 2
        )
    except struct.error:
        bad_spead_data = True
        spead_data = np.zeros(n).reshape(NSAMPS_PER_PACKET, NPOL, 2)
    return spead_data, bad_spead_data


def parse_packet(data):
    header = {}  # Create empty dictionary to hold a packet and data.
    stream = io.BytesIO(
        data
    )  # Create a buffered I/O implementation using in-memory bytes buffer (object).
    header.update(read_ethII_header(stream))
    header.update(read_ipv4_header(stream))
    header.update(read_udp_header(stream))
    header.update(
        read_spead_header(stream)
    )  # Read part of packet header adhering to the SPEAD definition.
    header.update(read_spead_descriptor(stream)),  # Read heap counter.
    header.update(read_spead_descriptor(stream)),  # Read packet length.
    header.update(read_spead_descriptor(stream)),  # Read scan ID.
    # Read parts of a packet deviating from SPEAD definition.
    header["csp_channel_info"] = struct.unpack("!H", stream.read(2))[
        0
    ]  # Read pointer to CSP channel info.
    header["logical_channel_id"] = struct.unpack("!H", stream.read(2))[0]
    header["subarray_beam_id"] = struct.unpack("!H", stream.read(2))[0]
    header["physical_frequency_channel"] = struct.unpack("!H", stream.read(2))[0]
    header["csp_antenna_info"] = struct.unpack("!H", stream.read(2))[
        0
    ]  # Read pointer to CSP antenna info.
    header["substation_id"] = struct.unpack("!b", stream.read(1))[0]
    header["subarray_id"] = struct.unpack("!b", stream.read(1))[0]
    header["station_id"] = struct.unpack("!H", stream.read(2))[0]
    header["reserved_aperture_id"] = struct.unpack("!H", stream.read(2))[0]
    header.update(read_spead_descriptor(stream))  # Read payload offset.
    header["data"], bad_spead_data = read_spead_data(stream, 8192)
    if bad_spead_data is True:
        print("\nPacket from {0:} was replaced with zeros.".format(header["src_addr"]))
    stream.close()  # Close the object.
    return header


DADA_HEADER = """HDR_VERSION 1.0
HDR_SIZE 4096
TELESCOPE {{telescope_name}}
STATION_NAME {{station_name}}
SUBARRAY_ID {{subarray_id}}
STATION_ID {{station_id}}
SUBSTATION_ID {{substation_id}}
RECEIVER {{receiver_name}}
INSTRUMENT {{instrument_name}}
FREQ {{frequency}}
COARSE_CHANNEL 1
LOGICAL_CHANNEL_ID {{logical_channel_id}}
PHYSICAL_FREQUENCY_CHANNEL {{physical_frequency_channel}}
BW {{bandwidth}}
BANDWIDTH_HZ 925925.937500
SAMPLE_RATE 925925.937500
FINE_CHAN_WIDTH_HZ 925925.937500
TSAMP 1.08
MODE PSR
SOURCE {{source_name}}
RA {{ra}}
DEC {{dec}}
UTC_START {{utc_start_time}}
PICOSECONDS {{picoseconds}}
UNIXTIME {{unix_time}}
UNIXTIME_MSEC {{unix_time_msec}}
OBS_OFFSET 0
HEAP_COUNTER {{heap_counter}}
SCAN_ID {{scan_id}}
SRC_MAC_ADDR {{src_mac_addr}}
SRC_ADDR {{src_addr}}
SRC_PORT {{src_port}}
DEST_MAC_ADDR {{dest_mac_addr}}
DEST_ADDR {{dest_addr}}
DEST_PORT {{dest_port}}
RESOLUTION 8192
NBIT {{nbit}}
NPOL {{npol}}
NCHAN {{nchan}}
NDIM 2
NTIMESAMPLES 1
NINPUTS 2
NINPUTS_XGPU 2
METADATA_BEAMS 2
BYTES_PER_SECOND {{bytes_per_second}}

"""


DADA_DEFAULTS = {
    "telescope_name": "SKA_Low",
    "station_name": "S8-1",
    "subarray_id": 1,
    "station_id": 1,
    "substation_id": 1,
    "receiver_name": "LFAASP",
    "instrument_name": "LFAASP",
    "frequency": 50.0,
    "nchan": 1,
    "npol": 2,
    "nbit": 8,
    "logical_channel_id": 0,
    "physical_frequency_channel": 64,
    "source_name": "J0000-0000",
    "ra": "00:00:00.00",
    "dec": "00:00:00.00",
    "utc_start_time": "1999-12-31-23:59:28",
    "picoseconds": 0,
    "unix_time": 946684768,
    "unix_time_msec": 0,
    "heap_counter": 0,
    "scan_id": 0,
    "src_mac_addr": "00:00:00:00:00:00",
    "src_addr": "0.0.0.0",
    "src_port": 0,
    "dest_mac_addr": "00:00:00:00:00:00",
    "dest_addr": "0.0.0.0",
    "dest_port": 0,
    "bandwidth": 0.9259259375,
    "bytes_per_second": 3703703.75,
}


STATION_NAME = {
    345: "S8-1",
    346: "S8-2",
    347: "S8-3",
    348: "S8-4",
    349: "S8-5",
    350: "S8-6",
    351: "S9-1",
    352: "S9-2",
    353: "S9-3",
    354: "S9-4",
    355: "S9-5",
    356: "S9-6",
    429: "S10-1",
    430: "S10-2",
    431: "S10-3",
    432: "S10-4",
    433: "S10-5",
    434: "S10-6",
}


def dada_defaults():
    return DADA_DEFAULTS.copy()


def render_dada_header(overrides):
    defaults = dada_defaults()
    defaults.update(overrides)
    return jinja2.Template(DADA_HEADER).render(**defaults)


class ChannelRingBuffer(
    object
):  # Instantiate object which will create and populate a single PSRDADA file.
    def __init__(self, physical_frequency_channel, args):
        print("\nNew ChannelRingBuffer instance created.")
        print("Physical Frequency Channel ID: {:d}".format(physical_frequency_channel))
        # This is a flag when the object is instantiated to first create PSRDADA header.
        self._first_pass = True
        # Create PSRDADA filename.
        if args.source_name:
            filename = "{:s}_{:s}_{:03d}.dada".format(
                args.source_name, args.prefix, physical_frequency_channel
            )
        else:
            filename = "{:s}_{:03d}.dada".format(
                args.prefix, physical_frequency_channel
            )
        self._file = open(filename, "wb")
        print("Output filename: {:s}".format(filename))

    def _write_dada_header(self, packet):
        dada_header = (
            dada_defaults()
        )  # Use PSRDADA dictionary with default header values.
        dada_header["src_addr"] = packet["src_addr"]
        dada_header["subarray_id"] = packet["subarray_id"]
        dada_header["station_id"] = packet["station_id"]
        dada_header["substation_id"] = packet["substation_id"]
        dada_header["logical_channel_id"] = packet["logical_channel_id"]
        dada_header["physical_frequency_channel"] = packet["physical_frequency_channel"]
        dada_header["heap_counter"] = packet["heap_counter"]
        dada_header["scan_id"] = packet["scan_id"]
        dada_header["src_mac_addr"] = packet["src_mac_addr"]
        dada_header["src_addr"] = packet["src_addr"]
        dada_header["src_port"] = packet["src_port"]
        dada_header["dest_mac_addr"] = packet["dest_mac_addr"]
        dada_header["dest_addr"] = packet["dest_addr"]
        dada_header["dest_port"] = packet["dest_port"]
        dada_header["station_name"] = STATION_NAME[packet["station_id"]]
        ADC_clock = 800
        ADC_samples = 1024
        dada_header["frequency"] = packet["physical_frequency_channel"] * (
            ADC_clock / ADC_samples
        )
        dada_header["ra"] = args.right_ascension
        dada_header["dec"] = args.declination
        dada_header["source_name"] = args.source_name
        oversampling_ratio = 32 / 27  # Frequency oversampling ratio.
        sample_time = 1 / (
            (ADC_clock / ADC_samples) * oversampling_ratio
        )  # Calculate sampling time.
        telescope_epoch = Time(
            "2000-01-01T00:00:00.000", format="isot", scale="tai"
        )  # Define telescope epoch as TAI 2000.0
        tai_start_time = telescope_epoch + (
            packet["heap_counter"] * sample_time * 2048 * u.microsecond
        )  # Calculate packet timestamp.
        utc_start_time = tai_start_time.utc
        utc_start_time.precision = 0
        dada_header["utc_start_time"] = utc_start_time.iso.replace(" ", "-")
        dada_header["unix_time"] = int(tai_start_time.utc.unix)
        dada_header["unix_time_msec"] = tai_start_time.utc.datetime.microsecond * 1e-3
        dada_header["picoseconds"] = tai_start_time.utc.datetime.microsecond * 1e6
        dada_header["nchan"] = 1
        dada_header["npol"] = NPOL
        dada_header["nbit"] = 8
        self._file.write(render_dada_header(dada_header).encode("ascii"))
        # 4096 bytes is a fixed PSRDADA header length.
        self._file.seek(4096)

    def add(self, packet):
        if self._first_pass:
            self._write_dada_header(packet)
            self._first_pass = (
                False  # Change flag value after PSRDADA file has been created.
            )
        self._data = np.zeros((NSAMPS_PER_PACKET, NPOL, 2), dtype=np.int8)
        self._data = packet["data"].astype(np.int8)  # Re-cast packet data to int8.
        self._data.tofile(self._file)

    def close(self):
        self._file.close()


class RingBufferManager(
    object
):  # Instantiate object which will instantiate another object (ChannelRingBuffer).
    def __init__(self, args):
        self._buffers = {}  # Create empty dictionary for ring buffers.
        self._args = args  # Pass arguments from options.

    def add(self, packet):  # Read and separate packets based on their frequency.
        physical_frequency_channel = packet["physical_frequency_channel"]
        # Check if physical_frequency_channel is not in _buffers dictionary (dictionary of objects), and create it if necessary.
        if physical_frequency_channel not in self._buffers:
            self._buffers[physical_frequency_channel] = ChannelRingBuffer(
                physical_frequency_channel, self._args
            )
        # Add a packet to a specific ChannelRingBuffer, basically sort packets per frequency.
        self._buffers[physical_frequency_channel].add(packet)

    def close_all(self):  # Close all objects.
        for key, rb in self._buffers.items():
            rb.close()


def stream_to_buffers(packet_generator, args):
    ring_buffer_manager = RingBufferManager(args)  # Create RingBufferManager object.
    max_packets = args.npackets
    packet_counter = 0
    # Start parsing packets.
    while True:
        packet = (
            packet_generator.next()
        )  # Might generate error if a pcap file is cut in the middle of a packet.
        # Check if packet header is empty (end of pcap file).
        if packet[0] is None:
            print("\n{:d} packets processed.".format(packet_counter))
            break
        else:
            packet_header = parse_packet(packet[1])
            if "physical_frequency_channel" in packet_header:
                ring_buffer_manager.add(
                    packet_header
                )  # Call parse_packet to read packet contents (reads entire packet).
            else:
                print(
                    '\nSkipping packet without "physical_frequency_channel" coming from {:s}'.format(
                        packet_header["src_addr"]
                    )
                )
        packet_counter += 1
        if (
            max_packets is not None
        ):  # Check if maximum requested number of packets was processed.
            if packet_counter >= max_packets:
                print(
                    "\nProcessed {:d} out of {:d} requested packets.".format(
                        packet_counter, max_packets
                    )
                )
                break
    # Call close_all function on ring_buffer_manager object (RingBufferManager class).
    ring_buffer_manager.close_all()
    return ring_buffer_manager


def main(args):
    # Open file and return the packet_generator object.
    packet_generator = pcapyplus.open_offline(args.pcap_file)
    # Call stream_to_buffers function.
    return stream_to_buffers(packet_generator, args)


if __name__ == "__main__":
    # This code will take care of negative values for declination.
    for i, arg in enumerate(sys.argv):
        if (arg[0] == "-") and arg[1].isdigit():
            sys.argv[i] = " " + arg
    parser = argparse.ArgumentParser(
        usage="lfaa_pcap_to_psrdada.py --file <pcap_file> [options]".format(
            prog=sys.argv[0]
        ),
        description="Convert LFAA (single SKA Low station) pcap data into PSRDADA files.",
        formatter_class=lambda prog: argparse.HelpFormatter(
            prog, max_help_position=100, width=250
        ),
        epilog="Copyright (C) 2025 by SKA Observatory, Maciej Serylak (maciej.serylak@skao.int)",
    )
    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "--file",
        dest="pcap_file",
        metavar="<pcap_file>",
        type=str,
        required=True,
        help="name of the pcap file to read",
    )
    parser.add_argument(
        "--source",
        dest="source_name",
        metavar="<source_name>",
        type=str,
        default=None,
        help="source name (default: None)",
    )
    parser.add_argument(
        "--ra",
        dest="right_ascension",
        metavar="<right_ascension>",
        type=str,
        default=None,
        help="source right ascension (default: None)",
    )
    parser.add_argument(
        "--dec",
        dest="declination",
        metavar="<declination>",
        type=str,
        default=None,
        help="source declination (default: None)",
    )
    parser.add_argument(
        "--prefix",
        dest="prefix",
        metavar="<prefix>",
        type=str,
        default="channel",
        help="prefix for output file names, (default: {source}_{prefix}_{physical_frequency_channel}.dada)",
    )
    parser.add_argument(
        "--npackets",
        dest="npackets",
        metavar="<npackets>",
        type=int,
        default=None,
        help="number of packets to read from file (default: read all packets)",
    )
    args = parser.parse_args()

    # Start timing the script.
    script_start_time = time.perf_counter()

    # Run the main function.
    main(args)

    # End timing the script and output running time.
    script_end_time = time.perf_counter()
    print(
        "\nScript running time: {:.1f} s.\n".format(script_end_time - script_start_time)
    )
