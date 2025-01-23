# speady: capture and reorder SKA SPEAD beamformed data

speady is a simple tool to capture SKA SPEAD UDP packets that contain station beam data. It can also be used to reprocess previously captured .pcap files (that were created with tcpdump).

The SKA SPEAD packet format is [here](https://confluence.skatelescope.org/display/SE/SP-3800+CBF+SPEAD+format+according+to+ICD+version).

By default speady will write to stdout. When reading file input it will ready from stdin. Debug and info messages go to stderr.

The SKA SPEAD packets contain 2048 samples of complex voltages (8+8 bit real/imag) for a single channel per packet, with subsequent packets containing samples for different channels. This has implicit ordering of data as `[time][freq][time][pol]` where the outer time index is over later packets.

speady offers several options for reordering data as required, using the -r command-line option.
By default, speady will reorder the data to standard streaming order, where the time jumps from packet boundaries are removed

    [time][freq][pol]

The order expected by dspsr is the native packet order, hence use -r 0 -n 1 to match the format of a .data file for a single coarse channel.

speady will also optionally insert data from missing packets via the -P option, where the inserted data is just a copy of the previous (not zeroes).

### Usage
type ./speady -h for a usage message and command-line arguments

### Dependencies
a C compiler, make

### Build
just type make. There is no "make install", but the executable is standalone, so you can copy it anywhere you like.

### Examples

Capture the default number of channels from UDP packets on the default port with default reordering, writing to a file:

speady > /tmp/myfile

Capture and send to next process in the a pipeline (e.g a spectrometer, which reads from stdin):

speady | my_spectrometer 

Reprocess a pcap file, extracting the first coarse channel, and write to a file that can have a dada header prepended:

speady -m 0 -n 1 -r 0 < raw_data.pcap > output_single_channel

