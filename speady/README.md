# speady: capture and reorder SKA SPEAD beamformed data

speady is a simple tool to capture SKA SPEAD UDP packets that contain station beam data. It can also be used to reprocess previously captured files (what were created with tcpdump).

The SKA SPEAD packet format is [here](https://confluence.skatelescope.org/display/SE/SP-3800+CBF+SPEAD+format+according+to+ICD+version).

By default speady will write to stdout. When reading file input it will ready from stdin. Debug and info messages go to stderr.

The SKA SPEAD packets contain 2048 samples of complex voltages (8+8 bit real/imag) for a single channel per packet, with subsequent packets containing samples for different channels. This has implicit ordering of data as `[time][freq][time][pol]` where the outer time index is over later packets.
speady will reorder the data to standard streaming order

    [time][freq][pol]

### Usage
type ./speady -h for a usage message and command-line arguments

### Dependencies
a C compiler, make

### Build
just type make. There is no "make install", but the executable is standalone, so you can copy it anywhere if you like.


