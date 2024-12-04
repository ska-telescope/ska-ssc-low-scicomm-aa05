#include <stdint.h>

#define SKA_SPEAD_PAYLOAD_LEN 8192

typedef struct {
  uint8_t	magic;			// fixed byte = 0x53
  uint8_t	version;		// "fixed" as 0x04;
  uint8_t	item_prt_width;		// "fixed" as 0x02;
  uint8_t	heap_addr_width;	// "fixed" as 0x06;
  uint32_t	num_items;		// possibly only 2 bytes of a 4 byte slot used.
  int16_t	heap_counter;		
  uint16_t	packets_since_pad;	// packets since is actually 6 bytes...
  uint32_t	packets_since;
  int16_t	packet_length;
  uint16_t	payload_length_pad;	// payload length fills 6 bytes...
  uint32_t	packet_payload_length;
  int16_t	scan_id_pad;		// scan ID is 6 bytes...
  uint32_t	scan_id;
  int16_t	csp_chan_info;
  int16_t	logical_chan_id;	// not sure if these can meaningfully be negative.
  int16_t	subarray_beam_id;
  int16_t	physical_freq_channel;
  int16_t	csp_antenna_info;
  uint8_t	substation_id;
  uint8_t	subarray_id;
  uint16_t	station_id;
  uint16_t	aperture_id;
  uint16_t	samples;
  uint16_t	offset_pad1;		// these don't seem to be used. 6 bytes.
  uint32_t	offset_pad2;
  //int8_t	data[8192];
} ska_spead_t;

