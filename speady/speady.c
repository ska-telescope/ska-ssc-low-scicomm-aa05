/*
  Standalone simple C code to catch SKA SPEAD UDP packets and optionally reorder
  and fill in missing data. This code also supports reading data captured with tcpdump, in which case
  there is extra info in the file compared to what is captured from an actual socket

  This is test/debugging code for SKA-Low commissioning

  The SPEAD packet format is defined in:
  https://confluence.skatelescope.org/display/SE/SP-3800+CBF+SPEAD+format+according+to+ICD+version

  Each packet contains (among other thigns):
  - timestamp like thing (packet counter)
  - coarse channel index
  - 8192 bytes of complex voltages where each sample is 8 bit signed int real/imag, for both pols

  Hence each packet contains 8192/2/2 = 2048 time samples of data for a single channel, both pols
  Packets are not sent in logical channel order, but all packets for a given timestamp are sent before
  packets for the next timestamp.

  Randall Wayth. Dec 2024. randall.wayth@skao.int
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <arpa/inet.h>		// this seems to be required for sockaddr_in
#include <assert.h>
#include <errno.h>
#include <signal.h>
#include "ska_spead.h" 

// define double buffers for receiving
# define N_RX_BUFS 2
#define  N_POL     2

typedef struct {
  int8_t  *dat;            // allocated at runtime depending on how many channels are caught
  uint64_t pkt_since_full; // this contains the full-sized packet counter, which is 6 bytes in the header
  int	n_added;
} rx_buf;

FILE *fpin=NULL,*fpd=NULL,*fpout=NULL;
int mode=1;	//0== file input, 1== UDP input
int sock_fd=-1, debug=0,skip_start=24,skip_prefix=58,done=0,firstread=1,n_chan=8,max_packets=0,reorder_mode=1;
in_port_t port=4660;
struct sockaddr_in from_addr;
rx_buf rx_bufs[N_RX_BUFS],out_buf;	// double receive buffers
volatile int curr_rx_buf=0;

void print_usage(char * const argv[]) {
    fprintf(stderr,"Usage:\n%s [options]\n",argv[0]);
    fprintf(stderr,"\t-p port\t\tPort to listen on. Default: %d\n",(int)port);
    fprintf(stderr,"\t-m mode\t\tMode. 0==stdin input. 1==UDP input. Default: %d\n",(int)mode);
    fprintf(stderr,"\t-n nchan\tNum coarse chans to capture. Default: %d\n",(int)n_chan);
    fprintf(stderr,"\t\t\tIf nchan < number of channels in the data, then upper channels in data are discarded\n");
    fprintf(stderr,"\t-s num  \tStop after capturing num packets. No default.\n");
    fprintf(stderr,"\t-r mode \tReorder mode: 0: no reorder. 1: time,freq,pol. 2: freq,pol,time. Default: %d\n",reorder_mode);
    fprintf(stderr,"\t\t        Default data order from TPMs is [freq][time][pol]\n");
    fprintf(stderr,"\t-d      \twrite debug and runtime info to stderr\n");
    exit(1);
}


void parse_cmdline(int argc, char * const argv[]) {
    int c;
    char optstring[]="dp:m:n:s:r:";

    while ((c=getopt(argc,argv,optstring)) != -1) {
        switch(c) {
            case 'p':
                port = atoi(optarg);
                if (port <=0 || port > 65535) {
                    fprintf(stderr,"bad port: %d\n", port);
                    print_usage(argv);
                }
                break;
            case 'm':
                mode = atoi(optarg);
                if (mode <0 || mode > 1) {
                    fprintf(stderr,"bad mode: %d\n", mode);
                    print_usage(argv);
                }
                break;
            case 's':
                max_packets = atoi(optarg);
                if (max_packets <0 ) {
                    fprintf(stderr,"bad max_packets: %d\n", max_packets);
                    print_usage(argv);
                }
                break;
            case 'n':
                n_chan = atoi(optarg);
                if (n_chan <0 || n_chan > 96) {
                    fprintf(stderr,"bad n_chan: %d\n", n_chan);
                    print_usage(argv);
                }
                break;
            case 'd':
                debug += 1;
                break;
            default:
                fprintf(stderr,"unknown option %c\n",c);
                print_usage(argv);
        }
    }
    // set constants based on mode
    if(mode==1) {
        
    }
}

void sig_handler(int sig) {
    fprintf(stderr,"Received signal %d\n",sig);
    done=1;
}

int open_UDP(int port) {
    int fd;
    struct sockaddr_in local_addr;

    fd = socket(PF_INET, SOCK_DGRAM, 0);
    assert(fd != -1);

    local_addr.sin_family = AF_INET;
    local_addr.sin_addr.s_addr = INADDR_ANY;
    local_addr.sin_port = htons(port);
    if(bind(fd, (struct sockaddr *)&local_addr, sizeof(local_addr)) < 0){
        perror("bind");
        exit(1);
    }

    return fd;
}

void print_spead_header(ska_spead_t *hdr,FILE *fp) {
    //fprintf(fp,"Magic:\t%x\tNumitm:\t%d\t",(int)hdr->magic,(int)hdr->num_items);
    fprintf(fp,"Pktsnc:\t%u\tsncpad:\t%d\tscanid:\t%d\t",(unsigned int)hdr->packets_since,(int)hdr->packets_since_pad,(int)hdr->scan_id);
    fprintf(fp,"Lgchan:\t%d\tphschan:\t%d\n",(int)hdr->logical_chan_id,(int)hdr->physical_freq_channel);
}

/*
   unpack the spead header into the data structure. The structure is designed to be a direct map of the received data,
   but multi-byte items will need to be converted from network order to host order
   (and more work might be required for the 6-byte items)
*/
int decode_spead_header(uint8_t *in, ska_spead_t *out) {
    uint32_t t32;

    memcpy(out,in,sizeof(ska_spead_t));
    out->num_items = ntohl(out->num_items);
    out->packets_since_pad = ntohs(out->packets_since_pad);
    out->packets_since = ntohl(out->packets_since);
    out->scan_id = ntohl(out->scan_id);
    out->logical_chan_id = ntohs(out->logical_chan_id);
    out->physical_freq_channel = ntohs(out->physical_freq_channel);
    // not everything is unpacked yet, but can be done as required...
    return 0;
}


/* the order of the data from the packets into the rx buffer is [freq][time][pol]
   where there are effectively two time dimensions, because each buffer is a new time.
   We need to reorder the data to be [time][freq][pol].
   Since the data are 2-byte complex samples, and the pol is the innermost dimension in both,
   we can move data in 4-byte chunks
*/
void reorder_buf_time_outermost(void *in, void *out) {
    uint32_t *inp,*outp;
    inp = (uint32_t *)in;	// cast to fake arrays with 4-byte item sizes
    outp= (uint32_t *)out;
    for (int t=0; t<SKA_SPEAD_PAYLOAD_LEN/4; t++) {
        for(int c=0; c<n_chan; c++) {
            outp[t*n_chan + c] = inp[c*(SKA_SPEAD_PAYLOAD_LEN/4) + t];
        }
    }
}


/* alternative reordering to match requirement for DSPSR (?) 
   We need to reorder to [freq][pol][time]
*/
void reorder_buf_time_innermost(void *in, void *out) {
    uint16_t *inp,*outp;
    inp = (uint16_t *)in;	// cast to fake arrays with 2-byte item sizes
    outp= (uint16_t *)out;
    for(int c=0; c<n_chan; c++) {
        for (int t=0; t<SKA_SPEAD_PAYLOAD_LEN/2; t++) {
            for(int p=0; p<N_POL; p++) {
                outp[t + p*N_POL + c*SKA_SPEAD_PAYLOAD_LEN/2] = inp[c*(SKA_SPEAD_PAYLOAD_LEN/2) + t*N_POL + p];
            }
        }
    }
}

/* reorder and copy to output */
void reorder_buf(void *in, void *out) {

    switch(reorder_mode) {
        case 0:
            // no reorder, just copy
            memcpy(out,in,SKA_SPEAD_PAYLOAD_LEN*n_chan);
            break;
        case 1:
            reorder_buf_time_outermost(in,out);
            break;
        case 2:
            reorder_buf_time_innermost(in,out);
            break;
        default:
            fprintf(stderr,"Bad reorder mode %d\n",reorder_mode);
            exit(1);
            break;
    }
}

/*
  unpack the header and copy the payload of the packet into a receiver buffer
  a receive buffer will contain all samples for all freq channels for a single timestamp.
*/
int process_packet(void *pkt) {
    ska_spead_t spead_header;
    uint64_t pkt_since_full;

    //memset(&spead_header,0,sizeof(ska_spead_t)); // this not necessary if memcpy used in decode_spead_header

    decode_spead_header(pkt,&spead_header);
    if (debug>1) {
        print_spead_header(&spead_header,stderr);
    }
    // extract full (6-byte) packet counter into a 64-bit int
    pkt_since_full = ((uint64_t)spead_header.packets_since_pad<<32) + spead_header.packets_since;

    // on startup, wait for the start of a new batch of packets
    if (rx_bufs[curr_rx_buf].pkt_since_full==0 && spead_header.logical_chan_id==0) {
        // set current buffer to be catching packets with this timestamp
        rx_bufs[curr_rx_buf].pkt_since_full = pkt_since_full;
    }

    if (pkt_since_full < rx_bufs[curr_rx_buf].pkt_since_full) {
        // a packet is late and has arrived after we've moved to a new buffer... just discard it.
        fprintf(stderr,"WARNING: late packet with time %lu. Current: %lu\n",(unsigned long)pkt_since_full,(unsigned long)rx_bufs[curr_rx_buf].pkt_since_full);
        return 0;
    }

    if (pkt_since_full > rx_bufs[curr_rx_buf].pkt_since_full) {
        // check if the buffer was full before moving to next one...
        if(rx_bufs[curr_rx_buf].n_added < n_chan) {
            fprintf(stderr,"Only read %d packets for timestamp %llu\n",rx_bufs[curr_rx_buf].n_added,rx_bufs[curr_rx_buf].pkt_since_full);
        }
        reorder_buf(rx_bufs[curr_rx_buf].dat, out_buf.dat);
	fwrite(out_buf.dat,n_chan*SKA_SPEAD_PAYLOAD_LEN,1,fpout);

        // new timestamp, so move to the next buffer
	curr_rx_buf = (curr_rx_buf+1) % N_RX_BUFS;
        if (debug) {
            fprintf(fpd,"New packet timestamp: %llu. Moving to buffer %d\n",(unsigned long long) pkt_since_full,(int)curr_rx_buf);
        }
	rx_bufs[curr_rx_buf].pkt_since_full = pkt_since_full;
	rx_bufs[curr_rx_buf].n_added =0;
    }
    // insert packet payload into buffer
    if (spead_header.logical_chan_id < n_chan) {
        // dicard packet if we are capturing fewer channels than are in the data
        uint8_t *ptmp = ((uint8_t *)pkt + sizeof(ska_spead_t));
        memcpy(rx_bufs[curr_rx_buf].dat+spead_header.logical_chan_id*SKA_SPEAD_PAYLOAD_LEN,ptmp,SKA_SPEAD_PAYLOAD_LEN);
    }
    rx_bufs[curr_rx_buf].n_added += 1;

    return 0;
}

int rx_packet_udp() {
    char message[SKA_SPEAD_PAYLOAD_LEN+sizeof(ska_spead_t)+58]; // extra 58 here is probably not necessary, just checking
    socklen_t from_addr_len=sizeof(from_addr);
    int res=0,n_to_read;

    n_to_read = SKA_SPEAD_PAYLOAD_LEN+sizeof(ska_spead_t);
    // Receive message:
    res = recvfrom(sock_fd, message, n_to_read, 0, (struct sockaddr*)&from_addr, &from_addr_len);
    if (res <0 ){
        fprintf(stderr,"Couldn't receive\n");
        return -1;
    }
    if (debug>1) {
        fprintf(fpd,"Received message from IP: %s and port: %d with size %d bytes\n",
           inet_ntoa(from_addr.sin_addr), ntohs(from_addr.sin_port),res);
    }
    return process_packet(message);
}

int rx_packet_file() {
    int res=0,n_to_read;
    char buf[24+58+56+SKA_SPEAD_PAYLOAD_LEN];

    n_to_read = SKA_SPEAD_PAYLOAD_LEN+skip_prefix+sizeof(ska_spead_t);
    assert(n_to_read <= sizeof(buf));
    // discard the first 24 (skip_start) bytes if from file created by tcpdump
    if (firstread) {
        fread(buf,skip_start,1,fpin);
        firstread=0;
    }

    res = fread(buf,n_to_read,1,fpin);
    if (res <1) return -1;
    res = process_packet(buf+skip_prefix);
    return res;
}

int rx_packet(int mode) {
    int res=0;

    if (mode==1) {
        res = rx_packet_udp();
    } else {
        res = rx_packet_file();
    }
    return res;
}

int main(int argc, char* argv[]) {
    int res=0;
    uint64_t n_read=0;

    fpd = stderr;
    fpout = stdout;
    fpin = stdin;

    parse_cmdline(argc,argv);

    //signal(SIGINT,&sig_handler);
    signal(SIGHUP,&sig_handler);

    if (debug) fprintf(fpd,"size of SPEAD header: %ld\n",sizeof(ska_spead_t));

    int rx_buf_size=n_chan*SKA_SPEAD_PAYLOAD_LEN;
    for (int i=0; i< N_RX_BUFS; i++) {

        rx_bufs[i].dat = malloc(rx_buf_size);
	assert(rx_bufs[i].dat != NULL);
	rx_bufs[i].pkt_since_full = 0;
	rx_bufs[i].n_added=0;
    }
    out_buf.dat = malloc(rx_buf_size);
    assert(out_buf.dat != NULL);

    if (mode==1) {
        // open UDP listen socket
        sock_fd = open_UDP(port);
    }

    while (!done) {
        if(debug>1) {
             fprintf(fpd,"Reading packet %ld\n",(long)n_read);
        }
        res = rx_packet(mode);
        if (res != 0) {
            done=1;
            fprintf(fpd,"Done after %ld reads\n",(long)n_read);
        }
        n_read += 1;
	if (max_packets>0 && n_read>max_packets) {
            fprintf(stderr,"Quitting after reading %d packets\n",max_packets);
            done=1;
        }
    }
    return 0;
}

