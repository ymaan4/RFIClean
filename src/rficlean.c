/*
 * rficlean -- the overall wrapper code.
 *
 * Yogesh Maan <maan@astron.nl>  2018.
 * 
 */
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "rficlean.h"
#include "version.h"

/*
  RFICLEAN -- clean some RFI from filterbank files
*/

/*
TODO: better whitening
TODO: Get 0-dm timeseries before and after cleaning, and make histograms.
TODO: Get 8-subband rms, 0-dm timeseries post cleaning, use that to identify outlier samples.
TODO: 2,4,8-sample averaged versions of 0-dm timeseries to identify outliers.
TODO: Keep track of min,max,mean,rms for both time and frequency, and use that to produce channel and sample masks
TODO: write out above sample and channel masks in rfifind-output format
TODO: write out all the plot data as well as reduced data products

TODO: ftclean in totpow
TODO: check kurtosis/skewness for each sample
TODO: parallelization
TODO: module for reading and writing data in (psr)fits format
TODO: fit and plot gaussians over pre/post-cleaning histograms
TODO: Gaussianity tests of pre/post-cleaning histograms
*/


int help_required(char *string)
{
  if (strings_equal(string,"--help")) return(1);
  if (strings_equal(string,"-h")) return(1);
  return(0);
}

void rficlean_help()
{
  puts("");
  puts("rficlean - clean some periodic and/or narrow-band RFIs in filterbank data\n");
  puts("usage: rficlean -{options} {input-filename} \n");
  puts("options:\n");
  puts("       <filename>   - filterbank data file (def=stdin)");
  puts("-t     <numsamps>   - number of time samples in a block (for periodicity search; def=4096)");
  puts("-ft    <thres.>     - specify threshold for FT-cleaning (def=6.0)");
  puts("-st    <thres.>     - specify threshold for timeseries-cleaning (def=10.0)");
  puts("-rt    <thres.>     - specify threshold for chan-diff cleaning (def=4.0)");
  puts("-n     <numbits>    - specify output number of bits (def=input)");
  puts("-white              - whitten each spectrum before RFI excision");
  puts("-zeordm             - remove 0-dm powers before writing out data");
  puts("-pcl                - assume that input data may have some parts replaced");
  puts("                      by 0s or mean/median by some other RFI excision program");
  puts("-nharm <nh>         - Number of fundamental+harmonics to be removed");
  puts("                      (def=1, i.e., only fundamental)");
  puts("-T     <nsamp>      - number of time samples to avg before output (def=0)");
  puts("-gm    <filename>   - For gmrt raw-data, name of the file that contains header information");
  puts("-gmtstamp <filename> - For gmrt data, name of the file that contains time-stamp information");
  puts("-ps    <filename>   - File-name of the output diagnostic plot (def=rficlean_output.ps)");
  puts("-o     <filename>   - specify output filename (def=stdout)");
  puts("-psrf  <F0>         - fundamental freq. (Hz) of pulsar to be protected (def=none)");
  puts("-psrfdf <dF>        - Delta-frequency (on either side) of the fundamental & harmonic");
  puts("-psrfbins <Nb>      - *Nbins (on either side) of the fundamental & harmonic");
  puts("                      (*specify either psrfbins or psrfdf!)");
  puts("                      frequencies to be protected (def=2)");
  puts("-headerless         - do not broadcast resulting header (def=broadcast)");
  puts("-bst   <bstart>     - Starting block number to be processed (def=1, i.e. start of file)");
  puts("-nbl   <nblocks>    - No. of blocks to be processed (def=till EoF)");
  puts("");
}

int file_exists(char *filename)
{
  if ((fopen(filename,"rb"))==NULL) { return(0);}
  else { return(1); }
}

FILE *open_file(char *filename, char *descriptor)
{
  FILE *fopen(), *fptr;
  if ((fptr=fopen(filename,descriptor)) == NULL) {
    fprintf(stderr,"Error in opening file: %s\n",filename);
    exit(1);
  }
  return fptr;
}

void error_message(char *message)
{
  fprintf(stderr,"ERROR: %s\n",message);
  exit(1);
}

void print_version(char *program, char *argument)
{
  if ( (strings_equal(argument,"-v")) ||
       (strings_equal(argument,"--version"))) {
    printf("PROGRAM: %s is part of rfiClean version: %.1f\n",program,RFICLEAN_VERSION);
    exit(0);
  }
}


void main (int argc, char *argv[])
{
  int i, ibits, nc, headersize, headerless=0,gm=0;
  char string[80];

  /* set up default global variables */
  obits=headerless=naddt=nsamp=0;
  ibits=0;
  fthresh = 6.0;
  forcefthresh = 1000.0;
  rthresh = 4.0;
  sthresh = 10.0;
  clipthresh = 5.0;
  chanfrac = 0.85;  // the fraction of bad channels that will mask a full block/sample
  sampfrac = 0.85;  // the fraction of bad samples that will mask a full channel in a block
  nharm = 1;  // Max-Number of fundamental+harmonics to be removed
  fnorm = 0;
  tnorm = 0;
  pcl   = 0;
  zerodm = 0;
  iwhite = 0;
  psrf = 1000000.0;
  psrfbins = -1;
  psrfdf = -1.0;
  bl_start = 1;
  nblocks  = 9999999;
  byte_offset = 0;
  
  strcpy(inpfile,"dummy98151");
  strcpy(outfile,"dummy98151");
  strcpy(gmhdrfile,"dummy98151");
  strcpy(psfile,"rficlean_output.ps");

  if (argc > 1) {
    /* check command-line parameters */ 
    print_version(argv[0],argv[1]);
    i=1;
    while (i<argc) {
      if (strings_equal(argv[i],"-t")) {
	i++;
	naddt=atoi(argv[i]);
      } else if (strings_equal(argv[i],"-normf")) {
	fnorm=1;
      } else if (strings_equal(argv[i],"-normt")) {
	tnorm=1;
      } else if (strings_equal(argv[i],"-pcl")) {
	pcl=1;
      } else if (strings_equal(argv[i],"-zerodm")) {
	zerodm=1;
      } else if (strings_equal(argv[i],"-force")) {
	forcefthresh=-100;
      } else if (strings_equal(argv[i],"-nharm")) {
	i++;
	nharm=atoi(argv[i]);
      } else if (strings_equal(argv[i],"-o")) {
	/* get and open file for output */
        strcpy(outfile,argv[++i]);
	output=fopen(outfile,"wb");
      } else if (strings_equal(argv[i],"-ibits")) {
        i++;
        ibits=atoi(argv[i]);
      } else if (strings_equal(argv[i],"-T")) {
	i++;
	nsamp=atoi(argv[i]); 
      } else if (strings_equal(argv[i],"-gmtstamp")) {
	i++;
	gm=1;
	strcpy(gmhdrfile,argv[i]);
      } else if (strings_equal(argv[i],"-gm")) {
	i++;
	gm=1;
	strcpy(gminfofile,argv[i]);
      } else if (strings_equal(argv[i],"-ft")) {
	i++;
	fthresh=atof(argv[i]);
      } else if (strings_equal(argv[i],"-rt")) {
	i++;
	rthresh=atof(argv[i]);
      } else if (strings_equal(argv[i],"-st")) {
	i++;
	sthresh=atof(argv[i]);
      } else if (strings_equal(argv[i],"-n")) {
	i++;
	obits=atoi(argv[i]);
      } else if (strings_equal(argv[i],"-ps")) {
	i++;
	strcpy(psfile,argv[i]);
      } else if (strings_equal(argv[i],"-headerless")) {
	headerless=1;
      } else if (strings_equal(argv[i],"-white")) {
	iwhite=1;
      } else if (strings_equal(argv[i],"-psrfbins")) {
	i++;
	psrfbins=atoi(argv[i]);
      } else if (strings_equal(argv[i],"-psrfdf")) {
        i++;
        psrfdf=atof(argv[i]);
      } else if (strings_equal(argv[i],"-psrf")) {
	i++;
	psrf=atof(argv[i]);
      } else if (strings_equal(argv[i],"-bst")) {
	i++;
	bl_start=atoi(argv[i]);
        if(bl_start<1) {
          printf("*** Start block number (-bst) not valid! ***\n\n");
	  rficlean_help();
          exit(0);
        }
      } else if (strings_equal(argv[i],"-nbl")) {
	i++;
	nblocks=atoi(argv[i]);
      } else if (help_required(argv[1])) {
	rficlean_help();
	exit(0);
      } else if (file_exists(argv[i])) {
	strcpy(inpfile,argv[i]);
	input=open_file(inpfile,"rb");
      } else {
	rficlean_help();
	sprintf(string,"unknown argument (%s) passed to rficlean",argv[i]);
	error_message(string);
      }
      i++;
    }
  }
  else {
    rficlean_help();
    exit(0);
  }
  if ((strings_equal(inpfile,"dummy98151")) || (strings_equal(outfile,"dummy98151"))) {
    rficlean_help();
    exit(0);
  }

  // some sanity checks
  if (psrf<1000000.0 && fthresh<4.0 && forcefthresh>0) fthresh=4.0;
  if (psrfdf>0.0 && psrfbins > 0){
    printf ("Both psrfdf and psrfbins are specified!!");
    printf ("Using the delta-F information from psrfdf.");
    psrfbins = -1;
  }
  if (psrfdf<0.0 && psrfbins<0){
    printf ("Nether psrfdf nor psrfbins is specified!!");
    printf ("Using psrfbins=8");
    psrfbins = 8;
  }

  /* read in the header to establish what the input data are... */
  if (nifs>1){
   printf ("Cannot process data with IFs more than 1.");
   printf ("Halting !!");
   exit(0);
  }
  if (gm>0) {
    if (strings_equal(gmhdrfile,"dummy98151")) {
      strcpy(gmhdrfile, inpfile);
      strcat(gmhdrfile,".hdr");
    }
    printf ("\n Reading the GMRT header & time-stamp file...\n");
    read_gmheader(gminfofile, gmhdrfile);
    if (naddt <= 1) naddt=4096;
    if (ibits > 0) nbits=ibits;
    if (obits == 0) obits=nbits;
    if( nsamp>0) tsamp = tsamp*nsamp ; 
    if (!headerless) {
      printf (" Broadcasting the sigproc-format header...\n");
      bcast_header();
      printf (" Done! \n");
    }
    data_type=1;
    totsamp = (long long) (long double) (sizeof_file(inpfile))/ (((long double) nbits) / 8.0)
                 /(long double) nchans;
    byte_offset = (long) ((long)(bl_start-1)*naddt*nchans*(nbits/8.0));
  }
  else if ((headersize=read_header(input))) {
    printf (" Reading the sigproc-file header... \n");
    if(machine_id==8) machine_id=14; // change backend name from "PULSAR2000" to unknown
    totsamp = nsamples(inpfile,headersize,nbits,nifs,nchans);
    //printf("total samp: %ld\n",totsamp);
    if( nsamp>0) tsamp = tsamp*nsamp ; 
    switch (data_type) {
    case 1:
      break;
    case 2:
      nchans=1;
      break;
    case 6:
      break;
    default:
      error_message("ERROR: Input data to rficlean is not in filterbank format");
      break;
    }
    /* check number of time samples to process in one chunk */
    if (naddt <= 1) naddt=4096;
    if (obits == 0) obits=nbits;
    if (obits==1) fprintf(stderr,"WARNING output of 1-bit data will result in vastly reduced S/N!\nselect a higher output bit size with the -n option\n");
    /* all ok - broadcast the new header */
    if (!headerless) {
      printf (" Broadcasting the sigproc-format header...\n");
      bcast_header();
      printf (" Done! \n");
    }
    byte_offset = (long)(headersize + ((long)(bl_start-1)*naddt*nchans*(nbits/8.0)));
  } else {
    error_message("input data file is of unknown origin!!!");
  }
  
  /* finally clean and output the data */
  rficlean_data(input,output);

  printf("\n");
  exit(0);
}
