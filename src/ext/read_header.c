/* read_header.c - general handling routines for SIGPROC headers */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "header.h"
#include "rficlean.h"
#include <string.h>
int nbins;
double period;
/* read a string from the input which looks like nchars-char[1-nchars] */
void get_string(FILE *inputfile, int *nbytes, char string[]) /* includefile */
{
  int nchar;
  strcpy(string,"ERROR");
  fread(&nchar, sizeof(int), 1, inputfile);
  if (feof(inputfile)) exit(0);
  if (nchar>80 || nchar<1){
    printf("\nnchar: %d ; more than allowed limit of 80.\n",nchar);
    *nbytes=sizeof(int);
    fread(string, 79, 1, inputfile);
    printf("Partial string: '%s'\n",string);
    return;
  }
  *nbytes=sizeof(int);
  fread(string, nchar, 1, inputfile);
  string[nchar]='\0';
  *nbytes+=nchar;
}

/* attempt to read in the general header info from a pulsar data file */
int read_header(FILE *inputfile) /* includefile */
{
  char string[80], message[80];
  int itmp,nbytes,totalbytes,expecting_rawdatafile=0,expecting_source_name=0; 
  int expecting_frequency_table=0,channel_index;


  /* try to read in the first line of the header */
  get_string(inputfile,&nbytes,string);
  if (!strings_equal(string,"HEADER_START")) {
	/* the data file is not in standard format, rewind and return */
	rewind(inputfile);
	return 0;
  }
  /* store total number of bytes read so far */
  totalbytes=nbytes;

  /* loop over and read remaining header lines until HEADER_END reached */
  while (1) {
    get_string(inputfile,&nbytes,string);
    if (strings_equal(string,"HEADER_END")) break;
    totalbytes+=nbytes;
    if (strings_equal(string,"rawdatafile")) {
      expecting_rawdatafile=1;
    } else if (strings_equal(string,"source_name")) {
      expecting_source_name=1;
    } else if (strings_equal(string,"FREQUENCY_START")) {
      expecting_frequency_table=1;
      channel_index=0;
    } else if (strings_equal(string,"FREQUENCY_END")) {
      expecting_frequency_table=0;
    } else if (strings_equal(string,"az_start")) {
      fread(&az_start,sizeof(az_start),1,inputfile);
      totalbytes+=sizeof(az_start);
    } else if (strings_equal(string,"za_start")) {
      fread(&za_start,sizeof(za_start),1,inputfile);
      totalbytes+=sizeof(za_start);
    } else if (strings_equal(string,"src_raj")) {
      fread(&src_raj,sizeof(src_raj),1,inputfile);
      totalbytes+=sizeof(src_raj);
    } else if (strings_equal(string,"src_dej")) {
      fread(&src_dej,sizeof(src_dej),1,inputfile);
      totalbytes+=sizeof(src_dej);
    } else if (strings_equal(string,"tstart")) {
      fread(&tstart,sizeof(tstart),1,inputfile);
      totalbytes+=sizeof(tstart);
    } else if (strings_equal(string,"tsamp")) {
      fread(&tsamp,sizeof(tsamp),1,inputfile);
      totalbytes+=sizeof(tsamp);
    } else if (strings_equal(string,"period")) {
      fread(&period,sizeof(period),1,inputfile);
      totalbytes+=sizeof(period);
    } else if (strings_equal(string,"fch1")) {
      fread(&fch1,sizeof(fch1),1,inputfile);
      totalbytes+=sizeof(fch1);
    } else if (strings_equal(string,"fchannel")) {
      fread(&frequency_table[channel_index++],sizeof(double),1,inputfile);
      totalbytes+=sizeof(double);
      fch1=foff=0.0; /* set to 0.0 to signify that a table is in use */
    } else if (strings_equal(string,"foff")) {
      fread(&foff,sizeof(foff),1,inputfile);
      totalbytes+=sizeof(foff);
    } else if (strings_equal(string,"nchans")) {
      fread(&nchans,sizeof(nchans),1,inputfile);
      totalbytes+=sizeof(nchans);
    } else if (strings_equal(string,"telescope_id")) {
      fread(&telescope_id,sizeof(telescope_id),1,inputfile);
      totalbytes+=sizeof(telescope_id);
    } else if (strings_equal(string,"machine_id")) {
      fread(&machine_id,sizeof(machine_id),1,inputfile);
      totalbytes+=sizeof(machine_id);
    } else if (strings_equal(string,"data_type")) {
      fread(&data_type,sizeof(data_type),1,inputfile);
      totalbytes+=sizeof(data_type);
    } else if (strings_equal(string,"ibeam")) {
      fread(&ibeam,sizeof(ibeam),1,inputfile);
      totalbytes+=sizeof(ibeam);
    } else if (strings_equal(string,"nbeams")) {
      fread(&nbeams,sizeof(nbeams),1,inputfile);
      totalbytes+=sizeof(nbeams);
    } else if (strings_equal(string,"nbits")) {
      fread(&nbits,sizeof(nbits),1,inputfile);
      totalbytes+=sizeof(nbits);
    } else if (strings_equal(string,"barycentric")) {
      fread(&barycentric,sizeof(barycentric),1,inputfile);
      totalbytes+=sizeof(barycentric);
    } else if (strings_equal(string,"pulsarcentric")) {
      fread(&pulsarcentric,sizeof(pulsarcentric),1,inputfile);
      totalbytes+=sizeof(pulsarcentric);
    } else if (strings_equal(string,"nbins")) {
      fread(&nbins,sizeof(nbins),1,inputfile);
      totalbytes+=sizeof(nbins);
    } else if (strings_equal(string,"nsamples")) {
      /* read this one only for backwards compatibility */
      fread(&itmp,sizeof(itmp),1,inputfile);
      totalbytes+=sizeof(itmp);
    } else if (strings_equal(string,"nifs")) {
      fread(&nifs,sizeof(nifs),1,inputfile);
      totalbytes+=sizeof(nifs);
    } else if (strings_equal(string,"npuls")) {
      fread(&npuls,sizeof(npuls),1,inputfile);
      totalbytes+=sizeof(npuls);
    } else if (strings_equal(string,"refdm")) {
      fread(&refdm,sizeof(refdm),1,inputfile);
      totalbytes+=sizeof(refdm);
    } else if (expecting_rawdatafile) {
      strcpy(rawdatafile,string);
      expecting_rawdatafile=0;
    } else if (expecting_source_name) {
      strcpy(source_name,string);
      expecting_source_name=0;
    } else {
      sprintf(message,"read_header - unknown parameter: %s\n",string);
      fprintf(stderr,"ERROR: %s\n",message);
      exit(1);
    } 
  } 

  /* add on last header string */
  totalbytes+=nbytes;

  /* return total number of bytes read */
  return totalbytes;
}



/* read in the general header info and timestamp for GMRT-obs from text files*/
int read_gmheader(char gminfofile[], char gmhdrfile[]) /* includefile */
{
  char string[80], message[80], str1[32],str2[32],str3[32];
  int expecting_source_name=0; 
  int hh,mm,day,month,year,nyear,nday,istat;
  double signed_bwidth,bwidth;
  double seconds,mjd_day,mjd_fracday;
  FILE *infofile,*timefile;

  infofile = fopen(gminfofile,"rb");
  fscanf(infofile, "%lf \n", &tsamp);
  tsamp = tsamp/1000.0 ; // now in seconds
  fscanf(infofile, "%lf \n", &fch1);
  fscanf(infofile, "%lf \n", &signed_bwidth);
  fscanf(infofile, "%d \n", &nchans);
  fscanf(infofile, "%s \n", source_name);
  if( fscanf(infofile, "%lf \n", &src_raj) == EOF) src_raj=000000.00;
  if( fscanf(infofile, "%lf \n", &src_dej) == EOF) src_dej=000000.00;
  fclose(infofile);

  iflip = 0;
  if (signed_bwidth > 0) iflip=1 ; //we need to flip the band before writing
  bwidth = fabs(signed_bwidth);
  foff = -bwidth/(double)nchans;
  expecting_source_name=1;
  nbits = 16;
  nifs = 1;
  nbeams = 1;
  ibeam = 1;
  telescope_id=7;   // for GMRT
  //machine_id=17; //unknown
  machine_id=7; //GMRTFB

  timefile = fopen(gmhdrfile,"rb");
  fscanf(timefile, "%s %s %s %s\n", string,string,string,string);
  fscanf(timefile, "%s %s %d:%d:%lf  \n", str1,str2,&hh,&mm,&seconds);
  fscanf(timefile, "%s  %d:%d:%d  \n", str1,&day,&month,&year);
  fclose(timefile);

  slaCldj(year,month,day, &mjd_day, &istat);
  mjd_fracday = (hh + (mm + (seconds / 60.0)) / 60.0) / 24.0;
  tstart = mjd_day + mjd_fracday - 5.5/24.0 ;  // correct for 5.5 hours diff between IST and UT
  //printf("mjd_day, mjd_fracday, tstart  %lf  %lf  %.12lf\n",mjd_day, mjd_fracday, tstart);

  return 1;
}

void slaCldj ( int iy, int im, int id, double *djm, int *j ) /*includefile*/
/*
**  - - - - - - - -
**   s l a C l d j
**  - - - - - - - -
**
**  Gregorian calendar to Modified Julian Date.
**
**  Given:
**     iy,im,id     int    year, month, day in Gregorian calendar
**
**  Returned:
**     *djm         double Modified Julian Date (JD-2400000.5) for 0 hrs
**     *j           int    status:
**                           0 = OK
**                           1 = bad year   (MJD not computed)
**                           2 = bad month  (MJD not computed)
**                           3 = bad day    (MJD computed)
**
**  The year must be -4699 (i.e. 4700BC) or later.
**
**  The algorithm is derived from that of Hatcher 1984 (QJRAS 25, 53-55).
**
**  Last revision:   29 August 1994
**
**  Copyright P.T.Wallace.  All rights reserved.
*/
{
   long iyL, imL;

/* Month lengths in days */
   static int mtab[12] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };



/* Validate year */
   if ( iy < -4699 ) { *j = 1; return; }

/* Validate month */
   if ( ( im < 1 ) || ( im > 12 ) ) { *j = 2; return; }

/* Allow for leap year */
   mtab[1] = ( ( ( iy % 4 ) == 0 ) &&
             ( ( ( iy % 100 ) != 0 ) || ( ( iy % 400 ) == 0 ) ) ) ?
             29 : 28;

/* Validate day */
   *j = ( id < 1 || id > mtab[im-1] ) ? 3 : 0;

/* Lengthen year and month numbers to avoid overflow */
   iyL = (long) iy;
   imL = (long) im;

/* Perform the conversion */
   *djm = (double)
        ( ( 1461L * ( iyL - ( 12L - imL ) / 10L + 4712L ) ) / 4L
        + ( 306L * ( ( imL + 9L ) % 12L ) + 5L ) / 10L
        - ( 3L * ( ( iyL - ( 12L - imL ) / 10L + 4900L ) / 100L ) ) / 4L
        + (long) id - 2399904L );
}

