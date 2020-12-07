/*
 *
 * rficlean_data --- handle reading/writing of data blocks, and
 *                   intitialization etc.
 *
 * Yogesh Maan <maan@astron.nl>   2018
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "rficlean.h"
#include <fftw3.h>
#include <omp.h>
#include <time.h>

void rficlean_data(FILE *input, FILE *output)
{
  char string[80],plotdevice[100];
  float *fblock,min,max,hpower,realtime,fsaved[2], *ts0dm, atemp;
  unsigned short *sblock;
  unsigned char  *cblock;
  int nsaved=0,opened=0, wpout=10;
  long int ns,nsblk,nout,iter,i,j,k, jt1, jt2, iblock;
  long int itemp, isum, nsize, istart, ii, jj, kk, n0;
  time_t s, start, cpu_time;
  clock_t t_global, t_function, t_function_sum;
  int i_function;
  i_function = 0;

  s = time(NULL);
  t_global = clock();

  last_tvar = -1.0;
  last_fvar = -1.0;
  nsect = 32;
  nvar = 0;
  tot_now = 0;
  maxhist = 256;
  numhist = 64;
  nsub = 4;
  nsub = 1;
  nsubchans = nchans/nsub ;
  strcpy(plotdevice, psfile);
  ////strcat(plotdevice,"/CPS");
  nsize = naddt;
  if(nsize < nchans) nsize = nchans;
  nsize = 2*nsize;

  printf ("\n Preparing all buffers...\n");
/***/
    // chandata = (double *) malloc(nsize*sizeof(double));
    vspec = (double *) malloc(nsize*sizeof(double));
    rspec = (double *) malloc(nsize*sizeof(double));
  coff = (long int *) malloc((nchans)*sizeof(long int));
  wt = (double *) malloc(nsize*sizeof(double));
  mspec = (double *) malloc(nsize*sizeof(double));
  wspec = (double *) malloc(nsize*sizeof(double));
/***/

  // in = fftw_malloc (sizeof(fftw_complex)*(naddt+5));
  // out = fftw_malloc (sizeof(fftw_complex)*(naddt+5));
  // fplan = fftw_plan_dft_1d( naddt, in, out, FFTW_FORWARD, FFTW_ESTIMATE );
  // bplan = fftw_plan_dft_1d( naddt, in, out, FFTW_BACKWARD, FFTW_ESTIMATE );
  // ai = (double *) malloc(sizeof(double)*(nsize));
  psrifs = (long int *) malloc(sizeof(long int)*(naddt));

  #ifdef _OPENMP
  if (fftw_init_threads()) fftw_plan_with_nthreads(omp_get_max_threads());
  #endif

  chandata = (double **) malloc(omp_get_max_threads() * sizeof(double *));
  in = (fftw_complex**) malloc (omp_get_max_threads() * sizeof(fftw_complex*));
  out = (fftw_complex**) malloc (omp_get_max_threads() * sizeof(fftw_complex*));
  fplan = (fftw_plan*) malloc (omp_get_max_threads() * sizeof(fftw_plan));
  bplan = (fftw_plan*) malloc (omp_get_max_threads() * sizeof(fftw_plan));
  ai = (double **) malloc(omp_get_max_threads() * sizeof(double *));
  for (i=0; i<omp_get_max_threads(); i++) {
    chandata[i] = (double *) malloc(nsize * sizeof(double));
    ai[i] = (double *) malloc(nsize * sizeof(double));
    in[i] = fftw_malloc ((naddt+5) * sizeof(fftw_complex));
    out[i] = fftw_malloc ((naddt+5) * sizeof(fftw_complex));
    fplan[i] = fftw_plan_dft_1d(naddt, in[i], out[i], FFTW_FORWARD, FFTW_ESTIMATE);
    bplan[i] = fftw_plan_dft_1d(naddt, in[i], out[i], FFTW_BACKWARD, FFTW_ESTIMATE);
  }

  printf (" Preparing the median-filtering plan ... ");
  // plan the median-filter window sizes for spectral whitening
  //itemp = 16; // starting length of median-filter
  itemp = 8; // starting length of median-filter
  iter = 0;
  redplan[iter] = itemp;
  isum = redplan[iter]+itemp;
  while(isum<naddt/2){    // will compute mean/rms only for half of ft-arrays
   iter = iter + 1;
   redplan[iter] = itemp;
   itemp = itemp*2;
   isum = isum + itemp;
  }
  isum = isum-itemp;
  iter++;
  redplan[iter] = naddt/2-isum-1;
  nred = iter+1;
  // wmean = (double *) malloc((nred)*sizeof(double));
  // wrms = (double *) malloc((nred)*sizeof(double));
  wmean = (double **) malloc(omp_get_max_threads() * sizeof(double *));
  wrms = (double **) malloc(omp_get_max_threads() * sizeof(double *));
  for (i=0; i<omp_get_max_threads(); i++) {
    wmean[i] = (double *) malloc(nred * sizeof(double));
    wrms[i] = (double *) malloc(nred * sizeof(double));
  }
  printf (" Done! \n");

  // identify the fft-bin corresponding to pulsar's fundamental to be protected
  printf (" Planning which fundamental and harmonics to be safeguarded ... ");
  atemp = 1.0/(tsamp*(float)naddt);  // delta-f
  psrif = floor((psrf+atemp*0.5)/atemp) ; // index in the fft-array
  if(psrif<(psrfbins+2)) {
     psrif=-100;
     printf("\n ERROR: Pulsar's fundamental too low, not possible to safeguard !");
     exit(1);
  }
  if(psrif>naddt/2) {
     psrif=-100;
  }
  if(psrif>0) psrnf = (naddt/2)/psrif;
  //printf("\n *** %ld ",psrfbins);
  for (i=0;i<naddt;i++) psrifs[i]=1;
  if (psrif>0){  // protect specified pulsar spin frequency
    for (i=0;i<psrnf;i++) {
      j = (i+1)*psrif;
      //printf("-- %ld %ld --",j,naddt-j);
      for (k=(j-psrfbins);k<(j+psrfbins+1);k++) {
        psrifs[k] = -1;
        psrifs[naddt-k] = -1;
      }
    }
  }
  printf (" Done! \n");



  nsblk=nchans*nifs*naddt;
  fblock=(float *) malloc(nsblk*sizeof(float));
  sblock=(unsigned short *) malloc(nsblk*sizeof(unsigned short));
  cblock=(unsigned char *) malloc(nsblk*sizeof(unsigned short));
  realtime=min=0.0;
  max=(float) pow(2.0,(double)obits) -1.0;
  hpower = (max+min)/2.0;

  ts0dm = (float *) malloc((totsamp+5)*sizeof(float));

  last_mspec = (double *) malloc((nchans)*sizeof(double));
  for (iter=0;iter<nchans;iter++) last_mspec[iter]=0.0;
  fftstat=(float *) malloc(nsblk*sizeof(float));
  chanstat=(float *) malloc(nchans*sizeof(float));
  for(iter=0;iter<nsblk;iter++) fftstat[iter]=0.0;
  for(iter=0;iter<nchans;iter++) chanstat[iter]=0.0;
  nints = ceil((double)totsamp/(double)naddt) ;
  tfvar = (double *) malloc((totsamp+5)*sizeof(float));
  tfmean = (double *) malloc((totsamp+5)*sizeof(float));
  predist = (float *) malloc(maxhist*nsub*sizeof(float));
  postdist = (float *) malloc(maxhist*nsub*sizeof(float));
  finaldist = (float *) malloc(maxhist*nsub*sizeof(float));
  xpredist = (float *) malloc(maxhist*sizeof(float));
  xpostdist = (float *) malloc(maxhist*sizeof(float));
  xfinaldist = (float *) malloc(maxhist*sizeof(float));
  for(iter=0;iter<maxhist*nsub;iter++){
     predist[iter]=0.0;
     postdist[iter]=0.0;
     finaldist[iter]=0.0;
  }
  for(iter=0;iter<maxhist;iter++){
     xpredist[iter]=0.0;
     xpostdist[iter]=0.0;
     xfinaldist[iter]=0.0;
  }
  printf (" All buffers prepared! \n\n");

  /* main loop */
  printf (" Now rfiClean-ing (and 0-DM cleaning & downsampling, if asked for), and writing out the data...\n \n");
  istart = 0;
  iblock = 0;
  while ((ns=read_block(input,nbits,fblock,nsblk,byte_offset))>0 && iblock<nblocks) {
    byte_offset = byte_offset + (long) (naddt*nchans*(nbits/8.0));
    iblock = iblock + 1;
    n0 = ns/nchans;
    if(n0<naddt){
      jj = naddt*nchans;
      for (j=ns;j<jj;j++) fblock[j] = 0.0;
    }

    if (i_function == 0) t_function_sum = t_function = clock();
    else t_function = clock();
    cleanit(fblock,nchans,naddt);
    t_function_sum += ((double) clock() - t_function);
    i_function += 1;
    //-------------------------------------------------
    // compute post-cleaning 0-DM tseries
    for (j=0;j<n0;j++) {
       atemp = 0.0;
       for (i=0; i<nchans; i++){
          k = j*nchans + i;
          atemp = atemp + fblock[k];
       }
       ts0dm[j+istart] = atemp/(float)nchans;
    }
    //-------------------------------------------------
    // subtract 0-DM tseries from data
    if(zerodm>0) {
      for (j=0;j<n0;j++) {
         for (i=0; i<nchans; i++){
            k = j*nchans + i;
            fblock[k] = fblock[k] - ts0dm[j+istart] + hpower;
         }
      }
    }
    //printf("istart, totsamp: %ld,  %ld\n", istart,totsamp);
    //-------------------------------------------------

    istart = istart+n0;


    if (!opened) {
      opened=1;
    }
    nout=ns;
    if(nsamp>0) {
      //printf("downsampling by %d \n",nsamp);
      nout = ns/nsamp;
      kk = n0/nsamp;
      for (j=0;j<kk;j++) {
         jt1 = j*nsamp;
         jt2 = j*nchans;
         for (i=0; i<nchans; i++){
            atemp = 0.0;
            for (jj=0;jj<nsamp;jj++){
              k = (jt1+jj)*nchans + i;
              atemp = atemp + fblock[k];
            }
            k = jt2 + i;
            fblock[k] = atemp/(float)nsamp;
         }
      }
    }
    switch (obits) {
    case 32:
      fwrite(fblock,sizeof(float),nout,output);
      break;
    case 16:
      float2short(fblock,nout,min,max,sblock);
      fwrite(sblock,sizeof(unsigned short),nout,output);
      break;
    case 8:
      float2char(fblock,nout,min,max,cblock);
      fwrite(cblock,sizeof(unsigned char),nout,output);
      break;
    case 4:
      if (nout==1) {
	/* must have at least two samples for four-bit packing save this one */
	fsaved[nsaved]=fblock[0];
	nsaved++;
	if (nsaved==2) {
	  /* we have 2 saved! write out */
	  float2four(fsaved,nsaved,min,max,cblock);
	  fwrite(cblock,sizeof(unsigned char),1,output);
	  nsaved=0;
	}
      } else {
	/* normal case */
	float2four(fblock,nout,min,max,cblock);
	fwrite(cblock,sizeof(unsigned char),nout/2,output);
      }
      break;
    }
    //realtime+=(float) tsamp * (float) ns/(float) nchans/(float) nifs;
    //sprintf(string,"time:%.1fs",realtime);
  }

  printf (" Data rfiClean-ed and written out! \n\n");
  printf (" Now making diagnostic plots ... ");

  plot_data(plotdevice,wpout);
  printf (" Done! \n ");

  fclose(input);
  fclose(output);

  // fftw_destroy_plan ( fplan );
  // fftw_destroy_plan ( bplan );
  // fftw_free ( in );
  // fftw_free ( out );
  // free (ai);
  // free (wmean);
  // free (wrms);
  free (tfvar);
  free (tfmean);
  free (fblock);
  free (sblock);
  free (cblock);
  free (fftstat);
  free (chanstat);
  free (predist);
  free (postdist);
  free (finaldist);
  free (xpredist);
  free (xpostdist);
  free (xfinaldist);
  free (last_mspec);
  free (ts0dm);
  free (psrifs);
  free(mspec);
  free(rspec);
  free(vspec);
  free(wspec);
  free(wt);
  free(coff);

  for (i=0; i<omp_get_max_threads(); i++) {
    fftw_free (in[i]);
    fftw_free (out[i]);
    fftw_destroy_plan ( fplan[i] );
    fftw_destroy_plan ( bplan[i] );
    free (ai[i]);
    free (wmean[i]);
    free (wrms[i]);
    free(chandata[i]);
  }
  free (fplan);
  free (bplan);
  free (in);
  free (out);
  free (ai);
  free (wmean);
  free (wrms);
  free (chandata);

  FILE *file_timing = fopen("timing.csv", "a+");
  fprintf(file_timing, "%d,%d,%d,%.3f,%.3f\n",
    naddt,
    omp_get_max_threads(),
    time(NULL) - s,
    ((double) clock() - t_global) / (CLOCKS_PER_SEC) / (double) omp_get_max_threads(),
    ((double) t_function_sum / (i_function * CLOCKS_PER_SEC)) / (double) omp_get_max_threads()
  );
  // fprintf(file_timing, "%s,%d,%d\n", time_section, -1, time(NULL) - s);
  fclose(file_timing);

  printf (" Freed the buffers. All done!\n\n ");
  printf ("%d seconds\n\n", time(NULL) - s);
}
