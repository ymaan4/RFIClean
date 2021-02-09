/*
 * cleanit  --- the actual cleaning/processing of data.
 * 
 * 
 * Yogesh Maan <maan@astron.nl>   2018.
 * 
 * 
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rficlean.h"
#include <fftw3.h>

//--------------------------------------------------------------
int all_samef(double *arr, long int np)
{
   int i;
   double atemp;
   atemp = arr[0];
   for (i=1; i<np; i++){
     if(atemp != arr[i]) return 0;
   }
   return 1;
}
int all_same(float *arr, int np)
{
   int i;
   float atemp;
   atemp = arr[0];
   for (i=1; i<np; i++){
     if(atemp != arr[i]) return 0;
   }
   return 1;
}

//============== Conventional mean rms ===============================
/* TO COMPUTE MEAN and RMS OF 'arr' */

void simple_meanrms(double *arr, long int np)
{
        long int i;
        double amean,rms,diff,avar;

        amean = 0.0;
        for (i=0; i<np; i++){
          amean=arr[i]+amean;
        }
        amean=amean/np;

        rms=0.0;
        for (i=0;i<np;i++){
          diff=arr[i]-amean;
          rms=rms+diff*diff;
        }
        avar=rms/(np-1.0);
        rms=sqrt(rms/(np-1.0));

        arr[np] = amean;
        arr[np+1] = rms;
        arr[np+2] = avar;
        return;
}
//============================================================================


//============== robust mean rms ===============================
/* TO COMPUTE MEAN and RMS OF 'arr' BY EXCLUDING WHAT MAY BE
SOME CONTRIBUTION FROM INTERFERENCE */

void robust_meanrms(double *arr, long int np)
{
        long int i,iter,maxiter;
        double an,amean,amean0,rms,rms0,diff;

        iter=0;
        maxiter=4;
        rms=0.0;
        arr[np] = 0.0;
        arr[np+1] = 1.0;

l101:   iter=iter+1;
        amean=0.0;
        an=0.0;

        for (i=0; i<np; i++){
          if(iter==1) goto l1;
          diff=fabs(arr[i]-amean0);
          if(diff<=(4.0*rms)) goto l1;
          goto l2;
l1:       amean=arr[i]+amean;
          an=an+1.0;
l2:       continue;
        }

        rms0=rms;
        if(an>0)amean=amean/an;

        rms=0.0;
        an=0.0;
        for (i=0;i<np;i++){
          diff=fabs(arr[i]-amean);
          if(iter==1) goto l11;
          if(diff <= (4.0*rms0)) goto l11;
          goto l12;
l11:      rms=rms+diff*diff;
          an=an+1.0;
l12:       continue;
        }

        if(an > 0.0)rms=sqrt(rms/an);
        amean0=amean;

        if(iter==1) goto l101;
        if(rms == 0.0) return;
        if(fabs((rms0/rms)-1.0) > 0.05 && iter < maxiter) goto l101;

        arr[np] = amean;
        arr[np+1] = rms;
        if(rms<=0.0) arr[np+1] = 100000.0;

        return;
}
//============================================================================

//============== FT cleaning =========================================
void fftclean(fftw_complex *in, fftw_complex *out, long int npts, int ioff)
{
  int i,j,k,nfft,ih,nfft2,ii,i1,i2,kk;
  short int ft_sign;
  double anfft, asign ;
  double ath, bth;

  nfft = npts ;
  nfft2 = nfft/2;
  fftw_execute ( fplan );
  anfft = (double)nfft;
  for (i=0;i<nfft;i++) {
      in[i][0] = out[i][0];
      in[i][1] = out[i][1];
  }

  if(iwhite==1){
    k = 1;  // exclude dc
    for (ii=0;ii<nred;ii++){
      kk = redplan[ii];
      for (j=0;j<kk;j++){
        ai[j] = out[k][1];  // imaginary part
        k++;
      }
      robust_meanrms(ai,kk);
      wmean[ii] = ai[kk];
      wrms[ii] = ai[kk+1];
    }

  }    // end of block for whitened estimates of mean,rms
  else {
    k = nfft-1;
    for (i=0;i<k;i++) ai[i]=out[i+1][1];  // exclude DC
    ai[k] = ai[k+1] = 0.0;
    robust_meanrms(ai,k);
    for (ii=0;ii<nred;ii++){
      wmean[ii] = ai[k];
      wrms[ii] = ai[k+1];
    }
  }  // end of no-whitening block

  // now the actual cleaning
  i1 = 1;  // exclude dc
  for (ii=0;ii<nred;ii++){
    ath = fthresh*wrms[ii] + wmean[ii];
    bth = ath;
    i2 = i1 + redplan[ii];
    for (i=i1;i<i2;i++) {
      if( (fabs(out[i][0])>ath || fabs(out[i][1])>bth) && psrifs[i]>0 ){
        fftstat[ioff+i] = fftstat[ioff+i] + 1.0;
        for (j=i,ih=1; ih<=nharm && j<nfft2; ih++,j=i*ih){
          in[j][0] = wmean[ii];
          in[j][1] = wmean[ii];
          in[nfft-j][0] = wmean[ii];
          in[nfft-j][1] = wmean[ii];
        }
     }
    }
    i1 = i2;
  }

/*  if (psrif>0){  // protect specified pulsar spin frequency
    for (i=0;i<psrnf;i++) {
      j = (i+1)*psrif;
      in[j][0] = out[j][0];
      in[j][1] = out[j][1];
      in[nfft-j][0] = out[nfft-j][0];
      in[nfft-j][1] = out[nfft-j][1];
      in[j+1][0] = out[j+1][0];
      in[j+1][1] = out[j+1][1];
      in[nfft-j-1][0] = out[nfft-j-1][0];
      in[nfft-j-1][1] = out[nfft-j-1][1];
      in[j-1][0] = out[j-1][0];
      in[j-1][1] = out[j-1][1]; 
      in[nfft-j+1][0] = out[nfft-j+1][0];
      in[nfft-j+1][1] = out[nfft-j+1][1]; 
    }
  }  */

  fftw_execute ( bplan );
  for (i=0;i<nfft;i++) {
    in[i][0] = out[i][0]/anfft;
  }

}
//============================================================================


//============== timeseries clipping =========================================
void tsclip(double *tdata, long int npts, float thresh)
{
  long int i,j;
  double anpt,ath;
  double m1,r1;

  for (i=0;i<npts;i++) ai[i]=tdata[i];
  robust_meanrms(ai,npts);
  m1 = ai[npts];
  r1 = ai[npts+1];
  if(r1<=0.0){
    r1 = 0.0;
    m1 = 0.0;
  }
  ath = m1+thresh*r1;
  for (i=0;i<npts;i++) {
    if( fabs(ai[i])>ath) ai[i] = (ai[i]/fabs(ai[i]))*ath;
  }
  for (i=0;i<npts;i++) tdata[i]=ai[i];
  tdata[npts] = m1;
  tdata[npts+1] = r1;
}
//============================================================================


//============== timeseries cleaning =========================================
void tsclean(double *tdata, long int npts, float thresh)
{
  long int i,j;
  double anpt,ath;
  double m1,r1;

  for (i=0;i<npts;i++) ai[i]=tdata[i];
  robust_meanrms(ai,npts);
  m1 = ai[npts];
  r1 = ai[npts+1];
  if(r1<=0.0){
    r1 = 0.0;
    m1 = 0.0;
  }
  ath = thresh*r1;
  for (i=0;i<npts;i++) {
    if( fabs((ai[i]-m1))>ath) ai[i] = m1;
  }
  for (i=0;i<npts;i++) tdata[i]=ai[i];
  tdata[npts] = m1;
  tdata[npts+1] = r1;
}
//============================================================================

//============== spectrum cleaning =========================================
void spfind(double *cdata, long int npts, float thresh, double *wt)
{
  long int i,j,k,n;
  double anpt,ath;
  double m1,r1,m2;

  n = npts-1;
  for (i=0;i<n;i++) ai[i]=fabs(cdata[i+1]-cdata[i]) ;
  ai[n] = 0.0;
  ai[n+1] = 0.0;
  robust_meanrms(ai,n);
  m1 = ai[n];
  r1 = ai[n+1];
  if(r1<=0.0) r1 = 10000.0;

  ath = thresh*r1;
  for (i=0;i<n;i++) {
    if( fabs(ai[i]-m1)>ath){
      wt[i] = -1.0;
      wt[i+1] = -1.0;
    }
  }

}
//============================================================================

//============== timeseries rfi finding =========================================
void tsfind(double *ttdata, long int npts, float thresh, double *wt)
{
  long int i,j,n;
  double anpt,ath;
  double m1,r1,m2;


  for (i=0;i<npts;i++) ai[i]=ttdata[i] ;
  ai[npts] = 0.0;
  ai[npts+1] = 0.0;
  robust_meanrms(ai,npts);
  m1 = ai[npts];
  r1 = ai[npts+1];
  if(r1<=0.0) r1 = 10000.0;

  ath = thresh*r1;
  for (i=0;i<npts;i++) {
    if( fabs((ai[i]-m1))>ath){
      wt[i] = -1.0;
    }
  }
  wt[npts] = m1;

}
//============================================================================


//============== spectral clipping =========================================
void spclip(double *ttdata, long int npts, float thresh)
{
  long int i,j,n;
  double ath;
  double m1,r1,m2;

  ttdata[npts] = 0.0;
  ttdata[npts+1] = 0.0;
  robust_meanrms(ttdata,npts);
  m1 = ttdata[npts];
  r1 = ttdata[npts+1];
  if(r1<=0.0) r1 = 10000.0;

  ath = m1+thresh*r1;
  for (i=0;i<npts;i++) {
    if( fabs(ttdata[i])>ath){
      ttdata[i] = (ttdata[i]/fabs(ttdata[i]))*ath;
    }
  }

}
//============================================================================



//============== main function that calls all cleaning functions ===================
void cleanit(float *data, int nchans, long int nadd)

{
  int inc,isame;
  int ii,jj,kk,channum;
  long int nxc, c, t, i;
  float thresh;
  double lg,lgr, an, atemp;


// get some pre-cleaning statistics

  for (channum=0; channum<nchans; channum++) coff[channum]=nadd*channum ;
  for (channum=0; channum<nchans; channum++) {
    inc = coff[channum];
  /* Select the correct channel */
    for (ii = 0, jj = channum; ii < nadd; ii++, jj += nchans) chandata[ii] = (double) data[jj];
    for (ii = 0; ii<nadd; ii++)
    { in[ii][0] = chandata[ii];
      in[ii][1] = 0.0; }
    if( pcl == 0 ){
      isame = 0;
    } else {
      isame = all_samef(chandata,nadd);
    }
    if( isame == 0){
      if(rfiFDx){
        fftclean(in,out,nadd,inc); // clean some periodic RFIs
      }
      for (ii = 0; ii<nadd; ii++) chandata[ii] = in[ii][0];
      if(rfiTx){
        tsclip(chandata,nadd,sthresh); // clean some spiky RFIs
        mspec[channum]=chandata[nadd];
        rspec[channum]=chandata[nadd+1];
      } else {
        robust_meanrms(chandata,nadd);
        mspec[channum]=chandata[nadd];
        rspec[channum]=chandata[nadd+1];
      }
    }
    else {
      mspec[channum]=chandata[0];
      rspec[channum]=0.0;
    }
    vspec[channum] = rspec[channum]*rspec[channum];
    for (ii = 0, jj = channum; ii < nadd; ii++, jj += nchans) data[jj]= (float) chandata[ii];
  }
// get some post- fft and timeseries cleaning statistics


  for (i=0; i<nchans; i++) wspec[i]=+1.0;
  if(rfiMSx){
    spfind(mspec,nchans,rthresh,wspec);
    tsfind(mspec,nchans,rthresh,wspec);
  }
  if(rfiVSx){
    spfind(vspec,nchans,rthresh,wspec);
    tsfind(vspec,nchans,rthresh,wspec);
  }
  if(pcl>0){
     kk = (int)(0.9*nchans);
     for (i=kk;i<nchans;i++) wspec[i]=0.0;
     for (i=0;i<kk;i++) wspec[i]=0.0;
  }
  for (i=0;i<nchans; i++){
    if(wspec[i] > 0.0) {
      lg = mspec[i];
      lgr = rspec[i];
      break;
    }
  }
  for (i=0;i<nchans;i++){
    if(wspec[i]<0.0) {
      mspec[i] = lg;
      rspec[i] = lgr;
    }
    if(wspec[i]>0.0) {
      lg = mspec[i];
      lgr = rspec[i];
    }
  }

  // log the information which channels got rejected
  for (i=0;i<nchans;i++) if(wspec[i]<0.0) chanstat[i]=chanstat[i]+1.0;
  

  for (channum=0; channum<nchans; channum++) {
    if(wspec[channum] < 0.0) {
       for (ii = 0, jj = channum; ii < nadd; ii++, jj += nchans) data[jj] = mspec[channum];
    }
  }
// see if the whole block needs to be rejected
  an = 0.0;
  for (ii=0; ii<nchans; ii++){ if(wspec[ii]<0.0) an=an+1.0; }
  if((an/nchans)>chanfrac){
    for (channum=0; channum<nchans; channum++) {
         for (ii = 0, jj = channum; ii < nadd; ii++, jj += nchans) data[jj] = last_mspec[channum]; }
  } else {
    for (ii=0; ii<nchans; ii++){ if(wspec[channum] > 0.0) last_mspec[ii] = mspec[ii];}
  }
  


// Try clipping some channels in individual samples
  if(rfiSclip){
    for (t=0; t<nadd; t++){
        nxc = t*nchans;
        for (c=0; c<nchans; c++) chandata[c] = data[nxc+c];
        spclip(chandata,nchans,clipthresh);
        for (c=0; c<nchans; c++) data[nxc+c] = chandata[c];
    }
  }
// Now some timeseries cleaning
  if(rfiTx){
    an = (double)nadd;
    for (t=0; t<nadd; t++){
        nxc = t*nchans;
        chandata[t]=0.0;
        for (c=0; c<nchans; c++) chandata[t] = chandata[t]+data[nxc+c];
        chandata[t]=chandata[t]/an;
    }
    for (i=0; i<nadd; i++) wt[i]=+1.0;
    tsfind(chandata,nadd,sthresh,wt);
    spfind(chandata,nadd,sthresh,wt);
    for (t=0; t<nadd; t++){
        nxc = t*nchans;
        if(wt[t] < 0.0){
          for (c=0; c<nchans; c++) data[nxc+c] = mspec[c];
        }
    }
  }

 // sometimes gpt results might need additional checks
  if(pcl>0){
    kk = 16;
    nvar = nchans/kk;
    for (t=0; t<nadd; t++){
      nxc = t*nchans;
      isame = 0;
      for(jj=0; jj<nvar; jj++){
        ii = nxc + jj*kk;
        for (c=0; c<kk; c++) chandata[c] = data[ii+c];
        if(all_samef(chandata,kk)>0) isame = isame + 1;
      }
      if(isame >= (int)(0.8*nvar)){
        for (c=0; c<nchans; c++) data[c+nxc] = 0.0;
      }
    }
  }

  // Try weighting temporal sections by variance
  if(tnorm>0){
      //------------------------------------------
    if(last_tvar<=0.0){
      for (t=0; t<nadd; t++){
          nxc = t*nchans;
          if(wt[t] < 0.0){
           rspec[t] = 0.0;
          }
          else{
           for (c=0; c<nchans; c++) chandata[c] = data[nxc+c];
           robust_meanrms(chandata,nchans);
           //simple_meanrms(chandata,nchans);
           rspec[t] = chandata[nchans+1]*chandata[nchans+1];
          }
      }
      t = nadd;
      robust_meanrms(rspec,t);
      //simple_meanrms(rspec,t);
      last_tvar = rspec[t];
    }
    nvar = nadd/nsect;
    for (jj=0; jj<nsect; jj++){
      for (c=0; c<nchans; c++) chandata[c] = 0.0;
      for (t=0; t<nvar; t++){
        nxc = (jj*nvar+t)*nchans;
        for (c=0; c<nchans; c++) chandata[c] = chandata[c]+data[nxc+c];
      }
      for (c=0; c<nchans; c++) chandata[c] = chandata[c]/nvar;
      //robust_meanrms(chandata,nchans);
      simple_meanrms(chandata,nchans);
      an = chandata[nchans+1]*chandata[nchans+1];
      if(an<=0.0) an=10000000.0;
      for (t=0; t<nvar; t++){
        nxc = (jj*nvar+t)*nchans;
        for (c=0; c<nchans; c++) data[nxc+c] = data[nxc+c]*last_tvar/an;
      }
    } 
  }

  // Try weighting channels by variance
  if(fnorm>0){
    if(last_fvar<=0.0){
      robust_meanrms(vspec,nchans);
      last_fvar = vspec[nchans];
    }
    for (channum=0; channum<nchans; channum++) {
      if(wspec[channum] > 0.0) {
         for (ii = 0, jj = channum; ii < nadd; ii++, jj += nchans) data[jj] = data[jj]*last_fvar/vspec[channum];
      }
      else {
         for (ii = 0, jj = channum; ii < nadd; ii++, jj += nchans) data[jj] = 0.0;
      }
    }   
  }


  /* ***Now the following part is in rficlean_data.c ***
//  if(iflip==1){  // flip the band
//    for (t=0; t<nadd; t++){
//      nxc = t*nchans;
//      for (c=0; c<nchans; c++) mspec[c] = data[nxc+c];
//      for (c=0; c<nchans; c++) data[nxc+c] = mspec[nchans-c-1];
//    }
//  }
  ***/

}
//============================================================================
