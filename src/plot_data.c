#include "cpgplot.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "rficlean.h"

/* returns the minimum value of a vector */
float vmin(float *vec, int n) /*includefile*/
{
  float minvalue; int i;
  minvalue=vec[0];
  for (i=1;i<n;i++) if (vec[i]<minvalue) minvalue=vec[i];
  return (minvalue);
}

/* returns the maximum value of a vector */
float vmax(float *vec, int n) 
{
  float maxvalue; int i;
  maxvalue=vec[0];
  for (i=1;i<n;i++) if (vec[i]>maxvalue) maxvalue=vec[i];
  return (maxvalue);
}



void plot_data(char outdev[])
{
  char plotdev[100];
  int i,j,k,npt, i1,i2,j1,j2;
  long int ntot;
  float *xarr, *yarr, *fcstat;
  double *tmparr;
  float dmin, dmax, tr[6], locut, hicut, fdel;
  float left, right, top, bottom;
  float xl, xh, yl, yh, x1,x2,y1,y2;
  float tt, ft, th, fh;   /* thin and fat thicknesses and heights */
  float lm, rm, tm, bm;   /* LRTB margins */
  int ii, mincol, maxcol, numcol;

  npt = maxhist*nsub;
  ntot = nchans*naddt;
  xarr=(float *) malloc(ntot*sizeof(float));
  yarr=(float *) malloc(ntot*sizeof(float));
  tmparr=(double *) malloc((ntot+5)*sizeof(double));
  fcstat = (float *) malloc(nchans*naddt*sizeof(float));
  //sprintf(plotdev, "%s/CPS", plotdev);
  strcpy(plotdev, outdev);


  /* Open and prep the device */
  cpgopen(plotdev);
  cpgpap(10.25, 8.5 / 11.0);
  cpgpage();
  //cpgiden();
  //cpgscf(2);
  cpgsch(0.7);
  cpgqcir(&mincol, &maxcol);
  numcol = maxcol - mincol + 1;
  for (ii = mincol; ii <= maxcol; ii++) {
      float color;
      color = (float) (maxcol - ii) / (float) numcol;
      cpgscr(ii, color, color, color);
  }

   /* Set thicknesses and margins */

   lm = 0.04;
   rm = 0.04;
   bm = 0.08;
   tm = 0.05;
   ft = 3.0;      /* This sets fat thickness = 3 x thin thickness */
   tt = 0.92 / (6.0 + 4.0 * ft);
   ft *= tt;
   fh = 0.55;
   th = tt * 11.0 / 8.5;

   /*=======================================================*/
   igetstat = 0;
   if(igetstat>0)
   { /* Power Histograms */
      float *theo, *hist, *hpows, *tpows,llimit;
      int numtheo = 200, ibin, numpows, ist, iend, itemp, ihist;
      double dtheo, dhist, spacing;

      for(ihist=0; ihist<3; ihist++){
        if(ihist==0){
          for(i=0; i<npt ; i++) yarr[i] = log10(predist[i]+1) ;
        } else if(ihist==1) {
          for(i=0; i<npt ; i++) yarr[i] = log10(postdist[i]+1) ;
        } else if(ihist==2) {
          for(i=0; i<npt ; i++) yarr[i] = log10(finaldist[i]+1) ;
        } 
        dmin = vmin(yarr,npt);
        dmax = vmax(yarr,npt);
        llimit =  (pow(10.0,dmin)+(pow(10.0,dmax)-pow(10.0,dmin))*0.1);
        if(llimit> (pow(10.0,dmin+1.0))) llimit=pow(10.0,dmin+1.0);
        //printf("dmin,dmax,llimit :  %f %f %f \n",dmin,dmax,llimit);
        ist=maxhist;
        for (j=0; j<nsub; j++){
          itemp=0;
          k = j*maxhist;
          for (i=0; i<maxhist/2; i++){
            itemp = i;
            if(ihist==0){ if(predist[i+k]>llimit) break;}
            else if(ihist==1){ if(postdist[i+k]>llimit) break;}
            else if(ihist==2){ if(finaldist[i+k]>llimit) break;}
          }
          if(ist>itemp) ist=itemp;
        }
        iend=0;
        for (j=0; j<nsub; j++){
          itemp=maxhist;
          k = j*maxhist;
          for (i=0; i<maxhist/2; i++){
            ibin = maxhist-i-1;
            itemp = ibin;
            if(ihist==0) { if(predist[ibin+k]>llimit) break;}
            else if(ihist==1) { if(postdist[ibin+k]>llimit) break;}
            else if(ihist==2) { if(finaldist[ibin+k]>llimit) break;}
          }
          if(iend<itemp) iend=itemp;
        }
        //printf("ist,iend,maxhist: %d %d   %d\n",ist,iend,maxhist);
        if(ihist==0){
          for(i=ist; i<iend ; i++) xarr[i-ist] = xpredist[i] ;
          left = lm;
          right = lm + ft + tt;
        } else if(ihist==1) {
          for(i=ist; i<iend ; i++) xarr[i-ist] = xpostdist[i] ;
          left = lm + ft + 2.0*tt;
          right = lm + 2.0*ft + 3.0*tt;
        } else if(ihist==2) {
          for(i=ist; i<iend ; i++) xarr[i-ist] = xfinaldist[i] ;
          left = lm + 2.0*ft + 4.0*tt;
          right = lm + 3.0*ft + 5.0*tt;
        } 
        bottom = 0.80;
        top = 0.96;
        cpgsvp(left, right, bottom, top);
        xl = xarr[0] - (xarr[iend-ist-1]-xarr[0])*0.1;
        xh = xarr[iend-ist-1] + (xarr[iend-ist-1]-xarr[0])*0.1;
        yl = fmax(log10(llimit*0.9),1.0);
        yh = dmax*1.1;
        cpgswin(xl, xh, yl, yh);
        cpgmtxt("L", 1.1, 0.5, 0.5, "Number");
        cpgmtxt("B", 2.1, 0.5, 0.5, "Power");
        k = iend-ist-1;
        for (j=0; j<nsub; j++){
          if(ihist==0){
            for(i=ist;i<iend;i++) yarr[i-ist]=log10(predist[i+maxhist*j]+1);}
          else if(ihist==1){
            for(i=ist;i<iend;i++) yarr[i-ist]=log10(postdist[i+maxhist*j]+1);}
          else if(ihist==2){
            for(i=ist;i<iend;i++) yarr[i-ist]=log10(finaldist[i+maxhist*j]+1);}
          cpgbin(k, xarr, yarr, 1);
          //if(ihist==2){ for(i=0;i<(iend-ist-1);i++) printf("%f  %f \n",xarr[i],yarr[i]);}
        }
        cpgscr(maxcol, 0.5, 0.5, 0.5);
        cpgsci(maxcol);     /* Grey */
        if(nsub==1){
          if(ihist==0){
             xarr[0] = xpredist[(int)((maxhist/2.0)+(numhist/2.0/5.0)*sthresh)];}
          else if(ihist==1){
             xarr[0] = xpostdist[(int)((maxhist/2.0)+(numhist/2.0/5.0)*sthresh)];}
          else if(ihist==2){
             xarr[0] = xfinaldist[(int)((maxhist/2.0)+(numhist/2.0/5.0)*sthresh)];}
          xarr[1] = xarr[0];
          yarr[0] = yl;
          yarr[1] = yh;
          cpgsls(4);   
          cpgscr(maxcol, 1.0, 0.0, 0.0);
          cpgsci(maxcol); 
          cpgline(2, xarr, yarr);
        }
        cpgsls(1);     
        cpgsci(1);    
        cpgbox("BCNST", 0.0, 0, "BCLST", 0.0, 0);
      }
   }
   /*=======================================================*/
      
   /*==== FFT birdies zap stats ==========================*/
        left = lm;
        right = lm + 3.0 * ft + 4.0 * tt;
        bottom = bm;
        top = bm + fh;
        fdel = (1.0/(tsamp*naddt));
        xl = log10(fdel/2.0);
        xh = log10(1.0/(2.0*tsamp));
        yl = fch1;
        yh = fch1 + nchans*foff;
        cpgsvp(left, right, bottom, top);
        cpgswin(xl, xh, yl, yh);
        //printf("%f %f %f %f\n",xl,xh,yl,yh);
        //printf("%f %f %f %f\n",pow(xl,10.0),pow(xh,10.0),yl,yh);
        cpgbox("BCNLST", 0.0, 0, "BNST", 0.0, 0);
        cpgscr(maxcol, 1.0, 0.0, 0.0);  /* Red */
        for(i=0;i<ntot;i++) tmparr[i]=fftstat[i];
        robust_meanrms(tmparr,ntot);
        locut = fmax(0.0,tmparr[ntot]-3.0*tmparr[ntot+1]);
        hicut = vmax(fftstat,ntot/2);
        if(hicut <= locut) hicut=locut+1.0;
        //hicut = tmparr[ntot]+1.5*tmparr[ntot+1];
        tr[2] = tr[4] = 0.0;
        tr[1] = (xh - xl) / (naddt/2);
        tr[0] = xl - (tr[1] / 2);
        tr[5] = (yh - yl) / nchans;
        tr[3] = yl - (tr[5] / 2);
        ii = naddt;
        //cpgsitf(1);
        //cpgimag(fftstat, ii, nchans, 1, naddt/2, 1, nchans, locut, hicut, tr);
        // try plotting over log scale in fluc. freq.
        xarr[0]=0.0;
        yarr[0] = log10(fdel/2.0/1000.0);
        for(i=1;i<naddt/2;i++){
          xarr[i] = log10(fdel*i);
          yarr[i] = log10(fdel*i+fdel/2.0)-log10(fdel*i-fdel/2.0);
        }
        for(i=1;i<naddt/2;i++){
          for(j=0;j<nchans;j++) fcstat[j]=fftstat[j*naddt+i];
          xl = xarr[i];
          xh = xl;
          tr[2] = tr[4] = 0.0;
          tr[1] = yarr[i];
          tr[0] = xl - (yarr[i]);
          tr[5] = (yh - yl) / nchans;
          tr[3] = yl - (tr[5] / 2);
          ii = 1;
          cpgimag(fcstat, ii, nchans, 1, 1, 1, nchans, locut, hicut, tr);
        }
        xl = log10(fdel/2.0);
        xh = log10(1.0/(2.0*tsamp));
        yl = fch1;
        yh = fch1 + nchans*foff;
        cpgswin(xl, xh, yl, yh);
        cpgbox("BCNLST", 0.0, 0, "BNST", 0.0, 0);
        //cpgbox("", 0.0, 0, "CST", 0.0, 0);
        cpgmtxt("B", 2.6, 0.5, 0.5, "Fourier Frequency (Hz)");
        cpgmtxt("L", 2.1, 0.5, 0.5, "Radio Frequency (MHz)");
        /* Label */
        left = lm + 3.0 * ft + 4.0 * tt;
        right = lm + 3.0 * ft + 4.0 * tt + tt;
        bottom = bm + fh;
        top = bm + fh + th;
        cpgsvp(left, right, bottom, top);
        cpgswin(0.0, 1.0, 0.0, 1.0);
        cpgscr(maxcol, 1.0, 0.0, 0.0);
        cpgsci(maxcol);      
        cpgsci(1);    
        cpgsch(0.65);
        cpgptxt(0.65, 0.7, 0.0, 0.5, "Fractional");
        cpgptxt(0.65, 0.35, 0.0, 0.5, "Number");
        cpgsci(1);    
        cpgsch(0.7);

        /*  Total number zapped versus radio Frequency */
        left = lm + 3.0 * ft + 4.0 * tt;
        right = lm + 3.0 * ft + 4.0 * tt + tt;
        bottom = bm;
        top = bm + fh;
        cpgsvp(left, right, bottom, top);
        for (i=0; i<nchans; i++){
            j = i*naddt;
            tmparr[i] = 0.0;
            for (k=0;k<naddt;k++) tmparr[i] = tmparr[i]+fftstat[k+j];
        }
        for (i=1; i<nchans; i++) yarr[i-1]=tmparr[i];
        dmax = vmax(yarr,nchans-2);
        if(dmax<=0) dmax=1.0;
        for (i=0; i<nchans; i++) yarr[i]=tmparr[i]/dmax;
        for (i=0; i<nchans; i++) tmparr[i]=yarr[i];
        for (i=0; i<nchans; i++) xarr[i]=(float) i;
        dmin = vmin(yarr,nchans);
        dmax = vmax(yarr,nchans);
        robust_meanrms(tmparr,nchans);
        ////xh = fmin(tmparr[nchans] + 8.0*tmparr[nchans+1],dmax);
        ////xh = fmax(tmparr[nchans] + 8.0*tmparr[nchans+1],0.55);
        //xh = fmin(dmax,1.05);
        //xl = fmin((dmin - (xh-dmin)*0.05),-0.05);
        xh = 1.05;    // now fix this range
        xl = -0.05;
        yl = 0.0;
        yh = nchans-1;
        cpgswin(xl, xh, yl, yh);
        cpgbox("BCNT", 0.0, 0, "BCMST", 0.0, 0);
        cpgscr(maxcol, 1.0, 0.0, 0.0);
        cpgsci(maxcol);         
        cpgsls(1);             
        cpgsci(1);
        cpgbin(nchans, yarr, xarr, 1);
        cpgmtxt("R", 2.5, 0.5, 0.5, "Channel Number");

        /*  Total number zapped versus Fluc. Frequency */
        for (i=0;i<naddt; i++) tmparr[i]=0.0;
        for (i=0; i<ntot; i++){
            j = (int) floor((float)i/naddt);
            k = i - j*naddt;
            tmparr[k]=tmparr[k]+fftstat[i];
        }
        for (i=0;i<naddt; i++) yarr[i]=tmparr[i+1]/(long long)(nints*nchans);
        dmin = vmin(yarr,naddt);
        dmax = vmax(yarr,naddt);  
        left = lm;
        right = lm + 3.0 * ft + 4.0 * tt;
        bottom = bm + fh;
        top = bm + fh + th;
        cpgsvp(left, right, bottom, top);
        xl = log10(fdel/2.0);
        xh = log10(1.0/(2.0*tsamp));
        //yl = dmin - (dmax-dmin)*0.05;
        //yh = dmax + (dmax-dmin)*0.05;
        yl = -0.1;  // fix the min-max range
        yh = 1.1;
        cpgswin(xl, xh, yl, yh);
        cpgbox("BLST", 0.0, 0, "BCNST", 0.0, 0);
        cpgscr(maxcol, 1.0, 0.0, 0.0);
        cpgsci(maxcol);  
        xarr[0]=0.0;
        for(i=0;i<naddt/2;i++) xarr[i] = log10(fdel*(i+1));
        cpgsls(1);     
        cpgsci(1);  
        cpgbin(naddt/2, xarr, yarr, 1);
        //for (i=0;i<naddt/2;i++) printf("%f  %f\n",xarr[i],yarr[i]);
        cpgbox("CMLST", 0.0, 0, "", 0.0, 0);
        cpgmtxt("T", 1.8, 0.5, 0.5, "Fourier Frequency (Hz)");
     /*=======================================================*/

     /*==== Chan-reject stats ==========================*/
        left = lm + 3.0 * ft + 6.0 * tt + 0.25*tt;
        right = lm + 3.0 * ft + 6.0 * tt + 1.25*tt;
        bottom = bm;
        top = bm + fh;
        xl = -0.1;
        xh = 1.1;
        yl = 0.0;
        yh = nchans-1;
        cpgsvp(left, right, bottom, top);
        cpgswin(xl, xh, yl, yh);
        cpgbox("BCNST", 0.0, 0, "BNST", 0.0, 0); 
        for(i=0; i<nchans; i++) xarr[i]=i;
        for(i=0; i<nchans-2; i++) yarr[i]=chanstat[i+1];
        dmax = vmax(yarr,nchans-2);
        if(dmax<=0) dmax=1.0;
        for(i=0; i<nchans; i++) yarr[i]=chanstat[i]/dmax;
        //for(i=0; i<nchans; i++) printf("%f %f %f \n",xarr[i],yarr[i],chanstat[i]);
        cpgbin(nchans, yarr, xarr, 1);
        yl = fch1;
        yh = fch1 + nchans*foff;
        cpgswin(xl, xh, yl, yh);
        cpgbox("", 0.0, 0, "CMST", 0.0, 0); 
        cpgmtxt("R", 2.1, 0.5, 0.5, "Radio Frequency (MHz)");
     /*=======================================================*/




  cpgclos();
  free (xarr);
  free (yarr);
  free (fcstat);
  free (tmparr);

}

