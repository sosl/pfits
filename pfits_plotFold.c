//  Copyright (C) 2015,2016 George Hobbs
// This file is part of the pfits software package
//

/* pfits is free software: you can redistribute it and/or modify 
 * it under the terms of the GNU General Public License as published by 
 * the Free Software Foundation, either version 3 of the License, or 
 * (at your option) any later version. 
 * pfits is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
 * GNU General Public License for more details. 
 * You should have received a copy of the GNU General Public License 
 * along with pfits.  If not, see <http://www.gnu.org/licenses/>. 
*/

//  gcc -lm -o pfits_plotFold pfits_plotFold.c pfits_setup.c pfits_loader.c -I/Users/hob044/hob044/software/cfitsio/include/ -L/Users/hob044/hob044/software/cfitsio/lib -lcfitsio -lcpgplot -lpgplot

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "pfits.h"
#include <cpgplot.h>

void drawColourMap(dSetStruct *dSet,int pol);

int main(int argc,char *argv[])
{
  dSetStruct *dSet;
  char fname[1024];
  int debug=0;
  int i;
  
  for (i=0;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	strcpy(fname,argv[++i]);
    }
  dSet = (dSetStruct *)malloc(sizeof(dSetStruct));
  initialise(&dSet,debug);
  setFilename(fname,dSet,debug);

  pfitsOpenFile(dSet,debug);
  printf("Loading header\n");
  pfitsLoadHeader(dSet,debug);


  drawColourMap(dSet,0);
  
  cpgend();
  

  //  pfitsCloseFile(dSet,debug);
  free(dSet);
}

void drawColourMap(dSetStruct *dSet,int pol)
{
  float tr[6];
  int nchan = dSet->head->nchan;
  int nbin = dSet->head->nbin;
  float heat_l[] = {0.0, 0.2, 0.4, 0.6, 1.0};
  float heat_r[] = {0.0, 0.5, 1.0, 1.0, 1.0};
  float heat_g[] = {0.0, 0.0, 0.5, 1.0, 1.0};
  float heat_b[] = {0.0, 0.0, 0.0, 0.3, 1.0};
  char title[128];
  float *plotArr;
  float *plotArr_scl;
  float *plotArr_meanScl;
  int i,j;
  int debug=0;
  int colnum;
  int subint=0;
  int status=0;
  short int nVal=0;
  short int *sVals;
  int initflag=0;
  float minV,maxV;
  float minV_scl,maxV_scl;
  float minV_meanScl,maxV_meanScl;
  float *dat_offs;
  float *dat_scl;
  float *mean_scl;
  float nfVal = 0;
  float *expectedX,*expectedY;
  
  float mx,my,mx2,my2;
  char key;
  int plot=1;
  int plotExpected=-1;
  int plotY=1;
  int zap;
  
  float plot_minX,plot_minY,plot_maxX,plot_maxY;
  
  printf("nchan = %d, nbin = %d\n",nchan,nbin);
  
  plotArr = (float *)malloc(sizeof(float)*nbin*nchan);
  plotArr_scl = (float *)malloc(sizeof(float)*nbin*nchan);
  plotArr_meanScl = (float *)malloc(sizeof(float)*nbin*nchan);
  expectedX = (float *)malloc(sizeof(float)*nchan);
  expectedY = (float *)malloc(sizeof(float)*nchan);
  dat_offs = (float *)malloc(sizeof(float)*nchan);
  dat_scl = (float *)malloc(sizeof(float)*nchan);
  mean_scl = (float *)malloc(sizeof(float)*nchan);
  
  sVals = (short int *)malloc(sizeof(short int)*nbin*nchan);
  fits_movnam_hdu(dSet->fp,BINARY_TBL,"SUBINT",1,&status);
  fits_get_colnum(dSet->fp,CASEINSEN,"DAT_OFFS",&colnum,&status);
  fits_read_col(dSet->fp,TFLOAT,colnum,subint+1,1,nchan,&nfVal,dat_offs,&initflag,&status);
  fits_get_colnum(dSet->fp,CASEINSEN,"DAT_SCL",&colnum,&status);
  fits_read_col(dSet->fp,TFLOAT,colnum,subint+1,1,nchan,&nfVal,dat_scl,&initflag,&status);
  
  fits_get_colnum(dSet->fp,CASEINSEN,"DATA",&colnum,&status);

  fits_read_col(dSet->fp,TSHORT,colnum,subint+1,1,nchan*nbin,&nVal,sVals,&initflag,&status);

  for (i=0;i<nchan;i++)
    {
      mean_scl[i] = 0;
      expectedX[i] = 0;
      expectedY[i] = i;
    }
  // Do some zapping
  for (i=0;i<nchan;i++)
    {
      zap = 0;
      if (i > 1793 && i <2050)	zap=1;
      if (i > 2050 && i < 2062) zap=1;
      if (i > 2170 && i < 2188) zap=1;
      if (i > 2296 && i < 2315) zap=1;
      if (i > 2336 && i < 2355) zap=1;
      if (i > 2601 && i < 2646) zap=1;
      if (i > 2740 && i < 2818) zap=1;
      if (i > 3066 && i < 3088) zap=1;
      if (i > 3200) zap=1;
      if (zap==1)	
	{
	  for (j=0;j<nbin;j++)
	    sVals[i*nbin+j] = 0;
	}
    }
  
  for (i=0;i<nchan;i++)
    {
      for (j=0;j<nbin;j++)
	{
	  plotArr[i*nbin+j] = (float)sVals[i*nbin+j];
	  plotArr_scl[i*nbin+j] = plotArr[i*nbin+j]*dat_scl[i]+dat_offs[i];
	  plotArr_meanScl[i*nbin+j] = plotArr_scl[i*nbin+j] - mean_scl[i];
	  if (i==0 && j==0)
	    {
	      minV = maxV = plotArr[i*nbin+j];
	      minV_scl = maxV_scl = plotArr_scl[i*nbin+j];
	      minV_meanScl = maxV_meanScl = plotArr_meanScl[i*nbin+j];
	    }
	  else
	    {
	      if (minV > plotArr[i*nbin+j]) minV = plotArr[i*nbin+j];
	      if (maxV < plotArr[i*nbin+j]) maxV = plotArr[i*nbin+j];
	      if (minV_scl > plotArr_scl[i*nbin+j]) minV_scl = plotArr_scl[i*nbin+j];
	      if (maxV_scl < plotArr_scl[i*nbin+j]) maxV_scl = plotArr_scl[i*nbin+j];
	      if (minV_meanScl > plotArr_meanScl[i*nbin+j]) minV_meanScl = plotArr_meanScl[i*nbin+j];
	      if (maxV_meanScl < plotArr_meanScl[i*nbin+j]) maxV_meanScl = plotArr_meanScl[i*nbin+j];
	    }
	}
    }

  printf("Min/max = %g %g (%g %g) (%g %g)\n",minV,maxV,minV_scl,maxV_scl,minV_meanScl,maxV_meanScl);
  tr[0] = 0;  tr[1] = 1;  tr[2] = 0;
  tr[3] = 0;  tr[4] = 0;  tr[5] = 1;

  cpgbeg(0,"/xs",1,1);
  cpgask(0);
  cpgctab(heat_l,heat_r,heat_g,heat_b,5,1.0,0.5);

  plot_minX=0;
  plot_minY=0;
  plot_maxX=dSet->head->nbin;
  plot_maxY=dSet->head->nchan;
  
  do {
    cpgenv(plot_minX,plot_maxX,plot_minY,plot_maxY,0,1);
    if (plotY==1)
      cpglab("Bin number","Channel number","");
    else if (plotY==2)
      cpglab("Bin number","Frequency (MHz)","");
    if (plot==1)
	cpgimag(plotArr,nbin,nchan,1,nbin,1,nchan,minV,maxV,tr);
    else if (plot==2)
	cpgimag(plotArr_scl,nbin,nchan,1,nbin,1,nchan,minV_scl,maxV_scl,tr);
    else if (plot==3)
	cpgimag(plotArr_meanScl,nbin,nchan,1,nbin,1,nchan,minV_meanScl,maxV_meanScl,tr);

    if (plotExpected==1)
      {
	cpgsci(3);
	cpgpt(nchan,expectedX,expectedY,20);
	cpgsci(1);
      }
    cpgcurs(&mx,&my,&key);
    if (key=='1') plot=1;
    else if (key=='2') plot=2;
    else if (key=='3') plot=3;
    else if (key=='z')
      {
	cpgband(2,0,mx,my,&mx2,&my2,&key);
	plot_minX = mx;
	plot_maxX = mx2;
	plot_minY = my;
	plot_maxY = my2;
      }
    else if (key=='#')
      {
	printf("Current min/max = %g %g\n",minV_meanScl,maxV_meanScl);
	printf("Enter new values ");
	scanf("%f %f",&minV_meanScl,&maxV_meanScl);
      }
    else if (key=='x') // Zap channel with highest value
      {
	int maxChan=0;
	float maxValThisChan;
	float maxValTot;
	float val;
	for (i=0;i<nchan;i++)
	  {
	    for (j=0;j<nbin;j++)
	      {
		if (plot==1) val=plotArr[i*nbin+j];
		else if (plot==2) val = plotArr_scl[i*nbin+j];
		else if (plot==3) val = plotArr_meanScl[i*nbin+j];
		if (j==0) maxValThisChan = val;
		else if (val > maxValThisChan) maxValThisChan = val;
	      }
	    if (i==0) {maxValTot = maxValThisChan; maxChan = i;}
	    else if (maxValThisChan > maxValTot) {maxValTot = maxValThisChan; maxChan = i;}
	  }
	printf("Maximium value is %g in channels\n",maxValTot);
	for (i=0;i<nchan;i++)
	  {
	    for (j=0;j<nbin;j++)
	      {
		if (plot==1) val=plotArr[i*nbin+j];
		else if (plot==2) val = plotArr_scl[i*nbin+j];
		else if (plot==3) val = plotArr_meanScl[i*nbin+j];
		if (j==0) maxValThisChan = val;
		else if (val > maxValThisChan) maxValThisChan = val;
	      }
	    if (maxValThisChan >= maxValTot-0.5)
	      {
		printf(" ... %d\n",i);
		for (j=0;j<nbin;j++)
		  plotArr[i*nbin+j] = plotArr_scl[i*nbin+j] = plotArr_meanScl[i*nbin+j] = 0.0;
	      }
	  }
	// Recalculate min/max
	for (i=0;i<nchan;i++)
	  {
	    for (j=0;j<nbin;j++)
	      {
		if (i==0 && j==0)
		  {
		    minV = maxV = plotArr[i*nbin+j];
		    minV_scl = maxV_scl = plotArr_scl[i*nbin+j];
		    minV_meanScl = maxV_meanScl = plotArr_meanScl[i*nbin+j];
		  }
		else
		  {
		    if (minV > plotArr[i*nbin+j]) minV = plotArr[i*nbin+j];
		    if (maxV < plotArr[i*nbin+j]) maxV = plotArr[i*nbin+j];
		    if (minV_scl > plotArr_scl[i*nbin+j]) minV_scl = plotArr_scl[i*nbin+j];
		    if (maxV_scl < plotArr_scl[i*nbin+j]) maxV_scl = plotArr_scl[i*nbin+j];
		    if (minV_meanScl > plotArr_meanScl[i*nbin+j]) minV_meanScl = plotArr_meanScl[i*nbin+j];
		    if (maxV_meanScl < plotArr_meanScl[i*nbin+j]) maxV_meanScl = plotArr_meanScl[i*nbin+j];
		  }
	      }
	  }

      }
    else if (key=='Z') // Set min-max from region
      {
	float min=0,max=0;
	float xval,yval;
	float tt;
	float val;
	int t=0;
	cpgband(2,0,mx,my,&mx2,&my2,&key);
	if (mx > mx2) {tt = mx; mx = mx2; mx2 = tt;}
	if (my > my2) {tt = my; my = my2; my2 = tt;}
	for (i=0;i<nbin;i++)
	  {
	    for (j=0;j<nchan;j++)
	      {
		xval = tr[0]+tr[1]*i+tr[2]*j;
		yval = tr[3]+tr[4]*i+tr[5]*j;
		if (plot==1) val=plotArr[j*nbin+i];
		else if (plot==2) val = plotArr_scl[j*nbin+i];
		else if (plot==3) val = plotArr_meanScl[j*nbin+i];
		//		printf("(%g %g) and (%g %g),(%g,%g)\n",xval,yval,mx,my,mx2,my2);
		if (xval > mx && xval <= mx2 && yval > my && yval < my2)
		  {
		    //		    printf("In here %g\n",val);
		    if (t==0)
		      {min = max = val; t=1;}
		    else
		      {
			if (min > val) min = val;
			if (max < val) max = val;
		      }
		  }
	      }
	  }


	if (plot==1)
	  {
	    minV = min;
	    maxV = max;
	  }
	else if (plot==2)
	  {
	    minV_scl = min;
	    maxV_scl = max;
	  }
	else if (plot==3)
	  {
	    minV_meanScl = min;
	    maxV_meanScl = max;
	  }
	printf("Have min/max = %g %g\n",min,max);
      }
    else if (key=='e')
      plotExpected*=-1;
    else if (key=='y')
      {
	plotY++;
	if (plotY==3) plotY=1;

	if (plotY==1)
	  {
	    tr[0] = 0;  tr[1] = 1;  tr[2] = 0;
	    tr[3] = 0;  tr[4] = 0;  tr[5] = 1;	      
	    plot_minY = 0; plot_maxY = nchan;
	    for (i=0;i<nchan;i++)
	      expectedY[i] = i;
	  }
	else if (plotY==2)
	  {
	    tr[0] = 0;  tr[1] = 1;  tr[2] = 0;
	    tr[3] = dSet->head->chanFreq[0];  tr[4] = 0;  tr[5] = (dSet->head->chanFreq[nchan-1]-dSet->head->chanFreq[0])/(double)nchan;
	    for (i=0;i<nchan;i++)
	      expectedY[i] = tr[3]+tr[5]*i;
	    printf("tr[5] = %g\n",tr[5]);
	    plot_minY = dSet->head->chanFreq[0]; plot_maxY = dSet->head->chanFreq[nchan-1];	    
	  }
      }
    else if (key=='D') // dedisperse
      {
	double dm=67.97;
	double freq;
	double fref = dSet->head->chanFreq[0];
	double period = 0.089403911483094;
	double tt;
	int offsetBin;
	float *temp;
	int newbin;
	float maxV;
	int maxJ;
	
	temp = (float *)malloc(sizeof(float)*nbin*nchan);
	printf("Got here 1\n");
	for (i=0;i<nchan;i++)
	  {
	    for (j=0;j<nbin;j++)
	      temp[i*nbin+j] = plotArr_meanScl[i*nbin+j];
	  }
	printf("Got here 2\n");

	for (i=0;i<nchan;i++)
	  {
	    freq = dSet->head->chanFreq[i];
	    offsetBin = -(4.13e-3*dm*(pow(fref/1000.0,-2)-pow(freq/1000.0,-2)))*nbin/period;

	    // Silly way to get peak in middle
	    {
	      printf("In here %d/%d\n",i,nchan);
	      maxV=temp[i*nbin];
	      for (j=0;j<nbin;j++)
		{
		  if (maxV < temp[i*nbin+j]) {maxV = temp[i*nbin+j]; maxJ = j;}
		}
	      offsetBin = maxJ-512;
	    }
	    for (j=0;j<nbin;j++)
	      {
		newbin = j+offsetBin;
		while (newbin < 0) newbin+=nbin;
		while (newbin >= nbin) newbin-=nbin;
		plotArr_meanScl[i*nbin+j] = temp[i*nbin+newbin];
	      }
	  }
	free(temp);
      }
    else if (key=='d') // Set DM
      {
	double dm;
	double freq;
	double fref = dSet->head->chanFreq[0];
	double period = 0.089403949593904;
	double offset = 0; 
	printf("SETTING PERIOD TO VELA!!!\n");
	printf("Enter DM ");
	scanf("%lf",&dm);
	printf("Enter offset ");
	scanf("%lf",&offset);
	for (i=0;i<nchan;i++)
	  {
	    freq = dSet->head->chanFreq[i];
	    expectedX[i] = -(4.13e-3*dm*(pow(fref/1000.0,-2)-pow(freq/1000.0,-2)))*nbin/period + offset;
	    while (expectedX[i] >= nbin) expectedX[i]-=nbin;
	    while (expectedX[i] < 0) expectedX[i]+=nbin;
	  }
      }
    else if (key=='u')
      {
	plot_minX=0;
	plot_maxX=dSet->head->nbin;
	if (plotY==1)
	  {
	    plot_minY=0;	
	    plot_maxY=dSet->head->nchan;	
	  }
	else if (plotY==2)
	  {
	    plot_minY = dSet->head->chanFreq[0]; plot_maxY = dSet->head->chanFreq[nchan-1];	    
	  }
      }
    else if (key=='b') // Set baseline
      {
	float b1,b2;
	printf("Enter baseline bin number X1 X2 ");
	scanf("%f %f",&b1,&b2);
	for (i=0;i<nchan;i++)
	  {
	    mean_scl[i] = 0;
	    for (j=b1;j<=b2;j++)
	      mean_scl[i]+=plotArr_scl[i*nbin+j];
	    mean_scl[i]/=(float)(b2-b1);
	  }
	for (i=0;i<nchan;i++)
	  {
	    for (j=0;j<nbin;j++)
	      {
		plotArr_meanScl[i*nbin+j] = plotArr_scl[i*nbin+j] - mean_scl[i];
		if (i==0 && j==0)
		  minV_meanScl = maxV_meanScl = plotArr_meanScl[i*nbin+j];
		else
		  {
		    if (minV_meanScl > plotArr_meanScl[i*nbin+j]) minV_meanScl = plotArr_meanScl[i*nbin+j];
		    if (maxV_meanScl < plotArr_meanScl[i*nbin+j]) maxV_meanScl = plotArr_meanScl[i*nbin+j];
		  }
	      }
	  }
	
      }
  } while (key!='q');

  free(plotArr);
  free(plotArr_scl);
  free(plotArr_meanScl);
  free(sVals);
  free(dat_offs);
  free(dat_scl);
  free(mean_scl);
  free(expectedX);
  free(expectedY);
}
