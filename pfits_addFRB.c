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

//gcc -lm -o pfits_addFRB pfits_addFRB.c pfits_setup.c pfits_loader.c T2toolkit.c -I/Users/hob044/hob044/software/cfitsio/include/ -L/Users/hob044/hob044/software/cfitsio/lib -lcfitsio 

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "pfits.h"
#include <cpgplot.h>
#include "fitsio.h"
#include "T2toolkit.h"


#define MAX_FRBS 1000
#define STRLEN 1024

typedef struct frbStruct {
  char fname[STRLEN];  
  float burstTime;
  float burstDM;
} frbStruct;

int convertTo2Bits(unsigned char *inArray, unsigned char *outArray, unsigned int n);
void readBurstHeader(FILE *finBurst,unsigned int *burstSamples);

int main(int argc,char *argv[])
{
  frbStruct *frb;
  dSetStruct *dSet;
  FILE *fin;
  int debug=0;
  int status=0;
  int i,ii,j,k,a;
  int setOut=0;
  int colnum_out;
  int colnum_datoffs;
  int colnum_datscl;
  char outFile[1024],oname[1024]; // Output filename
  fitsfile *outfptr;
  FILE *finBurst;
  char frbFile[1024];
  float burstTime,burstDM;
  long sub;
  int bpos;
  int nTimeSamples;
  long nSamples;
  int nFrequencySamples;
  float *loadData;
  char *out;
  char *outputVals;
  float fch1;
  float freq,dispersionDelay;
  double bt;
  long int pos;
  double t0,t1;
  unsigned int burstSamples;
  unsigned int changeVals[4][4];
  int ele;
  int totC;
  long randVal;
  float rnum;
  long int seed = TKsetSeed();
  int nFRB=0;
  int process=0;
  float tmax,tmin;
  long int ipos;
  
    // Initialise everything
  initialise(&dSet,debug);

  if (!(frb = (frbStruct *)malloc(sizeof(frbStruct)*MAX_FRBS)))
    {
      printf("Unable to allocate enough memory for %d FRBs\n",MAX_FRBS);
      exit(1);
    }

  
  // Read inputs
  for (i=0;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	setFilename(argv[++i],dSet,debug);
      else if (strcmp(argv[i],"-frb")==0)
	strcpy(frbFile,argv[++i]);
      else if (strcmp(argv[i],"-o")==0)
	{
	  strcpy(outFile,argv[++i]);
	  setOut=1;
	}
    }


  
  if (setOut==0)
    errorStop("Must provide an output filename using -o\n",dSet,debug);


  pfitsOpenFile(dSet,debug);
  pfitsLoadHeader(dSet,debug);

  fch1 = dSet->head->chanFreq[0];
  
    // Open the output file
  printf("Making a copy of the input file\n");
  sprintf(oname,"!%s",outFile);
  if (!fits_create_file(&outfptr,oname,&status))
    {
      // Copy the file
      ii=1;
      while( !fits_movabs_hdu(dSet->fp, ii++, NULL, &status) )
	{
	  fits_copy_hdu(dSet->fp, outfptr, 0, &status);
	  if (status)
	    {
	      fits_report_error(stderr,status);
	      exit(1);
	    }
	}
    }


  status=0;
  printf("Complete making a copy\n");
  fits_movnam_hdu(outfptr,BINARY_TBL,"SUBINT",0,&status);
  fits_get_colnum(outfptr,CASEINSEN,"DATA",&colnum_out,&status);
  fits_get_colnum(outfptr,CASEINSEN,"DAT_OFFS",&colnum_datoffs,&status);
  fits_get_colnum(outfptr,CASEINSEN,"DAT_SCL",&colnum_datscl,&status);

  // Add FRB
  if (!(fin = fopen(frbFile,"r")))
    {
      printf("Unable to open file %s (use -frb option)\n",frbFile);
      exit(1);
    }
  while (!feof(fin))
    {
      if (fscanf(fin,"%s %f %f",frb[nFRB].fname,&(frb[nFRB].burstTime),&(frb[nFRB].burstDM))==3)
	nFRB++;
    }
  fclose(fin);
  printf("Loaded %d FRBs\n",nFRB);



  
  loadData = (float *)malloc(sizeof(float)*dSet->head->nsblk*dSet->head->nchan);
  outputVals = (unsigned char *)malloc(sizeof(unsigned char)*dSet->head->nsblk*dSet->head->nchan);
  out = (unsigned char *)malloc(sizeof(unsigned char)*dSet->head->nsblk*dSet->head->nchan);
  printf("nbits = %d\n",dSet->head->nbits);
  // Go through each subintegration
  for (i=0;i<dSet->head->nsub;i++)
    {
      t0 = dSet->head->nsblk*i*dSet->head->tsamp;
      t1 = dSet->head->nsblk*(i+1)*dSet->head->tsamp;
      process=0;
      for (ii=0;ii<nFRB;ii++)
	{
	  freq = (fch1 + (dSet->head->nchan-1)*dSet->head->chanbw)/1000;
	  dispersionDelay = 4.15e-3*frb[ii].burstDM*(1.0/(pow(freq,2))-1.0/(pow(fch1/1000,2))); // Note different units 
	  if (frb[ii].burstTime + dispersionDelay > frb[ii].burstTime)
	    {
	      tmin = frb[ii].burstTime;
	      tmax = frb[ii].burstTime+dispersionDelay;
	    }
	  else
	    {
	      tmax = frb[ii].burstTime;
	      tmin = frb[ii].burstTime+dispersionDelay;
	    }
	  //	  printf("minimum/maximum = %g %g\n",tmin,tmax);
	  
	  //	  if ((tmin >= t0 && tmin <= t1) || (tmax >= t0 && tmax <= t1)) 
	  if ((tmin > t0 && tmin < t1) // Start in subint
	      || (tmax > t0 && tmax < t1) // Stop in subint
	      || (tmin < t0 && tmax > t1)) // Cross the subint
	    {
	      //	      if (tmin >= t0 && tmin <= t1)
	      //		printf("Processing: type1 %d\n",process);
	      //	      else
	      //		printf("Processing: type2 %d\n",process);
	      if (process==0)
		{
		  pfits_read1pol_float(loadData,0,dSet,i,i,1,&nSamples,&nTimeSamples,&nFrequencySamples,debug);
		  for (j=0;j<nSamples;j++)
		    out[j] = (unsigned char)loadData[j];
		}
	      process=1;
	      // Now change the data
	      if (!(finBurst = fopen(frb[ii].fname,"rb")))
		{
		  printf("ERROR: cannot open %s\n",frb[ii].fname);
		  exit(1);
		}
	      readBurstHeader(finBurst,&burstSamples);
	      printf("Read header %d\n",burstSamples);
	      for (k=0;k<dSet->head->nchan;k++)
		{
		  freq = (fch1 + k*dSet->head->chanbw)/1000.0;
		  dispersionDelay = 4.15e-3*frb[ii].burstDM*(1.0/(pow(freq,2))-1.0/(pow(fch1/1000,2))); // Note different units for freq and fcentral
		  for (bpos=0;bpos<burstSamples;bpos++)
		    {
		      bt = frb[ii].burstTime+bpos*dSet->head->tsamp;
		      ipos = (long int)(((bt+dispersionDelay)-t0)/(double)dSet->head->tsamp);
		      pos = ipos*dSet->head->nchan + k;
		      //		      printf("subint %d pos = %d %d\n",i,pos,dSet->head->nsblk*dSet->head->nchan);
		      for (a=0;a<4;a++)
			fread(changeVals[a],sizeof(unsigned int),4,finBurst);
			  
		      if (ipos >= 0 && ipos < dSet->head->nsblk)
			{
			  //		      printf("setting: %d %d %g %g %g\n",pos,k,bt+dispersionDelay,freq,fch1);
			  ele = out[pos];
			  totC = changeVals[ele][0] + changeVals[ele][1] + changeVals[ele][2] + changeVals[ele][3];
			  rnum = TKranDev(&seed);
			  randVal = (int)rnum*totC;
			  if (randVal < changeVals[ele][0])
			    out[pos] = 0;
			  else if (randVal < changeVals[ele][0] + changeVals[ele][1])
			    out[pos] = 1;
			  else if (randVal < changeVals[ele][0] + changeVals[ele][1] + changeVals[ele][2])
			    out[pos] = 2;
			  else if (randVal < changeVals[ele][0] + changeVals[ele][1] + changeVals[ele][2] + changeVals[ele][3])
			    out[pos] = 3;
			}
		    }

		}	      
	      fclose(finBurst);
	    }
	}
      if (process==1)
	{
	  // Now update the data in the output file
	  printf("Writing out\n");
	  convertTo2Bits(out, outputVals, dSet->head->nchan*dSet->head->nsblk);
	  fits_write_col(outfptr,TBYTE,colnum_out,i+1,1,dSet->head->nchan*dSet->head->nsblk*dSet->head->nbits/8.0,outputVals,&status);
	  if (status)
	    {
	      fits_report_error(stderr,status);
	      exit(1);
	    }
	  
	}
    }



  fits_close_file(outfptr,&status);
    // De-allocate the memory
  deallocateMemory(&dSet,debug);
  free(frb);
  
}

int convertTo2Bits(unsigned char *inArray, unsigned char *outArray, unsigned int n)
{
  /* take a list of chars in inArray which only have the bottom two */
  /* bits set and pack them into outArray with 4 * 2 bits in each   */
  /* char. n is the number of 2 bit chars in inArray.               */
  /* This routine returns 0 if all is OK, otherwise -1              */

  unsigned char x = 0;
  unsigned int  i;
  unsigned char *outPtr = outArray;
  unsigned char *inPtr = inArray;

  if (n < 0)
    {
      return -1;
    }

  for (i = 0; i < n; i++)
    {
      if ((*inPtr & 0x03) != *inPtr)
	{
	  return -1;
	}

      x = (x << 2) + *inPtr;
      inPtr++;
      if (((i + 1) % 4) == 0)
	{
	  //	  *outPtr = x;
	  //	  printf("Have %d %d\n",*inPtr,x);
	  *outPtr = x;
	  outPtr++;
	  x = 0;
	}
    }
  *outPtr = x;
  return 0;
}

void readBurstHeader(FILE *finBurst,unsigned int *burstSamples)
{
  char format[1024];
  unsigned int nchan;
  double tsamp;
  unsigned int bs;
  unsigned int nbits;
  double tsys;
  double gain;
  double f0;
  double chanbw;
  double peakFlux;
  double width;
  double dm;
  unsigned int s0;
  
  fread(format,sizeof(char),1024,finBurst);                  
  fread(&nchan,sizeof(unsigned int),1,finBurst);
  fread(&tsamp,sizeof(double),1,finBurst);
  fread(&bs,sizeof(unsigned int),1,finBurst);
  fread(&nbits,sizeof(unsigned int),1,finBurst);
  fread(&tsys,sizeof(double),1,finBurst);
  fread(&gain,sizeof(double),1,finBurst);
  fread(&f0,sizeof(double),1,finBurst);
  fread(&chanbw,sizeof(double),1,finBurst);
  fread(&peakFlux,sizeof(double),1,finBurst);
  fread(&width,sizeof(double),1,finBurst);
  fread(&dm,sizeof(double),1,finBurst);
  fread(&s0,sizeof(unsigned int),1,finBurst);

  *burstSamples = bs;
  printf("Read: %s, peakFlux = %g, dm = %g\n",format,peakFlux,dm);
}
