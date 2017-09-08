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

//gcc -lm -o pfits_simRealData pfits_simRealData.c pfits_loader.c T2toolkit.c pfits_setup.c -I/Users/hob044/hob044/software/cfitsio/include/ -L/Users/hob044/hob044/software/cfitsio/lib -lcfitsio 

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "pfits.h"
#include <cpgplot.h>
#include "fitsio.h"
#include "T2toolkit.h"

void process(dSetStruct *dSet,fitsfile *outfptr,int debug,float dm,float probOn,float burstTime,float si,float width);
void processDataSegment(dSetStruct *dSet,unsigned char *inArrChar_1,int n,int nsamp,float *fltArray);
void multMatrix(unsigned char p1,unsigned char p2,unsigned char p3,unsigned char p4,float mat[4][4],
		float *out1,float *out2,float *out3,float *out4);
void digitise1bit(unsigned char *data,unsigned char *cdata,int n);
int digitise2bit(unsigned char *data,unsigned char *cdata,int n_in);

void help()
{
  printf("-f <input filename>\n");
  printf("-o <output filename>\n");
  printf("-h this help\n");
  printf("-dm <dm of burst event>\n");
  printf("-probOn <the probability that the digitiser = 1 for every sample in the burst - value between 0 and 1>\n");
  printf("-burstTime <time in seconds of burst from start of observation>\n");
  printf("-si <spectral index, default = -2>\n");
  printf("-width <width in seconds of burst>\n");
  exit(1);
}

int main(int argc,char *argv[])
{
  dSetStruct *dSet;
  int debug=0;
  int status=0;
  int i,ii;
  char outFile[1024],oname[1024]; // Output filename
  int setOut=0;
  fitsfile *outfptr;
  float fref=-1;

  float dm=800;
  float probOn = 0.8;
  float burstTime =  1726;
  float si=-2;
  float width = 0.03;
  
  // Initialise everything
  initialise(&dSet,debug);
  
  // Read inputs
  for (i=0;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	setFilename(argv[++i],dSet,debug);
      else if (strcmp(argv[i],"-dm")==0)
	sscanf(argv[++i],"%f",&dm);
      else if (strcmp(argv[i],"-probOn")==0)
	sscanf(argv[++i],"%f",&probOn);
      else if (strcmp(argv[i],"-burstTime")==0)
	sscanf(argv[++i],"%f",&burstTime);
      else if (strcmp(argv[i],"-si")==0)
	sscanf(argv[++i],"%f",&si);
      else if (strcmp(argv[i],"-width")==0)
	sscanf(argv[++i],"%f",&width);
      else if (strcmp(argv[i],"-h")==0)
	help();
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

  // Open the output file
  printf("Making a copy of the input file\n");
  sprintf(oname,"!%s",outFile);
  if (!fits_create_file(&outfptr,oname,&status))
    {
      // Copy the file
      ii=1;
      while( !fits_movabs_hdu(dSet->fp, ii++, NULL, &status) )
	fits_copy_hdu(dSet->fp, outfptr, 0, &status);
    }
  status=0;
  printf("Complete making a copy\n");
  // Close the file
  //  pfitsCloseFile(dSet,debug);

  process(dSet,outfptr,debug,dm,probOn,burstTime,si,width);



  fits_close_file(outfptr,&status);
  if (status)
    {
      fits_report_error(stderr,status);
      exit(1);
    }
  
  // De-allocate the memory
  deallocateMemory(&dSet,debug);
}

void process(dSetStruct *dSet,fitsfile *outfptr,int debug,float dm,float probOn,float burstTime,float si,float width)
{
  int status=0;

  double t0,t1,tt;
  float ranNum;
  long seed=TKsetSeed();
  float fref;
  int sub0,sub1;
  int isamp0,isamp1;
  int sub;
  int i0,i1,i;
  int colnum_in;
  int colnum_out;
  char *inArrChar;
  float *invBits;
  char *maskChar;
  unsigned char nval = 0;
  int initflag=0;
  int k;
  int bitsinbyte=8/dSet->head->nbits;
  double timeSample;

  fits_get_colnum(dSet->fp,CASEINSEN,"DATA",&colnum_in,&status);
  fits_get_colnum(outfptr,CASEINSEN,"DATA",&colnum_out,&status);

  fref = 0.5*(dSet->head->chanFreq[0] + dSet->head->chanFreq[dSet->head->nchan-1]);
  printf("fref = %g\n",fref);
  printf("Width = %g\n",width);
  inArrChar = (unsigned char *)malloc(sizeof(unsigned char)*dSet->head->nchan/bitsinbyte);
  maskChar = (unsigned char *)malloc(sizeof(unsigned char)*dSet->head->nchan);
  invBits = (float *)malloc(sizeof(float)*dSet->head->nchan);

  // Add a burst at a given time
  t0 = burstTime - (4.15e-3*dm*(pow(fref/1000.0,si)-pow((dSet->head->chanFreq[0]+dSet->head->chanbw*dSet->head->nchan)/1000.0,si)));
  t1 = burstTime - (4.15e-3*dm*(pow(fref/1000.0,si)-pow((dSet->head->chanFreq[0])/1000.0,si)));
  if (t0 > t1)
    {
      tt = t0;
      t0 = t1;
      t1 = tt;
    }

  sub0 = (int)((t0-width/2.0)/(dSet->head->tsamp*dSet->head->nsblk));
  sub1 = (int)((t1+width/2.0)/(dSet->head->tsamp*dSet->head->nsblk));
  isamp0 = ((t0-width/2.0)-sub0*(dSet->head->tsamp*dSet->head->nsblk))/dSet->head->tsamp;
  isamp1 = ((t1+width/2.0)-sub1*(dSet->head->tsamp*dSet->head->nsblk))/dSet->head->tsamp;

  printf("time %g %g\n",t0,t1);
  printf("Subint %d %d\n",sub0,sub1);
  printf("Sample number = %d %d\n",isamp0,isamp1);
  for (sub=sub0;sub<=sub1;sub++)
    {
      if (sub==sub0)
	i0 = isamp0;
      else
	i0 = 0;

      if (sub==sub1)
	i1 = isamp1;
      else
	i1 = dSet->head->nsblk-1;

      printf("Trying sub = %d\n",sub);

      for (i=i0;i<=i1;i++)
	{
	  // Now we want to read this data, change it and write it
	  fits_read_col_byt(dSet->fp,colnum_in,sub+1,i*dSet->head->nchan/bitsinbyte+1,dSet->head->nchan/bitsinbyte,nval,inArrChar,&initflag,&status); 
	  pfits_bytesToFloats(bitsinbyte,dSet->head->nchan,inArrChar,invBits);
	  for (k=0;k<dSet->head->nchan;k++)
	    {
	      maskChar[k] = (unsigned char)invBits[k];
	      t0 = burstTime - (4.15e-3*dm*(pow(fref/1000.0,si)-pow((dSet->head->chanFreq[k])/1000.0,si)));	      
	      timeSample = sub*dSet->head->nsblk*dSet->head->tsamp+i*dSet->head->tsamp;
	      //	      printf("Timesample = %g %g\n",timeSample,width);
	      // Check if this time corresponds to this sample
	      //	      printf("Timesample = %g\n",timeSample);
	      //	      printf("Checking: %d %g %g %g %g %g timesample = %g\n",i,t0,timeSample + dSet->head->tsamp/2.0 + width/2.0,timeSample - dSet->head->tsamp/2 - width/2.0,dSet->head->tsamp,width,timeSample);
	      if (t0 < timeSample + dSet->head->tsamp/2 + width/2.0 && t0 > timeSample - dSet->head->tsamp/2 - width/2.0)
		{
		  //		  printf("In here %d %d %g\n",i,k,dSet->head->tsamp);
		  ranNum = TKranDev(&seed);
		  if (ranNum > 1-probOn)
		    {
		      if (bitsinbyte == 8)
			maskChar[k] = 0; // This seems to be upside down	      
		      else if (bitsinbyte == 4)
			maskChar[k] = 3; 	      
		    }
		}
	    }
	  if (bitsinbyte==8)
	    digitise1bit(maskChar,inArrChar,dSet->head->nchan/bitsinbyte);
	  else if (bitsinbyte==4)
	    {
	      if (digitise2bit(maskChar,inArrChar,dSet->head->nchan)==-1)
		{
		  printf("ERROR: 2 bit failed\n");
		  exit(1);
		}
	    }
	  else
	    {
	      printf("ERROR: Don't know how to deal with %d bits in byte\n",bitsinbyte);
	      exit(1);
	    }
	  fits_write_col_byt(outfptr,colnum_out,sub+1,i*dSet->head->nchan/bitsinbyte+1,dSet->head->nchan/bitsinbyte,inArrChar,&status);
	}
    }
  printf("Finishing\n");


  free(inArrChar);
  free(invBits);
  free(maskChar);
}


void digitise1bit(unsigned char *data,unsigned char *cdata,int n_in)
{
  int i,j;
  unsigned char tc;
  double bit_level=0;
  long n=0;
  for (i=0;i<n_in;i++)
    {
      tc=0;
      for (j=0;j<8;j++)
	{
	  if (data[n] == 0)
	    tc = tc | (1 << (7-j));
	  n++;
	}
      cdata[i] = tc;
    }

}

int digitise2bit(unsigned char *inArray,unsigned char *outArray,int n)
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
	  *outPtr = x;
	  outPtr++;
	  x = 0;
	}
    }
  *outPtr = x;
  return 0;
}
