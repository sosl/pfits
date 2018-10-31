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

//gcc -lm -o pfits_getZeroDM pfits_getZeroDM.c pfits_setup.c pfits_loader.c -I/Users/hob044/hob044/software/cfitsio/include/ -L/Users/hob044/hob044/software/cfitsio/lib -lcfitsio 

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "pfits.h"
#include <cpgplot.h>
#include "fitsio.h"

int main(int argc,char *argv[])
{
  dSetStruct *dSet;
  int debug=0;
  int status=0;
  int i,j;
  float *data;

    int subint;
  long nSamples;
  int nTime,nFreq;
  float sum1,sum2;
  long totCount;
  int sub0,sub1;
  int pol=0;
  int histogramTot[255];

  for (i=0;i<255;i++)
    histogramTot[i]=0;
 
  sub0 = 0;
  sub1 = 1;
  
  // Initialise everything
  initialise(&dSet,debug);
  
  // Read inputs
  for (i=0;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	setFilename(argv[++i],dSet,debug);
      else if (strcmp(argv[i],"-s1")==0)
	sscanf(argv[++i],"%d",&sub0);
      else if (strcmp(argv[i],"-s2")==0)
	sscanf(argv[++i],"%d",&sub1);
      else if (strcmp(argv[i],"-pol")==0)
	sscanf(argv[++i],"%d",&pol);
    }

  pfitsOpenFile(dSet,debug);
  pfitsLoadHeader(dSet,debug);

  data = (float *)malloc(sizeof(float)*dSet->head->nsblk*dSet->head->nchan);
  totCount=0;

  for (subint=sub0;subint<sub1;subint++)
    {
      pfits_read1pol_float(data,pol,dSet,subint,subint,1,&nSamples,&nTime,&nFreq,debug);
      for (j=0;j<nTime;j++)
	{
	  sum1=0;
	  sum2=0;
	  for (i=0;i<dSet->head->nchan;i++)
	    {
	      if ((i > 53 && i < 130)) 
		sum1 += data[j*nFreq+i];
	      if ((i > 192 && i < 252)) 
		sum2 += data[j*nFreq+i];
	    }
	  printf("result: %d %g %g %d\n",totCount,sum1,sum2,dSet->head->nchan);
	  totCount++;
	}
    }

  
  free(data);
  // De-allocate the memory
  deallocateMemory(&dSet,debug);
}
