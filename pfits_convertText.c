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

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "pfits.h"
#include <cpgplot.h>

int main(int argc,char *argv[])
{
  dSetStruct *dSet;
  int i,j,k;
  long nSamp;
  int nTime,nFreq;
  int debug=0;
  float *dataVals;
  FILE *fout;
  long s1=0,s2=0;
  
  // Allocate memory for these files
  dSet = (dSetStruct *)malloc(sizeof(dSetStruct));
  
  // Initialise everything
  initialise(&dSet,debug);
  
  // Read inputs
  for (i=0;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	setFilename(argv[++i],dSet,debug);
      else if (strcmp(argv[i],"-s1")==0)
	sscanf(argv[++i],"%ld",&s1);
      else if (strcmp(argv[i],"-s2")==0)
	sscanf(argv[++i],"%ld",&s2);
    }

  // Open the file and load header information
  pfitsOpenFile(dSet,debug);
  pfitsLoadHeader(dSet,debug);
  dataVals = (float *)malloc(sizeof(float)*dSet->head->nsblk*dSet->head->nchan);
  if (s2==0) s2 = dSet->head->nsub;
  // Do the plot
  fout = fopen("psrData.txt","w");
  for (i=s1;i<s2;i++)
    {
      printf("Subintegration %d\n",i);
      pfits_read1pol_float(dataVals,0,dSet,i,i,1,&nSamp,&nTime,&nFreq,debug);
      for (j=0;j<nTime;j++)
	{
	  for (k=0;k<nFreq;k++)
	    fprintf(fout,"%f\n",dataVals[j*nFreq+k]);
	  fprintf(fout,"\n");
	}
    }
  fclose(fout);
  //  doPlot(dSet,nFiles,plot,debug);
 
  
  // De-allocate the memory
  deallocateMemory(&dSet,debug);
  free(dataVals);
}

