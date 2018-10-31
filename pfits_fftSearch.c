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
// gcc -lm -o pfits_fftSearch pfits_fftSearch.c pfits_setup.c -I/Users/hob044/hob044/software/cfitsio/include/ -L/Users/hob044/hob044/software/cfitsio/lib -lcfitsio


#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "pfits.h"
#include <cpgplot.h>
#include "fitsio.h"

#define MAX_INPUT_FILES 20

int main(int argc,char *argv[])
{
  dSetStruct **dSet;
  char fname[MAX_INPUT_FILES][1024];
  int nFiles=0;
  int debug=0;
  float tsub;
  float requestTimeToProces;
  int i;

  for (i=0;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	strcpy(fname[nFiles++],argv[++i]);            
    }
  if (nFiles==0)
    {
      printf("Haven't loaded any data. Use -f to specify PSRFITS files to load\n");
      exit(1);
    }

  
  dSet = (dSetStruct **)malloc(sizeof(dSetStruct *)*nFiles);
  
  for (i=0;i<nFiles;i++)
    {
      initialise(&(dSet[i]),debug);
      setFilename(fname[i],dSet[i],debug);
      pfitsOpenFile(dSet[i],debug);
      pfitsLoadHeader(dSet[i],debug);
      printf("Should check the consistancy of this file with other files - not doing this\n");
    }  
  printf("Doing processing\n");

  // Select suitable number of subintegrations for processing
  tsub = dSet[0]->head->nsblk*dSet[0]->head->tsamp;
  printf("Subintegration time is %g s\n",tsub);

  
  for (i=0;i<nFiles;i++)
    deallocateMemory(&dSet[i],debug);

  free(dSet);
}
