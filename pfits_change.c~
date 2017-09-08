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

//gcc -lm -o pfits_change pfits_change.c pfits_setup.c -I/Users/hob044/hob044/software/cfitsio/include/ -L/Users/hob044/hob044/software/cfitsio/lib -lcfitsio 

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
  FILE *fin;
  int debug=0;
  int status=0;
  int i,ii;
  char outFile[1024],oname[1024]; // Output filename
  int setOut=0;
  fitsfile *outfptr;
  char replaceFile[1024];
  int setReplace=0;
  int colnum_out;
  int colnum_datoffs;
  int colnum_datscl;

  // Initialise everything
  initialise(&dSet,debug);
  
  // Read inputs
  for (i=0;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	setFilename(argv[++i],dSet,debug);
      else if (strcmp(argv[i],"-o")==0)
	{
	  strcpy(outFile,argv[++i]);
	  setOut=1;
	}
      else if (strcmp(argv[i],"-r")==0)
	{
	  strcpy(replaceFile,argv[++i]);
	  setReplace=1;
	}
    }

  if (setOut==0)
    errorStop("Must provide an output filename using -o\n",dSet,debug);

  if (setReplace==0)
    errorStop("Must provide a filename containing the replace file using -r\n",dSet,debug);

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
  fits_movnam_hdu(outfptr,BINARY_TBL,"SUBINT",0,&status);
  fits_get_colnum(outfptr,CASEINSEN,"DATA",&colnum_out,&status);
  fits_get_colnum(outfptr,CASEINSEN,"DAT_OFFS",&colnum_datoffs,&status);
  fits_get_colnum(outfptr,CASEINSEN,"DAT_SCL",&colnum_datscl,&status);
  printf("column = %d %d %d\n",colnum_out,colnum_datoffs,colnum_datscl);

  // Now replace the necessary data
  if (!(fin = fopen(replaceFile,"r")))
    {
      printf("Unable to open file %s\n",replaceFile);
    }
  else
    {
      int s,i;
      float vals[1024];
      float dat_scl, dat_offs;
      short fvals[dSet->head->nbin];
      float min,max;
      float maxShort = 32766;
      float minShort = -32767;

      for (s=0;s<dSet->head->nsub;s++)
	{
	  printf("Subintegration %d\n",s);
	  for (i=0;i<dSet->head->nbin;i++)
	    {
	      fscanf(fin,"%f",&vals[i]);
	      if (i==0)
		{
		  min = max = vals[i];
		}
	      else
		{
		  if (min > vals[i]) min = vals[i];
		  if (max < vals[i]) max = vals[i];
		}
	      //	      printf("bin %d %g\n",i,vals[i]);
	    }
	  // Calculate DAT_SCL and DAT_OFFS;
	  printf("min/max = %g %g\n",min,max);
	  dat_scl = (max-min)/(maxShort-minShort);
	  dat_offs = max -maxShort*dat_scl;
	  printf("dat_scl/dat_offs = %g %g\n",dat_scl,dat_offs);
	  for (i=0;i<dSet->head->nbin;i++)
	    fvals[i] = (short)((vals[i]-dat_offs)/dat_scl);
	  fits_write_col(outfptr,TSHORT,colnum_out,s+1,1,dSet->head->nbin,fvals,&status);
	  if (status)
	    {
	      printf("Error in writing columns\n");
	      exit(1);
	    }
	  // Now update dat_offs and dat_scl
	  fits_write_col(outfptr,TFLOAT,colnum_datoffs,s+1,1,1,&dat_offs,&status);
	  fits_write_col(outfptr,TFLOAT,colnum_datscl,s+1,1,1,&dat_scl,&status);
	}
      fclose(fin);
    }

  fits_close_file(outfptr,&status);
  if (status)
    {
      fits_report_error(stderr,status);
      exit(1);
    }

  
  // De-allocate the memory
  deallocateMemory(&dSet,debug);
}

