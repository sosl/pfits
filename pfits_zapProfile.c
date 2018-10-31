// Code to zap pulse profiles
//
//gcc -lm -o pfits_zapProfile pfits_zapProfile.c pfits_setup.c -I/Users/hob044/hob044/software/cfitsio/include/ -L/Users/hob044/hob044/software/cfitsio/lib -lcfitsio -lcpgplot -O3

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "pfits.h"
#include <cpgplot.h>
#include "fitsio.h"

void calcMinMax(float *timeAve,int nbin,int nchan,float *minVal,float *maxVal,int *imax);
void saveFile(float *datWts, int nsub, int nchan,dSetStruct *dSet);
void get_dm_tbin(dSetStruct *dSet,double *dm,double *tbin);
void formFrequencyTimeAverages(int nsub,int nchan,int nbin,double dm,double tbin,float fref,float *fChan,float *loadAll_pol1,float *loadAll_pol2,float *timeAve,float *freqAve,float *profile,float *datWts);
void help();
void displayZapCommand(int *zapFreq,int nchan,int *zapSub,int nsub,float *fChan);
int modOp(int a,int b);

void help()
{
  printf("pfits_zapProfile -- program to remove RFI from pulsar fold-mode data sets\n");
  printf("Command line arguments when starting software\n");
  printf("-f <filename>\n");
  printf("-tbin <tbin (sec)> - bin time in seconds - should be read from PSRFITS header, but sometimes missing\n");
  printf("\n");
  printf("Key presses:\n");
  printf("1: plot frequency-phase\n");
  printf("2: plot sub-integration-phase\n");
  printf("3: plot profile\n");
  printf("\n");
  printf("h: this help\n");
  printf("m: set the maximum value for the colour-scale in the plots\n");
  printf("p: display paz zap command\n");
  printf("q: quit - note that the file is not automatically saved - must press 's' to save first\n");
  printf("r: recalculate scalings\n");
  printf("s: save zapped file\n");
  printf("u: un-zoom\n");
  printf("X: (right mouse button) delete individual channel or sub-integration\n");
  printf("z: zoom into a region on the plot\n");
  printf("Z: zap a range of channels or subintegrations\n");  
}


int main(int argc,char *argv[])
{
  dSetStruct *dSet;
  int debug=0;
  int status=0;
  int i;

  float *loadAll_pol1;
  float *loadAll_pol2;
  
  float *freqAve;
  float *timeAve;
  float *profile;
  float *profileX;
  
  short int *sval;
  short int n_val =0;
  float n_fval=0;
  long nchan,nsub,nbin,npol;
  long ii,jj,kk;
  int colnum_data;
  int colnum_datOffs;
  int colnum_datScl;
  int colnum_datWts;
  int initflag = 0;
  int polNum = 0;
  float heat_l[] = {0.0, 0.2, 0.4, 0.6, 1.0};
  float heat_r[] = {0.0, 0.5, 1.0, 1.0, 1.0};
  float heat_g[] = {0.0, 0.0, 0.5, 1.0, 1.0};
  float heat_b[] = {0.0, 0.0, 0.0, 0.3, 1.0};
  float tr[6];
  int   imaxVal_freqAve,imaxVal_timeAve,imaxVal_profile;
  float minVal_freqAve,maxVal_freqAve;
  float minVal_timeAve,maxVal_timeAve;
  float minVal_profile,maxVal_profile;
  
  float *datScl;
  float *datOffs;
  float *datWts;
  float minx,maxx,miny,maxy;
  float *fChan,fref;
  
  double dm=0;
  double tbin=0;
  
  float mx,my;
  char key;
  int didDelete=0;
  
  float val;
  
  double mean1,mean2;
  int plotType=1;
  int dedispBin;
  double tdiff;
  int binOff;

  int *zapFreq;
  int *zapSub;

  // Initialise everything
  initialise(&dSet,debug);
  // Read inputs
  for (i=0;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	setFilename(argv[++i],dSet,debug);
      else if (strcmp(argv[i],"-tbin")==0)
	sscanf(argv[++i],"%lf",&tbin);
    }

  pfitsOpenFile(dSet,debug);
  pfitsLoadHeader(dSet,debug);

  // Now load the data to plot
  // Average in polarisation
  // Produce data files also time averaged and frequency averaged

  nchan = dSet->head->nchan;
  nbin = dSet->head->nbin;
  nsub = dSet->head->nsub;
  npol = dSet->head->npol;

  get_dm_tbin(dSet,&dm,&tbin);
  
  printf("Loaded header\n");
  printf("Number of channels = %d\n",nchan);
  printf("Number of bins = %d\n",nbin);
  printf("Number of sub-integrations = %d\n",nsub);
  printf("Number of polarisations = %d\n",npol);
  printf("Dispersion measure = %g\n",dm);

  loadAll_pol1 = (float *)malloc(sizeof(float)*nchan*nbin*nsub);
  loadAll_pol2 = (float *)malloc(sizeof(float)*nchan*nbin*nsub);

  timeAve = (float *)calloc(nchan*nbin,sizeof(float));
  freqAve = (float *)calloc(nsub*nbin,sizeof(float));
  profile = (float *)calloc(nbin,sizeof(float));
  profileX = (float *)calloc(nbin,sizeof(float));
  fChan  = (float *)malloc(sizeof(float)*nchan);
  sval   = (short int *)calloc(nchan*nbin,sizeof(short int));

  zapFreq = (int *)calloc(nchan,sizeof(int));
  zapSub = (int *)calloc(nsub,sizeof(int));
  
  
  

  
  for (ii=0;ii<nbin;ii++)
    profileX[ii] = ii;
  
  datScl = (float *)malloc(sizeof(float)*nchan*nsub);
  datOffs = (float *)malloc(sizeof(float)*nchan*nsub);
  datWts = (float *)malloc(sizeof(float)*nchan*nsub);
  
  fits_movnam_hdu(dSet->fp,BINARY_TBL,"SUBINT",1,&status);
  fits_get_colnum(dSet->fp,CASEINSEN,"DATA",&colnum_data,&status);
  fits_get_colnum(dSet->fp,CASEINSEN,"DAT_OFFS",&colnum_datOffs,&status);
  fits_get_colnum(dSet->fp,CASEINSEN,"DAT_SCL",&colnum_datScl,&status);
  fits_get_colnum(dSet->fp,CASEINSEN,"DAT_WTS",&colnum_datWts,&status);

  printf("Loading all data\n");
  for (ii=0;ii<nchan;ii++)
      fChan[ii] = dSet->head->chanFreq[ii];
  fref = dSet->head->freq;
  printf("Reference freq = %g\n",fref);
  
  for (ii=0;ii<nsub;ii++)
    {
      // Polarisation 1
      polNum=0;
      fits_read_col(dSet->fp,TFLOAT,colnum_datOffs,ii+1,1+polNum*nchan,nchan,&n_fval,datOffs+ii*nchan,&initflag,&status);
      fits_read_col(dSet->fp,TFLOAT,colnum_datScl,ii+1,1+polNum*nchan,nchan,&n_fval,datScl+ii*nchan,&initflag,&status);
      fits_read_col(dSet->fp,TFLOAT,colnum_datWts,ii+1,1,nchan,&n_fval,datWts+ii*nchan,&initflag,&status);
      fits_read_col(dSet->fp,TSHORT,colnum_data,ii+1,1+polNum*nchan*nbin,nchan*nbin,&n_val,sval,&initflag,&status);
	   
      for (kk=0;kk<nchan;kk++)
	{
	  mean1=0;
	  for (jj=0;jj<nbin;jj++)
	    {
 	      val = (sval[kk*nbin+jj]*datScl[ii*nchan+kk]+datOffs[ii*nchan+kk]);
	      loadAll_pol1[ii*nbin*nchan+kk*nbin+jj] = val;
	      mean1 += val;
	    }
	  for (jj=0;jj<nbin;jj++)
	    loadAll_pol1[ii*nbin*nchan+kk*nbin+jj] -= mean1/(double)nbin;
	}
      
      // Polarisation 2
      polNum=1;
      fits_read_col(dSet->fp,TFLOAT,colnum_datOffs,ii+1,1+polNum*nchan,nchan,&n_fval,datOffs+ii*nchan,&initflag,&status);
      fits_read_col(dSet->fp,TFLOAT,colnum_datScl,ii+1,1+polNum*nchan,nchan,&n_fval,datScl+ii*nchan,&initflag,&status);
      // Don't need to read the weights as not polarisation dependent
      fits_read_col(dSet->fp,TSHORT,colnum_data,ii+1,1+polNum*nchan*nbin,nchan*nbin,&n_val,sval,&initflag,&status);

      for (kk=0;kk<nchan;kk++)
	{
	  mean2=0;
	  for (jj=0;jj<nbin;jj++)
	    {
 	      val = (sval[kk*nbin+jj]*datScl[ii*nchan+kk]+datOffs[ii*nchan+kk]);
	      loadAll_pol2[ii*nbin*nchan+kk*nbin+jj] = val;
	      mean2+=val;
	    }	      
	  for (jj=0;jj<nbin;jj++)
	    loadAll_pol2[ii*nbin*nchan+kk*nbin+jj] -= mean2/(double)nbin;
	}
    }

  formFrequencyTimeAverages(nsub,nchan,nbin,dm,tbin,fref,fChan,loadAll_pol1,loadAll_pol2,timeAve,freqAve,profile,datWts);
  printf("Completed loading all data\n");
  
  // Get minimum and maximum values
  calcMinMax(freqAve,nbin,nsub,&minVal_freqAve,&maxVal_freqAve,&imaxVal_freqAve);
  calcMinMax(timeAve,nbin,nchan,&minVal_timeAve,&maxVal_timeAve,&imaxVal_timeAve);
  calcMinMax(profile,nbin,1,&minVal_profile,&maxVal_profile,&imaxVal_profile);
  printf("Frequency averaging: Minimum value = %g, maximum value = %g (subint for max = %d)\n",minVal_freqAve,maxVal_freqAve,imaxVal_freqAve);
  printf("Time averaging: Minimum value = %g, maximum value = %g (freq. channel for max = %d)\n",minVal_timeAve,maxVal_timeAve,imaxVal_timeAve);
  // Do the plot

  tr[0] = 0;  tr[1] = 1;  tr[2] = 0;
  tr[3] = 0;  tr[4] = 0;  tr[5] = 1;

  //  maxVal = 2000;
  
  cpgbeg(0,"/xs",1,1);
  cpgask(0);

  minx = 0;  maxx = nbin;  miny = 0;  maxy = nchan;

  help();

  
  do {
    cpgenv(minx,maxx,miny,maxy,0,1);
    cpgctab(heat_l,heat_r,heat_g,heat_b,5,1.0,0.5);
    if (plotType==1)
      {
	cpgimag(timeAve,nbin,nchan,1,nbin,1,nchan,minVal_timeAve,maxVal_timeAve,tr);
	cpglab("Bin number","Frequency channel","");
      }
    else if (plotType==2)
      {
	cpgimag(freqAve,nbin,nsub,1,nbin,1,nsub,minVal_freqAve,maxVal_freqAve,tr);
	cpglab("Bin number","Subint number","");
      }
    else if (plotType==3)
      {
	cpgline(nbin,profileX,profile);
	cpglab("Bin number","","");
      }
    cpgcurs(&mx,&my,&key);
    if (key=='z')
      {
	float mx2,my2;
	cpgband(3,0,mx,my,&mx2,&my2,&key);
	if (my > my2) {miny = my2; maxy = my;}
	if (my2 > my) {miny = my; maxy = my2;}
      }
    else if (key=='h')
      help();
    else if (key=='r')
      {
	calcMinMax(freqAve,nbin,nsub,&minVal_freqAve,&maxVal_freqAve,&imaxVal_freqAve);
	calcMinMax(timeAve,nbin,nchan,&minVal_timeAve,&maxVal_timeAve,&imaxVal_timeAve);
	calcMinMax(profile,nbin,1,&minVal_profile,&maxVal_profile,&imaxVal_profile);

      }
    else if (key=='m')
      {
	if (plotType==1)printf("Current maximum value is %g\n",maxVal_timeAve);
	else if (plotType==2)printf("Current maximum value is %g\n",maxVal_freqAve);
	printf("Please enter new maximum value ");
	if (plotType==1)scanf("%f",&maxVal_timeAve);
	else if (plotType==2)scanf("%f",&maxVal_freqAve);	
      }
    else if (key=='1') {plotType=1;   minx = 0;  maxx = nbin;  miny = 0;  maxy = nchan;
      if (didDelete==1)
	formFrequencyTimeAverages(nsub,nchan,nbin,dm,tbin,fref,fChan,loadAll_pol1,loadAll_pol2,timeAve,freqAve,profile,datWts);
      didDelete=0;
    }
    else if (key=='2') {
      plotType=2;   minx = 0;  maxx = nbin;  miny = 0;  maxy = nsub;
      if (didDelete==1)
	formFrequencyTimeAverages(nsub,nchan,nbin,dm,tbin,fref,fChan,loadAll_pol1,loadAll_pol2,timeAve,freqAve,profile,datWts);
      didDelete=0;
    }    
    else if (key=='3') {
      plotType=3;   minx = 0;  maxx = nbin;  miny = minVal_profile;  maxy = maxVal_profile;
      if (didDelete==1)
	formFrequencyTimeAverages(nsub,nchan,nbin,dm,tbin,fref,fChan,loadAll_pol1,loadAll_pol2,timeAve,freqAve,profile,datWts);
      didDelete=0;
    }
    else if (key=='p')
      {
	displayZapCommand(zapFreq,nchan,zapSub,nsub,fChan);
      }
    else if (key=='u')
      {
	if (plotType==1)
	  { miny = 0; maxy = nchan; minx = 0; maxx = nbin;}
	else if (plotType==2)
	  { miny = 0; maxy = nsub; minx = 0; maxx = nbin;}
	else if (plotType==3)
	  { miny = minVal_profile; maxy = maxVal_profile; minx = 0; maxx = nbin;}
      }
    else if (key=='X' && plotType==1)
      {
	int zapChannel;
	didDelete=1;
	zapChannel = (int)(my-0.5); // Check -0.5
	printf("Zapping channel %d\n",zapChannel);
	zapFreq[zapChannel] = 1;
	for (jj=0;jj<nsub;jj++)
	  datWts[jj*nchan+zapChannel] = 0;

	for (jj=0;jj<nbin;jj++)
	  timeAve[zapChannel*nbin+jj] = 0;

	calcMinMax(timeAve,nbin,nchan,&minVal_timeAve,&maxVal_timeAve,&imaxVal_timeAve);
	printf("Minimum value = %g, maximum value = %g, channel number for maximum value = %d\n",minVal_timeAve,maxVal_timeAve,imaxVal_timeAve);
      }
    else if (key=='X' && plotType==2)
      {
	int zap_sub;
	didDelete=1;
	zap_sub = (int)(my-0.5); // Check -0.5
	zapSub[zap_sub] = 1;
	printf("Zapping sub-int %d\n",zap_sub);
	for (jj=0;jj<nchan;jj++)
	  datWts[zap_sub*nchan+jj] = 0;
	
	for (jj=0;jj<nbin;jj++)
	  freqAve[zap_sub*nbin+jj] = 0;
	calcMinMax(freqAve,nbin,nsub,&minVal_freqAve,&maxVal_freqAve,&imaxVal_freqAve);
	printf("Minimum value = %g, maximum value = %g, subint number for maximum value = %d\n",minVal_freqAve,maxVal_freqAve,imaxVal_freqAve);
      }
    else if (key=='X' && plotType==3)
      printf("Sorry -- cannot delete bins in a profile\n");
    else if (key=='s') // Save new file
      {
	saveFile(datWts,nsub,nchan,dSet);
      }
    else if (key=='Z' && plotType==1)
      {
	float mx2,my2;
	int zapChannel;
	float min,max;
	didDelete=1;
	cpgband(3,0,mx,my,&mx2,&my2,&key);
	if (my > my2) {min = my2; max = my;}
	if (my2 > my) {min = my; max = my2;}

	for (zapChannel = (int)(min-0.5);zapChannel <= (int)(max-0.5);zapChannel++)
	  {
	    if (zapChannel > 0 && zapChannel < nchan)
	      {
		zapFreq[zapChannel] = 1;
		printf("Zapping channel %d\n",zapChannel);
		for (jj=0;jj<nsub;jj++)
		  datWts[jj*nchan+zapChannel] = 0;

		for (jj=0;jj<nbin;jj++)
		  timeAve[zapChannel*nbin+jj] = 0;
	      }
	  }
	calcMinMax(timeAve,nbin,nchan,&minVal_timeAve,&maxVal_timeAve,&imaxVal_timeAve);
	printf("Minimum value = %g, maximum value = %g\n",minVal_timeAve,maxVal_timeAve);
      }
    else if (key=='Z' && plotType==2)
      {
	float mx2,my2;
	int zap_sub;
	float min,max;
	didDelete=1;
	cpgband(3,0,mx,my,&mx2,&my2,&key);
	if (my > my2) {min = my2; max = my;}
	if (my2 > my) {min = my; max = my2;}

	for (zap_sub = (int)(min-0.5);zap_sub <= (int)(max-0.5);zap_sub++)
	  {
	    if (zap_sub > 0 && zap_sub < nsub)
	      {
		zapSub[zap_sub] = 1;
		printf("Zapping subint %d\n",zap_sub);
		for (jj=0;jj<nchan;jj++)
		  datWts[zap_sub*nchan+jj] = 0;

		for (jj=0;jj<nbin;jj++)
		  freqAve[zap_sub*nbin+jj] = 0;
	      }
	  }
	calcMinMax(freqAve,nbin,nsub,&minVal_freqAve,&maxVal_freqAve,&imaxVal_freqAve);
	printf("Minimum value = %g, maximum value = %g\n",minVal_freqAve,maxVal_freqAve);
      }
    else if (key!='q')
      printf("Unknown key press %c\n",key);
  } while (key != 'q');
  displayZapCommand(zapFreq,nchan,zapSub,nsub,fChan);
  cpgend();  
  if (status)
    {
      fits_report_error(stderr,status);
      exit(1);
    }
  
  // De-allocate the memory
  //  pfitsCloseFile(dSet,debug);
  deallocateMemory(&dSet,debug);


  free(zapFreq);
  free(zapSub);
  free(loadAll_pol1);
  free(loadAll_pol2);
  free(timeAve);
  free(freqAve);
  free(profile);
  free(profileX);
  free(sval);
  free(fChan);
  free(datScl);
  free(datOffs);
  free(datWts);

}

void calcMinMax(float *timeAve,int nbin,int nchan,float *minVal,float *maxVal,int *imax)
{
  int ii,jj;
  for (ii=0;ii<nchan;ii++)
    {
      for (jj=0;jj<nbin;jj++)
	{
	  if (jj==0 && ii==0)
	    {
	      *imax = 0;
	      *minVal = *maxVal = timeAve[ii*nbin+jj];
	    }
	  else if (*minVal > timeAve[ii*nbin+jj]) *minVal = timeAve[ii*nbin+jj];
	  else if (*maxVal < timeAve[ii*nbin+jj]) {*maxVal = timeAve[ii*nbin+jj]; *imax = ii;}
	}
    }
}

void saveFile(float *datWts, int nsub, int nchan, dSetStruct *dSet)
{
  char outname[1024];
  char temp[1024];
  fitsfile *outfptr;
  int status=0;
  int hdu=1;
  int colnum_datWts;
  long ii;
  
  printf("Enter output filename ");
  scanf("%s",temp);
  sprintf(outname,"!%s",temp);
  printf("Output file name is %s\n",outname);

  fits_create_file(&outfptr, outname, &status);
  if (status) {fits_report_error(stderr, status); printf("NOT WRITING FILE\n"); status=0; return;}
  
  /* Copy every HDU until we get an error */
  while( !fits_movabs_hdu(dSet->fp, hdu++, NULL, &status) )
    {
      fits_copy_hdu(dSet->fp, outfptr, 0, &status);
    }
  /* Reset status after normal error */
  if (status == END_OF_FILE) status = 0;

  // Update the DAT_WTS column
  fits_movnam_hdu(outfptr,BINARY_TBL,(char *)"SUBINT",0,&status);
  fits_get_colnum(outfptr,CASEINSEN,"DAT_WTS",&colnum_datWts,&status);
  for (ii=0;ii<nsub;ii++)
    fits_write_col(outfptr,TFLOAT,colnum_datWts,ii+1,1,nchan,datWts+ii*nchan,&status);
  
  fits_close_file(outfptr, &status);
  if (status) {fits_report_error(stderr, status); printf("FILE MAY BE CORRUPT\n"); status=0; return;}
  printf("Data file written\n");
  
}

void get_dm_tbin(dSetStruct *dSet,double *dm,double *tbin)
{
  int status=0;
  float fdm,ftbin;

  fits_movnam_hdu(dSet->fp,BINARY_TBL,"SUBINT",1,&status);
  fits_read_key(dSet->fp,TFLOAT,"DM",&fdm,NULL,&status);
  *dm = fdm;
  if (*tbin == 0)
    {
      fits_read_key(dSet->fp,TFLOAT,"TBIN",&ftbin,NULL,&status);
      *tbin = ftbin;
    }
  printf("Loaded DM = %g\n",*dm);
  printf("Loaded TBIN = %g\n",*tbin);
}

void formFrequencyTimeAverages(int nsub,int nchan,int nbin,double dm,double tbin,float fref,float *fChan,float *loadAll_pol1,float *loadAll_pol2,float *timeAve,float *freqAve,float *profile,float *datWts)
{
  long int ii,jj,kk;
  double tdiff;
  int binOff;
  int dedispBin;
  int t=0;

  printf("Forming freqTime\n");
  printf("Resetting arrays to zero\n");
  // Must be quicker way to reset all to 0
  for (ii=0;ii<nsub*nbin;ii++)
      freqAve[ii]=0;
  for (ii=0;ii<nchan*nbin;ii++)
      timeAve[ii]=0;
  for (ii=0;ii<nbin;ii++)
    profile[ii] = 0;

  printf("Forming frequency and time averaged data sets\n");
  // Form frequency and time averaged data sets
  for (ii=0;ii<nsub;ii++)
    {
      for (kk=0;kk<nchan;kk++)
	{
	  tdiff = 4.15e-3*dm*(pow(fref/1000.0,-2)-pow(fChan[kk]/1000.0,-2));
	  binOff = (int)(tdiff/tbin+0.5);
	  if (datWts[ii*nchan+kk] > 0)
	    {
	      for (jj=0;jj<nbin;jj++)
		{
		  //
		  //		  dedispBin = jj+binOff;
		  //		  if (dedispBin >= nbin)
		  dedispBin = modOp(jj+binOff,nbin);
		  
		    //		  else
		    //		    while (dedispBin < 0) dedispBin+=nbin;
		  //		  dedispBin = jj+binOff;
		  //
		  //		  while (dedispBin >= nbin) dedispBin-=nbin; // DO THIS BETTER
		  timeAve[kk*nbin+dedispBin] += (loadAll_pol1[ii*nbin*nchan+kk*nbin+jj]+loadAll_pol2[ii*nbin*nchan+kk*nbin+jj]);
		  freqAve[ii*nbin+dedispBin] += (loadAll_pol1[ii*nbin*nchan+kk*nbin+jj]+loadAll_pol2[ii*nbin*nchan+kk*nbin+jj]);
		  profile[dedispBin] += (loadAll_pol1[ii*nbin*nchan+kk*nbin+jj]+loadAll_pol2[ii*nbin*nchan+kk*nbin+jj]);	      
		}
	    }
	}
    }
  printf("Complete forming data sets\n");
}

void displayZapCommand(int *zapFreq,int nchan,int *zapSub,int nsub,float *fChan)
{
  int i;
  int haveZapChannel=0;
  int haveZapSub=0;
  int t=0;
  double chanbw = fabs(fChan[1]-fChan[0]);

  printf("\n");
  printf("paz command with individual subintegrations and channel numbers:\n\n");
  
  printf("paz");
  for (i=0;i<nchan;i++)
    {
      if (zapFreq[i] == 1)
	haveZapChannel=1;
    }
  for (i=0;i<nsub;i++)
    {
      if (zapSub[i] == 1)
	haveZapSub=1;
    }
  
  if (haveZapChannel==1)
    {
      printf(" -z \"");
      t=0;
      for (i=0;i<nchan;i++)
	{
	  if (zapFreq[i]==1)
	    {
	      if (t==0){printf("%d",i); t=1;}
	      else {printf(" %d",i);}
	    }
	}
      printf("\"");
    }

    if (haveZapSub==1)
    {
      printf(" -w \"");
      t=0;
      for (i=0;i<nsub;i++)
	{
	  if (zapSub[i]==1)
	    {
	      if (t==0){printf("%d",i); t=1;}
	      else {printf(" %d",i);}
	    }
	}
      printf("\"");
    }

  printf("\n");

  printf("------------\n");
  printf("paz commands with ranges and frequencies instead of channels\n\n");
  printf("paz ");
  if (haveZapChannel==1)
    {
      int s1,s2;
      int pos=0;

      for (i=0;i<nchan;i++)
	{
	  if (zapFreq[i]==1 && pos==0)
	    {
	      if (pos==0){s1 = i; pos = 1;}		
	    }
	 else if (zapFreq[i]!=1 && pos==1)
	    {
	      s2=i-1;
	      if (fChan[s1] < fChan[s2])
		{printf(" -F \"%g %g\"",fChan[s1]-chanbw/2.,fChan[s2]+chanbw/2.); pos=0;}
	      else
		{printf(" -F \"%g %g\"",fChan[s2]-chanbw/2.,fChan[s1]+chanbw/2.); pos=0;}
	    }
	   }
    }
  
  if (haveZapSub==1)
    {
      int s1,s2;
      int pos=0;
      
      for (i=0;i<nsub;i++)
	{
	  if (zapSub[i]==1 && pos==0)
	    {
	      if (pos==0){s1 = i; pos = 1;}		
	    }
	  else if (zapSub[i]!=1 && pos==1)
	    {
	      s2=i-1;
	      if (s2==s1)
		printf(" -w %d",s1);
	      else
		printf(" -W \"%d %d\"",s1,s2);
	      pos=0;
	    }
	}
    }
  printf("\n");

  
}

int modOp(int a,int b)
{
  int r=a%b;
  return r<0 ? r+b:r;
}
