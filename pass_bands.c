#include "pass_bands.h"
#include <stdio.h>
#include <stdlib.h>

//#define REDSHIFT2

int n_pb;      //number of samples on the curve
double *x_pb;  //wavelength
double *y_pb;  //transmissivity

//#define FIT_LAE

//load the pass band from a file
void LoadPassBand(int l)
{
  FILE *fp;
  char fname[200];
  int i;
  double xb, yb;
  double fac = 1.0;

  //select the pass band to read in
  switch(l)
  {

#ifndef REDSHIFT2
#ifndef FIT_LAE
    case 0:
      sprintf(fname,"pass_bands/f336w_trans.data");
      break;
    case 1:
      sprintf(fname,"pass_bands/u_trans.data");
      break;
    case 2:
      sprintf(fname,"pass_bands/B_trans.data");
      break;
    case 3:
      sprintf(fname,"pass_bands/NB497_trans.data");
      break;
    case 4:
      sprintf(fname,"pass_bands/V_trans.data");
      break;
    case 5:
      sprintf(fname,"pass_bands/R_trans.data");
      break;
    case 6:
      sprintf(fname,"pass_bands/i_trans.data");
      break;
    case 7:
      sprintf(fname,"pass_bands/z_trans.data");
      break;
    case 8:
      sprintf(fname,"pass_bands/J_wfcam_trans.data");
      break;
    case 9:
      sprintf(fname,"pass_bands/f160w_trans.data");
      break;
    case 10:
      sprintf(fname,"pass_bands/K_wfcam_trans.data");
      break;
    case 11:
      sprintf(fname,"pass_bands/ch1_trans.data");
      fac = 1.0e4;
      break;
    case 12:
      sprintf(fname,"pass_bands/ch2_trans.data");
      fac = 1.0e4;
      break;
    case 13:
      sprintf(fname,"pass_bands/ch3_trans.data");
      fac = 1.0e4;
      break;
    case 14:
      sprintf(fname,"pass_bands/ch4_trans.data");
      fac = 1.0e4;
      break;
#else  //FIT_LAE

    //This case just tests a single high-redshift LAE
    case 0:
      sprintf(fname,"pass_bands/B_trans.data");
      break;
    case 1:
      sprintf(fname,"pass_bands/V_trans.data");
      break;
    case 2:
      sprintf(fname,"pass_bands/R_trans.data");
      break;
    case 3:
      sprintf(fname,"pass_bands/z_trans.data");
      break;
    case 4:
      sprintf(fname,"pass_bands/K_wfcam_trans.data");
      break;
    case 5:
      sprintf(fname,"pass_bands/ch1_trans.data");
      fac = 1.0e4;
      break;
    case 6:
      sprintf(fname,"pass_bands/ch2_trans.data");
      fac = 1.0e4;
      break;

#endif //FIT_LAE

#else //REDSHIFT2
    case 0:
      sprintf(fname,"pass_bands/f275w_trans.dat");
      break;
    case 1:
      sprintf(fname,"pass_bands/f336w_trans.dat");
      break;
    case 2:
      sprintf(fname,"pass_bands/f435w_trans.dat");
      break;
    case 3:
      sprintf(fname,"pass_bands/f606w_trans.dat");
      break;
    case 4:
      sprintf(fname,"pass_bands/f775w_trans.dat");
      break;
    case 5:
      sprintf(fname,"pass_bands/f814w_trans.dat");
      break;
    case 6:
      sprintf(fname,"pass_bands/f850lp_trans.dat");
      break;
    case 7:
      sprintf(fname,"pass_bands/f105w_trans.dat");
      break;
    case 8:
      sprintf(fname,"pass_bands/f125w_trans.dat");
      break;
    case 9:
      sprintf(fname,"pass_bands/f140w_trans.dat");
      break;
    case 10:
      sprintf(fname,"pass_bands/f160w_trans.dat");
      break;
#endif //REDSHIFT2
  }

  printf("Fname = %s\n",fname);

  //load the passband transmissivity curve
  if(!(fp = fopen(fname,"r")))
  {
    printf("Error loading file %s.\n",fname);
    exit(0);
  }
  fscanf(fp,"%d\n",&n_pb);

  //allocate the pass band array
  x_pb = (double *) calloc(n_pb,sizeof(double));
  y_pb = (double *) calloc(n_pb,sizeof(double));

  //load the passband, with the wavelength
  //array stored in Angstroms in the restframe
  for(i=0;i<n_pb;i++)
  {
    fscanf(fp,"%lf %lf\n",&xb,&yb);
    x_pb[i] = xb*fac;
    y_pb[i] = yb;
  }

  //close the passband file
  fclose(fp);
}
void FreePassBand(void)
{
  free(x_pb);
  free(y_pb);
  n_pb = 0;
}
