#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "bpass_models.h"
#include "pass_bands.h"
#include "gsl/gsl_interp.h"

//#define REDSHIFT2
#define IGM_ABS

#define ALTERNATE_LINE_EMISSION


//#define PBINTERP //old pass band integration; comment out for new model


//#define FIT_LAE

#ifndef IGM_ABS
#include "madau_absorption.h"
#else //IGM_ABS
#include "igm-absorption.h"
#endif //IGM_ABS

int n_age;          //number of age bins
int n_z_met;        //number of metallicity samples
int n_f_bin;        //number of binarity samples
int n_lambda;       //number of wavelength samples in the SED model
double *age;        //array of age values
double **N_Lyc;      //array of N_Lyc values

double *z_met;      //array of metallicity values
double *f_bin;      //array of binarity values
double *lambda;     //array of wavelengths
double *f_nu;       //array of SED samples
int flag_binary;    //what binarity are we fitting?


//line model 
int n_lambda_line;
double *f_nu_line;      //array of line f_nu
double *lambda_line;    //array of wavelengths

#ifdef ALTERNATE_LINE_EMISSION
int n_lines;
double *f_lines;
double *lambda_lines;
#endif //ALTERNATE_LINE_EMISSION

//cont model 
int n_lambda_cont;
double *f_nu_cont;      //array of cont f_nu
double *lambda_cont;    //array of wavelengths

int ID;             //galaxy ID

//source properties
double z;
int     n_data;
double *data_x;
double *data_y;
double *data_ye;
int *fit_flag;

//sed model samples
double ****sed_model;  
double ****line_model;  
double ****cont_model;  

//absorption
double *absorption;

//distance modulus
double distance_modulus(double redshift)
{
#ifndef FIT_LAE

#ifndef REDSHIFT2
  int ndm = 14;
  double z_arr[14] = {2.5, 3.0, 3.05, 3.1, 3.15, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0};
  double DM_arr[14] = {4.518714e+01, 4.552045e+01, 4.555013e+01, 4.557923e+01, 4.560779e+01, 4.563582e+01, 4.569035e+01, 4.574297e+01, 4.579379e+01, 4.584292e+01,  4.589048e+01, 4.593654e+01, 4.598120e+01, 4.602453e+01};
#else //REDSHIFT2
  int ndm = 30;
  double z_arr[30] = {1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0};
  double DM_arr[30] = {43.55067698, 43.73397670, 43.90075223, 44.05349845, 44.19420196, 44.32446912, 44.44561612, 44.55873390, 44.66473595, 44.76439412, 44.85836591, 44.94721554, 45.03143042, 45.11143418, 45.18759710, 45.26024451, 45.32966360, 45.39610904, 45.45980759, 45.52096195, 45.57975397, 45.63634735, 45.69088992, 45.74351559, 45.79434600, 45.84349196, 45.89105464, 45.93712665, 45.98179298, 46.02513176};
#endif //REDSHIFT2
  double dz;
  double dm;
  int i = gsl_interp_bsearch(z_arr,redshift,0,ndm);
  dz = (redshift-z_arr[i])/(z_arr[i+1] - z_arr[i]);
  dm = (1.0-dz)*DM_arr[i] + dz*DM_arr[i+1];
  return pow(10., -0.4*dm);
#else //FIT_LAE
  //return pow(10.,-0.4*4.666091e+01);//5.7
  return pow(10.,-0.4*4.563294e+01); //3.1

#endif //FIT_LAE
}

//allocate the sed buffers
void allocate_sed_buffers(void)
{
  n_lambda = 100000;
  lambda = (double *) malloc(n_lambda*sizeof(double));
  f_nu   = (double *) malloc(n_lambda*sizeof(double));

  n_lambda_line = 34966;
  lambda_line = (double *) malloc(n_lambda_line*sizeof(double));
  f_nu_line   = (double *) malloc(n_lambda_line*sizeof(double));
  n_lambda_cont = 34966;
  lambda_cont = (double *) malloc(n_lambda_cont*sizeof(double));
  f_nu_cont   = (double *) malloc(n_lambda_cont*sizeof(double));
}
//free sed buffer memory
void deallocate_sed_buffers(void)
{
  free(lambda);
  free(f_nu);
  free(lambda_line);
  free(f_nu_line);
  free(lambda_cont);
  free(f_nu_cont);
}

void allocate_data(void)
{
  data_x      = (double *) malloc(n_data*sizeof(double));
  data_y      = (double *) malloc(n_data*sizeof(double));
  data_ye     = (double *) malloc(n_data*sizeof(double));
  fit_flag    = (int *) malloc(n_data*sizeof(int));
}
void free_data(void)
{
  free(data_x);
  free(data_y);
  free(data_ye);
  free(fit_flag);

}
void read_data(char fbase[])
{
  FILE *fp;
  char fdir[200];
  char fname[200];
  char fbuf[2];
  char fgbuf[10];

  double xb, yb, yeb;
  int ffb;
  int i, ib, jb;

  fbuf[0] = fbase[14];
  fbuf[1] = fbase[15];
  //sprintf(fbuf,"%s",&fbase[15]);
  //printf("%s\n",fbuf);

  //determine the source ID
  ID = atoi(fbuf);
  printf("ID %d\n",ID);
  //exit(0);


#ifndef FIT_LAE
  //first, read in the redshift
  sprintf(fname,"redshifts/redshifts.txt");
  if(!(fp=fopen(fname,"r")))
  {
    printf("Error opening %s.\n",fname);
    exit(-1);
  }
  for(i=0;i<64;i++)
  {
    fscanf(fp,"%s %lf %d %d\n",&fgbuf[0],&z,&ib,&jb);
    //printf("%s %f\n",fgbuf,z);
    if(i==(ID-1))
      break;
  }
  fclose(fp);
#else //FIT_LAE
  //z = 5.7;
  z = 3.1;

#endif //FIT_LAE
  printf("z = %e\n",z);

#ifndef REDSHIFT2
  if(z<3.0)
  {
    printf("Object lies at z<3.\n");
    printf("Aborting....\n");
    exit(-1);
  }
#endif //REDSHIFT2

  //next, read in the photometry

  //sprintf(fdir,"photometry/fnudata_v2/");
  sprintf(fdir,"photometry/");
  sprintf(fname,"%s%s",fdir,fbase);

  printf("Opening %s.\n",fname);

  if(!(fp=fopen(fname,"r")))
  {
    printf("Error opening %s.\n",fname);
    exit(-1);
  }

  //all data now has 14 data points
  //n_data = 14;

  //all data now has 15 data points
  //n_data = 15;


#ifndef REDSHIFT2

  //all data now has 13 data points
  n_data = 13;

#ifdef FIT_LAE
  //lae 5.7
  n_data = 4;

  //lae 3.1
  n_data = 7;

#endif //FIT_LAE
#else //REDSHIFT2

  n_data = 11;
#endif //REDSHIFT2

  //read in the number of SED samples
  printf("Number of data points = %d.\n",n_data);

  //allocate the data arrays
  allocate_data();

  //read the SED samples
  for(i=0;i<n_data;i++)
  {
    fscanf(fp,"%lf %lf %lf %d\n",&xb,&yb,&yeb,&ffb);
    printf("x % e y % e ye % e fit flag %d\n",xb,yb,yeb,ffb);
    data_x[i] = xb/(1+z);  //wavelength in restframe
    data_y[i] = yb;
    data_ye[i] = yeb;
    fit_flag[i] = ffb;
  }

  //close the data file
  fclose(fp);
}

//allocate the sed sample array
void allocate_sed_samples(void)
{
  int i,j,k,l;
  
  //n_bin x n_met x n_age x n_data 
  //array of sed values

  sed_model = (double ****) malloc(n_f_bin*sizeof(double ***));
  line_model = (double ****) malloc(n_f_bin*sizeof(double ***));
  cont_model = (double ****) malloc(n_f_bin*sizeof(double ***));

  for(i=0;i<n_f_bin;i++)
  {
    //allocate metallicity values
    sed_model[i] = (double ***) malloc(n_z_met*sizeof(double **));
    line_model[i] = (double ***) malloc(n_z_met*sizeof(double **));
    cont_model[i] = (double ***) malloc(n_z_met*sizeof(double **));

    for(j=0;j<n_z_met;j++)
    {
      sed_model[i][j] = (double **) malloc(n_age*sizeof(double *));
      line_model[i][j] = (double **) malloc(n_age*sizeof(double *));
      cont_model[i][j] = (double **) malloc(n_age*sizeof(double *));

      for(k=0;k<n_age;k++)
      {
        sed_model[i][j][k] = (double *) malloc(n_data*sizeof(double));
        line_model[i][j][k] = (double *) malloc(n_data*sizeof(double));
        cont_model[i][j][k] = (double *) malloc(n_data*sizeof(double));

        for(l=0;l<n_data;l++)
        {
          sed_model[i][j][k][l] = 0.0;  //initialize to zero
          line_model[i][j][k][l] = 0.0;  //initialize to zero
          cont_model[i][j][k][l] = 0.0;  //initialize to zero

        }
      }
    }
  }

  //Now Lyman continuum production rates
  N_Lyc = (double **) malloc(n_z_met*sizeof(double *));
  for(i=0;i<n_z_met;i++)
  {
    N_Lyc[i] = (double *) malloc(n_age*sizeof(double));
  }
}

//deallocate the sed sample array
void deallocate_sed_samples(void)
{
  int i,j,k,l;
  
  //n_bin x n_met x n_age x n_data 
  //array of sed values

  for(i=0;i<n_f_bin;i++)
  {
    for(j=0;j<n_z_met;j++)
    {
      for(k=0;k<n_age;k++)
      {
        free(sed_model[i][j][k]);
        free(line_model[i][j][k]);
        free(cont_model[i][j][k]);
      }
      free(sed_model[i][j]);
      free(line_model[i][j]);
      free(cont_model[i][j]);

    }
    free(sed_model[i]);
    free(line_model[i]);
    free(cont_model[i]);

  }
  free(sed_model);
  free(line_model);
  free(cont_model);
  for(j=0;j<n_z_met;j++)
  {
    free(N_Lyc[j]);
  }
  free(N_Lyc);
}

void initialize_bpass_models(void)
{
  int i;
  double age_in[50]  = {6.10,6.20,6.30,6.40,6.50,6.60,6.70,6.80,6.90,7.00,7.10,7.20,7.30,7.40,7.50,7.60,7.70,7.80,7.90,8.00,8.10,8.20,8.30,8.40,8.50,8.60,8.70,8.80,8.90,9.00,9.10,9.20,9.30,9.40,9.50,9.60,9.70,9.80,9.90,10.00,10.10,10.20,10.30,10.40,10.50,10.60,10.70,10.80,10.90,11.00};
  double z_met_in[1] = {0.01};
  double f_bin_in[1] = {0.0};

  ///////////////////////////
  // AGE
  ///////////////////////////

  //set n_age = 50;
  n_age = 50;

  //allocate age
  age = (double *) malloc(n_age*sizeof(double));

  //initialize age array
  for(i=0;i<n_age;i++)
  {
    age[i] = age_in[i];
    printf("age %e\n",age[i]);
  }


  ///////////////////////////
  // METALLICITY
  ///////////////////////////
  
  //set n_z_met = 1;
  n_z_met = 1;

  //allocate z_met
  z_met = (double *) malloc(n_z_met*sizeof(double));

  //initialize metallicity array
  for(i=0;i<n_z_met;i++)
    z_met[i] = z_met_in[i];

  ///////////////////////////
  // BINARITY
  ///////////////////////////
  
  //set n_f_bin = 1;
  n_f_bin = 1;

  //allocate f_bin
  f_bin = (double *) malloc(n_f_bin*sizeof(double));

  //initialize binarity array
  for(i=0;i<n_f_bin;i++)
    f_bin[i] = f_bin_in[i];

  //at this point, the range of model parameter
  //values we will consider have been recorded.
  //we now need to sample the models and create
  //a grid of model parameters at the necessary
  //wavelengths.

  //first, allocate the array containing
  //the model SED samples
  allocate_sed_samples();

  //next, load the SED samples
  //while correcting for IGM attenuation
  load_and_attenuate_sed_samples();

  //our model is now fully loaded
  //and ready to use as an interpolation
  //grid for finding a good fit to the SED data
}

//load the zmet filenames
char **load_z_met_filenames(void)
{
  int i;
  char ** fname_z_met = (char **) malloc(n_z_met*sizeof(char *));
  for(i=0;i<n_z_met;i++)
    fname_z_met[i] = (char *) malloc(6*sizeof(char));

  //load the metallicity file name components

  //only one for right now
  //sprintf(fname_z_met[0],"z010");
  //sprintf(fname_z_met[0],"zem5");
  sprintf(fname_z_met[0],"z002");

  return fname_z_met;
}
//free memory for metallicity file name components
void free_z_met_filenames(char **fname_z_met)
{
  int i;
  for(i=0;i<n_z_met;i++)
    free(fname_z_met[i]);
  free(fname_z_met);
}

//now we load the SED model
void load_sed_model(char fname[], double z)
{
  int i;
  double xb, yb;
  FILE *fp;
  double A_10pc;
  double A_DM = distance_modulus(z); //look up the distance modulus
  double pc_in_cm = 3.085678e18; //pc in cm
  double uJy = 1.0e-29; // microJansky in erg/s/Hz/cm^2

  //area of sphere of 10pc radius in cm^2
  A_10pc = 4.0*M_PI*pow(10*pc_in_cm,2);


  //print out the redshift and distance modulus
  printf("z %e DM %e %e\n",z,A_DM,-2.5*log10(A_DM));
  printf("Loading file = %s.\n",fname);
  if(!(fp = fopen(fname,"r")))
  {
    printf("Error opening %s.\n",fname);
    exit(0);
  }
  for(i=0;i<n_lambda;i++)
  {
    fscanf(fp,"%lf\t%lf\n",&xb,&yb);
    lambda[i] = xb; //in Angstrom
    f_nu[i]   = yb*A_DM/A_10pc; //in erg/s/Hz/cm^2
    f_nu[i]  /= uJy;  //convert to uJy

    //f_nu is now in uJy
  }

  //close the SED model file
  fclose(fp);
}

#ifdef ALTERNATE_LINE_EMISSION
void free_line_strength_buffers(void)
{
  free(lambda_lines);
  free(f_lines);
}

void rescale_line_strengths(double z, double lNlyc)
{
  int i, nll;
  double xb, yb;
  FILE *fp;
  double A_10pc;
  double A_DM = distance_modulus(z); //look up the distance modulus
  double pc_in_cm = 3.085678e18;     //pc in cm
  double uJy = 1.0e-29;              // microJansky in erg/s/Hz/cm^2
  double c = 2.99792458e18;          //speed of light in Angstroms / sec
  char fname[200];
  double A_z3 = 4.0*M_PI*pow(8.158102e+28,2); 

  //file containing list of emission lines and strengths
  sprintf(fname,"line_strengths.txt");

  //area of sphere of 10pc radius in cm^2
  A_10pc = 4.0*M_PI*pow(10*pc_in_cm,2); 

  printf("A_DM %e A_10pc %e\n",A_DM,A_10pc);

  printf("Loading file = %s lNlyc = %e.\n",fname,lNlyc);
  if(!(fp = fopen(fname,"r")))
  {
    printf("Error opening %s.\n",fname);
    exit(0);
  }
  fscanf(fp,"%d\n",&nll);
  printf("nll = %d fn %s\n",nll,fname);

  //allocate line buffers
  n_lines = nll;
  lambda_lines = (double *) malloc(n_lines*sizeof(double));
  f_lines   = (double *) malloc(n_lines*sizeof(double));

  for(i=0;i<n_lines;i++)
  {
    fscanf(fp,"%lf\t%lf\n",&xb,&yb);
    lambda_lines[i] = xb * 1.0e4; //in Ang in the rest frame

    //record the line flux
    f_lines[i]   = yb;              //input line strength in erg/s/NLyc
    f_lines[i]  *= pow(10.,lNlyc);  //rescale by Lyman continuum, to erg/s
    f_lines[i]  *= A_DM/A_10pc;     //in erg/s/cm^2
    f_lines[i]  /= (1+z);           //rescale as flux, not flux density
    f_lines[i]  /= uJy;             //convert to uJy * Hz in rest frame

    if((i==16)||(i==18)||(i==19))
      printf("i %d lambda %e yb %e lNlyc %e A_DM %e A_10pc %e (1+z) %e A_comb %e A_z3 %e uJy %e f_lines %e\n",i,lambda_lines[i],yb,lNlyc,A_DM,A_10pc,(1+z),A_10pc*(1+z)/A_DM,A_z3,uJy,f_lines[i]);

    //f_lines is now in uJy * Hz
    //printf("lambda_lines[%d] %e f_nu_lines %e [erg/s/cm^2]\n",i,lambda_lines[i],f_nu_lines[i]*uJy);
  }


  fclose(fp);
  //exit(0);
}
#endif //ALTERNATE_LINE_EMISSION

//now we load the line model
void load_line_model(char fname[], double z, double lNlyc)
{
  int i;
  double xb, yb;
  FILE *fp;
  double A_10pc;
  double A_DM = distance_modulus(z); //look up the distance modulus
  double pc_in_cm = 3.085678e18;     //pc in cm
  double uJy = 1.0e-29;              // microJansky in erg/s/Hz/cm^2
  double c = 2.99792458e18;          //speed of light in Angstroms / sec
  double y_Hbeta;
  double nu_hi, nu_lo;
  double y_hi, y_lo;
  int nll;

  //area of sphere of 10pc radius in cm^2
  A_10pc = 4.0*M_PI*pow(10*pc_in_cm,2);

  //The distance modulus to z=3.10 = 4.563294e+01

  printf("Loading file = %s lNlyc = %e.\n",fname,lNlyc);
  if(!(fp = fopen(fname,"r")))
  {
    printf("Error opening %s.\n",fname);
    exit(0);
  }
  fscanf(fp,"%d\n",&nll);
  printf("nll = %d fn %s\n",nll,fname);
  nu_lo = 0;
  y_lo = 0;
  y_Hbeta = 0;
  for(i=0;i<n_lambda_line;i++)
  {
    fscanf(fp,"%lf\t%lf\n",&xb,&yb);
    lambda_line[i] = xb * 1.0e4; //in Ang

    nu_hi = nu_lo;
    nu_lo = c/lambda_line[i]; //Hz
    y_hi = y_lo;
    y_lo = yb;

    //integrate the Hbeta luminosity as a check
    if( (lambda_line[i] > 4850)&(lambda_line[i] <4900) )
      y_Hbeta += 0.5*(y_hi + y_lo)*(nu_hi-nu_lo);

    //record the line flux density
    f_nu_line[i]   = yb*A_DM/A_10pc;  //in erg/s/Hz/cm^2
    f_nu_line[i]  /= uJy;             //convert to uJy
    f_nu_line[i]  *= pow(10.,lNlyc);  //scale by Lyman continuum

    //f_nu is now in uJy
  }

  //print a check on the Hbeta luminosity
  //in the line emission profile
  printf("y_Hbeta = %e A_DM %e A_10pc %e Nlyc %e\n",y_Hbeta,A_DM,A_10pc,lNlyc);

  //close the file
  fclose(fp);
}

//now we load the cont model
void load_cont_model(char fname[], double z, double lNlyc)
{
  int i;
  double xb, yb;
  FILE *fp;
  double A_10pc;
  double A_DM = distance_modulus(z); //look up the distance modulus
  double pc_in_cm = 3.085678e18;     //pc in cm
  double uJy = 1.0e-29;              // microJansky in erg/s/Hz/cm^2

  //area of sphere of 10pc radius in cm^2
  A_10pc = 4.0*M_PI*pow(10*pc_in_cm,2);


  //Load the continuum model
  printf("Loading file = %s.\n",fname);
  if(!(fp = fopen(fname,"r")))
  {
    printf("Error opening %s.\n",fname);
    exit(0);
  }
  fscanf(fp,"%d\n",&n_lambda_cont);
  for(i=0;i<n_lambda_cont;i++)
  {
    fscanf(fp,"%lf\t%lf\n",&xb,&yb);
    lambda_cont[i] = xb * 1.0e4;      //in Ang
    f_nu_cont[i]   = yb*A_DM/A_10pc;  //in erg/s/Hz/cm^2
    f_nu_cont[i]  /= uJy;             //convert to uJy
    f_nu_cont[i]  *= pow(10.,lNlyc);  //rescale by Lyman continuum
    //f_nu is now in uJy
  }

  //close the template file
  fclose(fp);
}

//load the sed samples and attuate them
//owing to the IGM absorption
void load_and_attenuate_sed_samples(void)
{
  //indices
  int i,j,k,l;

  //flag for switching binary fitting
  //flag_binary = 0 -- only fitting single
  //flag_binary = 1 -- only fitting binary
  //flag_binary = 2 -- fitting range of binarity
  flag_binary = 1;
  //flag_binary = 0;


  //do we need to load the IGM absorption?
  int flag_igm_abs = 1;

  //filename properties
  char fname_single[200];
  char fname_binary[200];
  char fname_ion[200];
  char fname_base[200];
  char fname_sfr[200];
  char fname_suff[200];
  char fname[200];
  char fname_line[200];
  char fname_cont[200];
  double delta_lam[13] = {271.340266941,453.875506634,479.114240852,44.5821543373,488.294997293,574.200210656,753.216532216,621.36788108,809.5358419
,1372.04049973,1796.04296323,3.77550926068e3,4.7976597986e3};


  char **fname_z_met;

  //load the metallicity file name components
  fname_z_met = load_z_met_filenames();

  //the directory containing the single
  //star SED models
  sprintf(fname_single,"single");

  //the directory containing the binary
  //star SED models
  if(flag_binary==0)
  {
    sprintf(fname_binary,"single"); 
  }else if(flag_binary==1){
    sprintf(fname_binary,"binary"); 
  }

  //the base filename
  sprintf(fname_base,"fnu_cgs");

  //the base sfr extension
  sprintf(fname_sfr,"1Msun_per_year");

  //the base file suffix
  sprintf(fname_suff,"txt");

  //allocate the sed buffers
  allocate_sed_buffers();


  //Let's do Lyman continuum first
  sprintf(fname_ion,"ionizing");

  //Load the Lyman continuum flux vs. age
  //of stellar population
  for(j=0;j<n_z_met;j++)
  {
    sprintf(fname,"%s/%s.%s.imf135_100.%s.age_%s.%s",fname_binary,fname_ion,fname_binary,fname_sfr,fname_z_met[j],fname_suff);
    printf("fname %s\n",fname);
    FILE *fp = fopen(fname,"r");
    double ab, nlb, cb, db, eb;
    for(k=0;k<n_age;k++)
    {
      fscanf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\n",&ab,&nlb,&cb,&db,&eb);
      N_Lyc[j][k] = nlb;
    }
    //close the file
    fclose(fp);
  }

  //begin the loop over binarity
  for(i=0;i<n_f_bin;i++)
  {
    //begin the loop over metallicity
    for(j=0;j<n_z_met;j++)
    {
      //begin the loop over age
      for(k=0;k<n_age;k++)
      {
        //make the filename
        if(flag_binary==0)
        {
          //fitting only single
          sprintf(fname,"%s/%s.%s.%s.age_%3.1f.%s.%s",fname_single,fname_base,fname_single,fname_sfr,age[k],fname_z_met[j],fname_suff);

        }else if(flag_binary==1){

          //fitting only binary
          sprintf(fname,"%s/%s.%s.%s.age_%3.1f.%s.%s",fname_binary,fname_base,fname_binary,fname_sfr,age[k],fname_z_met[j],fname_suff);

        }else if(flag_binary==2){


          //fitting both binary fractions
          if(i==0)
          {
            //single
            sprintf(fname,"%s/%s.%s.%s.age_%3.1f.%s.%s",fname_single,fname_base,fname_single,fname_sfr,age[k],fname_z_met[j],fname_suff);
          }else{
            //binary
            sprintf(fname,"%s/%s.%s.%s.age_%3.1f.%s.%s",fname_binary,fname_base,fname_binary,fname_sfr,age[k],fname_z_met[j],fname_suff);
          }
        }


        //nebular continuum and line templates
        sprintf(fname_cont,"nebular_continuum_emission.HI_and_HeI.txt");
#ifndef ALTERNATE_LINE_EMISSION
        sprintf(fname_line,"line_emission.%s.txt",fname_z_met[j]);
#endif  //ALTERNATE_LINE_EMISSION

        //at this point we should have the filename correct
        printf("i %d j %d k %d fname %s\n",i,j,k,fname);
#ifndef ALTERNATE_LINE_EMISSION
        printf("i %d j %d k %d fname_line %s\n",i,j,k,fname_line);
#endif //ALTERNATE_LINE_EMISSION
        printf("i %d j %d k %d fname_cont %s\n",i,j,k,fname_cont);

        //now we load the SED model, line model, and continuum model
        load_sed_model(fname, z);
#ifndef ALTERNATE_LINE_EMISSION
        load_line_model(fname_line, z, N_Lyc[j][k]);
#else  //ALTERNATE_LINE_EMISSION
        rescale_line_strengths(z, N_Lyc[j][k]);
#endif //ALTERNATE_LINE_EMISSION
        load_cont_model(fname_cont, z, N_Lyc[j][k]);



        //first we attenuate the line and continuum models
        printf("attenuating\n");
#ifndef ALTERNATE_LINE_EMISSION
        attenuate_line_model(&flag_igm_abs);

        //note the lines are at >3000 angstroms, so this isn't necessary

#endif //ALTERNATE_LINE_EMISSION
        attenuate_cont_model(&flag_igm_abs);

        //now we attenuate the SED model
        attenuate_sed_model(&flag_igm_abs);

        //now we sample the SED model
        //sample_sed_model(i,j,k);
        printf("sampling\n");

        //now we sample the SED model, averaging over the pass bands
        pass_band_average_sed_model(i,j,k);        
        pass_band_average_line_model(i,j,k);        
        pass_band_average_cont_model(i,j,k);        

        for(l=0;l<n_data;l++)
        {
          double EW = pow(data_x[l]*(1+z),2)*(line_model[i][j][k][l]/(sed_model[i][j][k][l]+cont_model[i][j][k][l]))/(1.+z);
          printf("i %d j %d k %d l %d lambda_obs %e data_y %e model_sed %e model_line %e model_cont %e\n",i,j,k,l,data_x[l]*(1+z),data_y[l],sed_model[i][j][k][l],line_model[i][j][k][l],cont_model[i][j][k][l]);
          printf("i %d j %d k %d l %d delta_lam %e EW %e delta_mag %e\n",i,j,k,l,delta_lam[l],EW,-2.5*log10(1.+EW*(1+z)/delta_lam[l]));
          printf("i %d j %d k %d l %d OIII %e %e [erg/s/cm^2, uJy*Hz]\n",i,j,k,l,(f_lines[18]+f_lines[19])*1e-29,f_lines[18]+f_lines[19]);
          printf("*****\n");
        }
#ifdef ALTERNATE_LINE_EMISSION
        free_line_strength_buffers();
#endif //ALTERNATE_LINE_EMISSION

      }
    }
  }

  //now our grid sed_model has the 
  //correct SED values to compare
  //against the data

  //free the z_met filenames
  free_z_met_filenames(fname_z_met);

  //deallocate the sed buffer memory
  deallocate_sed_buffers();

  //exit(0);

}
///now we attenuate the cont model
//owing to IGM absorption along the
//line of sight
void attenuate_cont_model(int *flag_igm_abs)
{
  int i;
  //double *absorption;

  if(*flag_igm_abs)
  {
    //we need to load the IGM absorption
#ifndef IGM_ABS
    absorption = madau_absorption(lambda_cont, n_lambda_cont, z);
#else //IGM_ABS
    absorption = igm_absorption(lambda_cont, n_lambda_cont, z);
#endif //IGM_ABS

    //print the absorption to a file
#ifndef IGM_ABS
    print_madau_absorption(lambda_cont, n_lambda_cont, z);
#else //IGM_ABS
    print_igm_absorption(lambda_cont, n_lambda_cont, z);
#endif //IGM_ABS

    //remember that we've loaded the attenuation model
    *flag_igm_abs = 0;
  }

  for(int i=0;i<n_lambda_cont;i++)
  {
    //adjust f_nu according to absorption
    f_nu_cont[i] *= absorption[i];
  }
  //free the absorption array
  free(absorption);
  *flag_igm_abs = 1;
}
///now we attenuate the line model
//owing to IGM absorption along the
//line of sight
void attenuate_line_model(int *flag_igm_abs)
{
  int i;
  //double *absorption;

  if(*flag_igm_abs)
  {
    //we need to load the IGM absorption
#ifndef IGM_ABS
    absorption = madau_absorption(lambda_line, n_lambda_line, z);
#else //IGM_ABS
    absorption = igm_absorption(lambda_line, n_lambda_line, z);
#endif //IGM_ABS

    //print the absorption to a file
#ifndef IGM_ABS
    print_madau_absorption(lambda_line, n_lambda_line, z);
#else //IGM_ABS
    print_igm_absorption(lambda_line, n_lambda_line, z);
#endif //IGM_ABS

    //remember that we've loaded the attenuation model
    *flag_igm_abs = 0;
  }

  for(int i=0;i<n_lambda_line;i++)
  {
    //adjust f_nu according to absorption
    f_nu_line[i] *= absorption[i];
  }
  //free the absorption array
  free(absorption);
  *flag_igm_abs = 1;
}

//now we attenuate the SED model
//owing to IGM absorption along the
//line of sight
void attenuate_sed_model(int *flag_igm_abs)
{
  int i;
  //double *absorption;

  if(*flag_igm_abs)
  {
    //we need to load the IGM absorption
#ifndef IGM_ABS
    absorption = madau_absorption(lambda, n_lambda, z);
#else //IGM_ABS
    absorption = igm_absorption(lambda, n_lambda, z);
#endif //IGM_ABS

    //print the absorption to a file
#ifndef IGM_ABS
    print_madau_absorption(lambda, n_lambda, z);
#else //IGM_ABS
    print_igm_absorption(lambda, n_lambda, z);
#endif //IGM_ABS

    //remember that we've loaded the attenuation model
    *flag_igm_abs = 0;
  }

  printf("rescaling\n");

  for(int i=0;i<n_lambda;i++)
  {
    //adjust f_nu according to absorption
    f_nu[i] *= absorption[i];
  }

  printf("done.\n");

  //free the absorption array
  //free(absorption);

}

//now we pass band average the SED model
void pass_band_average_sed_model(int i, int j, int k)
{
  int l, m;
  double f_nu_lo,f_nu_hi;
  double fT, T, y_T;
  int i_lo, i_hi;
  double dx_lo, dx_hi;
  double nu_lo, nu_hi;
  double nu_ave;

  //perform a check
  double c = 2.99792458e18; //speed of light in Ang/sec

  double oopz = 1./(1.+z);
#ifndef PBINTERP
  double p_lo, p_hi;
  double *x_pb_o;
#endif //PBINTERP

  //perform a linear interpolation at the values 
  //of the observed SED photometry

  //find the location of each sample in wavelength, and then
  //record the linear interpolation of this SED model at that
  //wavelength
  for(l=0;l<n_data;l++)
  {

    //note that each data point corresponds to a pass band
    //so we load the pass band

    LoadPassBand(l);


#ifndef PBINTERP
    x_pb_o = (double *) malloc(n_pb*sizeof(double));
    for(m=0;m<n_pb;m++)
      x_pb_o[m] = x_pb[m]*oopz;
#endif //PBINTERP

    //now we just need to average the SED over the pass band

#ifdef PBINTERP
    //loop over band curve points and integrate
    fT = 0;
    T  = 0;
    for(m=0;m<n_pb-1;m++)
    {
      //index where x[i_lo] <= x < x[i_lo+1]
      i_lo = gsl_interp_bsearch(lambda,x_pb[m]*oopz,0,n_lambda);
      i_hi = gsl_interp_bsearch(lambda,x_pb[m+1]*oopz,0,n_lambda);

      //fractional distance between model points in wavelength
      dx_lo = (x_pb[m]*oopz-lambda[i_lo])/(lambda[i_lo+1]-lambda[i_lo]);
      dx_hi = (x_pb[m+1]*oopz-lambda[i_hi])/(lambda[i_hi+1]-lambda[i_hi]);
      f_nu_lo  = (1.0-dx_lo)*f_nu[i_lo] + dx_lo*f_nu[i_lo+1];
      f_nu_hi  = (1.0-dx_hi)*f_nu[i_hi] + dx_hi*f_nu[i_hi+1];

      //\int f_nu * T dnu
      fT += 0.5*(f_nu_lo*y_pb[m] + f_nu_hi*y_pb[m+1])*(x_pb[m+1]-x_pb[m]);
      T  += 0.5*(y_pb[m] + y_pb[m+1])*(x_pb[m+1]-x_pb[m]);
    }
#else //PBINTERP

    //loop over band curve points and integrate
    fT = 0;
    T  = 0;
    for(m=0;m<n_lambda-1;m++)
    {
      //restrict to portion of the SED model covered by the
      //pass band -- this could be sped up by finding these
      //limits before entering the loop
      if( (lambda[m]>=x_pb_o[0])&&(lambda[m+1]<=x_pb_o[n_pb-1]))
      {
        //index where x[i_lo] <= x < x[i_lo+1]
        i_lo = gsl_interp_bsearch(x_pb_o,lambda[m],0,n_pb);
        i_hi = gsl_interp_bsearch(x_pb_o,lambda[m+1],0,n_pb);

        //fractional distance between model points in wavelength
        dx_lo = (lambda[m]-x_pb_o[i_lo])/(x_pb_o[i_lo+1]-x_pb_o[i_lo]);
        dx_hi = (lambda[m+1]-x_pb_o[i_hi])/(x_pb_o[i_hi+1]-x_pb_o[i_hi]);

        //pass band at the high and low frequency
        p_lo  = (1.0-dx_lo)*y_pb[i_lo] + dx_lo*y_pb[i_lo+1];
        p_hi  = (1.0-dx_hi)*y_pb[i_hi] + dx_hi*y_pb[i_hi+1];



        //\int f_nu * T dnu
        nu_lo = c/lambda[m+1];
        nu_hi = c/lambda[m];
        nu_ave = 0.5*(nu_hi + nu_lo);
        fT += 0.5*(p_lo*f_nu[m] + p_hi*f_nu[m+1])*(nu_hi-nu_lo)/nu_ave;
        T  += 0.5*(p_hi + p_lo)*(nu_hi-nu_lo)/nu_ave;      
      }
    }
#endif //PBINTERP

    //record linearly interpolated model value
    y_T = fT/T;
    sed_model[i][j][k][l] = y_T;

    //free the memory holding the pass band data
    FreePassBand();
#ifndef PBINTERP
    free(x_pb_o);
#endif //PBINTERP
  }
}


//now we pass band average the line model
void pass_band_average_line_model(int i, int j, int k)
{
  int l, m;
  double f_nu_lo,f_nu_hi;
  double fT, T, y_T;
  int i_lo, i_hi;
  double dx_lo, dx_hi;
  double nu_lo, nu_hi;
  double nu_ave;
  double y_Hbeta;

  double oopz = 1./(1.+z);  //1/(1+z)
  double c = 2.99792458e18; //speed of light in Ang/sec

#ifdef ALTERNATE_LINE_EMISSION
  double y_pb_max = 0;      //maximum transmissivity of pass band
  double f_pb_define = 0.5; //fraction of y_pb_max that defines width of pass band
  double l_pb_lo;    //minimum wavelength of pass band
  double l_pb_hi;    //maximum wavelength of pass band
  double l_pb_max = 0;
  //double dnu_pb;      //width of pass band in Hz
  double dlam_lo;
  double dlam_hi;
  double dlam_lo_sum;
  double dlam_hi_sum;
  double l_pb_ave;
  double T_pb_ave;
#endif //ALTERNATE_LINE_EMISSION

#ifndef PBINTERP
  double p_lo, p_hi;
  double *x_pb_o;
#endif //PBINTERP

  //perform a linear interpolation at the values 
  //of the observed SED photometry

  //find the location of each sample in wavelength, and then
  //record the linear interpolation of this SED model at that
  //wavelength
  for(l=0;l<n_data;l++)
  {
#ifdef ALTERNATE_LINE_EMISSION
    y_pb_max = 0;
    dlam_lo = 0;
    dlam_hi = 0;
    dlam_lo_sum = 0;
    dlam_hi_sum = 0;
    l_pb_ave = 0.0;
    T_pb_ave = 0.0;
#endif // ALTERNATE_LINE_EMISSION

    //note that each data point corresponds to a pass band
    //so we load the pass band

    LoadPassBand(l);

#ifndef PBINTERP
    x_pb_o = (double *) malloc(n_pb*sizeof(double));
    for(m=0;m<n_pb;m++)
    {
      x_pb_o[m] = x_pb[m]*oopz;
#ifdef ALTERNATE_LINE_EMISSION
      //find the maximum passband transmissivity
      //and pass band maximum in observed frame
      T_pb_ave += y_pb[m];
      l_pb_ave += y_pb[m]*x_pb[m];
      if(y_pb[m]>y_pb_max)
      {
        l_pb_max = x_pb[m];
        y_pb_max = y_pb[m];
      }
#endif //ALTERNATE_LINE_EMISSION
    }
#endif //PBINTERP


#ifdef ALTERNATE_LINE_EMISSION

    //define the bandwidth of the passband in the observed frame
    l_pb_lo = -1;
    l_pb_hi = -1;
    l_pb_ave /= T_pb_ave; //average wavelength in Angstrom in observed frame
    for(m=0;m<n_pb;m++)
    {
      if(x_pb[m]<l_pb_max)
      {
        dlam_lo += y_pb[m]*(l_pb_max-x_pb[m]);
        dlam_lo_sum += y_pb[m];
      }
      if(x_pb[m]>=l_pb_max)
      {
        dlam_hi += y_pb[m]*(x_pb[m]-l_pb_max);
        dlam_hi_sum += y_pb[m];
      }
      //find short wavelength side of band
      if( (y_pb[m]>f_pb_define*y_pb_max)&&(l_pb_lo<0) )
      {
        l_pb_lo = x_pb[m];
      }
      //find long wavelength side of band
      if( (y_pb[m]<f_pb_define*y_pb_max)&&(l_pb_hi<0)&&(l_pb_lo>0)&&(x_pb[m]>l_pb_max) )
      {
        l_pb_hi = x_pb[m];
      }
      //if(l==8)
        //printf("BANDWIDTH l %e y %e y_pb_max %e\n",x_pb_o[m],y_pb[m],y_pb_max);
    }
    //compute the pass band portions in observed Ang
    dlam_lo /= dlam_lo_sum;
    dlam_hi /= dlam_hi_sum;


    //printf("BANDWIDTH for %d l_ave %f; width %f %f l_pb_lo %e l_pb_max %e l_pb_hi %e y_pb_max %e dn_pb %e check %e %e %e\n",l,l_pb_ave,(dnu_hi + dnu_lo)/oopz,(l_pb_hi - l_pb_lo)/oopz,l_pb_lo/oopz,l_pb_max/oopz,l_pb_hi/oopz,y_pb_max,dnu_pb,(l_pb_max-dnu_lo)/oopz,(dnu_hi+l_pb_max)/oopz,c*oopz/(l_pb_max-dnu_lo) + c*oopz/(l_pb_max+dnu_hi));
    //printf("BANDWIDTH for %d l_ave %f; width %f %f l_pb_lo %e l_pb_max %e l_pb_hi %e y_pb_max %e\n",l,l_pb_ave,(dlam_hi + dlam_lo)/oopz,(l_pb_hi - l_pb_lo)/oopz,l_pb_lo/oopz,l_pb_max/oopz,l_pb_hi/oopz,y_pb_max);

#endif // ALTERNATE_LINE_EMISSION

    //exit(0);
    //now we just need to average the SED over the pass band

#ifndef  ALTERNATE_LINE_EMISSION


#ifdef PBINTERP
    //loop over band curve points and integrate
    fT = 0;
    T  = 0;
    for(m=0;m<n_pb-1;m++)
    {
      //index where x[i_lo] <= x < x[i_lo+1]
      i_lo = gsl_interp_bsearch(lambda_line,x_pb[m]*oopz,0,n_lambda_line);
      i_hi = gsl_interp_bsearch(lambda_line,x_pb[m+1]*oopz,0,n_lambda_line);

      //fractional distance between model points in wavelength
      dx_lo = (x_pb[m]*oopz-lambda_line[i_lo])/(lambda_line[i_lo+1]-lambda_line[i_lo]);
      dx_hi = (x_pb[m+1]*oopz-lambda_line[i_hi])/(lambda_line[i_hi+1]-lambda_line[i_hi]);
      f_nu_lo  = (1.0-dx_lo)*f_nu_line[i_lo] + dx_lo*f_nu_line[i_lo+1];
      f_nu_hi  = (1.0-dx_hi)*f_nu_line[i_hi] + dx_hi*f_nu_line[i_hi+1];

      //\int f_nu * T dnu
      fT += 0.5*(f_nu_lo*y_pb[m] + f_nu_hi*y_pb[m+1])*(x_pb[m+1]-x_pb[m]);
      T  += 0.5*(y_pb[m] + y_pb[m+1])*(x_pb[m+1]-x_pb[m]);
    }
#else //PBINTERP

    //loop over band curve points and integrate
    fT = 0;
    T  = 0;

    y_Hbeta = 0;
    for(m=0;m<n_lambda_line-1;m++)
    {
      //restrict to wavelengths covered by the passband
      if( (lambda_line[m]>=x_pb_o[0])&&(lambda_line[m+1]<=x_pb_o[n_pb-1]))
      {
        //index where x[i_lo] <= x < x[i_lo+1]
        i_lo = gsl_interp_bsearch(x_pb_o,lambda_line[m],0,n_pb);
        i_hi = gsl_interp_bsearch(x_pb_o,lambda_line[m+1],0,n_pb);

        //fractional distance between model points in wavelength
        dx_lo = (lambda_line[m]-x_pb_o[i_lo])/(x_pb_o[i_lo+1]-x_pb_o[i_lo]);
        dx_hi = (lambda_line[m+1]-x_pb_o[i_hi])/(x_pb_o[i_hi+1]-x_pb_o[i_hi]);

        //passband transmission curve sampled at high and low frequency
        p_lo  = (1.0-dx_lo)*y_pb[i_lo] + dx_lo*y_pb[i_lo+1];
        p_hi  = (1.0-dx_hi)*y_pb[i_hi] + dx_hi*y_pb[i_hi+1];

        //\int f_nu * T dnu
        nu_lo = c/lambda_line[m+1];
        nu_hi = c/lambda_line[m];
        nu_ave = 0.5*(nu_hi + nu_lo);
        fT += 0.5*(p_lo*f_nu_line[m] + p_hi*f_nu_line[m+1])*(nu_hi-nu_lo)/nu_ave;
        T  += 0.5*(p_hi + p_lo)*(nu_hi-nu_lo)/nu_ave;

        //if( (lambda_line[m]>4850)&(lambda_line[m]<4900) )
          //y_Hbeta += 0.5*(p_lo*f_nu_line[m] + p_hi*f_nu_line[m+1])*(nu_hi-nu_lo)/nu_ave;

          //y_Hbeta += 0.5*(p_lo + p_hi)*(nu_hi-nu_lo)/nu_ave;

          //y_Hbeta += 0.5*(f_nu_line[m] + f_nu_line[m+1])*(nu_hi-nu_lo);
          //y_Hbeta += 0.5*(p_lo + p_hi)*(nu_hi-nu_lo)/nu_ave;
        //if( (lambda_line[i] > 4850)&(lambda_line[i] <4900) )
          //y_Hbeta += 0.5*(y_hi + y_lo)*(nu_hi-nu_lo);

      }
    }

#endif //PBINTERP
    

    //record linearly interpolated model value
    y_T = fT/T;
    line_model[i][j][k][l] = y_T;

#else //ALTERNATE_LINE_EMISSION

    //find integral over pass band in observed frame
    T  = 0;
    double T_norm = 0;
    for(m=0;m<n_pb-1;m++)
    {
      nu_lo = c/x_pb[m+1];
      nu_hi = c/x_pb[m];
      nu_ave = 0.5*(nu_hi + nu_lo);

      //only worry about the transmissivity above 0.001
      //when estimating the average transmissivity
      if(y_pb[m]>1.0e-3)
      {
        T  += 0.5*(y_pb[m+1] + y_pb[m])*(nu_hi-nu_lo)/nu_ave;
        T_norm  += (nu_hi-nu_lo)/nu_ave;
      }

    }
    //printf("T = %e\n",T);
    //exit(0);

    //x_pb_o contains the wavelengths
    //of the passband reduced by 1./(1+z)
    //we can check and see whether each line contributes
    //to this passband.  First, we need to find the maximum
    //of the pass band

    //first compute f_lambda HERE
    double dlambda_pb;
    int il;

    //bandwidth based on transmission-averaged wavelengths
    dlambda_pb = (dlam_hi + dlam_lo); //band width in Ang, in observed frame


    //dlambda_pb = (dlam_hi + dlam_lo); //band width in Ang, in restframe

    //bandwidth based on a fractional power
    //dlambda_pb = (l_pb_hi - l_pb_lo); //band width in Ang, in restframe
    //dlambda_pb = (l_pb_hi - l_pb_lo)*(1+z)*(1+z); //band width in Ang, in restframe

    line_model[i][j][k][l] = 0;
    double lambda_line_obs;
    double fac_pb;
    for(il=0;il<n_lines;il++)
    {
      //find observed wavelength of line
      lambda_line_obs = lambda_lines[il] * (1+z);

      //check whether line falls in this pass band
      if((lambda_line_obs>=x_pb[0])&&(lambda_line_obs<=x_pb[n_pb-1]) )
      {

        //this line falls in the pass band

        //find f_lambda in uJy * Hz/Ang
        y_T = f_lines[il] / dlambda_pb;

        //moderate by the pass band
        //nu_lo = c/x_pb[i_lo+1];
        //nu_hi = c/x_pb[i_lo];
        //nu_ave = 0.5*(nu_hi + nu_lo);
        //fac_pb *= (nu_hi-nu_lo)/nu_ave;
        //fac_pb /= T;
        //y_T *= fac_pb;

        //at this point y_T is in dF/dlambda x T(lambda)
        //we need to convert to a flux density
        // dF/dnu = dF/dlambda * dlambda/dnu
        // lambda = c/nu dlambda/dnu = c/nu^2 = lambda^2 / c HERE

        y_T *= pow(lambda_line_obs,2) / c; // to uJy in observed frame

        //We have to compute the average transmission
        //and then take the ratio of the transmission at the
        //line to the average transmission (since we treat the line flux as an average over the pass band)


        //transmission at the line
        i_lo = gsl_interp_bsearch(x_pb,lambda_line_obs,0,n_pb);
        fac_pb = y_pb[i_lo] + (lambda_line_obs - x_pb[i_lo])*(y_pb[i_lo+1] - y_pb[i_lo])/(x_pb[i_lo+1] - x_pb[i_lo]);

        //take ratio of transmissivity at line to the average transmissivity
        fac_pb /= (T/T_norm);

        printf("BANDWIDTH lambda_lines[%d] %e dl %e f_lines %e fac_pb %e T/T_norm %e y_T %e\n",il,lambda_lines[il],dlambda_pb,f_lines[il],fac_pb,T/T_norm,y_T);

        line_model[i][j][k][l] += y_T * fac_pb;

        //#error check factors of 1+z
      }
  	}
  	//line_model[i][j][k][l] /= T; //divide by transmission integral, leave in uJy in observed frame

#endif //ALTERNATE_LINE_EMISSION

    //free the memory holding the pass band data
    FreePassBand();
#ifndef PBINTERP
    free(x_pb_o);
#endif //PBINTERP
  }

  //exit(0);
}


void pass_band_average_cont_model(int i, int j, int k)
{
  int l, m;
  double f_nu_lo,f_nu_hi;
  double fT, T, y_T;
  int i_lo, i_hi;
  double dx_lo, dx_hi;
  double nu_lo, nu_hi;
  double c = 2.99792458e18; //speed of light in Ang/sec

  double oopz = 1./(1.+z);  //1/(1+z) of the model

#ifndef PBINTERP
  double nu_ave;
  double p_lo, p_hi;
  double *x_pb_o;
#endif //PBINTERP


  //perform a linear interpolation at the values 
  //of the observed SED photometry

  //find the location of each sample in wavelength, and then
  //record the linear interpolation of this SED model at that
  //wavelength
  for(l=0;l<n_data;l++)
  {

    //note that each data point corresponds to a pass band
    //so we load the pass band

    LoadPassBand(l);

#ifndef PBINTERP
    x_pb_o = (double *) malloc(n_pb*sizeof(double));
    for(m=0;m<n_pb;m++)
      x_pb_o[m] = x_pb[m]*oopz;
#endif //PBINTERP

    //now we just need to average the SED over the pass band

#ifdef PBINTERP

    //loop over band curve points and integrate
    fT = 0;
    T  = 0;
    for(m=0;m<n_pb-1;m++)
    {
      //index where x[i_lo] <= x < x[i_lo+1]
      i_lo = gsl_interp_bsearch(lambda_cont,x_pb[m]*oopz,0,n_lambda_cont);
      i_hi = gsl_interp_bsearch(lambda_cont,x_pb[m+1]*oopz,0,n_lambda_cont);

      //fractional distance between model points in wavelength
      dx_lo = (x_pb[m]*oopz-lambda_cont[i_lo])/(lambda_cont[i_lo+1]-lambda_cont[i_lo]);
      dx_hi = (x_pb[m+1]*oopz-lambda_cont[i_hi])/(lambda_cont[i_hi+1]-lambda_cont[i_hi]);
      f_nu_lo  = (1.0-dx_lo)*f_nu_cont[i_lo] + dx_lo*f_nu_cont[i_lo+1];
      f_nu_hi  = (1.0-dx_hi)*f_nu_cont[i_hi] + dx_hi*f_nu_cont[i_hi+1];

      //\int f_nu * T dnu
      fT += 0.5*(f_nu_lo*y_pb[m] + f_nu_hi*y_pb[m+1])*(x_pb[m+1]-x_pb[m])*oopz;
      T  += 0.5*(y_pb[m] + y_pb[m+1])*(x_pb[m+1]-x_pb[m])*oopz;
    }
#else //PBINTERP

    //loop over band curve points and integrate
    fT = 0;
    T  = 0;
    for(m=0;m<n_lambda_cont-1;m++)
    {
      if( (lambda_cont[m]>=x_pb_o[0])&&(lambda_cont[m+1]<=x_pb_o[n_pb-1]))
      {
        //index where x[i_lo] <= x < x[i_lo+1]
        i_lo = gsl_interp_bsearch(x_pb_o,lambda_cont[m],0,n_pb);
        i_hi = gsl_interp_bsearch(x_pb_o,lambda_cont[m+1],0,n_pb);

        //fractional distance between model points in wavelength
        dx_lo = (lambda_cont[m]-x_pb_o[i_lo])/(x_pb_o[i_lo+1]-x_pb_o[i_lo]);
        dx_hi = (lambda_cont[m+1]-x_pb_o[i_hi])/(x_pb_o[i_hi+1]-x_pb_o[i_hi]);

        //pass band sampled at low and high frequencies
        p_lo  = (1.0-dx_lo)*y_pb[i_lo] + dx_lo*y_pb[i_lo+1];
        p_hi  = (1.0-dx_hi)*y_pb[i_hi] + dx_hi*y_pb[i_hi+1];

        //\int f_nu * T dnu
        nu_hi = c/lambda_cont[m];
        nu_lo = c/lambda_cont[m+1];
        nu_ave = 0.5*(nu_hi + nu_lo);

        //\int f_nu * T dnu/nu   
        fT += 0.5*(p_lo*f_nu_cont[m] + p_hi*f_nu_cont[m+1])*(nu_hi-nu_lo)/nu_ave;
        T  += 0.5*(p_hi + p_lo)*(nu_hi-nu_lo)/nu_ave;    


        //just double check that the passbands aren't interpolated
        //incorrectly in a way that leads to negative transmissivity
        if(p_hi<0 || p_lo<0)
        {
          printf("dx_lo %e %e %e\n",dx_lo,(lambda_cont[m]-x_pb_o[i_lo]),(x_pb_o[i_lo+1]-x_pb_o[i_lo]));
          printf("dx_hi %e %e %e\n",dx_hi,(lambda_cont[m+1]-x_pb_o[i_hi]),(x_pb_o[i_hi+1]-x_pb_o[i_hi]));
          printf("ERROR l %d m %d m+1 %d l[m] %e l[m+1] %e nl %d i_lo %d i_hi %d x_pb_o[i_lo] %e x_pb_o[i_lo+1] %e x_pb_o[i_hi] %e x_pb_o[i_hi+1] %e n_pb %d y_pb_ilo %e i_pb_ihi %e p_lo %e p_hi %e dx_lo %e dx_hi %e\n",l,m,m+1,lambda_cont[m],lambda_cont[m+1],n_lambda_cont,i_lo,i_hi,x_pb_o[i_lo],x_pb_o[i_lo+1],x_pb_o[i_hi],x_pb_o[i_hi+1],n_pb,y_pb[i_lo],y_pb[i_hi],p_lo,p_hi,dx_lo,dx_hi);
          exit(-1);
        }
      }
    }
#endif //PBINTERP

    //record linearly interpolated model value
    y_T = fT/T;

    cont_model[i][j][k][l] = y_T;

    //free the memory holding the pass band data
    FreePassBand();

#ifndef PBINTERP
    free(x_pb_o);
#endif //PBINTERP
  }
}
