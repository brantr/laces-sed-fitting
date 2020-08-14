#ifndef BPASS_MODELS_H
#define BPASS_MODELS_H

//source properties
extern double z;          //source redshift
extern int n_data;       //number of data points
extern double *data_x;   //data points
extern double *data_y;   //data points
extern double *data_ye;  //uncertainty on data points
extern int ID;           //galaxy ID


//model properties
extern int n_age;          //number of age bins
extern int n_z_met;        //number of metallicity samples
extern int n_f_bin;        //number of binarity samples
extern int n_lambda;       //number of wavelength samples in the SED model
extern double *age;        //array of age values
extern double *z_met;      //array of metallicity values
extern double *f_bin;      //array of binarity values
extern double *lambda;     //array of wavelengths
extern double *f_nu;       //array of SED samples
extern int flag_binary;    //what binary are we fitting?

extern int n_lambda_line;       //number of wavelength samples in the SED model
extern int n_lambda_cont;       //number of wavelength samples in the SED model
extern double *f_nu_line;       //array of SED samples
extern double *f_nu_cont;       //array of SED samples


//n_bin x n_met x n_age x n_data 
//array of sed values
extern double ****sed_model;  
extern double ****line_model;  
extern double ****cont_model;  

void load_sed_model(char fname[], double z);
void load_line_model(char fname[], double z, double lNlyc);
void load_cont_model(char fname[], double z, double lNlyc);


//initialize the bpass models
void initialize_bpass_models(void);

//load and attenuate the SED model samples
void load_and_attenuate_sed_samples(void);

//allocate sed samples
void allocate_sed_samples(void);

//deallocate sed samples
void deallocate_sed_samples(void);

//allocate the data for an object
void allocate_data(void);

//free the data associated with an object
void free_data(void);

//read in the data for an object
void read_data(char fname[]);

//now we attenuate the SED model
void attenuate_sed_model(int *flag_igm_abs);
void attenuate_line_model(int *flag_igm_abs);
void attenuate_cont_model(int *flag_igm_abs);

//now we sample and passband-average the SED model
void pass_band_average_sed_model(int i, int j, int k);
void pass_band_average_line_model(int i, int j, int k);
void pass_band_average_cont_model(int i, int j, int k);

#endif  //BPASS_MODELS_H
