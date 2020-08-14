/*! \file igm_absorption.h
 *  \brief Function declaration for IGM absorption following Inoue ea 2014.
 *
 *  See Inoue et al. 2014
 */
#ifndef IGM_ABSORPTION_H
#define IGM_ABSORPTION_H
/*! \fn double *ign_absorption(double *lambda, int n_lambda, double z);
 *  \brief Wavelength-dependent absorption owing to neutral H in the IGM.
 */
double *igm_absorption(double *lambda, int n_lambda, double z);

//print to a file
void print_igm_absorption(double *lambda, int n_lambda, double z);

double tau_IGM(double lobs, double z);
double tau_DLA_LC(double lobs, double z);
double tau_LAF_LC(double lobs, double z);
double tau_DLA_LS(double lobs, double z);
double tau_LAF_LS(double lobs, double z);

#endif //IGM_ABSORPTION_H