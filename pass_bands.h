#ifndef PASS_BANDS_H
#define PASS_BANDS_H

extern int n_pb;      //number of samples in the transmissivity curve
extern double *x_pb;  //wavelength
extern double *y_pb;  //transmissivity

void LoadPassBand(int l);
void FreePassBand(void);

#endif /*PASS_BANDS_H*/