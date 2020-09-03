#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include "multinest.h"
#include "bpass_models.h"
#include <gsl/gsl_interp.h>

//#define CHECK_PHOT

#define SMC_DUST
#define EBV_PRIOR

//#define SFR_PRIOR

//#define REVIEW_INFO
//#define SHOW_CHI2
//#define FIT_LOG_SFR




/*
0 F336W 3.360000e+03
1 U 3.850000e+03
2 B 4.500000e+03
3 NB 4.975000e+03
4 V 5.500000e+03
5 R 6.500000e+03
6 i 7.600000e+03
7 z 8.900000e+03
8 J 1.250000e+04
9 H 1.550000e+04
10 K 2.200000e+04
11 IRAC1 3.600000e+04
12 IRAC2 4.500000e+04
*/


//avoid fitting to F336W
//#define NO_FIT_F336W

//avoid fitting to U
//#define NO_FIT_U

//avoid fitting to B
//#define NO_FIT_B

//avoid fitting to V
//#define NO_FIT_V

//avoid fitting to R
//#define NO_FIT_R

//avoid fitting to I
//#define NO_FIT_I

//avoid fitting to Z
//#define NO_FIT_Z

//avoid fitting to J
//#define NO_FIT_J


//avoid fitting to K
//#define NO_FIT_K

//avoid fitting to F160W
//#define NO_FIT_F160W

//avoid fitting to IRAC
//#define NO_FIT_IRAC

//#define NO_FIT_OPTICAL


//force small escape fractions
//#define NO_ESCAPE

//fit to the line strengths individually
#define FIT_LINE_STRENGTHS

//fit to the continuum strength individually
//#define FIT_CONT_STRENGTH

//affix CONT to B_lines rather than B_cont
//#define COUPLE_LINES_AND_CONT

#define FIT_NON_DETECTIONS


#define FIT_F_ESC
#define FIT_LINES
//#define FIT_CONTINUUM
#define FIT_EBV

#ifndef FIT_LOG_SFR
#define A_MAX 500.0
#define A_MIN  0.0
#else //FIT_LOG_SFR
#define A_MAX 3.0
#define A_MIN -3.0
#endif //FIT_LOG_SFR

#define F_ESC_MIN 0.0
#ifndef NO_ESCAPE
#define F_ESC_MAX 1.0
#else //NO_ESCAPE
#define F_ESC_MAX 0.001
#endif //NO_ESCAPE

#define EBV_MIN 0.0
#define EBV_MAX 0.25
//#define EBV_MAX 0.5
//#define EBV_MAX 1.0

#define B_LINE_MIN 0.0
#define B_LINE_MAX 10.0

#define B_CONT_MIN 0.0
#define B_CONT_MAX 10.0

double klambda_smc(double lambda)
{
	double lambda_smc[31] = {0.09,0.116,0.119,0.123,0.127,0.131,0.136,0.140,0.145,0.151,0.157,0.163,0.170,0.178,0.186,0.195,0.205,0.216,0.229,0.242,0.258,0.276,0.296,0.370,0.440,0.550,0.660,0.810,1.250,1.650,2.198};
	double k_smc[31] = {9.5,6.992,6.436,6.297,6.074,5.795,5.575,5.272,5.000,4.776,4.472,4.243,4.013,3.866,3.637,3.489,3.293,3.161,2.947,2.661,2.428,2.220,2.000,1.672,1.374,1.000,0.801,0.567,0.131,0.169,0.016};
	double Rv_smc = 2.74; // gordon et al. 2003
	int nlsmc = 31;
	int    ki = gsl_interp_bsearch(lambda_smc,lambda,0,nlsmc);
	double ksmc =  k_smc[ki] + (k_smc[ki+1] - k_smc[ki])*(lambda - lambda_smc[ki])/(lambda_smc[ki+1] - lambda_smc[ki]);
	return ksmc * Rv_smc;
}




double gaussian_func(double x, double sigma, double mu)
{
	double A = 1./sqrt(2.0*M_PI*sigma*sigma);
	double xx = (x-mu)/sigma;
	return A*exp(-0.5*xx*xx);
}


/******************************************** loglikelihood routine ****************************************************/
// Input arguments
// ndim 						= dimensionality (total number of free parameters) of the problem
// npars 						= total number of free plus derived parameters
// context						void pointer, any additional information
//
// Input/Output arguments
// Cube[npars] 						= on entry has the ndim parameters in unit-hypercube
//	 						on exit, the physical parameters plus copy any derived parameters you want to store with the free parameters
//	 
// Output arguments
// lnew 						= loglikelihood

void LogLikeExample(double *Cube, int *ndim, int *npars, double *lnew, void *context)
{
	double chi = 1.0;
	int i;
	double sigma = 10.0*Cube[0];
	double mu    = 10.0*Cube[1] - 5.0;
	double l = 0;
	double x;
	double y;
	for(i=0;i<n_data;i++)
	{
		y = gaussian_func(data_x[i],sigma,mu);
		x = (y-data_y[i])/data_ye[i];
	}
	Cube[0] = sigma;
	Cube[1] = mu;
	*lnew = l;
}


void LogLike(double *Cube, int *ndim, int *npars, double *lnew, void *context)
{

	//we have the following parameters
	//A     == amplitude
	//f_bin == binarity
	//z_met == metallicity
	//age   == age

	//we perform an interpolation on our model grid
	//and then compare the model SED with the data

	int ip = 0;

	double A = A_MIN + Cube[ip]*(A_MAX-A_MIN);								 //SED amplitude
#ifdef FIT_LOG_SFR
	A = pow(10.0,A);
#endif //FIT_LOG_SFR
	ip++;

	//limit age to 9.3081373786380386
	//double a_age   = age[0] + (age[n_age-1]-age[0])*Cube[1]; //log10 age
	//double age_max = 9.3081373786380386; // z=3
	double age_max = 9.5206145218782368; // z=2
	//double age_min = 7.0; 
	double age_min = age[0];
	//double age_max = 8.9894498176666922; // z=5.7
	//double a_age   = age[0] + (age_max-age[0])*Cube[ip]; //log10 age
	double a_age   = age_min + (age_max-age_min)*Cube[ip]; //log10 age

	ip++;

#ifdef FIT_F_ESC
	double f_esc = F_ESC_MIN + Cube[ip]*(F_ESC_MAX-F_ESC_MIN);
	ip++;
#endif //FIT_F_ESC


#ifdef FIT_EBV
	double ebv = EBV_MIN + Cube[ip]*(EBV_MAX-EBV_MIN);
	ip++;
#endif //FIT_EBV




#ifdef FIT_LINE_STRENGTHS
	double B_line = B_LINE_MIN + Cube[ip]*(B_LINE_MAX-B_LINE_MIN);
	ip++;
#endif //FIT_LINE_STRENGTHS

#ifdef FIT_CONT_STRENGTH
	double B_cont = B_CONT_MIN + Cube[ip]*(B_CONT_MAX-B_CONT_MIN);
	ip++;
#endif //FIT_CONT_STRENGTH

#ifdef CHECK_PHOT
	A = 0.139265689050300772E+02;
	a_age = 0.649921138293866285E+01;
	f_esc = 0.0;
	//f_esc = 0.654566039812631617E+00;
	//f_esc = 1.0;
	ebv = 0.331014621373948090E-01;
	//ebv = 0.0;
	//f_esc = 0.5;
	f_esc = 0.25;
	B_line = 0.1;


	//(1-f_esc) * B_line
	A = 0.126595501423916801E+02;
	a_age = 0.656199302980435029E+01;
	f_esc = 0.847544147895375577E+00;
	ebv = 0.392088356147152345E-01;
#ifdef FIT_LINE_STRENGTHS
	B_line = 0.216287787816939758E+01;
#endif // FIT_LINE_STRENGTHS

	//Split B_line entirely




	A = 0.256283612365030855E+02;
	a_age = 0.631154324433783476E+01;
	f_esc = 0.203669973252986802E+00;
	ebv = 0.370709643383516987E-01;
#ifdef FIT_LINE_STRENGTHS
	B_line = 0.311503571759184694E+00;
#endif // FIT_LINE_STRENGTHS


//decouple lines and continuum
	A = 0.513693289556891060E+02;
	a_age = 0.616947367726794837E+01;
	f_esc = 0.146148059582514034E+00;
	ebv = 0.553256976460332478E-01;
#ifdef FIT_LINE_STRENGTHS
	B_line = 0.241206755469208001E+00;
#endif // FIT_LINE_STRENGTHS
#ifdef FIT_CONT_STRENGTH
	B_cont = 0.800726399751585174E-02;
#endif // FIT_CONT_STRENGTH

	printf("A      = %e\n",A);
	printf("age    = %e\n",a_age);
	printf("f_esc  = %e\n",f_esc);
	printf("ebv    = %e\n",ebv);
#ifdef FIT_LINE_STRENGTHS
	printf("B_line = %e\n",B_line);
#endif // FIT_LINE_STRENGTHS

#ifdef FIT_CONT_STRENGTH
	printf("B_cont = %e\n",B_cont);
#endif // FIT_CONT_STRENGTH

#endif // CHECK_PHOT



	double a_f_bin = 0;	//not using currently
	double a_z_met = 0;	//not using currently

	int m;
	double l = 0;
	double x;
	double y;
#ifdef FIT_LINES
	double yl;
#endif //FIT_LINES
#ifdef FIT_CONTINUUM
	double yc;
#endif //FIT_CONTINUUM
	//interpolation along age direction
	int     k_age = gsl_interp_bsearch(age,a_age,0,n_age);
	double dx_age = (a_age - age[k_age])/(age[k_age+1] - age[k_age]);

	double chi2 = 0;

#ifdef CHECK_PHOT
	printf("PHOT k_age %d dx_age %e\n",k_age,dx_age);
#endif //CHECK_PHOT
	
	for(m=0;m<n_data;m++)
	{
		//just interpolate along age direction to start
		y = (1.0-dx_age)*sed_model[0][0][k_age][m] + dx_age*sed_model[0][0][k_age+1][m];
#ifdef FIT_LINES
		yl = (1.0-dx_age)*line_model[0][0][k_age][m] + dx_age*line_model[0][0][k_age+1][m];
#endif //FIT_LINES
#ifdef FIT_CONTINUUM
		yc = (1.0-dx_age)*cont_model[0][0][k_age][m] + dx_age*cont_model[0][0][k_age+1][m];
#endif //FIT_CONTINUUM

#ifdef FIT_F_ESC

#ifdef FIT_LINES
#ifndef FIT_LINE_STRENGTHS
		y += (1.0-f_esc)*yl;
#else  //FIT_LINE_STRENGTHS
		//y += (1.0-f_esc)*B_line*yl;
		y += B_line*yl;
#endif //FIT_LINE_STRENGTHS
#endif //FIT_LINES



#ifdef FIT_CONTINUUM
#ifndef FIT_CONT_STRENGTH

#ifndef COUPLE_LINES_AND_CONT
		y += (1.0-f_esc)*yc;
#else //COUPLE_LINES_AND_CONT
		y += B_line*yc;
#endif //COUPLE_LINES_AND_CONT
#else //FIT_CONT_STRENGTH
		y += B_cont*yc;
#endif //FIT_CONT_STRENGTH
#endif //FIT_CONTINUUM


		if(m==0)
			y*=f_esc;

#endif //FIT_F_ESC

#ifdef FIT_EBV
		double AEBV = 0;
		double kp;

#ifndef SMC_DUST
		double Rvp = 4.05;
		double xx = data_x[m] * 1.0e-4;	//note data_x is in rest frame

		kp = 2.659*(-1.857 + 1.040/xx) + Rvp;
		//printf("m %d data_x %e xx %e\n",m,data_x[m],xx);
		if(xx<0.63)
			kp = 2.659*(-2.156 + 1.509/xx - 0.198/(xx*xx) + 0.011/pow(xx,3))+ Rvp;
		if(xx<0.0912)
			kp = 0.0;
#else	//SMC_DUST
		kp = klambda_smc(data_x[m] * 1.0e-4);
		if(data_x[m] * 1.0e-4 < 0.0912)
			kp = 0.0;
#endif 	//SMC_DUST
		AEBV = -0.4*ebv*kp;
		if(AEBV>0)
		{
			printf("ERROR AEBV %e\n",AEBV);
			exit(0);
		}

//#ifdef CHECK_PHOT
//		printf("BEFORE DUST AEBV %e ebv %e y prev %e A*y_prev %e\n",AEBV,ebv,y,A*y);
//#endif //CHECK_PHOT

//		printf("rescaling y %e by %e\n",y,pow(10.0,AEBV));
		y *= pow(10.0,AEBV);
#endif //FIT_EBV

//#ifdef CHECK_PHOT
//		printf("AFTER DUST AEBV %e ebv %e y prev %e A*y_prev %e\n",AEBV,ebv,y,A*y);
//#endif //CHECK_PHOT


		x = (A*y-data_y[m])/data_ye[m];

#ifdef CHECK_PHOT
	printf("PHOT m %d lambda %e y %e data_y %e data_ye %e x %e\n",m,data_x[m]*(1+z),A*y,data_y[m],data_ye[m],x);
#endif //CHECK_PHOT

// OLD FITTING
/*
#ifndef FIT_NON_DETECTIONS
		if(((data_y[m]>1.0e-10)&&(fabs(data_x[m]-1215.67)>25.))) //avoid non-detections for the time being, and Lya
#else	//FIT_NON_DETECTIONS

		//fit non detections
		if(data_y[m]<1.0e-10)
			x = A*y/data_ye[m]; //zero mean gaussian, 1-sigma errors


		//if(fabs(data_x[m]-1215.67)>25.) //avoid Lya
		if( (fabs(data_x[m]-1215.67)>25.) && (data_y[m]>0) ) //avoid Lya and negative values
#endif //FIT_NON_DETECTIONS
		{
			//add multiplicative contribution to - log likelihood

			int flag_fit = 1;



#ifdef NO_FIT_F336W
			if(m!=0)
#endif //NO_FIT_F336W
#ifdef NO_FIT_IRAC
			if(m<n_data-2)
#endif //NO_FIT_F336W
#ifdef NO_FIT_OPTICAL
			if(m<n_data-3)
#endif //NO_FIT_F336W


			//if this data point contributes
			//then add the contribution to chi2

			if(flag_fit)
				l += -0.5*x*x;

		}
*/

// NEW FITTING
#ifndef FIT_NON_DETECTIONS
		if((data_y[m]>1.0e-10))
#else	//FIT_NON_DETECTIONS

		//fit non detections
		if(data_y[m]<1.0e-10)
			x = A*y/data_ye[m]; //zero mean gaussian

		if(1)
#endif //FIT_NON_DETECTIONS
		{
			//add multiplicative contribution to - log likelihood
			int flag_fit = 1;

			//avoid lyman alpha, never fit narrow band
			if(m==3)
				flag_fit = 0;


#ifdef NO_FIT_F336W
			if(m==0)
				flag_fit = 0;
#endif //NO_FIT_F336W

#ifdef NO_FIT_U
			if(m==1)
				flag_fit = 0;
#endif //NO_FIT_U

#ifdef NO_FIT_B
			if(m==2)
				flag_fit = 0;
#endif //NO_FIT_B


#ifdef NO_FIT_V
			if(m==4)
				flag_fit = 0;
#endif //NO_FIT_V


#ifdef NO_FIT_R
			if(m==5)
				flag_fit = 0;
#endif //NO_FIT_R

#ifdef NO_FIT_I
			if(m==6)
				flag_fit = 0;
#endif //NO_FIT_I

#ifdef NO_FIT_Z
			if(m==7)
				flag_fit = 0;
#endif //NO_FIT_Z

#ifdef NO_FIT_J
			if(m==8)
				flag_fit = 0;
#endif //NO_FIT_J

#ifdef NO_FIT_F160W
			if(m==9)
				flag_fit = 0;
#endif //NO_FIT_F160W

#ifdef NO_FIT_K
			if(m==10)
				flag_fit = 0;
#endif //NO_FIT_K

#ifdef NO_FIT_IRAC
			if(m>=n_data-2)
				flag_fit = 0;
#endif //NO_FIT_F336W
#ifdef NO_FIT_OPTICAL
			if(m>=n_data-3)
				flag_fit = 0;
#endif //NO_FIT_F336W


			//if this data point contributes
			//then add the contribution to chi2

			if(flag_fit)
				l += -0.5*x*x;

		}


		//save chi2 for this object
		chi2 += 0.5*x*x;
	}

	//these prior constraints go outside
	//and apply only once


	//add a prior on EBV?

#ifdef EBV_PRIOR
	double x_ebv = (ebv - 0.0)/0.05;
	//double x_ebv = (ebv - 0.0)/0.1;
	//double x_ebv = (ebv - 0.0)/0.2;

	l += -0.5*x_ebv*x_ebv;
#endif //EBV_PRIOR


	//add a prior on SFR?
#ifdef SFR_PRIOR
#ifdef FIT_LOG_SFR
	double x_sfr = (log10(A) - 0.0)/3.0;
#else //FIT_LOG_SFR
	double x_sfr = (A - 1.0)/10.;
#endif //FIT_LOG_SFR
	l += -x_sfr*x_sfr;
#endif //SFR_PRIOR

	//Rescale the parameters
	ip = 0;

#ifdef FIT_LOG_SFR
	A = log10(A);
#endif //FIT_LOG_SFR
	Cube[ip] = A;		//SFR in Msun/yr
	ip++;
	Cube[ip] = a_age;	//log_10 Age of stellar population
	ip++;
#ifdef FIT_F_ESC
	Cube[ip] = f_esc;	//escape fraction at <912 Angstroms
	ip++;
#endif //FIT_F_ESC
#ifdef FIT_EBV
	Cube[ip] = ebv;		//E(B-V)
	ip++;
#endif
#ifdef FIT_LINE_STRENGTHS
	Cube[ip] = B_line;  //multiplicative line strength
	ip++;
#endif// FIT_LINE_STRENGTHS
#ifdef FIT_CONT_STRENGTH
	Cube[ip] = B_cont;  //multiplicative continuum strength
	ip++;
#endif// FIT_CONT_STRENGTH

	*lnew = l;			//save the likelihood

#ifdef SHOW_CHI2
	printf("CHI2: A %e age %e f_esc %e ebv %e chi2 %e\n",A,a_age,f_esc,ebv,chi2);
#endif //SHOW_CHI2

#ifdef CHECK_PHOT
	printf("CHI2: A %e age %e f_esc %e ebv %e chi2 %e\n",A,a_age,f_esc,ebv,chi2);
	exit(0);
#endif //CHECK_PHOT
}

/***********************************************************************************************************************/




/************************************************* dumper routine ******************************************************/

// The dumper routine will be called every updInt*10 iterations
// MultiNest doesn not need to the user to do anything. User can use the arguments in whichever way he/she wants
//
//
// Arguments:
//
// nSamples 						= total number of samples in posterior distribution
// nlive 						= total number of live points
// nPar 						= total number of parameters (free + derived)
// physLive[1][nlive * (nPar + 1)] 			= 2D array containing the last set of live points (physical parameters plus derived parameters) along with their loglikelihood values
// posterior[1][nSamples * (nPar + 2)] 			= posterior distribution containing nSamples points. Each sample has nPar parameters (physical + derived) along with the their loglike value & posterior probability
// paramConstr[1][4*nPar]:
// paramConstr[0][0] to paramConstr[0][nPar - 1] 	= mean values of the parameters
// paramConstr[0][nPar] to paramConstr[0][2*nPar - 1] 	= standard deviation of the parameters
// paramConstr[0][nPar*2] to paramConstr[0][3*nPar - 1] = best-fit (maxlike) parameters
// paramConstr[0][nPar*4] to paramConstr[0][4*nPar - 1] = MAP (maximum-a-posteriori) parameters
// maxLogLike						= maximum loglikelihood value
// logZ							= log evidence value
// logZerr						= error on log evidence value
// context						void pointer, any additional information

void dumper(int *nSamples, int *nlive, int *nPar, double **physLive, double **posterior, double **paramConstr, double *maxLogLike, double *logZ, double *logZerr, void *context)
{
	// convert the 2D Fortran arrays to C arrays
	
	
	// the posterior distribution
	// postdist will have nPar parameters in the first nPar columns & loglike value & the posterior probability in the last two columns
	
	int i, j;
	
	double postdist[*nSamples][*nPar + 2];
	for( i = 0; i < *nPar + 2; i++ )
		for( j = 0; j < *nSamples; j++ )
			postdist[j][i] = posterior[0][i * (*nSamples) + j];
	
	
	
	// last set of live points
	// pLivePts will have nPar parameters in the first nPar columns & loglike value in the last column
	
	double pLivePts[*nlive][*nPar + 1];
	for( i = 0; i < *nPar + 1; i++ )
		for( j = 0; j < *nlive; j++ )
			pLivePts[j][i] = physLive[0][i * (*nlive) + j];
}

/***********************************************************************************************************************/




/************************************************** Main program *******************************************************/



int main(int argc, char *argv[])
{
	int i;
	
	// set the MultiNest sampling parameters
	
	
	//int mmodal = 1;					// do mode separation?
	int mmodal = 0;					// do mode separation?
	
	int ceff = 0;					// run in constant efficiency mode?
	
	//int nlive = 10000;				// number of live points
	int nlive = 2000;				// number of live points
	//int nlive = 50000;				// number of live points

	double efr = 1.0;				// set the required efficiency
	
	double tol = 0.5;				// tol, defines the stopping criteria
	
	int ndims = 2;					// dimensionality (no. of free parameters)
	
	int nPar = 2;					// total no. of parameters including free & derived parameters
	
	int nClsPar = 2;				// no. of parameters to do mode separation on

#ifdef FIT_F_ESC
	ndims++;					// dimensionality (no. of free parameters)
	nPar++;					// dimensionality (no. of free parameters)
	nClsPar++;					// dimensionality (no. of free parameters)
#endif //FIT_F_ESC
#ifdef FIT_EBV
	ndims++;					// dimensionality (no. of free parameters)
	nPar++;					// dimensionality (no. of free parameters)
	nClsPar++;					// dimensionality (no. of free parameters)
#endif //FIT_EBV
#ifdef FIT_LINE_STRENGTHS
	ndims++;					// dimensionality (no. of free parameters)
	nPar++;					// dimensionality (no. of free parameters)
	nClsPar++;					// dimensionality (no. of free parameters)
#endif //FIT_LINE_STRENGTHS
#ifdef FIT_CONT_STRENGTH
	ndims++;					// dimensionality (no. of free parameters)
	nPar++;					// dimensionality (no. of free parameters)
	nClsPar++;					// dimensionality (no. of free parameters)
#endif //FIT_CONT_STRENGTH
		

	int updInt = 100;				// after how many iterations feedback is required & the output files should be updated
							// note: posterior files are updated & dumper routine is called after every updInt*10 iterations
	
	double Ztol = -1E90;				// all the modes with logZ < Ztol are ignored
	
	int maxModes = 10;				// expected max no. of modes (used only for memory allocation)
	
	int pWrap[ndims];				// which parameters to have periodic boundary conditions?
	for(i = 0; i < ndims; i++) pWrap[i] = 0;
	
	char root[100] = "chains/test_gaussian-";		// root for output files
	
	int seed = -1;					// random no. generator seed, if < 0 then take the seed from system clock
	
	int fb = 1;					    // need feedback on standard output?
	
	int resume = 0;					// resume from a previous job?
	
	int outfile = 1;				// write output files?
	
	int initMPI = 1;				// initialize MPI routines?, relevant only if compiling with MPI
							// set it to F if you want your main program to handle MPI initialization
	
	double logZero = -DBL_MAX;			// points with loglike < logZero will be ignored by MultiNest
	
	int maxiter = 0;				// max no. of iterations, a non-positive value means infinity. MultiNest will terminate if either it 
							// has done max no. of iterations or convergence criterion (defined through tol) has been satisfied
	
	void *context = 0;				// not required by MultiNest, any additional information user wants to pass

	//load the data first
	char fname[200];
	sprintf(fname,"fnudata_total_02_90340");
	if(argc>1)
	{
		sprintf(fname,"%s",argv[1]);
	}

	printf("Data file name = %s.\n",fname);

	//read data
	read_data(fname);

	//initialize the bpass models
	//1) read the bpass files
	//2) based on the source redshift, apply IGM attenuation
	//3) remember only the sampled locations on the SED
	initialize_bpass_models();

#ifdef REVIEW_INFO
	//stop here to review diagnostic
	//information about the galaxies
	exit(0);
#endif //REVIEW_INFO

	// calling MultiNest
	
	run(mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, 
	logZero, maxiter, LogLike, dumper, context);
}

/***********************************************************************************************************************/
