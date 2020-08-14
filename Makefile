NESTLIBDIR = /Users/brant/code/multinest/lib
#LIBS =  -L$(NESTLIBDIR) -lnest3 -L/usr/local/gfortran/lib $(LAPACKLIB) -lblas -llapack -lgfortran -lgsl
LIBS =  -L$(NESTLIBDIR) -lnest3 -L/usr/local/gfortran/lib $(LAPACKLIB) -lgfortran -lgsl  -lblas -llapack -lm
OBJFILES = fit_bpass_seds.o bpass_models.o madau_absorption.o pass_bands.o igm-absorption.o

CFLAGS = -DREDSHIFT2

CC = cc
all: fit_bpass_seds

%.o: %.c
	$(CC) $(CFLAGS) -c $*.c
 
fit_bpass_seds: $(OBJFILES)
	$(CC) $(FFLAGS) -m64 -o fit_bpass_seds $(OBJFILES) $(LIBS)

clean:
	rm -f *.o *.mod fit_bpass_seds
