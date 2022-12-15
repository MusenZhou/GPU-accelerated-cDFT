# Makefile for none
# Copyright (C) 2012 
# $Header$
# $Log$


FC = ifort -I "/home/jiafu/intel/mkl/include"
LIB =  -w -L "/home/jiafu/intel/mkl/lib/intel64" -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm  -liomp5  -ldl 


OBJS95 = ADFT.o fftw_fit.f90 temps.o\
	Compute_Torrence.o sysdef2.o Compute_initial_rho.o  \
        Compute_Vext2.o Compute_DCFs.o Compute_new_rho.o Compute_Vint.o \
        Compute_bridge.o Compute_Intra_Factor.o Get_cos_table.o \
        Compute_ele_pot_3D.o final_output.o


SRCS90 = $(OBJS90:.o=.f90)
SRCS95 = $(OBJS95:.o=.f)

all: $(OBJS90) $(OBJS95)
	$(FC) $(OBJS90) $(OBJS95) $(LIB) -o dft 

$(OBJS90):$(SRCS90)
	ifort -c $(SRCS90)

$(OBJS95):$(SRCS95)
	ifort -c $(SRCS95)

clean:
	rm -f $(OBJS95) $(OBJS90) dft

