SHELL=/bin/sh
#
# This makefile was produced by /usr/bin/ftnmgen
# at 05:07:18 PM on 02/29/00
#
# If it is invoked by the command line
#	make -f makefile
# it will compile the fortran modules indicated by SRCS into the object
# modules indicated by OBJS and produce an executable named mis.
#
# To remove all the objects but leave the executables use the command line
#	make -f makefile clean
#
# To remove everything but the source files use the command line
#	make -f makefile clobber
#
# To remove the source files created by /usr/bin/ftnmgen
# and this makefile use the command line
#	make -f makefile void
#
# The parameters SRCS and OBJS should not need to be changed.  If, however,
# you need to add a new source file add the name of the source file to the
# SRCS parameter and add the name of the resulting object file to the OBJS
# parameter.
#
#!SRCS=	modules.f90 flow2d_i.f90 del-flux.f90 stab.f90 cell.f90 ausmdv1_flux.f90 bc.f90 im_fl5.f90 muscl_prim.f90 bc_vis.f90 flux_vis.f90

#OBJS=	modules.o flow2d_i.o del-flux.o stab.o cell.o ausmdv1_flux.o bc.o im_fl5.o muscl_prim.o bc_vis.o flux_vis.o
SRCS=	modules.f90 flow2d_i.f90 del-flux.f90 stab.f90 cell.f90 inviscid_flux.f90 bc.f90 im_fl5.f90 muscl_prim.f90 bc_vis.f90 flux_vis.f90 results.f90 read_input.f90 turbulencemodels.f90 

OBJS=	modules.o flow2d_i.o del-flux.o stab.o cell.o inviscid_flux.o bc.o im_fl5.o muscl_prim.o bc_vis.o flux_vis.o  results.o  read_input.o turbulencemodels.o

# Tunable parameters
#
# FC		Name of the fortran compiling system to use
# FFLAGS	Flags to the compiler
#
#  Origin 200: -64 -extend_source -mips4 -r10000 -O3
#  O2        : -extend_source -mips4 -r5000 -O3
#  Challenge : -64 -extend_source -mips4 -SWP:=ON -O3 -OPT:fold_arith_limit=2599
#
# LDFLAGS	Flags to the loader
# LIBS		List of libraries
# CMD		Name of the executable
#

#FC  =		f95
FC  =		ifort
#----------------------
#FC  =	gfortran
# gfortran Flag
#FFLAGS =  -Wall -O5
FFLAGS =   -O5
#----------------------
#FFLAGS =  -CB
#FFLAGS =  -openmp  -parallel
#FFLAGS =  -C -openmp 
#FFLAGS =  -ipo -O3
#FFLAGS =   -g -C
#FFLAGS =  -e -C
#FFLAGS = -e -fast -pad -xipo -stackvar
#FFLAGS =  -ipo 
#FFLAGS = -extend_source -O3 
#LDFLAGS =	
#LIBS =		
CMD =		mis

# Rules to build mis.
#
all:		$(CMD)
$(CMD):		$(OBJS)
	   $(FC) $(LDFLAGS) $(FFLAGS) -o $(CMD) $(OBJS) $(LIBS)

modules.o:     modules.f90
	   $(FC) -c $(LDFLAGS) $(FFLAGS)  modules.f90 $(LIBS)

flow2d_i.o:     flow2d_i.f90
	   $(FC) -c $(LDFLAGS) $(FFLAGS) flow2d_i.f90 $(LIBS)

del-flux.o:     del-flux.f90
	   $(FC) -c $(LDFLAGS) $(FFLAGS) del-flux.f90 $(LIBS)

stab.o:     stab.f90
	   $(FC) -c $(LDFLAGS) $(FFLAGS) stab.f90 $(LIBS)

cell.o:   cell.f90
	   $(FC) -c $(LDFLAGS) $(FFLAGS) cell.f90 $(LIBS)

inviscid_flux.o:  inviscid_flux.f90 
	   $(FC) -c $(LDFLAGS) $(FFLAGS) inviscid_flux.f90 $(LIBS)

#ausmdv1_flux.o:   ausmdv1_flux.f90
#	   $(FC) -c $(LDFLAGS) $(FFLAGS) ausmdv1_flux.f90 $(LIBS)

bc.o:   bc.f90
	   $(FC) -c $(LDFLAGS) $(FFLAGS) bc.f90 $(LIBS)

im_fl5.o:   im_fl5.f90
	   $(FC) -c $(LDFLAGS) $(FFLAGS) im_fl5.f90 $(LIBS)

muscl_prim.o:   muscl_prim.f90
	   $(FC) -c $(LDFLAGS) $(FFLAGS) muscl_prim.f90 $(LIBS)

bc_vis.o:   bc_vis.f90
	   $(FC) -c $(LDFLAGS) $(FFLAGS) bc_vis.f90 $(LIBS)

flux_vis.o:   flux_vis.f90
	   $(FC) -c $(LDFLAGS) $(FFLAGS) flux_vis.f90 $(LIBS)

read_input.o:  read_input.f90 
	   $(FC) -c $(LDFLAGS) $(FFLAGS) read_input.f90 $(LIBS)

results.o: results.f90 
	   $(FC) -c $(LDFLAGS) $(FFLAGS) results.f90 $(LIBS)

turbulencemodels.o:  turbulencemodels.f90 
	   $(FC) -c $(LDFLAGS) $(FFLAGS) turbulencemodels.f90 $(LIBS)

clean:
	-rm -f $(OBJS) *.mod

clobber:	clean
	-rm -f $(CMD)

void:	clobber
	-rm -f $(SRCS) makefile
