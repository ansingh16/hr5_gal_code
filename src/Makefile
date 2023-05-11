#############################################################################
# Fortran compiler options and directives




OPTIONS = -cpp -fopenmp
OPTIONS += -DDISPLAY_PROGRESS_PERCENT



F90         = /opt/ohpc/pub/intel/2020up4/compilers_and_libraries_2020.4.304/linux/bin/intel64/ifort
#FFLAGS      = -g -traceback -fpp -check all -debug -warn all -ftrapuv 
##FFLAGS      = -O3 -fpp -ftz  -xavx -sox -ipo -fno-alias -fno-fnalias -no-prec-div -no-prec-sqrt -align all
FFLAGS      = -O3 -fpp -ftz -ipo -fno-alias -fno-fnalias -no-prec-div -no-prec-sqrt -align all
##FFLAGS      = -p -g
FFLAGS += $(OPTIONS)
LDFLAGS     = 



# 3/ pure HI
COMPOBJS = module_HI_1216_model.o module_gas_composition_HI.o


GASOBJS = module_constants.o module_random.o module_utils.o module_voigt.o module_domain.o coolrates_module.o module_ramses.o module_uparallel.o $(COMPOBJS)

OBJSDOM = $(GASOBJS) module_select.o module_mesh.o


#############################################################################
#.SUFFIXES: .o .f90
#.f90.o:
#	$(F90) -c  $(FFLAGS) $<

%.o:%.f90
	$(F90) $(FFLAGS) -c $^ -o $@

#############################################################################
all: CreateDomDump get_time


CreateDomDump: $(OBJSDOM) CreateDomDump.o
	$(F90) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJSDOM) CreateDomDump.o

get_time:
	f2py -c -m time_ramses time_ramses.f90 --f90exec=${F90} --f77exec=${F90} 


clean:
	rm -f $(OBJSDOM)  *.o *.mod
#############################################################################
