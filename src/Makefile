BINPATH = ../bin/
#
# MORSE POTENTIAL
tise_1D_morse_SRC = nrtype.f90 numerov_global_vars.f90 phys_constants.f90 pot_param.f90 numerov_algorithm.f90 \
f_morse.f90 x_sum_rule.f90 x2_sum_rule.f90 ew_sum_rules.f90\
tise_1D_numerov.f90
#
tise_1D_morse_OBJ = $(tise_1D_morse_SRC:.f90=.o)
#
#
# PTELLER POTENTIAL
tise_1D_pteller_SRC = nrtype.f90 numerov_global_vars.f90 phys_constants.f90 pot_param.f90 numerov_algorithm.f90 \
f_pteller.f90 x_sum_rule.f90 x2_sum_rule.f90 ew_sum_rules.f90\
tise_1D_numerov.f90
#
tise_1D_pteller_OBJ = $(tise_1D_pteller_SRC:.f90=.o)
#
#
# WSAXON POTENTIAL
tise_1D_wsaxon_SRC = nrtype.f90 numerov_global_vars.f90 phys_constants.f90 pot_param.f90 numerov_algorithm.f90 \
f_wsaxon.f90 x_sum_rule.f90 x2_sum_rule.f90 ew_sum_rules.f90\
tise_1D_numerov.f90
#
tise_1D_wsaxon_OBJ = $(tise_1D_wsaxon_SRC:.f90=.o)
#
#
FC = gfortran
FOPT	= -c 
FOPT   += -Wall
#MODLIB  = -I/usr/local/lib/minpack
#FOPT   += $(MODLIB)
#
LIB_PATH = -L/usr/local/lib
LIB_MINPACK = -lminpack 
LIB_INTLIB = -lintlib 
# NAG
NAGLIB = -L/usr/local/lib/NAG -lnag
#
#
LIBS     = $(LIB_PATH) $(LIB_MINPACK) $(LIB_INTLIB) $(NAGLIB)
#
all: morse pteller wsaxon
morse: tise_1D_morse 
pteller: tise_1D_pteller 
wsaxon: tise_1D_wsaxon 
.PHONY : all
#
# Standard compilation rules
#
%.o: %.f90
	$(FC) $(FOPT) -c $<
#
%.o: %.F90
	$(FC) $(FOPT) -c $<
#
%.o: %.f
	$(FC) $(FOPT) $<
#
#######################
tise_1D_morse: $(tise_1D_morse_OBJ)
	$(FC) -o $(BINPATH)/tise_1D_morse  $(tise_1D_morse_OBJ) $(LIBS) 
#######################
tise_1D_pteller: $(tise_1D_pteller_OBJ)
	$(FC) -o $(BINPATH)/tise_1D_pteller  $(tise_1D_pteller_OBJ) $(LIBS) 
#######################
tise_1D_wsaxon: $(tise_1D_wsaxon_OBJ)
	$(FC) -o $(BINPATH)/tise_1D_wsaxon  $(tise_1D_wsaxon_OBJ) $(LIBS) 
#######################
clean:
	@rm -f $(tise_1D_morse_OBJ) $(tise_1D_pteller_OBJ) $(tise_1D_wsaxon_OBJ) *.mod




