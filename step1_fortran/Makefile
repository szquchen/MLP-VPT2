FC := ifx -r8 -O
LIBS := -qmkl=sequential

OBJ := \
pes_basis.o  pes_shell.o \
util_opt.o  util_paxis.o \
util_nma.o  util_force.o \
vpt2.o 

# Targets:
.PHONY : clean 

%.o : %.f90
	$(FC) -c $(FFLAGS) $< 

vpt2.x :  $(OBJ)
	$(FC) -o $@  $^ $(FFLAGS) $(LIBS)

clean :
	rm *.o *.mod 
