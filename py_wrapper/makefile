#compile the emulator
FLAGS = -c -O3 -fPIC
LIBS = -llapack -lblas
FCOMPILER = gfortran
emulator: service_functions.o lin_mod.o
	f2py -c -m emulib $(LIBS) sig.pyf emulator_f.f90 $^
service_functions.o: ../core/service_functions.f90
	$(FCOMPILER) $(FLAGS) ../core/service_functions.f90 $(LIBS)
lin_mod.o: ../core/lin_mod.f90
	$(FCOMPILER) $(FLAGS) ../core/lin_mod.f90 $(LIBS)
emulator_f.o: emulator_f.f90
	$(FCOMPILER) $(FLAGS) emulator_f.f90 $(LIBS)
clean:
	rm *.mod *.o
signature:
	f2py emulator_f.f90 -m emulib -h sig.pyf --overwrite-signature
#end of makefile
