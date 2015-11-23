#compile the emulator
FLAGS = -c -O3 -fPIC
LIBS = -llapack -lblas
FCOMPILER = gfortran
F2PYCOMPILER= f2py
emulator: service_functions.o lin_mod_changed.o symmlq.o
	$(F2PYCOMPILER) -c -m emulib $(LIBS) sig.pyf emulator_f.f90 $^
service_functions.o: ../core/service_functions.f90
	$(FCOMPILER) $(FLAGS) ../core/service_functions.f90 $(LIBS)
lin_mod_changed.o: ../core/lin_mod_changed.f90
	python parseConfig.py $(CONFIG_FILE)
	$(FCOMPILER) $(FLAGS) $^ $(LIBS)
emulator_f.o: emulator_f.f90
	$(FCOMPILER) $(FLAGS) $^ $(LIBS)
symmlq.o: ../core/symmlq.f90
	$(FCOMPILER) $(FLAGS) $^ $(LIBS)
clean:
	rm *.mod *.o
signature:
	$(F2PYCOMPILER) emulator_f.f90 -m emulib -h sig.pyf --overwrite-signature
#end of makefile
