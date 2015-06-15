#compile the emulator
LIBS = -llapack -lblas -L/usr/lib64/R/lib 
# switch the 2 comments to switch between linux and windows
# FCOMPILER = gfortran
# FLAGS = -c -O3 -fPIC
FLAGS = -c -O3
FCOMPILER = i686-w64-mingw32-gfortran
emulator: service_functions.o lin_mod.o emulator_f.o
	$(FCOMPILER) -shared -L/usr/lib64/R/lib -Wl,-O3 -o emulator_f.o $^ -o emulator_f.dll $(LIBS)
service_functions.o: ../adli_common/service_functions.f90
	$(FCOMPILER) $(FLAGS) ../adli_common/service_functions.f90
lin_mod.o: ../adli_common/lin_mod.f90
	$(FCOMPILER) $(FLAGS) ../adli_common/lin_mod.f90
emulator_f.o: emulator_f.f90
	$(FCOMPILER) $(FLAGS) emulator_f.f90
clean:
	rm *.mod *.o
#end of makefile
