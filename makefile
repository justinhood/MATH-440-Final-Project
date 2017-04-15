FC90 = mpif90

main_file = initialpar.f90 initialMod.f90


all: main_exe 

main_exe:$(main_file)
	$(FC90)	$(main_file)	-o	$@

clean:
	rm *_exe *.mod
