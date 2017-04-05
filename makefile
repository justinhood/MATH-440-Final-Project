FC90 = gfortran

main_file = initialserial.f90 initialMod.f90


all: main_exe 

main_exe:$(main_file)
	$(FC90)	$(main_file)	-o	$@

clean:
	rm *_exe *.mod
