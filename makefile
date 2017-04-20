FC90 = mpif90

main_file = par2.f90 Mod2.f90


all: par2_exe 

par2_exe:$(main_file)
	$(FC90)	$(main_file)	-o	$@

clean:
	rm *_exe *.mod
