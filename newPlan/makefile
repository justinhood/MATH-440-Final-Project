FC90 = mpif90

main_file = noScatter.f90 Mod2.f90

all:	noScatter_exe 

noScatter_exe:$(main_file)
	$(FC90)	$(main_file)	-o	$@

clean:
	rm *_exe *.mod
