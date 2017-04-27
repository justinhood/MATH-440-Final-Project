FC90 = mpif90

main_file = par2.f90 Mod2.f90
main_file_implicit = par3.f90 Mod2.f90

all:par3_exe # par2_exe par3_exe

#par2_exe:$(main_file)
#	$(FC90)	$(main_file)	-o	$@
par3_exe:$(main_file_implicit)
	$(FC90)	$(main_file_implicit)	-o	$@

clean:
	rm *_exe *.mod
