FC = gfortran

.PHONY: run
run: mansaft.exe
	./mansaft.exe ../NAG/dmm_methanol.in | tee output.txt

mansaft.exe: mansaft.f90 minpack.o Global_mod.o Types_mod.o GL_mod.o Input_mod.o Setup_mod.o Zeff_mod.o
	${FC} $^ -o $@

minpack.o: minpack.f dogleg.f dpmpar.f enorm.f fdjac1.f hybrd1.f hybrd.f qform.f qrfac.f r1mpyq.f r1updt.f
	${FC} -c minpack.f

Types_mod.o: Types_mod.f95
	${FC} -c $<

Global_mod.o: Global_mod.f95 Types_mod.o
	${FC} -c $<

GL_mod.o: GL_mod.f95 Global_mod.o Types_mod.o
	${FC} -c $<

Input_mod.o: Input_mod.f95 Global_mod.o Types_mod.o GL_mod.o
	${FC} -c $<

Setup_mod.o: Setup_mod.f95 Global_mod.o Types_mod.o
	${FC} -c $<

Zeff_mod.o: Zeff_mod.f95 Global_mod.o Types_mod.o
	${FC} -c $<
	
.PHONY: clean
clean:
	rm -f *.mod *.o
