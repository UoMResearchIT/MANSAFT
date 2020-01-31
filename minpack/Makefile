.PHONY: run
run: mansaft.exe
	./mansaft.exe ../NAG/dmm_methanol.in | tee output.txt

mansaft.exe: mansaft.f90 minpack.o
	gfortran $< minpack.o  -o $@

minpack.o: minpack.f dogleg.f dpmpar.f enorm.f fdjac1.f hybrd1.f hybrd.f qform.f qrfac.f r1mpyq.f r1updt.f
	gfortran -c minpack.f

.PHONY: clean
clean:
	rm -f *.mod *.o
