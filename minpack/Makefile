INTERFACE_BLOCKS=/home/mbexegc2/NAG/fll6i26dfl/nag_interface_blocks
NAG_LIBRARY=/home/mbexegc2/NAG/fll6i26dfl/lib/libnag_nag.a

mansaft.exe: mansaft.f90 minpack.o
	gfortran -I${INTERFACE_BLOCKS} $< minpack.o ${NAG_LIBRARY} -lstdc++ -o $@

minpack.o: minpack.f dogleg.f dpmpar.f enorm.f fdjac1.f hybrd1.f hybrd.f qform.f qrfac.f r1mpyq.f r1updt.f
	gfortran -c minpack.f

.PHONY: run
run: mansaft.exe
	./mansaft.exe ../NAG/dmm_methanol.in

.PHONY: clean
clean:
	rm -f *.mod *.o
