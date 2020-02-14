INTERFACE_BLOCKS=/home/mbexegc2/NAG/fll6i26dfl/nag_interface_blocks
NAG_LIBRARY=/home/mbexegc2/NAG/fll6i26dfl/lib/libnag_nag.a

optimiser.exe: opt_vle.f90
	gfortran -I${INTERFACE_BLOCKS} $< ${NAG_LIBRARY} -o $@

.PHONY: run
run: optimiser.exe
	./optimiser.exe dmm_methanol_opt_assoc.in opt_dmm_methanol_assoc.in

.PHONY: clean
clean:
	rm -f *.mod
	rm -f *.o
