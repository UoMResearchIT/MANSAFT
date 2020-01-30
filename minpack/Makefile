INTERFACE_BLOCKS=/home/mbexegc2/NAG/fll6i26dfl/nag_interface_blocks
NAG_LIBRARY=/home/mbexegc2/NAG/fll6i26dfl/lib/libnag_nag.a

mansaft.exe: mansaft.f90 dpmpar.f
	gfortran -I${INTERFACE_BLOCKS} $< dpmpar.f ${NAG_LIBRARY} -lstdc++ -o $@

.PHONY: run
run: 
	./mansaft.exe ../NAG/dmm_methanol.in

.PHONY: clean
clean:
	rm -f *.mod *.o
