FC = gfortran
INTERFACE_BLOCKS = /home/mbexegc2/NAG/fll6i26dfl/nag_interface_blocks
NAG_LIBRARY = /home/mbexegc2/NAG/fll6i26dfl/lib/libnag_nag.a
EXE = optimiser.exe
OFILES = Types.o

.PHONY: run
run: ${EXE}
	./$< dmm_methanol_opt_assoc.in opt_dmm_methanol_assoc.in | tee output.txt

optimiser.exe: opt_vle.f90 ${OFILES}
	${FC} -I${INTERFACE_BLOCKS} $^ ${NAG_LIBRARY} -o $@

Types.o: Types.f90
	${FC} -c $<

.PHONY: clean
clean:
	rm -f *.mod
	rm -f *.o
	rm -f ${EXE}
