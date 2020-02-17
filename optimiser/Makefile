FC = gfortran
INTERFACE_BLOCKS = /home/mbexegc2/NAG/fll6i26dfl/nag_interface_blocks
NAG_LIBRARY = /home/mbexegc2/NAG/fll6i26dfl/lib/libnag_nag.a
EXE = optimiser.exe
OFILES = Types.o Global.o GL.o Input.o Zig.o Input_opt.o Setup.o Z_eff.o Ideal.o \
         Mono.o Chain.o Assoc.o Ion.o

.PHONY: run
run: ${EXE}
	./$< dmm_methanol_opt_assoc.in opt_dmm_methanol_assoc.in | tee output.txt

optimiser.exe: opt_vle.f90 ${OFILES}
	${FC} -I${INTERFACE_BLOCKS} $^ ${NAG_LIBRARY} -o $@

Types.o: Types.f90
	${FC} -c $<

Global.o: Global.f90 Types.o
	${FC} -c $<

GL.o: GL.f90 Types.o Global.o
	${FC} -c $<

Input.o: Input.f90 Types.o Global.o
	${FC} -c $<

Zig.o: Zig.f90 Types.o
	${FC} -c $<

Input_opt.o: Input_opt.f90 Types.o Global.o Zig.o
	${FC} -c $<

Setup.o: Setup.f90 Global.o Types.o
	${FC} -c $<

Z_eff.o: Z_eff.f90 Global.o Types.o
	${FC} -c $<

Ideal.o: Ideal.f90 Global.o Types.o
	${FC} -c $<

Mono.o: Mono.f90 Types.o Global.o Z_eff.o
	${FC} -c $<

Chain.o: Chain.f90 Types.o Global.o Z_eff.o
	${FC} -c $<

Assoc.o: Assoc.f90 Types.o Global.o
	${FC} -c $<

Ion.o: Ion.f90 Types.o Global.o
	${FC} -c $<

.PHONY: clean
clean:
	rm -f *.mod
	rm -f *.o
	rm -f ${EXE}
