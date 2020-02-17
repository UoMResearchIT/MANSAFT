FC = gfortran
INTERFACE_BLOCKS = /home/mbexegc2/NAG/fll6i26dfl/nag_interface_blocks
NAG_LIBRARY = /home/mbexegc2/NAG/fll6i26dfl/lib/libnag_nag.a
EXE = optimiser.exe
OFILES = Types.o Global.o \
	${TG_OBJS} ${Z_OBJS} Zig.o Input_opt.o Z_eff.o  \
	Diff.o Pressure.o ChemPot.o Vol.o
TG_OBJS = GL.o Input.o Setup.o Z_eff.o Ideal.o Assoc.o Ion.o
Z_OBJS = Mono.o Chain.o
COMMON_DEPS = Types.o Global.o Setup.o Ideal.o Mono.o Assoc.o Ion.o Diff.o

.PHONY: run
run: ${EXE}
	./$< dmm_methanol_opt_assoc.in opt_dmm_methanol_assoc.in | tee output.txt

optimiser.exe: opt_vle.f90 ${OFILES}
	${FC} -I${INTERFACE_BLOCKS} $^ ${NAG_LIBRARY} -o $@

Types.o: Types.f90
	${FC} -c $<

Global.o: Global.f90 Types.o
	${FC} -c $<

Zig.o: Zig.f90 Types.o
	${FC} -c $<

Input_opt.o: Input_opt.f90 Types.o Global.o Zig.o
	${FC} -c $<

Diff.o: Diff.f90 Types.o Global.o Setup.o
	${FC} -c $<

Pressure.o: Pressure.f90 ${COMMON_DEPS}
	${FC} -c $<

ChemPot.o: ChemPot.f90 ${COMMON_DEPS} Chain.o
	${FC} -c $<

Vol.o: Vol.f90 ${COMMON_DEPS} Pressure.o Chain.o ChemPot.o
	${FC} -c $<

# Rules for objects with common dependency lists
${TG_OBJS}: %.o: %.f90 Types.o Global.o
	${FC} -c $<

${Z_OBJS}: %.o : %.f90 Types.o Global.o Z_eff.o
	${FC} -c $<

.PHONY: clean
clean:
	rm -f *.mod
	rm -f *.o
	rm -f ${EXE}
