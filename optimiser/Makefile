FC = gfortran
NAG_LIB_DIR = /home/mbexegc2/NAG/fll6i26dfl
INTERFACE_BLOCKS = ${NAG_LIB_DIR}/nag_interface_blocks
NAG_LIBRARY = ${NAG_LIB_DIR}/lib/libnag_nag.a
EXE = optimiser.exe
OBJECT_FILES = Types.o Global.o \
	${TG_OBJS} ${Z_OBJS} Zig.o Input_opt.o Z_eff.o  \
	Diff.o Pressure.o ChemPot.o Vol.o Quote.o Input.o
TG_OBJS = GL.o Setup.o Z_eff.o Ideal.o Assoc.o Ion.o
Z_OBJS = Mono.o Chain.o
COMMON_DEPS = Types.o Global.o Setup.o Ideal.o Mono.o Chain.o Assoc.o Ion.o Diff.o
MP_DIR = ../minpack

.PHONY: run
run: ${EXE}
	./$< dmm_methanol_opt_assoc.in opt_dmm_methanol_assoc.in | tee output.txt

optimiser.exe: optimiser.f90 ${OBJECT_FILES}
	${FC} -I${INTERFACE_BLOCKS} $^ ${NAG_LIBRARY} -o $@

Types.o: Types.f90
	${FC} -c $<

Global.o: Global.f90 Types.o
	${FC} -c $<

Zig.o: Zig.f90 Types.o
	${FC} -c $<

Input_opt.o: Input_opt.f90 Types.o Global.o Zig.o
	${FC} -c $<

Quote.o: Quote.f90
	${FC} -c $<

Input.o: Input.f90 Types.o Global.o
	${FC} -c $<


# Reuse module source code from minpack directory
Diff.o: ${MP_DIR}/Diff.f90 Types.o Global.o Setup.o
	${FC} -c $<

Pressure.o: ${MP_DIR}/Pressure.f90 ${COMMON_DEPS}
	${FC} -c $<

ChemPot.o: ${MP_DIR}/ChemPot.f90 ${COMMON_DEPS}
	${FC} -c $<

Vol.o: ${MP_DIR}/Vol.f90 ${COMMON_DEPS} Pressure.o ChemPot.o
	${FC} -c $<

# Rules for objects with common dependency lists
${TG_OBJS}: %.o: ${MP_DIR}/%.f90 Types.o Global.o
	${FC} -c $<

${Z_OBJS}: %.o : ${MP_DIR}/%.f90 Types.o Global.o Z_eff.o
	${FC} -c $<

.PHONY: clean
clean:
	rm -f *.mod
	rm -f *.o
	rm -f ${EXE}
