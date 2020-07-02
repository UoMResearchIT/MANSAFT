FC = gfortran
MANSAFT_OBJECT_FILES = minpack.o Global.o Types.o  Input.o  \
                       ${TG_OBJECT_FILES} ${Z_OBJECT_FILES} \
                       Diff.o Pressure.o ChemPot.o Vol.o solver.o pure_phase.o
MINPACK_DEPENDENCIES = dogleg.f dpmpar.f enorm.f fdjac1.f hybrd1.f hybrd.f qform.f qrfac.f r1mpyq.f r1updt.f
TG_OBJECT_FILES = GL.o Setup.o Z_eff.o Ideal.o Assoc.o Ion.o
Z_OBJECT_FILES = Mono.o Chain.o
COMMON_DEPENDENCIES = Types.o Global.o Setup.o Ideal.o Mono.o Assoc.o Ion.o Diff.o
NL_DIR = ${HOME}/Downloads/nlopt
NLOPT = -I${NL_DIR}/include -L${NL_DIR}/lib -lnlopt -lm 

.PHONY: all
all: mansaft.exe optimiser.exe

.PHONY: run_optimiser_new_example
run_optimiser_new_example: optimiser.exe 
	./$< dmm_methanol_opt_assoc.in opt_dmm_methanol_assoc.in 2>&1 | tee optimiser_output_new_example.txt

mansaft.exe: mansaft.f90 ${MANSAFT_OBJECT_FILES}
	${FC}  $^ -o $@ 

optimiser.exe: optimiser.f90 ${MANSAFT_OBJECT_FILES} Input_opt.o nlopt.o
	${FC} $^ ${NLOPT} -o $@ 

minpack.o: minpack.f ${MINPACK_DEPENDENCIES}
	${FC} -c $<

Types.o: Types.f90
	${FC} -c $<

Global.o: Global.f90 Types.o
	${FC} -c $<

Input.o: Input.f90 Types.o Global.o GL.o
	${FC} -c $<

GL.o: GL.f90 Types.o Global.o
	${FC} -c $<

Setup.o: Setup.f90 Types.o Global.o
	${FC} -c $<

Z_eff.o: Z_eff.f90 Types.o Global.o
	${FC} -c $<

Ideal.o: Ideal.f90 Types.o Global.o
	${FC} -c $<
	
Mono.o: Mono.f90 Types.o Global.o Z_eff.o
	${FC} -c $<
	
Chain.o: Chain.f90 Types.o Global.o Z_eff.o
	${FC} -c $<

Assoc.o: Assoc.f90 Types.o Global.o
	${FC} -c $<

Ion.o: Ion.f90 Types.o Global.o
	${FC} -c $<

Diff.o: Diff.f90 Types.o Global.o Setup.o
	${FC} -c $<

Zig.o: Zig.f90 Types.o
	${FC} -c $<

Pressure.o: Pressure.f90 Types.o Global.o Setup.o Ideal.o Mono.o Chain.o Assoc.o Ion.o Diff.o
	${FC} -c $<

ChemPot.o: ChemPot.f90 Types.o Global.o Setup.o Ideal.o Mono.o Chain.o Assoc.o Ion.o Diff.o
	${FC} -c $<

Vol.o: Vol.f90 Types.o Global.o Setup.o Pressure.o Ideal.o Mono.o Chain.o Assoc.o Ion.o ChemPot.o
	${FC} -c $<

solver.o: solver.f90 Types.o Global.o Pressure.o ChemPot.o Vol.o
	${FC} -c $<
	
pure_phase.o: pure_phase.f90 Types.o Global.o Pressure.o ChemPot.o Vol.o
	${FC} -c $<

Input_opt.o: Input_opt.f90 Types.o Global.o Zig.o
	${FC} -c $<

nlopt.o: nlopt.f90 Types.o Global.o Pressure.o ChemPot.o Vol.o Input.o Input_opt.o solver.o
	${FC} -c $< ${NLOPT}

.PHONY: clean
clean:
	rm -f *.mod *.o *.exe
