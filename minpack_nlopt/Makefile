FC = gfortran
OBJECT_FILES = minpack.o Global.o Types.o GL.o Input.o Setup.o \
               Z_eff.o Ideal.o Mono.o Chain.o Assoc.o Ion.o \
               Diff.o Pressure.o ChemPot.o Vol.o solver.o pure_phase.o
MINPACK_DEPENDENCIES = dogleg.f dpmpar.f enorm.f fdjac1.f hybrd1.f hybrd.f qform.f qrfac.f r1mpyq.f r1updt.f
STD_O_FILES = Types.o Global.o
NL_DIR = ${HOME}/Downloads/nlopt
NLOPT = -I${NL_DIR}/include -L${NL_DIR}/lib -lnlopt -lm 

all: mansaft.exe optimiser.exe

mansaft.exe: mansaft.f90 ${OBJECT_FILES}
	${FC}  $^ -o $@ 

optimiser.exe: optimiser.f90 ${OBJECT_FILES} Input_opt.o nlopt.o
	${FC} ${NLOPT} $^ -o $@ 

minpack.o: minpack.f ${MINPACK_DEPENDENCIES}
	${FC} -c $<

Types.o: Types.f90
	${FC} -c $<

Global.o: Global.f90 Types.o
	${FC} -c $<

GL.o: GL.f90 ${STD_O_FILES}
	${FC} -c $<

Input.o: Input.f90 ${STD_O_FILES} GL.o
	${FC} -c $<

Setup.o: Setup.f90 ${STD_O_FILES}
	${FC} -c $<

Z_eff.o: Z_eff.f90 ${STD_O_FILES}
	${FC} -c $<

Ideal.o: Ideal.f90 ${STD_O_FILES}
	${FC} -c $<
	
Mono.o: Mono.f90 ${STD_O_FILES} Z_eff.o
	${FC} -c $<
	
Chain.o: Chain.f90 ${STD_O_FILES} Z_eff.o
	${FC} -c $<

Assoc.o: Assoc.f90 ${STD_O_FILES}
	${FC} -c $<

Ion.o: Ion.f90 ${STD_O_FILES}
	${FC} -c $<

Diff.o: Diff.f90 ${STD_O_FILES} Setup.o
	${FC} -c $<

Zig.o: Zig.f90 Types.o
	${FC} -c $<

Pressure.o: Pressure.f90 ${STD_O_FILES} Setup.o Ideal.o Mono.o Chain.o Assoc.o Ion.o Diff.o
	${FC} -c $<

ChemPot.o: ChemPot.f90 ${STD_O_FILES} Setup.o Ideal.o Mono.o Chain.o Assoc.o Ion.o Diff.o
	${FC} -c $<

Vol.o: Vol.f90 ${STD_O_FILES} Setup.o Pressure.o Ideal.o Mono.o Chain.o Assoc.o Ion.o ChemPot.o
	${FC} -c $<

solver.o: solver.f90 ${STD_O_FILES} Pressure.o ChemPot.o Vol.o
	${FC} -c $<
	
pure_phase.o: pure_phase.f90 ${STD_O_FILES} Pressure.o ChemPot.o Vol.o
	${FC} -c $<

Input_opt.o: Input_opt.f90 ${STD_O_FILES} Zig.o
	${FC} -c $<

nlopt.o: nlopt.f90 ${STD_O_FILES} Pressure.o ChemPot.o Vol.o Input.o Input_opt.o solver.o
	${FC} ${NLOPT} -c $< 

.PHONY: clean
clean:
	rm -f *.mod *.o *.exe
