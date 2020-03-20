FC = gfortran
OBJECT_FILES = minpack.o Global.o Types.o GL.o Input.o Setup.o \
               Z_eff.o Ideal.o Mono.o Chain.o Assoc.o Ion.o \
               Diff.o Pressure.o ChemPot.o Vol.o solver.o
MINPACK_DEPENDENCIES = dogleg.f dpmpar.f enorm.f fdjac1.f hybrd1.f hybrd.f qform.f qrfac.f r1mpyq.f r1updt.f
STD_O_FILES = Types.o Global.o 

ifeq (${FC}, nagfor)
  OPT = -C=all -u -C=undefined -g -gline -mtrace=on,line,paranoia -Warn=allocation -Warn=subnormal -ieee=stop
endif

.PHONY: run
run: mansaft.exe
	./mansaft.exe ../NAG/dmm_methanol.in | tee output.txt 

.PHONY: diff
diff: output.txt
	@echo "Diff against output from NAG code:"
	@diff output.txt ../NAG/output.txt

output.txt: run

mansaft.exe: mansaft.f90 ${OBJECT_FILES}
	${FC} ${OPT} $^ -o $@

minpack.o: minpack.f ${MINPACK_DEPENDENCIES}
	${FC} ${OPT} -c $<

Types.o: Types.f90
	${FC} ${OPT} -c $<

Global.o: Global.f90 Types.o
	${FC} ${OPT} -c $<

GL.o: GL.f90 ${STD_O_FILES}
	${FC} ${OPT} -c $<

Input.o: Input.f90 ${STD_O_FILES} GL.o
	${FC} ${OPT} -c $<

Setup.o: Setup.f90 ${STD_O_FILES}
	${FC} ${OPT} -c $<

Z_eff.o: Z_eff.f90 ${STD_O_FILES}
	${FC} ${OPT} -c $<

Ideal.o: Ideal.f90 ${STD_O_FILES}
	${FC} ${OPT} -c $<
	
Mono.o: Mono.f90 ${STD_O_FILES} Z_eff.o
	${FC} ${OPT} -c $<
	
Chain.o: Chain.f90 ${STD_O_FILES} Z_eff.o
	${FC} ${OPT} -c $<

Assoc.o: Assoc.f90 ${STD_O_FILES}
	${FC} ${OPT} -c $<

Ion.o: Ion.f90 ${STD_O_FILES}
	${FC} ${OPT} -c $<

Diff.o: Diff.f90 ${STD_O_FILES} Setup.o
	${FC} ${OPT} -c $<

Pressure.o: Pressure.f90 ${STD_O_FILES} Setup.o Ideal.o Mono.o Chain.o Assoc.o Ion.o Diff.o
	${FC} ${OPT} -c $<

ChemPot.o: ChemPot.f90 ${STD_O_FILES} Setup.o Ideal.o Mono.o Chain.o Assoc.o Ion.o Diff.o
	${FC} ${OPT} -c $<

Vol.o: Vol.f90 ${STD_O_FILES} Setup.o Pressure.o Ideal.o Mono.o Chain.o Assoc.o Ion.o ChemPot.o
	${FC} ${OPT} -c $<

solver.o: solver.f90 ${STD_O_FILES} Pressure.o ChemPot.o Vol.o
	${FC} ${OPT} -c $<

.PHONY: clean
clean:
	rm -f *.mod *.o
	rm -f output.txt
