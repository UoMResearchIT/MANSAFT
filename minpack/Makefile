FC = gfortran
OBJECT_FILES = minpack.o Global_mod.o Types_mod.o GL_mod.o Input_mod.o Setup_mod.o \
               Zeff_mod.o Ideal_mod.o Mono_mod.o Chain_mod.o Assoc_mod.o Ion_mod.o \
               Diff_mod.o
MINPACK_DEPENDENCIES = dogleg.f dpmpar.f enorm.f fdjac1.f hybrd1.f hybrd.f qform.f qrfac.f r1mpyq.f r1updt.f
STD_MOD_FILES = Types_mod.o Global_mod.o 

.PHONY: run
run: mansaft.exe
	./mansaft.exe ../NAG/dmm_methanol.in | tee output.txt

mansaft.exe: mansaft.f90 ${OBJECT_FILES}
	${FC} $^ -o $@

minpack.o: minpack.f ${MINPACK_DEPENDENCIES}
	${FC} -c $<

Types_mod.o: Types_mod.f95
	${FC} -c $<

Global_mod.o: Global_mod.f95 Types_mod.o
	${FC} -c $<

GL_mod.o: GL_mod.f95 ${STD_MOD_FILES}
	${FC} -c $<

Input_mod.o: Input_mod.f95 ${STD_MOD_FILES} GL_mod.o
	${FC} -c $<

Setup_mod.o: Setup_mod.f95 ${STD_MOD_FILES}
	${FC} -c $<

Zeff_mod.o: Zeff_mod.f95 ${STD_MOD_FILES}
	${FC} -c $<

Ideal_mod.o: Ideal_mod.f95 ${STD_MOD_FILES}
	${FC} -c $<
	
Mono_mod.o: Mono_mod.f95 ${STD_MOD_FILES} Zeff_mod.o
	${FC} -c $<
	
Chain_mod.o: Chain_mod.f95 ${STD_MOD_FILES} Zeff_mod.o
	${FC} -c $<

Assoc_mod.o: Assoc_mod.f95 ${STD_MOD_FILES}
	${FC} -c $<

Ion_mod.o: Ion_mod.f95 ${STD_MOD_FILES}
	${FC} -c $<

Diff_mod.o: Diff_mod.f95 ${STD_MOD_FILES} Setup_mod.o
	${FC} -c $<

.PHONY: clean
clean:
	rm -f *.mod *.o
