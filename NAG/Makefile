interface_blocks=/home/mbexegc2/NAG/fll6i26dfl/nag_interface_blocks 
nag_library=/home/mbexegc2/NAG/fll6i26dfl/lib/libnag_nag.a
exe=mansaft.exe

.PHONY: run
run: ${exe}
	./${exe} dmm_methanol.in 2>&1 | tee output.txt

.PHONY: compile
compile: ${exe}

${exe}: NAG_for_binary_system.f90
	gfortran -I${interface_blocks} $^ ${nag_library} -lstdc++ -o $@

.PHONY: clean
clean:
	rm -f ${exe}
	rm -f *.mod
