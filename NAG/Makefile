interface_blocks=/home/mbexegc2/NAG/fll6i26dfl/nag_interface_blocks 
src=NAG_for_binary_system
nag_library=/home/mbexegc2/NAG/fll6i26dfl/lib/libnag_nag.a 

${src}.exe: ${src}.f90
	gfortran -I${interface_blocks} ${src}.f90 ${nag_library} -lstdc++ -o ${src}.exe
