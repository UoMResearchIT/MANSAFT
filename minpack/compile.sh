gfortran -c minpack.f
gfortran -I/home/mbexegc2/NAG/fll6i26dfl/nag_interface_blocks minpack.o mansaft.f90 /home/mbexegc2/NAG/fll6i26dfl/lib/libnag_nag.a -lstdc++ -o mansaft.exe
./mansaft.exe ../NAG/dmm_methanol.in
