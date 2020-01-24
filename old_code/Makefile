FF		= 	gfortran 
FFLAGS  =   
LDFLAGS = 
OPT		=   -o3
DEBUG	= 	-Og -g -Wall -Wextra -pedantic -fbounds-check -fbacktrace -fcheck=all 

all: mansaft.exe optimiser.exe

mansaft.exe: Input_mod.o Global_mod.o Setup_mod.o Types_mod.o Press_mod.o Vol_mod.o Mu_mod.o Ideal_mod.o Mono_mod.o Chain_mod.o Assoc_mod.o Ion_mod.o Zeff_mod.o Act_ion_mod.o Pure_phase_mod.o GL_mod.o Zig_mod.o Simplex_mod.o Error_mod.o Quote_mod.o Simplex_mod.o Diff_mod.o Therm_mod.o Phase_simplex_mod.o Main.o
	$(FF) $(FFLAGS) $(DEBUG) $^ $(LDFLAGS) -o $@

optimiser.exe: Input_mod.o Input_opt_mod.o Global_mod.o Zig_mod.o Simplex_mod.o Error_mod.o Setup_mod.o Types_mod.o Press_mod.o Vol_mod.o Mu_mod.o Ideal_mod.o Mono_mod.o Chain_mod.o Assoc_mod.o Ion_mod.o Zeff_mod.o GL_mod.o Zig_mod.o Quote_mod.o Pure_phase_mod.o Diff_mod.o Therm_mod.o Act_ion_mod.o Optimiser.o
	$(FF) $(FFLAGS) $(OPT)  $^ $(LDFLAGS) -o $@
	
Main.o: Input_mod.o Global_mod.o Setup_mod.o Types_mod.o Press_mod.o Vol_mod.o Act_ion_mod.o Pure_phase_mod.o Quote_mod.o Simplex_mod.o Therm_mod.o Phase_simplex_mod.o Main.f95 $^
	$(FF) -c Main.f95

Optimiser.o: Input_mod.o Input_opt_mod.o Global_mod.o Simplex_mod.o Setup_mod.o Types_mod.o Press_mod.o Vol_mod.o Simplex_mod.o Zig_mod.o Therm_mod.o Quote_mod.o Optimiser.f95	
	$(FF) -c Optimiser.f95
	
Press_mod.o: Global_mod.o Types_mod.o Setup_mod.o Ideal_mod.o Mono_mod.o Chain_mod.o Assoc_mod.o Ion_mod.o Diff_mod.o Press_mod.f95 $^
	$(FF) -c Press_mod.f95

Mu_mod.o: Global_mod.o Types_mod.o Setup_mod.o Ideal_mod.o Mono_mod.o Chain_mod.o Assoc_mod.o Ion_mod.o Diff_mod.o Mu_mod.f95 $^
	$(FF) -c Mu_mod.f95
		
Vol_mod.o: Global_mod.o Types_mod.o Setup_mod.o Ideal_mod.o Mono_mod.o Chain_mod.o Assoc_mod.o Ion_mod.o Press_mod.o Mu_mod.o Vol_mod.f95 $^
	$(FF) -c Vol_mod.f95
	
Types_mod.o: Types_mod.f95 $^
	$(FF) -c Types_mod.f95
	
Global_mod.o: Types_mod.o Global_mod.f95 $^
	$(FF) -c Global_mod.f95

Zeff_mod.o: Types_mod.o Global_mod.o Zeff_mod.f95 $^
	$(FF) -c Zeff_mod.f95
	
Input_mod.o: Global_mod.o Types_mod.o GL_mod.o Input_mod.f95 $^
	$(FF) -c Input_mod.f95

Input_opt_mod.o: Global_mod.o Types_mod.o Zig_mod.o GL_mod.o Input_opt_mod.f95 $^
	$(FF) -c Input_opt_mod.f95
	
Setup_mod.o: Global_mod.o Types_mod.o Setup_mod.f95 $^
	$(FF) -c Setup_mod.f95
	
Ideal_mod.o: Global_mod.o Types_mod.o Ideal_mod.f95 $^
	$(FF) -c Ideal_mod.f95
	
Mono_mod.o: Global_mod.o Types_mod.o Zeff_mod.o Mono_mod.f95 $^
	$(FF) -c Mono_mod.f95

Chain_mod.o: Global_mod.o Types_mod.o Zeff_mod.o Chain_mod.f95 $^
	$(FF) -c Chain_mod.f95

Assoc_mod.o: Global_mod.o Types_mod.o Assoc_mod.f95 $^
	$(FF) -c Assoc_mod.f95

Ion_mod.o: Global_mod.o Types_mod.o Ion_mod.f95 $^
	$(FF) -c Ion_mod.f95
	
GL_mod.o: Global_mod.o Types_mod.o GL_mod.f95 $^
	$(FF) -c GL_mod.f95	
	
Zig_mod.o: Types_mod.o Zig_mod.f95
	$(FF) -c Zig_mod.f95

Simplex_mod.o:	Types_mod.o Error_mod.o Simplex_mod.f95
	$(FF) -c Simplex_mod.f95

Error_mod.o:	Global_mod.o Setup_mod.o Types_mod.o Press_mod.o Vol_mod.o Mu_mod.o Ideal_mod.o Mono_mod.o Chain_mod.o Assoc_mod.o Ion_mod.o Zeff_mod.o Pure_phase_mod.o Therm_mod.o Act_ion_mod.o Error_mod.f95
	$(FF) -c Error_mod.f95

Act_ion_mod.o:	Global_mod.o Setup_mod.o Types_mod.o Vol_mod.o Ideal_mod.o Mono_mod.o Chain_mod.o Assoc_mod.o Ion_mod.o Mu_mod.o Act_ion_mod.f95
	$(FF) -c Act_ion_mod.f95
		
Pure_phase_mod.o: Global_mod.o Setup_mod.o Types_mod.o Ideal_mod.o Mono_mod.o Chain_mod.o Assoc_mod.o Ion_mod.o Press_mod.o Mu_mod.o Pure_phase_mod.f95
	$(FF) -c Pure_phase_mod.f95		
	
Therm_mod.o:	Global_mod.o Setup_mod.o Types_mod.o Vol_mod.o Press_mod.o Ideal_mod.o Mono_mod.o Chain_mod.o Assoc_mod.o Ion_mod.o Therm_mod.f95
	$(FF) -c Therm_mod.f95
		
Diff_mod.o: Global_mod.o Types_mod.o Setup_mod.o Diff_mod.f95 $^
	$(FF) -c Diff_mod.f95		

Phase_simplex_mod.o: Global_mod.o Types_mod.o Zig_mod.o Mu_mod.o Press_mod.o Vol_mod.o Ideal_mod.o Mono_mod.o Chain_mod.o Assoc_mod.o Ion_mod.o Phase_simplex_mod.f95
	$(FF) -c Phase_simplex_mod.f95			
	
Quote_mod.o: Zig_mod.o Quote_mod.f95 $^
	$(FF) -c Quote_mod.f95	
	
.PHONY: clean
clean:
	rm -rf *.o *.exe *.mod *.stackdump
