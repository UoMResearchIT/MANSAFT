!***************************************************************************************************
!
!       MANSAFT MAIN
!       ============
!
!       Version 2.7
!***************************************************************************************************
!   SAFT gamma mie program to:
!       1. Read general parameter list form input and
!       2. Calculate system parameters
!       3. Optimise parameters according to SIMPLEX method
!***************************************************************************************************
!   
!   Authors:
!       SJ Halstead (simon.halstead@manchester.ac.uk)
!       AJ Masters
!   
!       Contributions & advice from Durham & Imperial
!
!***************************************************************************************************
!   Created by SJ Halstead Oct 2015
!       Update 2  Jan 2016
!           - tidied and standardised units
!       Update 3  Feb 2016
!           - added simple P / V calculator
!       Update 4  Mar 2016
!           - standardised MONO units
!       Update 5  Mar 2016
!           - reorganised to focus on propety calculation
!       Update 6  May 2016
!           - property calculations enabled for:
!               - P/V (able to select correct V using lowest G)
!               - Chem pot
!               - Phase equilibrium
!       Update 7  June 2016
!           - tided and efficiency improved (finger crossed)
!           - chemical potential working correclty, 4 point derivative
!       Update 8  September 2016
!           - Reorganised modules and program according to best practice
!       Update 9  October 2016
!           - Ion working w.r.t. density
!           - Activity coefficients added
!       Update 10 March 2017
!           - Numerical Mu values insufficient, updating to all analytical
!           - Analytical pressures
!       Update 23 April 2-17
!           - Mu correct, checked all intrinsics are DP
!       Update 22 Feb 2018
!           - Adding flash calculations for phase equilibrium
!       Update 2 Apr 2018
!           - Simple binary phase equilibrium algorithm added
!***************************************************************************************************
!       REFERENCES
!       ==========
!
!       [1] Accurate statistical associating fluid theory for chain molecules formed from Mie segments
!               J Chem Phys 193, 154504 (2013)
!               T. Lafitte, A. Apostolakou, C. Avendano, A. Galindo, C. Adjiman, E. Muller, G. Jackson
!
!               (Main reference)
!
!
!       [2] Prediction of thermodynamic propertied and phase behavior of fluids and mixtures with the
!           SAFT-gamma Mie group-contribution equation of state
!               J. Chem. Eng. Data, 59, 3272-3288 (2014)
!               S. Dufal, V. Papaioannou, M. Sadeqzadeh, T. Pogiatzis, A. Chremos, C.S. Adjiman, 
!               G. Jackson, A. Galindo
!
!               (Associations)
!
!       
!       [3] The A in SAFT: developing the contribution of association to the Helmholtz free energy
!           within a Wertheim TPT1 treatment of generic Mie fluids
!               Mol. Phys. 113, 948-984 (2015)
!               S. Dufal, T. Lafitte, A.J. Haslam, A. Galindo, G.N.I. Clark, C.Vega, G. Jackson            
!
!               (Association coefficients)            
!
!
!***************************************************************************************************
program Main
!***************************************************************************************************
!Modules
!=======
    use Types_mod           ! Definitions of types and double precision
    use Input_mod           ! Read the input
    use Global_mod          ! Important global parameters
    use Press_mod           ! Pressure function
    use Vol_mod             ! Volume function
    use Act_ion_mod         ! Activity coefficient
    use Pure_phase_mod      ! Pure phase equilibrium
    use Therm_mod           ! Various thermodynamic state values
    use Phase_simplex_mod   ! Phase equilibrium simplex
    use Quote_mod           ! Random quote
!***************************************************************************************************
!Variables
!=========
    implicit none
    
    !Clock
    integer             ::  start, finish, rate
    !Property calculation
    character(len=1)    ::  prop
    real(kind=DP)       ::  vl_out, vv_out, g_out, h_out
    real(kind=DP), allocatable  ::  t_mu_out(:), mu_total_out(:)
    !Counter
    integer             ::  i, j   
!***************************************************************************************************
!Program header
!==============
    print*,
    print*, "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    print*, "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    print*, "++                                                                             ++"
    print*, "++                                                                             ++"
    print*, "++            M     M    A    N   N   SSSS    A    FFFFF  TTTTT                ++"
    print*, "++            MM   MM   A A   NN  N  SS      A A   F        T                  ++"
    print*, "++            M M M M  A   A  N N N   SSS   A   A  FFF      T                  ++"
    print*, "++            M  M  M  A A A  N  NN     SS  A A A  F        T                  ++"
    print*, "++            M     M  A   A  N   N  SSSS   A   A  F        T                  ++"
    print*, "++                                                                             ++"
    print*, "++                                                                             ++"
    print*, "++                              MANSAFT                                        ++"
    print*, "++                     Manchester SAFT Gamma Mie                               ++"
    print*, "++                                                                             ++"
    print*, "++                VERSION 2.7 - The Ant on the Elephant                        ++"
    print*, "++                                2018                                         ++"
    print*, "++                                                                             ++"
    print*, "++                                                                             ++"
    print*, "++                       S.J. Halstead & A.J.Masters                           ++"
    print*, "++                                                                             ++"
    print*, "++         School of Chemical Engineering and Analytical Science               ++"
    print*, "++                     The University of Manchester                            ++"
    print*, "++                                                                             ++"
    print*, "++               contact:  simon.halstead@manchester.ac.uk                     ++"
    print*, "++                         simon.halstead@gmail.com                            ++"
    print*, "++                                                                             ++"
    print*, "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    print*, "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    print*,
    
    call flush( )
!***************************************************************************************************
!Start clock
!===========
    call System_clock( start, rate )
!***************************************************************************************************
!Read input
!==========
    call Read_input( )
!***************************************************************************************************    
!Calculate properties
!====================
        select case(properties%type)
        
            case('p','P')
                
                print*, "Pressure calculation"
                print*, "===================="
                print*,
                print*, " T(K)  V(m^3 mol^-1)  P(MPa)"  
                print*,
                
                do i = 1, properties%n                    
                    t = properties%t(i)
                    v = properties%v(i)                                   
                    p = Press( )
                    write(*, '(f8.2, 2e16.6)') t, v, p / 1.e6_DP
                    call flush( )
                end do
!---------------------------------------------------------------------------------------------------     
            case('px','PX')
                
                print*, "Pressure calculation"
                print*, "===================="
                print*,
                print*, " T(K)  V(m^3 mol^-1)  P(MPa)"  
                print*,
                
                do i = 1, properties%n                    
                    t = properties%t(i)
                    v = properties%v(i)    
                    Comp_array(1:nctypes)%xi = properties%xi(i, 1:nctypes)                               
                    p = Press( )
                    write(*, '(f8.2, 2e16.6)') Comp_array(1:nctypes)%xi, t, v, p / 1.e6_DP
                    call flush( )
                end do
!---------------------------------------------------------------------------------------------------            
            case('opt','OPT','Opt')
                
                print*, 
                print*, "Optimisation requested"
                print*, " - please use Optimiser.exe, not Mansaft.exe "
                print*,   
                print*,               
!---------------------------------------------------------------------------------------------------
            case('eopt','EOPT')
                
                print*, 
                print*, "Electrolyte optimisation requested"
                print*, " - please use Optimiser.exe, not Mansaft.exe "
                print*,   
                print*,               
!---------------------------------------------------------------------------------------------------
            case('ilopt','ILOPT')
                
                print*, 
                print*, "Ionic liquid optimisation requested"
                print*, " - please use Optimiser.exe, not Mansaft.exe "
                print*,   
                print*,               
!---------------------------------------------------------------------------------------------------
            case('v','V')                
                print*, "Volume calculation"
                print*, "=================="
                print*,
                print*, " T(K)  P(MPa)  V(m^3 mol^-1)"  
                print*,
                
                do i = 1, properties%n                    
                    t = properties%t(i)
                    p = properties%p(i)                 
                    v = Vol_dens_g( )                  
                    write(*, '(f8.2, 2e16.6)') t, v, p / 1.e6_DP
                    call flush( )
                end do            
!---------------------------------------------------------------------------------------------------     
            case('vx','VX')                
                print*, "Volume calculation"
                print*, "=================="
                print*,
                print*, " T(K)  P(MPa)  V(m^3 mol^-1)"  
                print*,
                
                do i = 1, properties%n 
                    Comp_array(1:nctypes)%xi = properties%xi(i, 1:nctypes)                   
                    t = properties%t(i)
                    p = properties%p(i)                 

                    v = Vol_dens_g( )
                                    
                    write(*, '(f8.2, 10e16.6)') t, v, p / 1.e6_DP, Comp_array(:)%xi
                    call flush( )
                end do            
!--------------------------------------------------------------------------------------------------- 
            case('act','ACT','Act')
                
                print*, "Activity coefficient"
                print*, "===================="
                print*,
                print*, " T(K)  P(MPa)"  
                print*, t, p
                print*,
                print*, "molality    activity coeff"
                print*,
              
              
                do i = 1, properties%n                    
                    Comp_array(1:nctypes)%xi = properties%xi(i, 1:nctypes)
                    
                    if(ion_switch) then
                        Call Ion_act( )
                    end if

                    call flush( )
                end do             
!---------------------------------------------------------------------------------------------------
            case('pha1','PHA1','Pha1')
                
                print*, "Phase equilibrium"
                print*, "================="
                print*,
                print*, " T(K)  V liq , vap(m^3 mol^-1)  P(MPa)"  
                print*,

                do i = 1, properties%n                    
                    t = properties%t(i)
           
                    call Pure_phase(vl_out, vv_out)
              
                    write(*, '(f8.2, 3f16.10)') t, vl_out, vv_out, p / 1.e6_DP                  
                    call flush( )
                end do             
!---------------------------------------------------------------------------------------------------  
            case('vle','VLE')

                print*, "V-L Phase equilibrium"
                print*, "====================="
                print*,
                print*, " T(K)  P(MPa)  x1  y1  error"  
                print*,
                
                call Bin_simp()             
!---------------------------------------------------------------------------------------------------
            case('cp','Cp','CP')                
                print*, "Cp calculation"
                print*, "=============="
                print*,
                print*, " xi  T(K)  P(MPa) Cp(kJ mol^-1 K^-1"  
                print*,
                
                do i = 1, properties%n                    
                    t = properties%t(i)
                    p = properties%p(i)                 
                    Comp_array(1:nctypes)%xi = properties%xi(i, 1:nctypes)
                    g_out=Cp( )
                    write(*, '(2f20.2, 3e16.6)') Comp_array(1:nctypes)%xi, t, p / 1.e6_DP,v, g_out
                    call flush( )
                end do            
!---------------------------------------------------------------------------------------------------    
            case('Therm','therm','THERM')                
                print*, "Thermodynamic state variable calculations"
                print*, "========================================="
                print*,
                print*, " xi  T(K)  P(MPa) V A U H S G CP (all in kJ / mol / K as appropriate)"  
                print*,
                
                do i = 1, properties%n                    
                    t = properties%t(i)
                    p = properties%p(i)                 
                    Comp_array(1:nctypes)%xi = properties%xi(i, 1:nctypes)
                    call Therm_all( i )
                    call flush( )
                end do               
!---------------------------------------------------------------------------------------------------     
            case('Mu','mu','MU')
            print*, "v(m^3/mol),       mu_excess(i) (j/mol),   Mu(i),          G,       H"
                
                allocate(t_mu_out(1:nctypes), mu_total_out(1:nctypes))
                
                do i = 1, properties%n
				    Comp_array(1:nctypes)%xi = properties%xi(i, 1:nctypes)                   
                    t = properties%t(i)
                    p = properties%p(i)                 

                    v = Vol_dens_g( )				
                                        
                    h_out = EnthV( )
                    g_out = ThermGV( )
                    
                    do j = 1,nctypes                   
                        t_mu_out(j) = Mu_res(j)
						mu_total_out(j) = Mu(j)
                    end do
                    
					
                    print*, v, t_mu_out(:), mu_total_out(:), g_out, h_out
                   
                    call flush( )
                end do   
!---------------------------------------------------------------------------------------------------                                            
            case default
                print*, "Unrecognised property ",prop
                print*, " - Should have been caught in input read?!"
                stop
        end select
!***************************************************************************************************
!End time
!========
        call System_clock(finish)
       
        print*,
        print*, "#############################################################################################"
        print*, "                    ALL DONE  "
        write(*,'(a15, f14.2, a8)') "    Time taken: ",REAL(finish-start)/real(rate)," seconds"
        print*, "#############################################################################################"
        print*,
        
        call Quote( )
        print*,
        print*, "#############################################################################################"
        print*,
        print*,
        print*,
        print*,
        print*,
        
        call flush( )
!*************************************************************************************************** 

!***************************************************************************************************
    stop
end program Main
!***************************************************************************************************