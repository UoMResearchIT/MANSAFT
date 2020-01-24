!***************************************************************************************************
!
!       MANSAFT OPTIMISER
!       =================
!
!       Version 2.7
!***************************************************************************************************
!   SAFT gamma mie program to:
!       1. Read general parameter list form input and
!       2. Calculate system parameters
!       3. Optimise parameters according to simplex method
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
!       Update 7  Sept 2016
!           - restructured with FORTRAN best practice
!       Update 7
!           - added electrolyte optimiser
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
!               S. Dufal, V. Papaioannou, M. Sadeqzadeh, T> Pogiatzis, A. Chremos, C.S. Adjiman, 
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
program optimiser
!***************************************************************************************************
!Modules
!=======
    use Types_mod       ! Definitions of types and double precision
    use Input_mod       ! Read the input
    use Input_opt_mod   ! Read optimisation parameters
    use Global_mod      ! Important global parameters
    use Simplex_mod     ! Simplex optimiser
    use Quote_mod       ! Random quote
!***************************************************************************************************
!Variables
!=========
    implicit none
    
    !Clock
    integer             ::  start, finish, rate
    !Property calculation
    character(len=1)    ::  prop
    !Counter
    integer             ::  i
    !Optimiser
    integer                         ::  n_opt, opt_max
    real(kind=DP),allocatable       ::  min_opt(:), max_opt(:), paramin(:,:), paramout(:)
    real(kind=DP)                   ::  errout
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
    print*, "++                                                                             ++"
    print*, "++                      ****     *****     ********                            ++"
    print*, "++                     ******    ******    ********                            ++"
    print*, "++                    **    **   **  ***      **                               ++"
    print*, "++                    **    **   **  ***      **                               ++"
    print*, "++                    **    **   ******       **                               ++"
    print*, "++                    **    **   ***          **                               ++"
    print*, "++                     ******    ***          **                               ++"
    print*, "++                      ****     ***          **                               ++"
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
    print*,"optimiser"    
!***************************************************************************************************
!Read input
!==========
    call Read_input( )    
      
    if((properties%type/='opt').and.(properties%type/='OPT').and.(properties%type/='Opt')   &
    & .and.(properties%type/='eopt').and.(properties%type/='EOPT')                          &
    & .and.(properties%type/='ilopt').and.(properties%type/='ILOPT')) then
        stop "For optimisation, input run type must be OPT, EOPT or ILOPT"
    end if

    call Read_opt( n_opt, max_opt, min_opt, paramin )    
!***************************************************************************************************    
!Do simplex
!==========  
    opt_max = 1000
    allocate( paramout(1:n_opt) )   
    call Simp( n_opt, paramin(:,:), errout, paramout(:), opt_max, max_opt(:), min_opt(:))
    
    call flush( )  
    print*,
    print*, paramout, errout        
!***************************************************************************************************
!End time
!========
    call System_clock(finish)
   
    print*, 
    print*,  "##################################################################"
    print*,  "                    ALL DONE  "
    print*,  "    Time taken: ",REAL(finish-start)/real(rate)," seconds"
    print*,  "##################################################################"
    print*,
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
    stop
!***************************************************************************************************

!***************************************************************************************************
end program optimiser
!***************************************************************************************************