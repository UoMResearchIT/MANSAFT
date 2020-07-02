program optimiser
!***************************************************************************************************
!Modules
!=======
    use Types       ! Definitions of types and double precision
    use Input       ! Read the input
    use Input_opt   ! Read optimisation parameters
    use Global      ! Important global parameters
    !use Simplex_mod     ! Simplex optimiser
    !use Quote       ! Random quote
    use nlopt 
	
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
    integer                         ::  n_opt !, opt_max
    real(kind=DP),allocatable       ::  min_opt(:), max_opt(:), paramin(:,:), paramout(:), values_opt(:)!, &
									!&	min_init(:), max_init(:)
    real(kind=DP)                   ::  errout!, values_opt(opt_num)
	! real(kind=DP),allocatable,public::  init_values(:)
!***************************************************************************************************
!Program header
!==============
    print*, " "
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
    print*, " "
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
		!print*, opt_num, min_num(1:opt_num), max_num(1:opt_num)
!***************************************************************************************************    
		!allocate( init_values(1:opt_num))!, min_init(1:opt_num), max_init(1:opt_num))
		
		!do i = 1, opt_num
			!init_values(i) = (min_num(i) + max_num(i))/2
			!print *, init_values(i)
		!end do
	 !print *, init_values(1:opt_num)
	
	allocate(values_opt(1:opt_num), comp_x(1:nctypes))
	allocate(Mu_L(1:nctypes), Mu_V(1:nctypes), x_init(1:nctypes+1))

	call opt_vle(values_opt)
	do i = 1, opt_num
			if(param_key(i)==1) then
				values_opt(i) = values_opt(i)/1000  !1000 is to make nag get converged
        !print*, Seg_array(param_index(i,1))%sig, x(i)
            else if(param_key(i)==2) then
				values_opt(i) = values_opt(i)/10     !10 is to make nag get converged 
            else if(param_key(i)==3) then
				values_opt(i) = values_opt(i)/100   !100 is to make nag get converged
            else if(param_key(i)==5) then
				values_opt(i) = values_opt(i)/10000 !10000 is to make nag get converged
			else if(param_key(i)==7)then
				values_opt(i) = values_opt(i)/10  !10 is to make nag get converged				
            else if(param_key(i)==8)then
				values_opt(i) = values_opt(i)/100 !100 is to make nag get converged
            else if(param_key(i)==10)then
				values_opt(i) = values_opt(i)/10 !10 is to make nag get converged
			end if
		end do
	write(*,100)  values_opt(1:opt_num)

100 Format (10(2X,F9.4))
!***************************************************************************************************
!End time
!========
    call System_clock(finish)
   
    print*,  " "
    print*,  "##################################################################"
    print*,  "                    ALL DONE  "
    print*,  "    Time taken: ",REAL(finish-start)/real(rate)," seconds"
    print*,  "##################################################################"
    print*,  " "
    print*,  " "
    
!    call PrintQuote( )
    print*, " "
    print*, "#############################################################################################"
    print*, " "
    print*, " "
    print*, " "
    print*, " "
    print*, " "

!***************************************************************************************************            
    stop
!***************************************************************************************************

!***************************************************************************************************
end program optimiser
!
