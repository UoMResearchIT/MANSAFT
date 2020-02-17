
!***************************************************************************************************

    Module c05qbfe_mod

!     C05QBF Example Program Module:
!            Parameters and User-defined Routines

!     .. Use Statements ..
      Use nag_library, Only: nag_wp
	  use Types           ! Definitions of types and double precision
	  use Global          ! Important global parameters
	  use Pressure
	  use ChemPot
    use Vol
        !     .. Implicit None Statement ..
      Implicit None
!     .. Accessibility Statements ..
      Private
      Public                           :: fcn
!     .. Parameters ..
      Integer, Parameter, Public       :: n_vle = 3, nout_vle = 6
	  Real (Kind=nag_wp), Public			   :: P_L, P_V, Mu_L_1, Mu_L_2, Mu_V_1, Mu_V_2, x_1, x_init_1, x_init_2, x_init_3, x_init_4
    Contains
      Subroutine fcn(n_vle,x,fvec,iuser,ruser,iflag)

!       .. Scalar Arguments ..
        Integer, Intent (Inout)        :: iflag
        Integer, Intent (In)           :: n_vle
!       .. Array Arguments ..
        Real (Kind=nag_wp), Intent (Out) :: fvec(n_vle)
        Real (Kind=nag_wp), Intent (Inout) :: ruser(*)
        Real (Kind=nag_wp), Intent (In) :: x(n_vle)
        Integer, Intent (Inout)        :: iuser(*)
!       .. Executable Statements ..
		v   = x(1)				! Liquid volume is x1
		Comp_array(1)%xi = x_1 !x1
		Comp_array(2)%xi = 1.0_nag_wp - x_1 !x2
		p = Press()
		P_L = p
		Mu_L_1 = Mu(1)
		Mu_L_2 = Mu(2)
		v   = x(2)				! Vapor volume is x2
		Comp_array(1)%xi = x(3) ! y1
		Comp_array(2)%xi =1.0_nag_wp - x(3) !y2
		p = Press()
		P_V = p
		Mu_V_1 = Mu(1)
		Mu_V_2 = Mu(2)
		
		
		
        fvec(1) = P_L - P_V
        fvec(2) = Mu_L_1 - Mu_V_1
		fvec(3) = Mu_L_2 - Mu_V_2
		
        !       Set iflag negative to terminate execution for any reason.
        iflag = 0
        Return
      End Subroutine fcn
    End Module c05qbfe_mod
!***************************************************************************************************    
subroutine c05qbfe(output_vle)

!     C05QBF Example Main Program

!     .. Use Statements ..
      Use c05qbfe_mod, Only: fcn, n_vle, nout_vle, x_init_1, x_init_2, x_init_3, x_init_4

      Use nag_library, Only: c05qbf, dnrm2, nag_wp, x02ajf
!     .. Implicit None Statement ..
      Implicit None
!     .. Local Scalars ..
      Real (Kind=nag_wp)               :: fnorm, xtol
      Integer                          :: i, ifail
!     .. Local Arrays ..
      Real (Kind=nag_wp), Allocatable  :: fvec(:), x(:)
      Real (Kind=nag_wp)               :: ruser(1)
      Integer                          :: iuser(1)
      !     .. Intrinsic Procedures ..
      Intrinsic                        :: sqrt
	  Real  (Kind=nag_wp)              :: output_vle(3)
!     .. Executable Statements ..

      Allocate (fvec(n_vle),x(n_vle))

!     The following starting values provide a rough solution.

      x(1) = x_init_1
      x(2) = x_init_2
      x(3) = x_init_3
      
      xtol = sqrt(x02ajf())

      ifail = -1
      Call c05qbf(fcn,n_vle,x,fvec,xtol,iuser,ruser,ifail)

      If (ifail==0 .Or. ifail==2 .Or. ifail==3 .Or. ifail==4) Then
        If (ifail==0) Then
!         The NAG name equivalent of dnrm2 is f06ejf
          fnorm = dnrm2(n_vle,fvec,1)
        Else
          Write (nout_vle,*)
          Write (nout_vle,*) 'Approximate solution'
		  x(3)=10
        End If
        !Write (nout,99998)(x(i),i=1,n)
		Do i = 1, n_vle
           output_vle(i) = x(i)
        End Do
      End If

99999 Format (1X,A,E12.4)
99998 Format (1X,3E12.4)
  return
End

!   E04JCF Example Program Text
!   Mark 26.1 Release. NAG Copyright 2017.
    Module e04jcfe_mod

!     E04JCF Example Program Module:
!            Parameters and User-defined Routines

!     .. Use Statements ..
      Use nag_library, Only: nag_wp
	  use Types           ! Definitions of types and double precision
	  use Global          ! Important global parameters
	  use Pressure
	  use Input
	  use Input_opt   ! Read optimisation parameters
	  use Vol
	  use c05qbfe_mod
	  use ChemPot
!     .. Implicit None Statement ..
      Implicit None
!     .. Accessibility Statements ..
      Private
      Public                           :: monfun, objfun
!     .. Parameters ..
      Integer, Parameter, Public       :: nout = 6 
      ! integer, public                         ::  n_opt, opt_max
      !real(kind=DP),allocatable, public       ::  min_bound(:), max_bound(:)
      
    Contains
      Subroutine objfun(n,x,f,iuser,ruser,inform)

!       .. Parameters ..
        ! Real (Kind=nag_wp), Parameter  :: five = 5.0_nag_wp
        ! Real (Kind=nag_wp), Parameter  :: ten = 1.0E1_nag_wp
        ! Real (Kind=nag_wp), Parameter  :: two = 2.0_nag_wp
!       .. Scalar Arguments ..
        Real (Kind=nag_wp), Intent (Out) :: f
        Integer, Intent (Out)          :: inform
        Integer, Intent (In)           :: n
!       .. Array Arguments ..
        Real (Kind=nag_wp), Intent (Inout) :: ruser(*)
        Real (Kind=nag_wp), Intent (In) :: x(n)
        Integer, Intent (Inout)        :: iuser(*)
		integer 						:: i
		real(Kind=nag_wp)				::  p_oj, y_oj(2),values(3), mu_1, mu_2, &
										& p_squ, y_squ(2)
!       .. Executable Statements ..
        inform = 0
		
		!allocate(min_bound(1:opt_num),max_bound(1:opt_num))
		do i = 1, opt_num
			if(param_key(i)==1) then               
				!min_bound(i)=min_num(i) * ANG
				!max_bound(i)=max_num(i) * ANG
				sig(param_index(i,1),param_index(i,1)) = x(i)* ANG/1000  !1000 is to make nag get converged
        !print*, Seg_array(param_index(i,1))%sig, x(i)
            else if(param_key(i)==2) then
				!min_bound(i)=min_num(i)
				!max_bound(i)=max_num(i)
				eps(param_index(i,1),param_index(i,1)) = x(i)/10     !10 is to make nag get converged 
            else if(param_key(i)==3) then
				!min_bound(i)=min_num(i)
				!max_bound(i)=max_num(i)
				lr(param_index(i,1),param_index(i,1)) = x(i)/100   !100 is to make nag get converged
            else if(param_key(i)==4) then
				!min_bound(i)=min_num(i)
				!max_bound(i)=max_num(i)
				la(param_index(i,1),param_index(i,1)) = x(i)
            else if(param_key(i)==5) then
				!min_bound(i)=min_num(i)
				!max_bound(i)=max_num(i)
				Seg_array(param_index(i,1))%sf = x(i)/10000 !10000 is to make nag get converged
            else if(param_key(i)==6) then
				!min_bound(i)=min_num(i)
				!max_bound(i)=max_num(i)
				Seg_array(param_index(i,1))%nseg = x(i)
			else if(param_key(i)==7)then
				!min_bound(i)=min_num(i)
				!max_bound(i)=max_num(i)
				eps(param_index(i,1), param_index(i,2))=x(i)/10  !10 is to make nag get converged
				eps(param_index(i,2), param_index(i,1))=x(i)/10 !10 is to make nag get converged				
            else if(param_key(i)==8)then
				!min_bound(i)=min_num(i)
				!max_bound(i)=max_num(i)
				lr(param_index(i,1), param_index(i,2))=x(i)/100 !100 is to make nag get converged
				lr(param_index(i,2), param_index(i,1))=x(i)/100 !100 is to make nag get converged
            else if(param_key(i)==9)then
				!min_bound(i)=min_num(i)
				!max_bound(i)=max_num(i)
				ehb(param_index(i,1),param_index(i,2), &
				& param_index2(i,1),param_index2(i,2)) = x(i)
				ehb(param_index2(i,1),param_index2(i,2), &
				& param_index(i,1),param_index(i,2)) = x(i)
            else if(param_key(i)==10)then
				!min_bound(i)=min_num(i) * ANG**3.0_DP * NA
				!max_bound(i)=max_num(i) * ANG**3.0_DP * NA
				khb(param_index(i,1),param_index(i,2), &
				& param_index2(i,1),param_index2(i,2)) = x(i)* ANG**3.0_DP * NA/10 !10 is to make nag get converged
				khb(param_index2(i,1),param_index2(i,2), &
				& param_index(i,1),param_index(i,2)) = x(i)* ANG**3.0_DP * NA/10 !10 is to make nag get converged
			end if
		end do
		
		p_oj = 0
		y_oj(1)=0
		y_oj(2)=0
		p_squ=0
		y_squ(1)=0
		y_squ(2)=0
		


	if (properties%opt_l(6)) then
		do i = 1, properties%nmu 
          Comp_array(1:nctypes)%xi = properties%ximu(i,1:nctypes)
		  t = properties%t_opt(i)
		  p = properties%p_opt(i)
          x_1 = properties%ximu(i,1)
		  P = 2*P
		  v = Vol_dens_g( )
		  mu_1 = Mu_res(1)
		  mu_2 = Mu_res(2)
		  x_init_1 = v
		  x_init_2 = 8.314*t/(properties%p_opt(i))
		  x_init_3 = exp(mu_1/(8.314*t))*x_1/(exp(mu_1/(8.314*t))*x_1+exp(mu_2/(8.314*t))*(1-x_1))
		  p = properties%p_opt(i)

          call c05qbfe(values )
			p_oj=p_oj+ABS(((P_L+P_V)/2-properties%p_opt(i))/properties%p_opt(i))
			p_squ=p_squ+(((P_L+P_V)/2-properties%p_opt(i))/properties%p_opt(i))**2
			y_oj(1)=y_oj(1)+ABS((values(3)-properties%yimu(i,1))/properties%yimu(i,1))
			y_squ(1)=y_squ(1)+((values(3)-properties%yimu(i,1))/properties%yimu(i,1))**2
			y_oj(2)=y_oj(2)+ABS((1.0_nag_wp-values(3)-properties%yimu(i,2))/properties%yimu(i,2))
			y_squ(2)=y_squ(2)+((1.0_nag_wp-values(3)-properties%yimu(i,2))/properties%yimu(i,2))**2
		end do
			p_oj=p_oj/properties%nmu
			p_squ=p_squ/properties%nmu
			y_oj(1)=y_oj(1)/properties%nmu
			y_squ(1)=y_squ(1)/properties%nmu
			y_oj(2)=y_oj(2)/properties%nmu
			y_squ(2)=y_squ(2)/properties%nmu
		end if		
        
		
		f = p_squ + y_squ(1) + y_squ(2)
   !print*, "f          x(i)"
   write(*, 101)"f=",f, "p_oj=", p_oj, "y_oj1=", y_oj(1), "y_oj2=", y_oj(2),"x(i)=",x(1:opt_num)
101 Format ((2X,A2),(2X,F6.4),(2X,A5),(2X,F6.4),(2X,A6),(2X,F6.4),(2X,A6),(2X,F6.4),(2X,A5),*(2X,F10.4))

        Return

      End Subroutine objfun
      Subroutine monfun(n,nf,x,f,rho,iuser,ruser,inform)

!       .. Scalar Arguments ..
		!use Types       ! Definitions of types and double precision
		!use Global,only: opt_num, min_num, max_num, init_values      ! Important global parameters
		
        Real (Kind=nag_wp), Intent (In) :: f, rho
        Integer, Intent (Out)          :: inform
        Integer, Intent (In)           :: n, nf
!       .. Array Arguments ..
        Real (Kind=nag_wp), Intent (Inout) :: ruser(*)
        Real (Kind=nag_wp), Intent (In) :: x(n)
        Integer, Intent (Inout)        :: iuser(*)
!       .. Local Scalars ..
        Logical                        :: verbose_output
!       .. Executable Statements ..
        inform = 0

        Write (nout,Fmt=99999) 'Monitoring: new trust region radius =', rho

!       Set this to .True. to get more detailed output
        verbose_output = .false.

        If (verbose_output) Then
          Write (nout,Fmt=99998) 'Number of function calls =', nf
          Write (nout,Fmt=99997) 'Current function value =', f
          Write (nout,Fmt=99996) 'The corresponding X is:', x(1:n)
        End If

        Return
99999   Format (/,4X,A,1P,E13.3)
99998   Format (4X,A,I16)
99997   Format (4X,A,1P,E12.4)
99996   Format (4X,A,/,(4X,5E12.4))
      End Subroutine monfun
    End Module e04jcfe_mod

module nag_e04jcfe_mod
contains    
	subroutine e04jcfe(output)

!     Example problem for E04JCF.

!     .. Use Statements ..
      Use e04jcfe_mod, Only: monfun, nout, objfun!, min_bound, max_bound
      Use nag_library, Only: e04jcf, nag_wp, x02alf
	  !use Types       ! Definitions of types and double precision
	  use Global,only: opt_num, min_num, max_num, init_values, param_key, ANG, NA       ! Important global parameters
!     .. Implicit None Statement ..
      Implicit None
!     .. Local Scalars ..
      Real (Kind=nag_wp)               :: f, infbnd, rhobeg, rhoend
      Integer                          :: ifail, maxcal, n, nf, npt
!     .. Local Arrays ..
      Real (Kind=nag_wp), Allocatable  :: bl(:), bu(:), x(:), output(:)
      Real (Kind=nag_wp)               :: ruser(1)!, output(opt_num)
      Integer                          :: iuser(1), i
!     .. Executable Statements ..
      Write (nout,*) 'E04JCF Example Program Results'

      maxcal = 5000
      rhobeg = 100
      rhoend = 1.0E-6_nag_wp
      n = opt_num !4
      npt = ((n+1)*(n+2))/2   !2*n + 1

!     x(3) is unconstrained, so we're going to set bl(3) to a large
!     negative number and bu(3) to a large positive number.

      infbnd = x02alf()**0.25_nag_wp
	!print *, infbnd
      Allocate (bl(n),bu(n),x(n))
		do i=1, opt_num
			bl(i) = min_num(i) !(/1.0_nag_wp,-2.0_nag_wp,-infbnd,1.0_nag_wp/)
			bu(i) = max_num(i) !(/3.0_nag_wp,0.0_nag_wp,infbnd,3.0_nag_wp/)
			x(i) = init_values(i) !(/3.0_nag_wp,-1.0_nag_wp,0.0_nag_wp,1.0_nag_wp/)
			!print *, bl(i),bu(i),x(i)
		end do
		!print*, bl(1:n),bu(1:n),x(1:n),output(1:n)
   print*, "initial values"
   write(*,102) x(1:n)
102 Format (*(2X,F9.4))
      ifail = -1
      Call e04jcf(objfun,n,npt,x,bl,bu,rhobeg,rhoend,monfun,maxcal,f,nf,iuser, &
        ruser,ifail)
	!print*, "test"
      Select Case (ifail)
      Case (0,2:5)
	  !print*, "test"
        If (ifail==0) Then
          Write (nout,Fmt=99999) 'Successful exit from E04JCF.',               &
            'Function value at lowest point found =', f
        Else
          Write (nout,Fmt=99998)                                               &
            'On exit from E04JCF, function value at lowest point found =', f
        End If
		!Write (nout,Fmt=99997) 'The corresponding X is:', x(1:n)
		!print*, "test"
		!Allocate (output(n))	
		do i = 1, opt_num
			
			!print*, "test", x(i)
				output(i) = x(i)
			!print*, "test2", x(i)	
		end do
	!print*, "test"
      End Select

99999 Format (2(/,1X,A),1P,E13.3)
99998 Format (/,1X,A,1P,E13.3)
99997 Format (1X,A,/,(2X,5E13.3))
	return
 End
 end 
 
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
    use Types       ! Definitions of types and double precision
    use Input       ! Read the input
    use Input_opt   ! Read optimisation parameters
    use Global      ! Important global parameters
    !use Simplex_mod     ! Simplex optimiser
    use Quote       ! Random quote
    use nag_e04jcfe_mod 
	
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
    real(kind=DP),allocatable       ::  min_opt(:), max_opt(:), paramin(:,:), paramout(:), values(:)!, &
									!&	min_init(:), max_init(:)
    real(kind=DP)                   ::  errout!, values(opt_num)
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
	
	allocate(values(1:opt_num))
	call e04jcfe(values)
	do i = 1, opt_num
			if(param_key(i)==1) then
				values(i) = values(i)/1000  !1000 is to make nag get converged
        !print*, Seg_array(param_index(i,1))%sig, x(i)
            else if(param_key(i)==2) then
				values(i) = values(i)/10     !10 is to make nag get converged 
            else if(param_key(i)==3) then
				values(i) = values(i)/100   !100 is to make nag get converged
            else if(param_key(i)==5) then
				values(i) = values(i)/10000 !10000 is to make nag get converged
			else if(param_key(i)==7)then
				values(i) = values(i)/10  !10 is to make nag get converged				
            else if(param_key(i)==8)then
				values(i) = values(i)/100 !100 is to make nag get converged
            else if(param_key(i)==10)then
				values(i) = values(i)/10 !10 is to make nag get converged
			end if
		end do
	write(*,100)  values(1:opt_num)

100 Format (2X,F9.4,$)
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
    
    call PrintQuote( )
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
