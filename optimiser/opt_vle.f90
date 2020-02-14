!***************************************************************************************************
!   SAFT Module to:
!       1. Define double precision
!       2. Define types
!***************************************************************************************************
!
!***************************************************************************************************
module Types_mod
!***************************************************************************************************
!Modules
!=======    
    use,intrinsic   ::  iso_fortran_env
!***************************************************************************************************
    implicit none
!***************************************************************************************************
!   DOUBLE PRECISION DEFINITION
!   ===========================
    integer,parameter       ::  DP=REAL64
!***************************************************************************************************
!   SEGMENTS
!   ~~~~~~~~   
!
!   defintion array
!   ~~~~~~~~~~~~~~~ 
!   lab     -   segment name
!   sig     -   mie sigma in / m
!   eps     -   mie eps in K
!   sf      -   shape factor
!   m       -   molar mass in kg
!   nseg    -   number of identical spherical segments compsosing segment
!   xsi     -   fraction of segments in mixture
!   dxsin   -   differential of xsi wrt N
!   nassoc  -   number of types of association sites
!   nassoctype  -   number of each type of association site present in segment
!   charged -   is this group charged?
!   q       -   charge
!   sigb    -   born sigma
!
    type SEGMENT
        character(len=10)               ::  lab
        real(kind=DP)                   ::  sig, eps, la, lr, sigb
        real(kind=DP)                   ::  sf, m, nseg
        real(kind=DP)                   ::  xsi,dxsin
        integer                         ::  nassoc
        integer,allocatable             ::  nassoctyp(:)
        logical                         ::  charged
        real(kind=DP)                   ::  q
    end type SEGMENT
!***************************************************************************************************
!   COMPONENTS / MOLECULES
!   ======================
!   
!   definition array
!   ~~~~~~~~~~~~~~~~
!
!   lab     -   name
!   xi      -   mole fraction
!   dxin    -   differential of mole fraction wrt N
!   m       -   mass in kg
!   db3     -   thermal debroglie wavelength / m
!   rho     -   number density moles / m3
!   ms      -   number of segments in molecule
!   comp    -   composition array - number of each segment type
!   nm      -   number of moles - effectively mole fraction at present
!   solv    -   is this a solvent particle?
!   dv / dt -   solvent parameters for dielectric
!
   type COMPONENT
       character(len=10)                ::  lab
       real(kind=DP)                    ::  xi,m,nm,dxin
       real(kind=DP)                    ::  db3,rho      
       integer                          ::  ms
       integer,allocatable              ::  comp(:)
       logical                          ::  solv
       real(kind=DP)                    ::  dv, dt
    end type COMPONENT
!***************************************************************************************************
!   SYSTEM
!   ======
!   
!   definition array
!   ~~~~~~~~~~~~~~~~
!
!   type    -   property to calculate
!   n       -   number of properties
!   v       -   volume array
!   t       -   temperaure array
!   p       -   pressure array
!   phase   -   phase array
!   xi      -   mole fraction arrays
!   A,B,C,D,E   -   ideal heat capacity (Cp) coefficients
!   nmu / nliq  -   no. opt points
!   v_l / v_v   - liquid and vapour volumes for opt
!   p_opt   -   p for optimisations
!   p_liq / t_liq / v_liq /cp   -   liquid parameters
!   anion / cation  - integer listing component for activity calculations
!   an_stoich, cat_stoich   -   stoichiometries
!   p_l     -   phase equib pressures in data limited optimisations
!
    type SYSTEM
        character(len=5)                    ::  type
        integer                             ::  n, nmu, nliqv, nliqh, nliqc, gi
        integer                             ::  anion, cation, cat_stoich, an_stoich
        real(kind=DP),allocatable           ::  v(:), t(:), v_liq(:), p(:),     &
                                            &   xi(:,:), v_l(:), v_v(:), p_opt(:), t_opt(:), cp(:), h(:), &
                                            &   xiv(:,:), xih(:,:), xic(:,:), xip(:,:),ximu(:,:), yimu(:,:), &				!yimu   yichun
                                            &    t_liqv(:), p_liqv(:), t_liqh(:), p_liqh(:), t_liqc(:), p_liqc(:)      
        character(len=1),allocatable        ::  phase(:)
        character(len=1)                    ::  opt(1:6)		!from 5 to 6  yichun
        logical                             ::  opt_l(1:6)		!from 5 to 6  yichun
        real(kind=DP)                       ::  p1, p2 !binary vle calculations
        real(kind=DP)                       ::  A, B, C, D, E
        
        !for electrolyte opt
        integer                             ::  nliqv_e, nliqv_vle, nliqp_vle, nliqm, ndhmix
        real(kind=DP),allocatable           ::  xi_liqve(:,:), t_liqve(:), p_liqve(:), v_liqve(:)
        real(kind=DP),allocatable           ::  xiv_liqvle(:,:), t_liqvle(:), v_liqvle(:)
        real(kind=DP),allocatable           ::  xidhmixi(:,:), xidhmixf(:,:), tdhmix(:), pdhmix(:), dhmix(:)
        real(kind=DP),allocatable           ::  xip_liqvle(:,:), tp_liqvle(:), p_liqvle(:)
        real(kind=DP),allocatable           ::  xim(:,:), tm(:), pm(:), mu(:)
    end type SYSTEM
!***************************************************************************************************
end module Types_mod
!***************************************************************************************************

!***************************************************************************************************
!   SAFT Module to:
!       1. Hold key global variables
!       2. Define constants
!***************************************************************************************************
!
!***************************************************************************************************
module Global_mod
!***************************************************************************************************
!Modules
!=======
    use Types_mod       ! Definitions of types and double precision
!***************************************************************************************************
    implicit none
!***************************************************************************************************
!   CONSTANTS
!   =========
!
!   Fundamental constants from http://physics.nist.gov/cuu/Constants/
!
!   ANG     -   convert ANG to m
!   AMU     -   convert AMU to kg
!   NA      -   Avogadro's number
!   KB      -   Boltzmann constant
!   H       -   Planck's constant
!   PI / TWO PI - as name suggests...
!   QE      -   Electron charge
!   CBORN   -   Constant for Born calculation
!
    real(kind=DP),parameter  ::  ANG=1.e-10_DP, AMU=1.660539040e-27_DP
    real(kind=DP),parameter  ::  NA=6.022140857e23_DP, KB=1.38064852e-23_DP
    real(kind=DP),parameter  ::  H=6.626070040e-34_DP, PI=3.14159265359e0_DP
    real(kind=DP),parameter  ::  TWOPI=2.e0_DP*PI
    real(kind=DP),parameter  ::  QE=1.6021766208e-19_DP, E0=8.854187817e-12_DP
!***************************************************************************************************
!   Segments
!   ========
!
!   nstypes         -   number of segment types
!   Seg_array       -   array of segment details
!   eps,sig,la,lr   -   2 D interaction arrays for Mie epsilon, sigma, lambda attractive and repulsive
!   cij     -   Mie coefficients    [1] A3
!   dij     -   effective diameters [1] A9 / A46
!   un_eps  -   unlike eps values
!   un_lr   -   unlike lr values
!
    integer                         ::  nstypes
    type(SEGMENT),allocatable       ::  Seg_array(:)
    real(kind=DP),allocatable       ::  eps(:,:), sig(:,:), la(:,:), lr(:,:)
    real(kind=DP),allocatable       ::  cij(:,:), dij(:,:)
    real(kind=DP),allocatable       ::  un_eps(:), un_lr(:)
    integer, allocatable            ::  un_ieps(:), un_jeps(:), un_ilr(:), un_jlr(:)
    integer                         ::  un_neps, un_nlr
!***************************************************************************************************
!   Components (molecules)
!   ==========
!
!   nctypes         -   number of componenets (molecule types)
!   Comp_array      -   array of component details
    integer                         ::  nctypes
    type(COMPONENT),allocatable     ::  Comp_array(:)
!***************************************************************************************************
!   Property
!   ========
!
!   properties      -   properties to calculae
!   t, v, p         -   current temperature, volume, pressure
!   rho             -   density in mole per m^3
!   rhos            -   sphecial segment density
!   beta            -   1 / kB / T
!
    type(SYSTEM)        ::  properties
    real(kind=DP)       ::  t, v, p, rho, rhos
    real(kind=DP)       ::  beta, beta_K
!***************************************************************************************************
!   GAUSSIAN LEGENDRE PARAMETERS
!   ============================
!
!   N_GL    -   number of point - set to 64
!   X_GL    -   abssica - read from gl_x.txt
!   W_GL    -   weights - read from gl_w.txt
!
    integer,parameter       ::  N_GL=330
    real(kind=DP)           ::  X_GL(1:N_GL), W_GL(1:N_GL)
!***************************************************************************************************
!   CONSTANTS ARRAYS
!   ================
!
!   c_array     -   coefficients in [1] A17 (effective packing fraction)
!   phik        -   coefficients in [1] A26 (part of monomer contribution)
!   c_x         -   coefficients in [2] 24  (part of association, parameters originally from [3])
!
    real(kind=DP)           ::  c_array(1:4, 1:4), phik(0:7, 0:7), c_x(0:10, 0:10)  
!***************************************************************************************************
!   CALCULATION SWITCHES
!   ====================
!
    logical     ::  chain_switch, assoc_switch, ion_switch
!***************************************************************************************************
!   USED IN A MONO CALCULATION
!   ==========================
!
!   zl       - moments of number density (0,1,2,3)                  [1] A7
!   zetax    - packing fraction of pure fluid with diameter dij     [1] A13
!   KHS      - hard sphere isothermal compressibility               [1] A21
!   zetabar  - packing fraction of pure fluid with diameter sigma   [1] A23
!   fkij     - functions                                            [1] A26
!   alpha    -                                                      [1] A24
!
    real(kind=DP)                ::  zl(0:3), zetax, khs, zetabar
    real(kind=DP),allocatable    ::  fkij(:,:,:), alpha(:,:)
!***************************************************************************************************
!   USED IN A CHAIN CALCULATION
!   ===========================
!
!   sig3ij      -   average molecular segment diameter  [2] 41
!   d3ij        -   reference hard sphere diameter      [2] 42
!   epschij     -   molecular eps                       [2] 44
!   lachij      -   molecular la                        [2] 45
!   lrchij      -   molecular lr                        [2] 45
!   zetax_sum   -   sum of all zetax, used in gradient calculations
!
    real(kind=DP),allocatable    ::  sig3ij(:,:), d3ij(:), epschij(:,:), lachij(:), lrchij(:)
    real(kind=DP)                ::  zetax_sum
!***************************************************************************************************
!   USED IN A ASSOC CALCULATION
!   ===========================
!
!   x       - mass action array                 [2] 22
!   delx    - association strengths             [2] 23
!   sigx3   -                                   [2] 25 
!   ehb     - associaiton parameters            [1] A41
!   khb     - input parameter
!
    real(kind=DP),allocatable    ::  x(:,:,:), delx(:,:,:,:,:,:), ehb(:,:,:,:), khb(:,:,:,:)
    real(kind=DP)                ::  sigx3  
!***************************************************************************************************
!   DIFFERENTIAL BITS
!   =================
!
!   all labels are as before plus 
!       d____n for Ni differentials
!       d____v for V differentials
!

!Chemical potential calculation
    real(kind=DP)               ::  nsum    !sum of n particles
    real(kind=DP)               ::  dsigx3n, di_xn
    real(kind=DP),allocatable   ::  ddelxn(:,:,:,:,:,:), dxn(:,:,:)
    real(kind=DP)               ::  drhon, drhosn, dzln(0:3)
    real(kind=DP)               ::  dzetaxn, dzetabarn, dKHSn, dzetax_sumn

!Pressure calculation
    real(kind=DP)               ::  drhov, drhosv, dzlv(0:3), dzetaxv, dzetabarv, dKHSv
    real(kind=DP),allocatable   ::  ddelxv(:,:,:,:,:,:), dxv(:,:,:)
!***************************************************************************************************
!   OPTIMISER
!   =========
!   param_key - which parameters being optimised
!   param_index - which beads relate to param_key
!       1 sig, 2 eps, 3 lr, 4 la, 5 sf, 6 nseg, 7 eps ij, 8, lr ij, 9 ehb, 10 khb
    
    integer, allocatable        ::  param_key(:), param_index(:,:), param_index2(:,:)
	real(kind=DP),allocatable	:: init_values(:), min_num(:), max_num(:)
	integer 					:: opt_num
    
!***************************************************************************************************
!   PHASE EQUILIBRIUM
!   =================
!
    !inputs 
    real(kind=DP), allocatable  ::  fl_p_crit_in(:), fl_t_crit_in(:),fl_w_crit_in(:)  
!***************************************************************************************************
!   OTHER
!   =====
!   info_on - logical switch to turn on vebose data printing for debugging
!
    logical :: info_on
!***************************************************************************************************
END MODULE Global_mod
!***************************************************************************************************
!***************************************************************************************************
! Marsaglia & Tsang generator for random normals & random dexponentials.
! Translated from C by Alan Miller (amiller@bigpond.net.au)

! Marsaglia, G. & Tsang, W.W. (2000) `The ziggurat method for generating
! random variables', J. Statist. Software, v5(8).

! This is an electronic journal which can be downloaded from:
! http://www.jstatsoft.org/v05/i08

! N.B. It is assumed that all integers are 32-bit.
! N.B. The value of M2 has been halved to compensate for the lack of
!      unsigned integers in Fortran.

! Latest version - 1 January 2001

! Modified by SJH to bring inline with FORTRAN standards
!***************************************************************************************************
module Zig_mod
    use Types_mod
    implicit none
    
    public ::   zigset, uni
    real(kind=DP), parameter    ::  m1=2147483648.0e0_DP,   m2=2147483648.0e0_DP,half=0.5e0_DP                         
    real(kind=DP)               ::  dn=3.442619855899e0_DP, tn=3.442619855899e0_DP, vn=0.00991256303526217e0_DP,           &
                                    q, de=7.697117470131487e0_DP, te=7.697117470131487e0_DP,            &
                                    ve=0.003949659822581572e0_DP
    integer, save               ::  iz, jz, jsr=123456789, kn(0:127), ke(0:255), hz                       
    real(kind=DP), save         ::  wn(0:127), fn(0:127), we(0:255), fe(0:255)
    logical, save               ::  initialized=.false.
!***************************************************************************************************
contains
!***************************************************************************************************
subroutine zigset( jsrseed )

   integer, intent(in)  :: jsrseed
   integer  :: i

   !  Set the seed
   jsr = jsrseed

   !  Tables for RNOR
   q = vn*dexp(half*dn*dn)
   kn(0) = (dn/q)*m1
   kn(1) = 0
   wn(0) = q/m1
   wn(127) = dn/m1
   fn(0) = 1.0_DP
   fn(127) = dexp( -half*dn*dn )
   DO  i = 126, 1, -1
      dn = dsqrt( -2.0_DP * dlog( vn/dn + dexp( -half*dn*dn ) ) )
      kn(i+1) = (dn/tn)*m1
      tn = dn
      fn(i) = dexp(-half*dn*dn)
      wn(i) = dn/m1
   END DO

   !  Tables for Rdexp
   q = ve*dexp( de )
   ke(0) = (de/q)*m2
   ke(1) = 0
   we(0) = q/m2
   we(255) = de/m2
   fe(0) = 1.0_DP
   fe(255) = dexp( -de )
   
   do  i = 254, 1, -1
      de = -dlog( ve/de + dexp( -de ) )
      ke(i+1) = m2 * (de/te)
      te = de
      fe(i) = dexp( -de )
      we(i) = de/m2
   end do
   initialized = .true.
   return
end subroutine zigset
!***************************************************************************************************
!  Generate uniformly distributed random numbers
function uni( ) result( fn_val )
   real(kind=DP)  ::  fn_val

   fn_val = half + 0.2328306e-9_DP * shr3( )
   return
end function uni
!***************************************************************************************************
function shr3( ) result( ival )
   integer  ::  ival

   jz = jsr
   jsr = IEOR( jsr, ISHFT( jsr,  13 ) )
   jsr = IEOR( jsr, ISHFT( jsr, -17 ) )
   jsr = IEOR( jsr, ISHFT( jsr,   5 ) )
   ival = jz + jsr
   
    return
end function shr3
!***************************************************************************************************
end module Zig_mod
!***************************************************************************************************
!***************************************************************************************************
!   SAFT Module to:
!       1. SETUP COMMON DATA
!***************************************************************************************************
!   Created by SJ Halstead Nov 2015
!       Update Jan 2016
!           -   Making more efficient / functions
!       Update Mar 2016
!           -   Removed units nm -> m for clarity
!           -   Made all arrays to deallocate if previously allocated
!       Update June 2016
!           -   Tidied and efficiency approved
!***************************************************************************************************
!
!***************************************************************************************************
!***************************************************************************************************
!   SAFT Module to:
!       1. Calculate GaussianLegednre quadrature abscissas and weights
!***************************************************************************************************
!
!   This module is slightly modified from:
!       http://web.aeromech.usyd.edu.au//wwwcomp/subrout.html
!
!       Gaussm3.f90
!
!***************************************************************************************************
module GL_mod
!***************************************************************************************************
!Modules
!=======    
    use Types_mod           ! Definitions of types and double precision
    use Global_mod          ! Important global parameters 
!***************************************************************************************************
    implicit none
!***************************************************************************************************    
    contains
!***************************************************************************************************
	subroutine GL_arrays( )

        implicit none
        integer                     :: i, j, m
        real(kind=DP)               :: p1, p2, p3, pp, z, z1
        real(kind=DP), parameter    :: EPS=3.0e-15_DP       

	    m = (N_GL + 1) / 2
	    
	    do i=1, m				
     		z = dcos( PI * (i-0.25e0_DP) / (N_GL + 0.5e0_DP) )

100     	p1 = 1.0e0_DP
        	p2 = 0.0e0_DP

        	do j = 1, N_GL
                p3 = p2
                p2 = p1
                p1 = ((2.0e0_DP * j - 1.0e0_DP) * z * p2 - (j - 1.0e0_DP) * p3) / j
        	enddo

        	pp = N_GL * (z * p1 - p2) / (z * z - 1.0e0_DP)
        	z1 = z
        	z  = z1 - p1 / pp            

        	if (dabs(z - z1) > EPS) goto  100

      	    X_GL(i) =  - z                    	
            X_GL(N_GL+1-i) =  + z               
            W_GL(i) = 2.0d0/((1.0d0-z*z)*pp*pp) 
      	    W_GL(N_GL+1-i) = W_GL(i)               
        end do     
    end subroutine GL_arrays
!***************************************************************************************************
end module GL_mod
!***************************************************************************************************
!   SAFT Module to:
!       1. READ INPUT AND ASSIGN ARRAYS
!***************************************************************************************************
!   Created by SJ Halstead Oct 2015
!       1. Tidied 23/6/2016
!***************************************************************************************************
!
!***************************************************************************************************
module Input_mod
!***************************************************************************************************
!MODULES
use Global_mod     !CONTAINS GLOBAL VARIABLES
use Types_mod
use GL_mod         !Calculate G-L dabscissas ans weights
!***************************************************************************************************
    implicit none
!***************************************************************************************************
contains
!***************************************************************************************************
    subroutine Read_input( testinput )
!--------------------------------------------------------------------------
        implicit none
!--------------------------------------------------------------------------        
!Variables        
        !local do counters
        integer         ::  i1,i2,i3,i4,ierr
        integer         ::  j1,j2,j3,j4
        !input file
        character(len=80)               ::  filein
        character(len=80)               ::  line
        character(len=100)               ::  longline
        character(len=*),optional       ::  testinput
        !temp values
        real(kind=DP)   ::  doub,doub2
        !gl input
        integer                         ::  cfin        
        character(len=80)               ::  path,path1        
!--------------------------------------------------------------------------        
!Read file        
        if(present(testinput)) then
            filein = testinput
        else  
            call getarg(1,filein)        
        end if
   
        open(11, file=filein, status="old", action="read", iostat=ierr)
     
        if(ierr/=0) then
            print*, "Failed to open ",filein
            stop
        end if
    
        read(11,'(a80)') line
        read(11,*) 
      
        
        print*, "Reading file:"
        write(*,'(a5, a80)') "     ", line       
!--------------------------------------------------------------------------             
!Segments      
        chain_switch = .false.       
        assoc_switch = .false.
        ion_switch   = .false.
  
        read(11,*) nstypes       

        allocate(Seg_array(1:nstypes), stat=ierr)
      
        if(ierr/=0) then
            print*, "Failed to allocate segment array, nstypes = ",nstypes
            stop
        end if
      
        do i1=1, nstypes
            read(11,*, iostat=ierr)                                                 &
            &   Seg_array(i1)%lab, Seg_array(i1)%sig, Seg_array(i1)%eps,            &
            &   Seg_array(i1)%lr, Seg_array(i1)%la, Seg_array(i1)%sf,               &
            &   Seg_array(i1)%nseg, Seg_array(i1)%m, Seg_array(i1)%q,               & 
            &   Seg_array(i1)%sigb, Seg_array(i1)%nassoc                                

            !Set ion_switch
            if(.not.ion_switch) then
                if(dabs(Seg_array(i1)%q)>1.0e-6_DP) ion_switch=.true.
            end if
            
            if(ierr/=0) then
                print*, "Problem reading component ",i1
                stop
            end if
            
            Seg_array(i1)%charged=.true. 
            if(dabs(Seg_array(i1)%q) <=1.0e-6_DP) Seg_array(i1)%charged=.false.  

            !Number of associations types
            i3 = Seg_array(i1)%nassoc
            allocate(Seg_array(i1)%nassoctyp(1:i3), stat=ierr)
            
            if(ierr/=0) then
                print*, "Failed to allocate associations array for component ",i1
                stop
            end if
            
            !Number of sites for each association type on a segment
            do i2=1,i3
                read(11,*) Seg_array(i1)%nassoctyp(i2)               
            end do
            
            !Adjust units to SI
            Seg_array(i1)%sig  = Seg_array(i1)%sig  * ANG
            Seg_array(i1)%sigb = Seg_array(i1)%sigb * ANG
            Seg_array(i1)%m    = Seg_array(i1)%m    * AMU
        end do
        
        read(11,*)      
!--------------------------------------------------------------------------      
!Form interaction arrays        
        allocate(eps(1:nstypes, 1:nstypes), stat=ierr)
        if(ierr/=0) then
            print*, "Failed to allocate eps array "
            stop
        end if
            
        allocate(sig(1:nstypes, 1:nstypes), stat=ierr)
        if(ierr/=0) then
            print*, "Failed to allocate sig array "
            stop
        end if
        
        allocate(la(1:nstypes, 1:nstypes), stat=ierr)
        if(ierr/=0) then
            print*, "Failed to allocate lambda association array "
            stop
        end if
        
        allocate(lr(1:nstypes, 1:nstypes), stat=ierr)
        if(ierr/=0) then
            print*, "Failed to allocate lambda repulsion array "
            stop
        end if
       
        do i2=1,nstypes                                                                 !Full 2D array
        do i1=i2,nstypes
            sig(i1, i2)   =  (Seg_array(i1)%sig + Seg_array(i2)%sig) / 2.0_DP           ![1] A45
            eps(i1, i2)  =  dsqrt(Seg_array(i1)%sig**3.0_DP * Seg_array(i2)%sig**3.0) / sig(i1,i2)**3.0_DP &
            &   * dsqrt(Seg_array(i1)%eps * Seg_array(i2)%eps)                           ![1] A51
            
            la(i1, i2) = dsqrt((Seg_array(i1)%la - 3.0_DP) * (Seg_array(i2)%la - 3.0_DP)) + 3.0_DP     ![1] A56
            lr(i1, i2) = dsqrt((Seg_array(i1)%lr - 3.0_DP) * (Seg_array(i2)%lr - 3.0_DP)) + 3.0_DP     ![1] A56
            
            !Reflect 
            sig(i2, i1) = sig(i1, i2)
            eps(i2, i1) = eps(i1, i2)
            la(i2, i1)  = la(i1, i2)
            lr(i2, i1)  = lr(i1, i2)
        end do
        end do
!--------------------------------------------------------------------------     
!Read unlikes
        !unlike with specific eps       
        read(11,*) i2
        un_neps = i2
        allocate(un_eps(1:un_neps), un_ieps(1:un_neps), un_jeps(1:un_neps) )
      
        do i1=1,i2
            read(11,*) j1, j2, doub
          
            eps(j1,j2) = doub
            eps(j2,j1) = doub
            
            un_eps(i1)  = doub
            un_ieps(i1) = j1
            un_jeps(i1) = j2
        end do
       
        read(11,*)
        read(11,*) i2
        un_nlr = i2
        allocate(un_lr(1:un_nlr), un_ilr(1:un_nlr), un_jlr(1:un_nlr) )
       
        do i1=1,i2
            read(11,*) j1, j2, doub
         
            lr(j1,j2) = doub
            lr(j2,j1) = doub
            
            un_lr(i1)  = doub
            un_ilr(i1) = j1
            un_jlr(i1) = j2
        end do
              
        read(11,*)   
!--------------------------------------------------------------------------      
!Read associations
        allocate(ehb(1:nstypes, 1:3, 1:nstypes, 1:3), stat=ierr)
        if(ierr/=0) then
            print*, "Failed to allocate ehb association array "
            stop
        end if
        
        allocate(khb(1:nstypes, 1:3, 1:nstypes, 1:3), stat=ierr)
        if(ierr/=0) then
            print*, "Failed to allocate khb association array "
            stop
        end if
        
        ehb=0.0_DP   !energy term
        khb=0.0_DP   !volume term
        
        read(11,*) i3
        
        !Set assoc switch
        if(i3/=0) assoc_switch=.true.

        do i1=1,i3
            read(11,*) j1,j2,j3,j4,doub,doub2        
        
            ehb(j1,j2,j3,j4) = doub
            khb(j1,j2,j3,j4) = doub2 * ANG**3.0_DP * NA  !Scale to m^3 mol^-1
            
            ehb(j3,j4,j1,j2) = doub
            khb(j3,j4,j1,j2) = doub2 * ANG**3.0_DP * NA  !Scale to m^3 mol^-1
        end do                                             
!--------------------------------------------------------------------------     
!Read molecule descriptions
        read(11,*) nctypes   !number of molecules / components
        allocate(Comp_array(1:nctypes),stat=ierr)
       
        if(ierr/=0) then
            print*, "Failed to allocate comp array "
            stop
        end if    

        do i1=1,nctypes
            !Molecule description
            read(11,*, iostat=ierr)  Comp_array(i1)%lab, Comp_array(i1)%xi, Comp_array(i1)%m,     &
            &                        Comp_array(i1)%dv, Comp_array(i1)%dt,                        &   
            &                        Comp_array(i1)%ms   !name, mole fraction, mass, number of segments

            if(ierr/=0) then
                print*, "Problem reading component array, component ",i1
                stop
            end if
        
            Comp_array(i1)%solv = .true.
            if((dabs(Comp_array(i1)%dv)<1.0e-6_DP).and.(dabs(Comp_array(i1)%dv)<1.0e-6_DP))       &
            &   Comp_array(i1)%solv=.false.
        
            Comp_array(i1)%m = Comp_array(i1)%m * AMU  !to kg per particle
            
            allocate(Comp_array(i1)%comp(1:nstypes), stat=ierr)
            if(ierr/=0) then
                print*, "Failed to allocate comp molecule list array "
                stop
            end if
        
            Comp_array(i1)%comp = 0
            
            !List number of each type of segment in molecule
            do i2 = 1, Comp_array(i1)%ms
                read(11,*) j1, j2
                Comp_array(i1)%comp(j1) = j2                              
            end do
            !set no. segments to total no.
            Comp_array(i1)%ms = sum(Comp_array(i1)%comp(1:Comp_array(i1)%ms))
            !Set chain switch         
            if(.not.chain_switch) then
                if(Comp_array(i1)%ms>1) chain_switch=.true.
                
                do i2=1, nstypes
                    if(Comp_array(i1)%comp(i2)/=0) then
                        if(Seg_array(i2)%nseg>1.0e0_DP) then
                            chain_switch=.true.
                        end if
                    end if
                end do
                
            end if

        end do
    
        read(11,*)     
!--------------------------------------------------------------------------     
!Calculation type
        read(11,*) properties%type, properties%n

        if(allocated(properties%v)) deallocate(properties%v)
        allocate(properties%v(1:properties%n), stat=ierr)
        if(ierr/=0) then
            print*, "Failed to allocate property v array "
            stop
        end if 
        
        if(allocated(properties%t)) deallocate(properties%t)
        allocate(properties%t(1:properties%n), stat=ierr)
        if(ierr/=0) then
            print*, "Failed to allocate property t array "
            stop
        end if 
        
        if(allocated(properties%p)) deallocate(properties%p)
        allocate(properties%p(1:properties%n), stat=ierr)
        if(ierr/=0) then
            print*, "Failed to allocate property p array "
            stop
        end if 
        
        if(allocated(properties%phase)) deallocate(properties%phase)
        allocate(properties%phase(1:properties%n), stat=ierr)
        if(ierr/=0) then
            print*, "Failed to allocate property phase array "
            stop
        end if 
        
        if(allocated(properties%xi)) deallocate(properties%xi)
        allocate(properties%xi(1:properties%n, 1:nctypes), stat=ierr)
        if(ierr/=0) then
            print*, "Failed to allocate property xi array "
            stop
        end if 
        
        select case(properties%type)
            !Pressure
            case('p','P')  
                do i1=1,properties%n
                    read(11,*) properties%t(i1), properties%v(i1)
                end do
            !Pressure for varying conc
            case('px','PX')  
                do i1=1,properties%n
                    read(11,*) properties%t(i1), properties%v(i1), properties%xi(i1, 1:nctypes)
                end do
            !Volume    
            case('v','V')   
                do i1=1,properties%n
                    read(11,*) properties%t(i1), properties%p(i1)
                end do
            !Volume for varying mole fractions   
            case('vx','VX')   
                do i1=1,properties%n                    
                    read(11,*) properties%t(i1), properties%p(i1), properties%xi(i1, 1:nctypes)
                end do
            !Pure phase equilibrium
            case('pha1','PHA1','Pha1')
                if(nctypes>1) stop "Method  PHA1  is only for single component (pure) phase equilibrium"
                
                do i1=1,properties%n
                    read(11,*) properties%t(i1)
                end do
            !Binary VLE
            case('vle','VLE')
                if(nctypes/=2) print*, "Method  VLE  is only for binary liquid-vapour equilibrium"
               
                do i1=1,properties%n
                    read(11,*) properties%t(i1), properties%p(i1),properties%xi(i1, 1)                  
                end do
               
            !Activity coefficient
            case('act','ACT','Act')
                read(11,*) t, p
                read(11,*) properties%cation, properties%cat_stoich
                read(11,*) properties%anion, properties%an_stoich
                     
                do i1=1,properties%n
                    read(11,*) properties%xi(i1, 1:nctypes)
                end do   
            !Cp heat capacity at constant pressure
            case('cp','Cp','CP')  
                read(11,*) properties%A, properties%B, properties%C, properties%D, properties%E
            
                do i1=1,properties%n
                    read(11,*) properties%xi(i1, 1:nctypes),properties%t(i1), properties%p(i1)
                end do 
            !Themodynamic state data
            case('Therm','therm','THERM')  
                read(11,*) properties%A, properties%B, properties%C, properties%D, properties%E
                
                do i1=1,properties%n
                    read(11,*) properties%xi(i1, 1:nctypes),properties%t(i1), properties%p(i1)
                end do 
            !Gibbs free energy
            case('Gval','gval','GVAL')  
                do i1=1,properties%n
                    read(11,*) properties%xi(i1, 1:nctypes),properties%t(i1), properties%p(i1)
                end do
            !Excess chemical potential
            case('Mu','mu','MU')  
                do i1=1,properties%n
                    read(11,*) properties%xi(i1, 1:nctypes), properties%v(i1), properties%t(i1)
                end do                 
!********************************************************************************************************************************
            !special - temporary testng input choice
            case('spec')
!                do i1=1,properties%n
!                    read(11,*) properties%t(i1),properties%p(i1)!,properties%xi(i1,1:nctypes)
!                end do
!********************************************************************************************************************************                                               
            !Optimiser
            case('OPT','Opt','opt')   
                read(11,*) properties%opt(:) 
                properties%opt_l(:) = .false.
                
                if((properties%opt(1)=='y').or.(properties%opt(1)=='Y')) properties%opt_l(1)=.true. !Liq V
                if((properties%opt(2)=='y').or.(properties%opt(2)=='Y')) properties%opt_l(2)=.true. !Liq Cp
                if((properties%opt(3)=='y').or.(properties%opt(3)=='Y')) properties%opt_l(3)=.true. !Liq H
                if((properties%opt(4)=='y').or.(properties%opt(4)=='Y')) properties%opt_l(4)=.true. !VLE Psat
                if((properties%opt(5)=='y').or.(properties%opt(5)=='Y')) properties%opt_l(5)=.true. !VLE V
				if((properties%opt(6)=='y').or.(properties%opt(6)=='Y')) properties%opt_l(6)=.true. !VLE binary  yichun
                       
                read(11,*)
                
                !Liquid volumes
                if(properties%opt_l(1)) then
                    read(11,*) properties%nliqv
                    
                    allocate(properties%t_liqv(1:properties%nliqv), stat=ierr)               
                    allocate(properties%p_liqv(1:properties%nliqv), stat=ierr)
                    allocate(properties%v_liq(1:properties%nliqv), stat=ierr)
                    allocate(properties%xiv(1:properties%nliqv,1:nctypes), stat=ierr)               
                    
                    do i1=1, properties%nliqv
                        read(11,*) properties%xiv(i1,1:nctypes),properties%t_liqv(i1),properties%p_liqv(i1), &
                        &   properties%v_liq(i1)
                    end do
                    
                    read(11,*)
                end if
                
                !Liquid heat capacities
                if(properties%opt_l(2)) then
                    read(11,*) properties%nliqc
                    read(11,*) properties%A, properties%B, properties%C, properties%D, properties%E
                    
                    allocate(properties%t_liqc(1:properties%nliqc), stat=ierr)               
                    allocate(properties%p_liqc(1:properties%nliqc), stat=ierr)
                    allocate(properties%xic(1:properties%nliqc,1:nctypes), stat=ierr)               
                    allocate(properties%cp(1:properties%nliqc), stat=ierr)
                    
                    do i1=1, properties%nliqc
                        read(11,*) properties%xic(i1,1:nctypes),properties%t_liqc(i1),properties%p_liqc(i1), &
                        &   properties%cp(i1)
                    end do
                    
                    read(11,*)
                end if
                
                !Liquid total enthalpies
                if(properties%opt_l(3)) then
                    read(11,*) properties%nliqh
                    
                    allocate(properties%t_liqh(1:properties%nliqh), stat=ierr)               
                    allocate(properties%p_liqh(1:properties%nliqh), stat=ierr)
                    allocate(properties%xih(1:properties%nliqh,1:nctypes), stat=ierr)               
                    allocate(properties%h(1:properties%nliqh), stat=ierr)
                    
                    do i1=1, properties%nliqh
                        read(11,*) properties%xih(i1,1:nctypes),properties%t_liqh(i1),properties%p_liqh(i1), &
                        &   properties%h(i1)
                    end do
                    
                    read(11,*)
                end if
                 
                if((properties%opt_l(4)).or.(properties%opt_l(5))) then
                    if(nctypes/=1) stop "Can only optimise pure VLE data at present!"
                
                    read(11,*) properties%nmu
                    allocate(properties%v_l(1:properties%nmu))
                    allocate(properties%v_v(1:properties%nmu))
                    allocate(properties%p_opt(1:properties%nmu))
                    allocate(properties%t_opt(1:properties%nmu))
                    allocate(properties%ximu(1:properties%nmu,1:nctypes))
                    
                    properties%v_l(:)   = 0.0e0_DP
                    properties%v_v(:)   = 0.0e0_DP
                    properties%p_opt(:) = 0.0e0_DP
                    properties%t_opt(:) = 0.0e0_DP
                    
                    do i1=1, properties%nmu
                        read(11,'(a100)') longline  
                        
                        if((properties%opt_l(4)).and.(properties%opt_l(5))) then
                            read(longline,*) properties%ximu(i1,1:nctypes), properties%t_opt(i1), properties%p_opt(i1), &
                                &   properties%v_l(i1), properties%v_v(i1)
                        else if(properties%opt_l(5)) then
                            read(longline,*)  properties%ximu(i1,1:nctypes), properties%t_opt(i1), properties%p_opt(i1),  &
                            &   properties%v_l(i1), properties%v_v(i1)
                        else if(properties%opt_l(4)) then    
                            read(longline,*)  properties%ximu(i1,1:nctypes), properties%t_opt(i1),  &
                            &   properties%v_l(i1), properties%v_v(i1)
                        else
                            print*, "Why am I in this loop?"
                        end if
                    end do
                end if
				
				!The followings are added by yichun
				if(properties%opt_l(6)) then
                    read(11,*) properties%nmu
                    
                    allocate(properties%t_opt(1:properties%nmu), stat=ierr)               
                    allocate(properties%p_opt(1:properties%nmu), stat=ierr)
                    allocate(properties%ximu(1:properties%nmu,1:nctypes), stat=ierr)               
                    allocate(properties%yimu(1:properties%nmu,1:nctypes), stat=ierr)
                    
                    do i1=1, properties%nmu
                        read(11,*) properties%t_opt(i1),properties%p_opt(i1),properties%ximu(i1,1:nctypes), &
                        &   properties%yimu(i1,1:nctypes)
                    end do
                    
                end if
!********************************************************************************************************************************                                               
            !Electrolyte optimiser
            case('eopt','EOPT')   
                
                read(11,*) properties%nliqv_e
                
                allocate(properties%xi_liqve(1:properties%nliqv_e,1:nctypes),    &
                &   properties%t_liqve(1:properties%nliqv_e),  &
                &   properties%p_liqve(1:properties%nliqv_e),  &
                &   properties%v_liqve(1:properties%nliqv_e))
                
                do i1=1, properties%nliqv_e
                    read(11,*) properties%t_liqve(i1),properties%p_liqve(i1),properties%xi_liqve(i1,1:nctypes), &
                    &   properties%v_liqve(i1)
                end do
                
                read(11,*)
                read(11,*) properties%ndhmix
                
                allocate(properties%xidhmixi(1:properties%ndhmix,1:nctypes), &
                &   properties%xidhmixf(1:properties%ndhmix,1:nctypes), &
                &   properties%tdhmix(1:properties%ndhmix), properties%pdhmix(1:properties%ndhmix), &
                &   properties%dhmix(1:properties%ndhmix))
                
                do i1=1, properties%ndhmix
                    read(11,*) properties%tdhmix(i1), properties%pdhmix(i1),  &
                    &   properties%xidhmixi(i1,1:nctypes), properties%xidhmixf(i1,1:nctypes), properties%dhmix(i1)
                end do
                
!********************************************************************************************************************************
            !Ionic liquid optimiser
            case('ilopt','ILOPT')   
                
                read(11,*) properties%nliqv
                
                allocate(properties%xiv(1:properties%nliqv,1:nctypes),    &
                &   properties%t_liqv(1:properties%nliqv),  &
                &   properties%p_liqv(1:properties%nliqv),  &
                &   properties%v_liq(1:properties%nliqv))
                
                do i1=1, properties%nliqv
                    read(11,*) properties%t_liqv(i1),properties%p_liqv(i1),properties%xiv(i1,1:nctypes), &
                    &   properties%v_liq(i1)
                end do            
                
                !Liquid heat capacities
                read(11,*)
                read(11,*) properties%nliqc
                    
                allocate(properties%t_liqc(1:properties%nliqc), stat=ierr)               
                allocate(properties%p_liqc(1:properties%nliqc), stat=ierr)
                allocate(properties%xic(1:properties%nliqc,1:nctypes), stat=ierr)               
                allocate(properties%cp(1:properties%nliqc), stat=ierr)
                
                do i1=1, properties%nliqc
                    read(11,*) properties%t_liqc(i1),properties%p_liqc(i1),properties%xic(i1,1:nctypes), &
                    &   properties%cp(i1)
                end do
!********************************************************************************************************************************             
            case default
                print*, "Unrecognised property ", properties%type
                stop
        end select
!*************************************************************************************************** 
!Read in / setup up constant arrays
        !coefficients in [1] A17 (effective packing fraction)
        c_array(1, 1) = 0.81096e0_DP
        c_array(1, 2) = 1.7888e0_DP
        c_array(1, 3) = -37.578e0_DP
        c_array(1, 4) = 92.284e0_DP
        
        c_array(2, 1) = 1.0205e0_DP
        c_array(2, 2) = -19.341e0_DP
        c_array(2, 3) = 151.26e0_DP
        c_array(2, 4) = -463.50e0_DP
        
        c_array(3, 1) = -1.9057e0_DP
        c_array(3, 2) = 22.845e0_DP
        c_array(3, 3) = -228.14e0_DP
        c_array(3, 4) = 973.92e0_DP
        
        c_array(4, 1) = 1.0885e0_DP
        c_array(4, 2) = -6.1962e0_DP
        c_array(4, 3) = 106.98e0_DP
        c_array(4, 4) = -677.64e0_DP
        
        !coefficients in [1] A26 (part of monomer contribution)
        phik      = 0.e0_DP
        
        phik(1, 0) = 7.5365557e0_DP
        phik(1, 1) = -37.60463e0_DP
        phik(1, 2) = 71.745953e0_DP
        phik(1, 3) = -46.83552e0_DP
        phik(1, 4) = -2.467982e0_DP
        phik(1, 5) = -0.50272e0_DP
        phik(1, 6) = 8.0956883e0_DP
        
        phik(2, 0) = -359.44e0_DP
        phik(2, 1) = 1825.6e0_DP
        phik(2, 2) = -3168.0e0_DP
        phik(2, 3) = 1884.2e0_DP
        phik(2, 4) = -0.82376e0_DP
        phik(2, 5) = -3.1935e0_DP
        phik(2, 6) = 3.7090e0_DP
        
        phik(3, 0) = 1550.9e0_DP
        phik(3, 1) = -5070.1e0_DP
        phik(3, 2) = 6534.6e0_DP
        phik(3, 3) = -3288.7e0_DP
        phik(3, 4) = -2.7171e0_DP
        phik(3, 5) = 2.0883e0_DP
        phik(3, 6) = 0.0000e0_DP
        
        phik(4, 0) = -1.19932e0_DP
        phik(4, 1) = 9.063632e0_DP
        phik(4, 2) = -17.9482e0_DP
        phik(4, 3) = 11.34027e0_DP
        phik(4, 4) = 20.52142e0_DP
        phik(4, 5) = -56.6377e0_DP
        phik(4, 6) = 40.53683e0_DP
        
        phik(5, 0) = -1911.28e0_DP
        phik(5, 1) = 21390.175e0_DP
        phik(5, 2) = -51320.7e0_DP
        phik(5, 3) = 37064.54e0_DP
        phik(5, 4) = 1103.742e0_DP
        phik(5, 5) = -3264.61e0_DP
        phik(5, 6) = 2556.181e0_DP
        
        phik(6, 0) = 9236.9e0_DP
        phik(6, 1) = -129430e0_DP
        phik(6, 2) = 357230e0_DP
        phik(6, 3) = -315530e0_DP
        phik(6, 4) = 1390.2e0_DP
        phik(6, 5) = -4518.2e0_DP
        phik(6, 6) = 4241.6e0_DP
        
        phik(7, 0) = 10.e0_DP
        phik(7, 1) = 10.e0_DP
        phik(7, 2) = 0.57e0_DP
        phik(7, 3) = -6.7e0_DP
        phik(7, 4) = -8.e0_DP
        phik(7, 5) = 0.e0_DP
        phik(7, 6) = 0.e0_DP
        
        !coefficients in [2] 24  (part of association,  parameters originally from [3])
        c_x=0.e0_DP
        
        c_x(0, 0) = 7.56425183020431e-02_DP
        c_x(0, 1) = -1.28667137050961e-01_DP
        c_x(0, 2) = 1.28350632316055e-01_DP
        c_x(0, 3) = -7.25321780970292e-02_DP
        c_x(0, 4) = 2.57782547511452e-02_DP
        c_x(0, 5) = -6.01170055221687e-03_DP
        c_x(0, 6) = 9.33363147191978e-04_DP
        c_x(0, 7) = -9.55607377143667e-05_DP
        c_x(0, 8) = 6.19576039900837e-06_DP
        c_x(0, 9) = -2.30466608213628e-07_DP
        c_x(0, 10) = 3.74605718435540e-09_DP
        c_x(1, 0) = 1.34228218276565e-01_DP
        c_x(1, 1) = -1.82682168504886e-01_DP
        c_x(1, 2) = 7.71662412959262e-02_DP
        c_x(1, 3) = -7.17458641164565e-04_DP
        c_x(1, 4) = -8.72427344283170e-03_DP
        c_x(1, 5) = 2.97971836051287e-03_DP
        c_x(1, 6) = -4.84863997651451e-04_DP
        c_x(1, 7) = 4.35262491516424e-05_DP
        c_x(1, 8) = -2.07789181640066e-06_DP
        c_x(1, 9) = 4.13749349344802e-08_DP
        c_x(2, 0) = -5.65116428942893e-01_DP
        c_x(2, 1) = 1.00930692226792e+00_DP
        c_x(2, 2) = -6.60166945915607e-01_DP
        c_x(2, 3) = 2.14492212294301e-01_DP
        c_x(2, 4) = -3.88462990166792e-02_DP
        c_x(2, 5) = 4.06016982985030e-03_DP
        c_x(2, 6) = -2.39515566373142e-04_DP
        c_x(2, 7) = 7.25488368831468e-06_DP
        c_x(2, 8) = -8.58904640281928e-08_DP
        c_x(3, 0) = -3.87336382687019e-01_DP
        c_x(3, 1) = -2.11614570109503e-01_DP
        c_x(3, 2) = 4.50442894490509e-01_DP
        c_x(3, 3) = -1.76931752538907e-01_DP
        c_x(3, 4) = 3.17171522104923e-02_DP
        c_x(3, 5) = -2.91368915845693e-03_DP
        c_x(3, 6) = 1.30193710011706e-04_DP
        c_x(3, 7) = -2.14505500786531e-06_DP
        c_x(4, 0) = 2.13713180911797e+00_DP
        c_x(4, 1) = -2.02798460133021e+00_DP
        c_x(4, 2) = 3.36709255682693e-01_DP
        c_x(4, 3) = 1.18106507393722e-03_DP
        c_x(4, 4) = -6.00058423301506e-03_DP
        c_x(4, 5) = 6.26343952584415e-04_DP
        c_x(4, 6) = -2.03636395699819e-05_DP
        c_x(5, 0) = -3.00527494795524e-01_DP
        c_x(5, 1) = 2.89920714512243e+00_DP
        c_x(5, 2) = -5.67134839686498e-01_DP
        c_x(5, 3) = 5.18085125423494e-02_DP
        c_x(5, 4) = -2.39326776760414e-03_DP
        c_x(5, 5) = 4.15107362643844e-05_DP
        c_x(6, 0) = -6.21028065719194e+00_DP
        c_x(6, 1) = -1.92883360342573e+00_DP
        c_x(6, 2) = 2.84109761066570e-01_DP
        c_x(6, 3) = -1.57606767372364e-02_DP
        c_x(6, 4) = 3.68599073256615e-04_DP
        c_x(7, 0) = 1.16083532818029e+01_DP
        c_x(7, 1) = 7.42215544511197e-01_DP
        c_x(7, 2) = -8.23976531246117e-02_DP
        c_x(7, 3) = 1.86167650098254e-03_DP
        c_x(8, 0) = -1.02632535542427e+01_DP
        c_x(8, 1) = -1.25035689035085e-01_DP
        c_x(8, 2) = 1.14299144831867e-02_DP
        c_x(9, 0) = 4.65297446837297e+00_DP
        c_x(9, 1) = -1.92518067137033e-03_DP
        c_x(10, 0) = -8.67296219639940e-01_DP
        
!Gaussian Legendre setup for [1] A9 integral
        call GL_arrays( )
!***************************************************************************************************

!*************************************************************************************************** 
!Housekeeping
        close(11)   
!***************************************************************************************************  
        return
    end subroutine Read_input        
!***************************************************************************************************
end module Input_mod
!***************************************************************************************************
!   SAFT Module to:
!       1. READ OPTIMISATION INPUT
!***************************************************************************************************
!   Created by SJ Halstead Sept 2016
!
!
!***************************************************************************************************
!
!***************************************************************************************************
module Input_opt_mod
!***************************************************************************************************
!MODULES
use Global_mod      ! CONTAINS GLOBAL VARIABLES
use Types_mod       ! Definition of types and numerical precision
use Zig_mod         ! RNG
!***************************************************************************************************
    implicit none
!***************************************************************************************************
contains
!***************************************************************************************************
    subroutine Read_opt( n_param, max_param, min_param, in_vals )
!--------------------------------------------------------------------------
        implicit none
!--------------------------------------------------------------------------        
!Variables
        integer, intent(out)                    ::  n_param
        real(kind=DP), allocatable, intent(out) ::  max_param(:), min_param(:)    
        real(kind=DP), allocatable, intent(out) ::  in_vals(:,:) 
        
        integer                                 ::  cinter1, cinter2
        character(len=80)                       ::  line_l
        character(len=50)                       ::  fileopt 

        integer                                 ::  ierr
!seed
        integer                     ::  seed
        integer, dimension(1:8)     ::  info            
!--------------------------------------------------------------------------
!Open optimise input        
        call getarg(2,fileopt)
        
        open(11, file=fileopt, status="old", action="read", iostat=ierr)
    
        if(ierr/=0) then
            print*, "Failed to open ",fileopt
            stop
        end if
!--------------------------------------------------------------------------        
!Nparams to fit
        read(11,*) n_param
        opt_num = n_param
        allocate( max_param(1:n_param), min_param(1:n_param), &
        &   param_key(1:n_param), param_index(1:n_param,1:2), param_index2(1:n_param,1:2))
		allocate(min_num(1:opt_num), max_num(1:opt_num), init_values(1:opt_num))
!--------------------------------------------------------------------------       
!Read params
!key:
! 1 sig ii; 2 eps ii; 3 lr ii; 4 la ii; 5 sf; 6 nseg; 7 eps ij; 8 lr ij; 9 eHB; 10 kHB  
        do cinter1=1, n_param       
            read(11,'(a80)') line_l
            read(line_l,*) cinter2
        
            if((cinter2==1)) then
                read(line_l,*) param_key(cinter1), param_index(cinter1,1), &
                &   min_param(cinter1), max_param(cinter1), init_values(cinter1)
				min_param(cinter1)=min_param(cinter1)*1000					!1000 is to make NAG get converged
				max_param(cinter1)=max_param(cinter1)*1000					!1000 is to make NAG get converged
				init_values(cinter1)=init_values(cinter1)*1000
			else if ((cinter2==2)) then
                read(line_l,*) param_key(cinter1), param_index(cinter1,1), &
                &   min_param(cinter1), max_param(cinter1), init_values(cinter1)
				min_param(cinter1)=min_param(cinter1)*10					!10 is to make NAG get converged
				max_param(cinter1)=max_param(cinter1)*10					!10 is to make NAG get converged
				init_values(cinter1)=init_values(cinter1)*10
			else if ((cinter2==3)) then
                read(line_l,*) param_key(cinter1), param_index(cinter1,1), &
                &   min_param(cinter1), max_param(cinter1), init_values(cinter1)
				min_param(cinter1)=min_param(cinter1)*100					!1000 is to make NAG get converged
				max_param(cinter1)=max_param(cinter1)*100					!1000 is to make NAG get converged
				init_values(cinter1)=init_values(cinter1)*100
			else if ((cinter2==4)) then
                read(line_l,*) param_key(cinter1), param_index(cinter1,1), &
                &   min_param(cinter1), max_param(cinter1), init_values(cinter1)
			else if ((cinter2==5)) then
                read(line_l,*) param_key(cinter1), param_index(cinter1,1), &
                &   min_param(cinter1), max_param(cinter1), init_values(cinter1)
				min_param(cinter1)=min_param(cinter1)*10000					!10000 is to make NAG get converged
				max_param(cinter1)=max_param(cinter1)*10000					!10000 is to make NAG get converged
				init_values(cinter1)=init_values(cinter1)*10000
			else if ((cinter2==6)) then
                read(line_l,*) param_key(cinter1), param_index(cinter1,1), &
                &   min_param(cinter1), max_param(cinter1), init_values(cinter1)
            else if((cinter2==7))then
                read(line_l,*) param_key(cinter1), param_index(cinter1,1:2),&
                &   min_param(cinter1), max_param(cinter1), init_values(cinter1)
				min_param(cinter1)=min_param(cinter1)*10					!10 is to make NAG get converged
				max_param(cinter1)=max_param(cinter1)*10					!10 is to make NAG get converged
				init_values(cinter1)=init_values(cinter1)*10
            else if((cinter2==8))then
                read(line_l,*) param_key(cinter1), param_index(cinter1,1:2),&
                &   min_param(cinter1), max_param(cinter1), init_values(cinter1)
				min_param(cinter1)=min_param(cinter1)*100					!1000 is to make NAG get converged
				max_param(cinter1)=max_param(cinter1)*100					!1000 is to make NAG get converged
				init_values(cinter1)=init_values(cinter1)*100
            else if((cinter2==9))then
                read(line_l,*) param_key(cinter1), param_index(cinter1,1:2),&
                &   param_index2(cinter1,1:2), min_param(cinter1), max_param(cinter1), init_values(cinter1)
			else if((cinter2==10))then
                read(line_l,*) param_key(cinter1), param_index(cinter1,1:2),&
                &   param_index2(cinter1,1:2), min_param(cinter1), max_param(cinter1), init_values(cinter1)	
				min_param(cinter1)=min_param(cinter1)*10					!10 is to make NAG get converged
				max_param(cinter1)=max_param(cinter1)*10					!10 is to make NAG get converged
				init_values(cinter1)=init_values(cinter1)*10
            else
                stop "parameter not recognised"
            end if
			min_num(cinter1)=min_param(cinter1)
			max_num(cinter1)=max_param(cinter1)
			
        end do

        close(11)
!--------------------------------------------------------------------------
!Set rng
        call date_and_time( values=info )
        seed=mod(((info(3) + info(7) + info(1)) * info(6) * info(5) +   &
        &   info(8) + info(2)), info(8) * info(8))      
        call zigset(seed) !SEED ZIG RNG
!--------------------------------------------------------------------------
!Set simplex input
        allocate( in_vals(1:n_param + 1, 1:n_param))

        do cinter1=1, n_param+1
        do cinter2=1, n_param
            in_vals(cinter1, cinter2) = (max_param(cinter2) -        &
            &       min_param(cinter2)) * uni() + min_param(cinter2)
        end do
        end do              
!--------------------------------------------------------------------------  
        return
    end subroutine Read_opt        
!--------------------------------------------------------------------------    

!***************************************************************************************************
end module Input_opt_mod
!***************************************************************************************************

module Setup_mod
!***************************************************************************************************
!MODULES
use Global_mod     !CONTAINS GLOBAL VARIABLES
use Types_mod
!***************************************************************************************************
    implicit none
    
    integer             ::  i_set1, i_set2, i_set3, i_set4
!***************************************************************************************************
    public  ::  Sys_setup
    private ::  Debroglie_set, Nm_set, Rho_set, Xsi_set, Cij_set,    &
    &   Dij_set, Zetal_set, Zetax_set, Alpha_set, F_set, Chain_set
!***************************************************************************************************
contains
!***************************************************************************************************
    subroutine Sys_setup( is_mu_in )
        implicit none
        
        logical, optional, intent(in)   :: is_mu_in
        logical                         :: is_mu
        
        is_mu = .false.
        if(present(is_mu_in)) is_mu = is_mu_in
        
        call Debroglie_set()         !thermal debroglie wavelengths and beta     
        if(is_mu) then
            call Xi_set()            !determine xi of each molecule
        else
            call Nm_set()                !determine nm of each molecule       
        end if
        call Rho_set()               !Number density and spherical segment density         
        call Xsi_set()               !mole fractions of components
        call Cij_set()               !Cij values  and on
        call Dij_set()               !Dij effective diameters 
        call Zetal_set               !Moments of number density        
        call Zetax_set()             !Packing fractions   
        call Alpha_set()             !for mono
        call F_set()                 !for mono
        call Chain_set()             !for chains
        
        return
    end subroutine Sys_setup
!***************************************************************************************************
    subroutine Debroglie_set() !DEGROGLIE WAVELENGTH
        implicit none
        
        real(kind=DP)       ::  db_local
        
        !Set beta values, depending on T (2nd is because energy values in units of K)
        beta    = 1/(KB * t)
        beta_K  = 1/t
        
        do i_set1=1,nctypes
            db_local = H / dsqrt( TWOPI * Comp_array(i_set1)%m * KB * t)    
            Comp_array(i_set1)%db3 = db_local**3.e0_DP                     
        end do

        return
    end subroutine Debroglie_set
!*************************************************************************************************** 
    subroutine  Xi_set() !Calculates mole fraction - for Mu calc only
        implicit none   
        
        real(kind=DP)       ::  nm_tot

        nm_tot = 0.e0_DP
                       
        do i_set1=1, nctypes
            nm_tot = nm_tot + Comp_array(i_set1)%nm
        end do
        
        do i_set1=1, nctypes
            Comp_array(i_set1)%xi = Comp_array(i_set1)%nm / nm_tot
        end do
 
        nsum=nm_tot
    
        return
    end subroutine Xi_set
!***************************************************************************************************  
    subroutine  Nm_set() !Equates number of molecules to mole fraction
        implicit none   
                       
        do i_set1=1,nctypes    
            Comp_array(i_set1)%nm = Comp_array(i_set1)%xi * NA
        end do
        
        return
    end subroutine Nm_set
!*************************************************************************************************** 
    subroutine  Rho_set()     
        implicit none

        real(kind=DP)    ::  trho, tsum
      
        tsum = sum(Comp_array(:)%nm)

        !Density moles per m^3
        rho = 0.0e0_DP       
   
        do i_set1=1, nctypes
            trho = Comp_array(i_set1)%nm / v / NA

            rho = rho + trho
            Comp_array(i_set1)%rho = trho 
        end do
  
        !Density of spherical segments moles per m^3
        rhos = 0.e0_DP 

        do i_set1=1,nctypes
        do i_set2=1,nstypes
            rhos = rhos + Comp_array(i_set1)%xi * Comp_array(i_set1)%comp(i_set2) &  
            &   * Seg_array(i_set2)%sf * Seg_array(i_set2)%nseg          
        end do
        end do

        rhos = rhos * rho   !moles per m^3

        return
    end subroutine Rho_set
!***************************************************************************************************
    subroutine Xsi_set()     ![1] A8  
        implicit none
        
        real(kind=DP)    ::  xsi_tot, xsi_sum
        
        xsi_tot=0.0e0_DP
        
        do i_set1=1, nctypes
        do i_set2=1, nstypes
            xsi_tot = xsi_tot + Comp_array(i_set1)%comp(i_set2) * Seg_array(i_set2)%nseg *Comp_array(i_set1)%xi &
            &   * Seg_array(i_set2)%sf
        end do
        end do
        
        do i_set1=1, nstypes
        
            xsi_sum=0.0e0_DP
            
            do i_set2=1, nctypes
                xsi_sum = xsi_sum +  (Comp_array(i_set2)%comp(i_set1) * Seg_array(i_set1)%sf *      &  !Including shape factor!
                &                  Seg_array(i_set1)%nseg * Comp_array(i_set2)%xi) / xsi_tot
            end do
            
        !mole fraction for each segments
            Seg_array(i_set1)%xsi = xsi_sum           
        end do   
  
        return
    end subroutine Xsi_set
!***************************************************************************************************
    subroutine Cij_set() ![1] A3 Mie correction scaling factor
        implicit none
        
        if(allocated(cij)) deallocate(cij)
        allocate(cij(1:nstypes, 1:nstypes))
        
        do i_set1=1, nstypes
        do i_set2=1, nstypes
            cij(i_set1, i_set2) = ( lr(i_set1, i_set2) / (lr(i_set1, i_set2) - la(i_set1, i_set2)) ) * (lr(i_set1, i_set2) /   &
            &   la(i_set1, i_set2)) ** (la(i_set1, i_set2) / (lr(i_set1, i_set2) - la(i_set1, i_set2)))
        end do
        end do
      
        return 
    end subroutine Cij_set    
!***************************************************************************************************
    subroutine Dij_set()     ![1] A9 effective diameter of Mie segments 
        implicit none

        real(kind=DP)           ::  lg_sum,mie_x,ba2
        real(kind=DP)           ::  t_eps,t_sig,t_cij,t_lr,t_la
        
        if(allocated(dij)) deallocate(dij)
        allocate(dij(1:nstypes, 1:nstypes))

        !Gaussian-Legendre sum
        do i_set1=1, nstypes
            t_eps = eps(i_set1, i_set1)
            t_sig = sig(i_set1, i_set1)
            
            t_cij = cij(i_set1, i_set1)
            t_lr  = lr(i_set1, i_set1) 
            t_la  = la(i_set1, i_set1)
            
            ba2 = t_sig / 2.0e0_DP     !G-L limit
            lg_sum = 0.0e0_DP
     
            do i_set2=1, N_GL
                  mie_x = ba2 * (X_GL(i_set2) + 1.0e0_DP)              
                  lg_sum = lg_sum + W_GL(i_set2) * (1.0e0_DP - &
                  & dexp( -beta_K * Mieii( mie_x, t_sig, t_eps, t_cij, t_lr, t_la)))
            end do

            dij(i_set1, i_set1) = ba2 * lg_sum      !In metres               
        end do
  
        do i_set1=1, nstypes
        do i_set2=1, nstypes      
            dij(i_set1, i_set2) = ( dij(i_set1, i_set1) + dij(i_set2, i_set2)) / 2.0e0_DP    !1 A46           
        end do
        end do

        return
    end subroutine Dij_set
!***************************************************************************************************
    subroutine Zetal_set()     ![1]  A7 Moments of number density
        implicit none
        
        real(kind=DP)           ::  t_dii

        zl = 0.0e0_DP
        
        do i_set1=1,nstypes
            t_dii = dij(i_set1,i_set1)     
           
            zl(0) = zl(0) + Seg_array(i_set1)%xsi
            zl(1) = zl(1) + Seg_array(i_set1)%xsi * t_dii
            zl(2) = zl(2) + Seg_array(i_set1)%xsi * t_dii**2.0e0_DP
            zl(3) = zl(3) + Seg_array(i_set1)%xsi * t_dii**3.0e0_DP          
        end do

        zl(0) = zl(0) * PI * rhos / 6.0e0_DP * NA              
        zl(1) = zl(1) * PI * rhos / 6.0e0_DP * NA              
        zl(2) = zl(2) * PI * rhos / 6.0e0_DP * NA            
        zl(3) = zl(3) * PI * rhos / 6.0e0_DP * NA            

        return
    end subroutine Zetal_set
!***************************************************************************************************
    subroutine Zetax_set()         ![1] A13 packing fractions
        implicit none

        zetax = 0.0e0_DP
        
        do i_set1=1, nstypes
        do i_set2=1, nstypes
            zetax = zetax + Seg_array(i_set1)%xsi * Seg_array(i_set2)%xsi * dij(i_set1,i_set2)**3.0e0_DP
        end do
        end do

        zetax_sum = zetax      

        zetax = zetax * PI * rhos / 6.0e0_DP * NA
      
        KHS = ((1.d0 - zetax)**4.0e0_DP) / (1.d0 + 4.0e0_DP * zetax + 4.0e0_DP * zetax**2 - 4.0e0_DP *  &
        &       zetax**3.0e0_DP + zetax**4.0e0_DP)  ![1] A21 isothermal compressibility

        zetabar = 0.0e0_DP

        do i_set1=1, nstypes
        do i_set2=1, nstypes
            zetabar = zetabar + Seg_array(i_set1)%xsi * Seg_array(i_set2)%xsi * sig(i_set1,i_set2)**3.0e0_DP 
        end do
        end do
        !Averaged packing fractions
        zetabar = NA * zetabar * PI * rhos / 6.0e0_DP

    end subroutine Zetax_set         
!***************************************************************************************************
    subroutine Alpha_set()  ![1] A24
        implicit none
        
        if(allocated(alpha)) deallocate(alpha)
        allocate(alpha(1:nstypes, 1:nstypes))
        
        do i_set1=1, nstypes
        do i_set2=1, nstypes
            alpha(i_set1, i_set2) = cij(i_set1, i_set2) * (1.0e0_DP /( la(i_set1, i_set2)   &
            &   - 3.0e0_DP) - 1.0e0_DP / (lr(i_set1, i_set2) - 3.0e0_DP)) 
        end do
        end do
        
        return
    end subroutine  Alpha_set
!**************************************************************************************************
    subroutine F_set()    ![1] A26
        implicit none
        
        double precision    ::  alpha_l,f_numer,f_denom
           
        if(allocated(fkij)) deallocate(fkij)           
        allocate(fkij(1:6, 1:nstypes, 1:nstypes))
        
        do i_set1=1, nstypes
        
            do i_set2=1, nstypes
                alpha_l = alpha(i_set1, i_set2)
                
                do i_set3=1, 6    
                
                    f_numer = 0.0e0_DP
                    f_denom = 0.0e0_DP
                    
                do i_set4 = 0, 3
                    f_numer = f_numer + phik(i_set3, i_set4) * alpha_l**i_set4
                end do
                    
                    do i_set4=4, 6  
                        f_denom = f_denom + phik(i_set3, i_set4) * alpha_l**(i_set4-3)
                    end do
                    
                    fkij(i_set3, i_set1, i_set2) = f_numer / (1.0e0_DP + f_denom)
                
                end do
            end do
        end do
         
        return
    end subroutine F_set
!**************************************************************************************************    
    subroutine Chain_set() !Set up values for chain calculation
        implicit none

        real(kind=DP)    ::  zki(1:nctypes, 1:nstypes), zsum
        
        zki = 0.0e0_DP ![2] 13
        
        do i_set1=1, nctypes       
            zsum = 0.0e0_DP
        
            do i_set2=1, nstypes
                zsum = zsum + Comp_array(i_set1)%comp(i_set2) * Seg_array(i_set2)%nseg * Seg_array(i_set2)%sf
            end do
            
            do i_set2=1, nstypes
                zki(i_set1, i_set2) =  Comp_array(i_set1)%comp(i_set2) * Seg_array(i_set2)%nseg *       &
                    &        Seg_array(i_set2)%sf / zsum   
            end do
        end do
        
        if(allocated(sig3ij)) deallocate(sig3ij,d3ij,epschij,lachij,lrchij)  
        allocate(sig3ij(1:nctypes, 1:nctypes), d3ij(1:nctypes) ,                  &
        &       epschij(1:nctypes, 1:nctypes), lachij(1:nctypes), lrchij(1:nctypes))
        
        !All averaged parameters for molecule
        sig3ij  = 0.0e0_DP     ![2] 14  
        d3ij    = 0.0e0_DP     ![2] 15
        epschij = 0.0e0_DP     ![2] 16
        lachij  = 0.0e0_DP     ![2] 17
        lrchij  = 0.0e0_DP
        sigx3   = 0.0e0_DP     ![2] 25 
        
        do i_set1=1, nctypes   !Form molecular averages
            do i_set2=1, nstypes
            do i_set3=1, nstypes
                sig3ij(i_set1, i_set1)= sig3ij(i_set1, i_set1) + zki(i_set1, i_set2) *      &
                &                       zki(i_set1, i_set3) * sig(i_set2, i_set3)**3.0e0_DP

                d3ij(i_set1)          = d3ij(i_set1) + zki(i_set1, i_set2) * zki(i_set1, i_set3) *  &
                &                       dij(i_set2, i_set3)**3.0e0_DP            
                epschij(i_set1, i_set1) = epschij(i_set1,i_set1) + zki(i_set1, i_set2) * zki(i_set1, i_set3) *  &
                &                         eps(i_set2, i_set3)
                lachij(i_set1)          = lachij(i_set1) + zki(i_set1, i_set2) * zki(i_set1, i_set3) *      &
                &                         la(i_set2, i_set3)
                lrchij(i_set1)          = lrchij(i_set1) + zki(i_set1, i_set2) * zki(i_set1, i_set3) *      &
                &                         lr(i_set2, i_set3)              
            end do
            end do           
        end do

        do i_set1=1, nctypes   !Estimate intermolecular averages (for Aassoc)
        do i_set2=1, nctypes
            if(i_set1/=i_set2) then
                sig3ij(i_set1, i_set2)   =   (sig3ij(i_set2, i_set2) + sig3ij(i_set1, i_set1)) /2.0e0_DP
                epschij(i_set1, i_set2)  =   dsqrt(sig3ij(i_set2, i_set2) * sig3ij(i_set1, i_set1)) / &
                    &   sig3ij(i_set1, i_set2)  * dsqrt(epschij(i_set2, i_set2) * epschij(i_set1, i_set1))
            end if
        end do
        end do
    
        do i_set1=1, nstypes
        do i_set2=1, nstypes
            sigx3 = sigx3 + Seg_array(i_set1)%xsi * Seg_array(i_set2)%xsi * sig(i_set1, i_set2)**3.0e0_DP
        end do
        end do

        return
    end subroutine Chain_set 
!***************************************************************************************************
    function Mieii(xm_in, sm_in, em_in, cm_in, lrm_in, lam_in) result(mie_out) !A2 for self-interaction  
        implicit none
        
        real(kind=DP),intent(in)    ::  xm_in,sm_in,em_in,cm_in,lrm_in,lam_in
        real(kind=DP)               ::  mie_out
        real(kind=DP)               ::  t_mie

        t_mie = sm_in / xm_in        
        mie_out = cm_in * em_in * (t_mie**lrm_in - t_mie**lam_in)

        return
    end function MIEII
!***************************************************************************************************

!***************************************************************************************************
end module Setup_mod
!***************************************************************************************************
!***************************************************************************************************
!   SAFT Module to:
!       1. Calculate A ideal
!***************************************************************************************************
!
!***************************************************************************************************
!***************************************************************************************************
!   SAFT Module to:
!       1. Calculate all zeff terms
!***************************************************************************************************
!
!***************************************************************************************************
module Zeff_mod
!***************************************************************************************************
!Modules
!=======    
    use Types_mod           ! Definitions of types and double precision
    use Global_mod          ! Important global parameters
!***************************************************************************************************
    implicit none
!***************************************************************************************************
    contains
!***************************************************************************************************
    pure function zeff( lam_in ) result(zeff_out)
        implicit none
        
        integer                         ::  i_zeff
        real(kind=DP),intent(in)        ::  lam_in
        real(kind=DP)                   ::  zeff_out
        real(kind=DP)                   ::  ca(1:4)
        
        do i_zeff=1, 4
            ca(i_zeff) = c_array(i_zeff, 1) * 1.0e0_DP + c_array(i_zeff, 2) * 1.0e0_DP /  &
            &   lam_in + c_array(i_zeff, 3) * 1.0e0_DP / lam_in**2.0e0_DP                 &
            &   + c_array(i_zeff, 4) * 1.0e0_DP / lam_in**3.0e0_DP
        end do        

        zeff_out = ca(1) * zetax + ca(2) * zetax**2.0e0_DP + ca(3) * zetax**3.0e0_DP +  &
        &   ca(4) * zetax**4.0e0_DP      

        return
    end function zeff  
!***************************************************************************************************    
    pure function dzeffn( lam_in ) result(zeff_out)
        implicit none
        
        integer                         ::  i_zeff
        real(kind=DP),intent(in)        ::  lam_in
        real(kind=DP)                   ::  zeff_out
        real(kind=DP)                   ::  ca(1:4)
        
        do i_zeff=1, 4
            ca(i_zeff) = c_array(i_zeff, 1) * 1.0e0_DP + c_array(i_zeff, 2) * 1.0e0_DP /  &
            &   lam_in + c_array(i_zeff, 3) * 1.0e0_DP / lam_in**2.0e0_DP                 &
            &   + c_array(i_zeff, 4) * 1.0e0_DP / lam_in**3.0e0_DP
        end do        

        zeff_out = ca(1) * dzetaxn + ca(2) * 2.0e0_DP*dzetaxn*zetax + ca(3) * 3.0e0_DP*dzetaxn*zetax**2.0e0_DP +  &
        &   ca(4) * 4.0e0_DP*dzetaxn*zetax**3.0e0_DP      

        return
    end function dzeffn
!*************************************************************************************************** 
    pure function dzeffv( lam_in ) result(zeff_out)
        implicit none
        
        integer                         ::  i_zeff
        real(kind=DP),intent(in)        ::  lam_in
        real(kind=DP)                   ::  zeff_out
        real(kind=DP)                   ::  ca(1:4)
        
        do i_zeff=1, 4
            ca(i_zeff) = c_array(i_zeff, 1) * 1.0e0_DP + c_array(i_zeff, 2) * 1.0e0_DP /  &
            &   lam_in + c_array(i_zeff, 3) * 1.0e0_DP / lam_in**2.0e0_DP                 &
            &   + c_array(i_zeff, 4) * 1.0e0_DP / lam_in**3.0e0_DP
        end do        

        zeff_out = ca(1) * dzetaxv + ca(2) * 2.0e0_DP*dzetaxv*zetax + ca(3) * 3.0e0_DP*dzetaxv*zetax**2.0e0_DP +  &
        &   ca(4) * 4.0e0_DP*dzetaxv*zetax**3.0e0_DP       

        return
    end function dzeffv  
!***************************************************************************************************
    pure function zeffprime( lamin ) result(zeffprime_out)
        implicit none
        
        integer                     ::  i_zeffp
        real(kind=DP), intent(in)   ::  lamin
        real(kind=DP)               ::  ca(1:4), l_zetax
        real(kind=DP)               ::  zeffprime_out
        
        l_zetax = zetax / rhos
      
        do i_zeffp=1, 4
            ca(i_zeffp) = c_array(i_zeffp, 1) * 1.0e0_DP + c_array(i_zeffp, 2) * 1.0e0_DP / lamin +     &
            &   c_array(i_zeffp, 3) * 1.0e0_DP / lamin**2.0e0_DP                                        &
            &   + c_array(i_zeffp, 4) * 1.0e0_DP / lamin**3.0e0_DP
        end do        

        zeffprime_out = ca(1) * l_zetax + 2.0e0_DP * zetax * ca(2) * l_zetax           +       &
        &   3.0e0_DP * zetax**2.0e0_DP * ca(3) * l_zetax + 4.0e0_DP * zetax**3.0e0_DP *       &
        &   ca(4) * l_zetax    

        return
    end function zeffprime  
!***************************************************************************************************
    pure function dzeffprimen( lamin ) result(zeffprime_out)
        implicit none
        
        integer                     ::  i_zeffp
        real(kind=DP), intent(in)   ::  lamin
        real(kind=DP)               ::  ca(1:4), l_zetax, dl_zetaxn
        real(kind=DP)               ::  zeffprime_out
        
        l_zetax   = zetax / rhos
        dl_zetaxn = (rhos * dzetaxn - drhosn * zetax) / rhos**2.0e0_DP
      
        do i_zeffp=1, 4
            ca(i_zeffp) = c_array(i_zeffp, 1) * 1.0e0_DP + c_array(i_zeffp, 2) * 1.0e0_DP / lamin +     &
            &   c_array(i_zeffp, 3) * 1.0e0_DP / lamin**2.0e0_DP                                        &
            &   + c_array(i_zeffp, 4) * 1.0e0_DP / lamin**3.0e0_DP
        end do        

        zeffprime_out = ca(1) * dl_zetaxn + 2.0e0_DP * ca(2) * (zetax*dl_zetaxn + dzetaxn*l_zetax)          + &
        &   3.0e0_DP * ca(3) * (zetax**2.0e0_DP*dl_zetaxn + 2.0e0_DP*zetax*dzetaxn*l_zetax)                 + &
        &   4.0e0_DP * ca(4) * (zetax**3.0e0_DP*dl_zetaxn + 3.0e0_DP*zetax**2.0e0_DP*dzetaxn*l_zetax)  

        return
    end function dzeffprimen 
!***************************************************************************************************
    pure function dzeffprimev( lamin ) result(zeffprime_out)
        implicit none
        
        integer                     ::  i_zeffp
        real(kind=DP), intent(in)   ::  lamin
        real(kind=DP)               ::  ca(1:4), l_zetax, dl_zetaxv
        real(kind=DP)               ::  zeffprime_out
        
        l_zetax   = zetax / rhos
        dl_zetaxv = (rhos * dzetaxv - drhosv * zetax) / rhos**2.0e0_DP
      
        do i_zeffp=1, 4
            ca(i_zeffp) = c_array(i_zeffp, 1) * 1.0e0_DP + c_array(i_zeffp, 2) * 1.0e0_DP / lamin +     &
            &   c_array(i_zeffp, 3) * 1.0e0_DP / lamin**2.0e0_DP                                        &
            &   + c_array(i_zeffp, 4) * 1.0e0_DP / lamin**3.0e0_DP
        end do        

        zeffprime_out = ca(1) * dl_zetaxv + 2.0e0_DP * ca(2) * (zetax*dl_zetaxv + dzetaxv*l_zetax)          + &
        &   3.0e0_DP * ca(3) * (zetax**2.0e0_DP*dl_zetaxv + 2.0e0_DP*zetax*dzetaxv*l_zetax)                 + &
        &   4.0e0_DP * ca(4) * (zetax**3.0e0_DP*dl_zetaxv + 3.0e0_DP*zetax**2.0e0_DP*dzetaxv*l_zetax)  

        return
    end function dzeffprimev 
!***************************************************************************************************
!***************************************************************************************************
end module Zeff_mod
!***************************************************************************************************
module Ideal_mod
!***************************************************************************************************
!Modules
!=======    
    use Types_mod           ! Definitions of types and double precision
    use Global_mod          ! Important global parameters 
!***************************************************************************************************
    implicit none
!***************************************************************************************************
    public  ::  A_ideal, A_ideal_dv, A_ideal_dn
!***************************************************************************************************
    contains
!***************************************************************************************************
    function A_ideal() result(a_res)     ![1]  A1      [2] 5
        implicit none
        
        integer         ::  i_a
        real(kind=DP)   ::  a_res
       
        a_res = -1.0e0_DP
        
        do i_a=1,nctypes  
            if(Comp_array(i_a)%xi/=0.0e0_DP) a_res = a_res + Comp_array(i_a)%xi * &
            & dlog(Comp_array(i_a)%rho * NA * Comp_array(i_a)%db3)   !NOTE [1] A1 misprint - xi outside LN                       
        end do

        return
    end function A_ideal
!***************************************************************************************************   
    function A_ideal_dn( i_part ) result(a_res)
        implicit none
        
        integer                 ::  i_a
        integer, intent(in)     ::  i_part
        real(kind=DP)           ::  a_res

        a_res = 0.0e0_DP

        do i_a=1,nctypes
            if(Comp_array(i_a)%xi/=0.0e0_DP) a_res = a_res + Comp_array(i_a)%dxin*dlog(Comp_array(i_a)%nm/v * Comp_array(i_a)%db3)     
        end do
        
        a_res = nsum *( a_res +  Comp_array(i_part)%xi * 1.0e0_DP/Comp_array(i_part)%nm ) +  A_ideal()   
 
        return
    end function A_ideal_dn       
!***************************************************************************************************         
    function A_ideal_dv() result(a_res)     ![1]  A1      [2] 5
        implicit none
        
        integer         ::  i_a
        real(kind=DP)   ::  a_res
       
        a_res = 0.0e0_DP
        
        do i_a=1,nctypes  
            if(Comp_array(i_a)%xi/=0.0e0_DP) a_res = a_res + Comp_array(i_a)%xi / v                                 
        end do    

        return
    end function A_ideal_dv
!***************************************************************************************************  
end module Ideal_mod
!***************************************************************************************************
!***************************************************************************************************
!   SAFT Module to:
!       1. Calculate A mono
!
!       Analytical derivatives added March 2017
!
!***************************************************************************************************
!
!***************************************************************************************************
module Mono_mod
!***************************************************************************************************
!Modules
!=======    
    use Types_mod           ! Definitions of types and double precision
    use Global_mod          ! Important global parameters 
    use Zeff_mod            ! Zeff calculations
!***************************************************************************************************
    implicit none
!***************************************************************************************************
    public  ::  A_mono, A_mono_dn, A_mono_dv
!***************************************************************************************************    
    contains
!***************************************************************************************************    
    function A_mono() result(a_res)     ![1]  A4   [2] 6 ish
        implicit none
        
        integer         ::  i_a1, i_a2
        real(kind=DP)   ::  a_res
   
        a_res = 0.0e0_DP
        
        do i_a1=1,nctypes
        do i_a2=1,nstypes    
            a_res = a_res + Comp_array(i_a1)%xi * Comp_array(i_a1)%comp(i_a2) * &
            &   Seg_array(i_a2)%nseg * Seg_array(i_a2)%sf
        end do
        end do

        a_res = a_res * am()

        return
    end function A_mono 
!***************************************************************************************************
    !hacked version to separate hs and dispersion contributions
    subroutine A_mono2(a_res, amono_hs, amono_disp)   ![1]  A4   [2] 6 ish
        implicit none
        
        integer                     ::  i_a1, i_a2
        real(kind=DP),intent(out)   ::  a_res, amono_hs, amono_disp
        real(kind=DP)               ::  am_hs, am_disp, am_full
        
   
        a_res = 0.0e0_DP
        
        do i_a1=1,nctypes
        do i_a2=1,nstypes    
            a_res = a_res + Comp_array(i_a1)%xi * Comp_array(i_a1)%comp(i_a2) * &
            &   Seg_array(i_a2)%nseg * Seg_array(i_a2)%sf
        end do
        end do

        am_full     = am(am_hs, am_disp)
        
        amono_hs    = a_res * am_hs
        amono_disp  = a_res * am_disp
        
        a_res = a_res * am_full

        return
    end subroutine A_mono2 
!***************************************************************************************************
    function A_mono_dn() result(a_res)     ![1]  A4   [2] 6 ish
        implicit none
        
        integer         ::  i_a1, i_a2
        real(kind=DP)   ::  a_res, mono_u, mono_up
  
        mono_u   = 0.0e0_DP
        mono_up  = 0.0e0_DP
        
        do i_a1=1,nctypes       
            do i_a2=1,nstypes    
                mono_u = mono_u + Comp_array(i_a1)%xi * Comp_array(i_a1)%comp(i_a2) * &
                &   Seg_array(i_a2)%nseg * Seg_array(i_a2)%sf
                
                mono_up = mono_up + Comp_array(i_a1)%dxin * Comp_array(i_a1)%comp(i_a2) * &
                &   Seg_array(i_a2)%nseg * Seg_array(i_a2)%sf
            end do       
        end do
        
        a_res = mono_u * damn() * nsum + mono_up * am() * nsum + mono_u * am()

        return
    end function A_mono_dn 
!***************************************************************************************************
    !hacked version to separate hs and dispersion contributions
    subroutine A_mono_dn2(a_res, amono_hs, amono_disp) ![1]  A4   [2] 6 ish
        implicit none
        
        integer         ::  i_a1, i_a2
        real(kind=DP)   ::  a_res, amono_hs, amono_disp, mono_u, mono_up
        real(kind=DP)   ::  am_hs, am_disp, am_full
        real(kind=DP)   ::  damn_hs, damn_disp, damn_full
  
        mono_u   = 0.0e0_DP
        mono_up  = 0.0e0_DP
        
        do i_a1=1,nctypes              
            do i_a2=1,nstypes    
                mono_u = mono_u + Comp_array(i_a1)%xi * Comp_array(i_a1)%comp(i_a2) * &
                &   Seg_array(i_a2)%nseg * Seg_array(i_a2)%sf
                
                mono_up = mono_up + Comp_array(i_a1)%dxin * Comp_array(i_a1)%comp(i_a2) * &
                &   Seg_array(i_a2)%nseg * Seg_array(i_a2)%sf
              
            end do       
        end do
      
        am_full   = am(am_hs, am_disp)
        damn_full = damn(damn_hs, damn_disp)
        
        amono_hs    = mono_u * damn_hs * nsum + mono_up * am_hs * nsum + mono_u * am_hs  
        amono_disp  = mono_u * damn_disp * nsum + mono_up * am_disp * nsum + mono_u * am_disp  
  
        a_res = mono_u * damn_full * nsum + mono_up * am_full * nsum + mono_u * am_full

        return
    end subroutine A_mono_dn2 
!***************************************************************************************************
    function A_mono_dv() result(a_res)     ![1]  A4   [2] 6 ish
        implicit none
        
        integer         ::  i_a1, i_a2
        real(kind=DP)   ::  a_res, mono_u, mono_up
  
        mono_u   = 0.0e0_DP
        mono_up  = 0.0e0_DP
        
        do i_a1=1,nctypes
        if(Comp_array(i_a1)%xi/=0.0e0_DP) then        
            do i_a2=1,nstypes    
                mono_u = mono_u + Comp_array(i_a1)%xi * Comp_array(i_a1)%comp(i_a2) * &
                &   Seg_array(i_a2)%nseg * Seg_array(i_a2)%sf
            end do
        end if        
        end do

        a_res = mono_u * damv() 

        return
    end function A_mono_dv
!***************************************************************************************************
    function am(hs_out, disp_out) result(am_out)      !1 A5
        implicit none
 
        real(kind=DP)           ::  am_out
        real(kind=DP),optional  ::  hs_out, disp_out  
        real(kind=DP)           ::  local_hs, local_disp     
      
        local_hs   = ahs()
        local_disp = beta_K * a1() * NA + beta_K**2.0e0_DP * a2()*NA + beta_K**3.0e0_DP * a3() 
        am_out = local_hs + local_disp

        if(present(hs_out)) hs_out=local_hs
        if(present(disp_out)) disp_out=local_disp

        return
    end function am
!***************************************************************************************************
    function damn(hs_out, disp_out) result(am_out)      !1 A5
        implicit none
 
        real(kind=DP)   ::  am_out
        real(kind=DP),optional  ::  hs_out, disp_out  
        real(kind=DP)           ::  local_hs, local_disp  
        
        local_hs   = dahsn()
        local_disp = beta_K * da1n() * NA + beta_K**2.0e0_DP * da2n()*NA + beta_K**3.0e0_DP * da3n()

        am_out = local_hs + local_disp
        
        if(present(hs_out)) hs_out=local_hs
        if(present(disp_out)) disp_out=local_disp

        return
    end function damn
!***************************************************************************************************
    function damv() result(am_out)      !1 A5
        implicit none
 
        real(kind=DP)   ::  am_out

        am_out = dahsv() + beta_K * da1v() * NA + beta_K**2.0e0_DP * da2v()*NA + beta_K**3.0e0_DP * da3v()

        return
    end function damv
!***************************************************************************************************
    function ahs() result(ahs_out)     !1 A6   2 Eq. 7
        implicit none
        
        real(kind=DP)   ::  ahs_out
        
        ahs_out = 6.0e0_DP / PI / (rhos * NA) * (                                           &
        &       (( zl(2)**3.0e0_DP / zl(3)**2.0e0_DP) - zl(0)) * dlog(1.0e0_DP - zl(3)) +   &
        &       3.0e0_DP * zl(1) * zl(2) / (1.0e0_DP - zl(3)) + zl(2)**3.0e0_DP /          &
        &       (zl(3) * (1.0e0_DP - zl(3))**2.0e0_DP)        )  

        return
    end function ahs
!***************************************************************************************************
    function dahsn() result(dahsn_out)     !1 A6   2 Eq. 7
        implicit none
        
        real(kind=DP)       ::  dahsn_out
        real(kind=DP)       ::  ltemp1,ltemp2,ltemp3,ltemp4,ltemp5,ltemp6
        
        ltemp1 = -6.0e0_DP*drhosn/PI/rhos**2.0e0_DP / NA
        ltemp2 = (3.0e0_DP*zl(3)*zl(2)**2.0e0_DP*dzln(2) - 2.0e0_DP*zl(2)**3.0e0_DP*dzln(3))     / &
        &       zl(3)**3.0e0_DP - dzln(0)   
        ltemp3 = -dzln(3) / (1.0e0_DP-zl(3))
        ltemp4 = 3.0e0_DP * ( (1.0e0_DP-zl(3)) * (zl(1)*dzln(2) + dzln(1)*zl(2)) + zl(1)*zl(2)*dzln(3))     / &
        &       (1.0e0_DP - zl(3))**2.0e0_DP
        ltemp6 = zl(3)*(-2.0e0_DP * (1.0e0_DP-zl(3)) * dzln(3)) + dzln(3)*(1.0e0_DP-zl(3))**2.0e0_DP
        ltemp5 = (3.0e0_DP * zl(3) *(1.0e0_DP-zl(3))**2.0e0_DP * zl(2)**2.0e0_DP * dzln(2)   - &
        &       zl(2)**3.0e0_DP * ltemp6) / (zl(3) * (1.0e0_DP - zl(3))**2.0e0_DP )**2.0e0_DP 
        
        dahsn_out = ltemp1 * (                                                                  &
        &       (( zl(2)**3.0e0_DP / zl(3)**2.0e0_DP) - zl(0)) * dlog(1.0e0_DP - zl(3))      +  &
        &       3.0e0_DP * zl(1) * zl(2) / (1.0e0_DP - zl(3)) + zl(2)**3.0e0_DP             /   &
        &       (zl(3) * (1.0e0_DP - zl(3))**2.0e0_DP)        )                             +   &  
        &       6.0e0_DP / PI / (rhos * NA) * (                                                 &
        &       ltemp3 * (( zl(2)**3.0e0_DP / zl(3)**2.0e0_DP) - zl(0))                     +   &
        &       dlog(1.0e0_DP - zl(3)) * ltemp2 + ltemp4 + ltemp5)
      
        return
    end function dahsn
!***************************************************************************************************
    function dahsv() result(dahsv_out)     !1 A6   2 Eq. 7
        implicit none
        
        real(kind=DP)       ::  dahsv_out
        real(kind=DP)       ::  ltemp1,ltemp2,ltemp3,ltemp4,ltemp5,ltemp6
        
        ltemp1 = -6.0e0_DP*drhosv/PI/rhos**2.0e0_DP / NA
        ltemp2 = (3.0e0_DP*zl(3)*zl(2)**2.0e0_DP*dzlv(2) - 2.0e0_DP*zl(2)**3.0e0_DP*dzlv(3))     / &
        &       zl(3)**3.0e0_DP - dzlv(0)   
        ltemp3 = -dzlv(3) / (1.0e0_DP-zl(3))
        ltemp4 = 3.0e0_DP * ( (1.0e0_DP-zl(3)) * (zl(1)*dzlv(2) + dzlv(1)*zl(2)) + zl(1)*zl(2)*dzlv(3))     / &
        &       (1.0e0_DP - zl(3))**2.0e0_DP
        ltemp6 = zl(3)*(-2.0e0_DP * (1.0e0_DP-zl(3)) * dzlv(3)) + dzlv(3)*(1.0e0_DP-zl(3))**2.0e0_DP
        ltemp5 = (3.0e0_DP * zl(3) *(1.0e0_DP-zl(3))**2.0e0_DP * zl(2)**2.0e0_DP * dzlv(2)   - &
        &       zl(2)**3.0e0_DP * ltemp6) / (zl(3) * (1.0e0_DP - zl(3))**2.0e0_DP )**2.0e0_DP 
        
        dahsv_out = ltemp1 * (                                                                  &
        &       (( zl(2)**3.0e0_DP / zl(3)**2.0e0_DP) - zl(0)) * dlog(1.0e0_DP - zl(3))      +  &
        &       3.0e0_DP * zl(1) * zl(2) / (1.0e0_DP - zl(3)) + zl(2)**3.0e0_DP             /   &
        &       (zl(3) * (1.0e0_DP - zl(3))**2.0e0_DP)        )                             +   &  
        &       6.0e0_DP / PI / (rhos * NA) * (                                                 &
        &       ltemp3 * (( zl(2)**3.0e0_DP / zl(3)**2.0e0_DP) - zl(0))                     +   &
        &       dlog(1.0e0_DP - zl(3)) * ltemp2 + ltemp4 + ltemp5)
      
        return
    end function dahsv
!***************************************************************************************************
    function a1() result(a1_out)
        implicit none
        
        real(kind=DP)   ::  a1_out
        integer         ::  i_a1_1, i_a1_2
        
        a1_out=0.0e0_DP
      
        do i_a1_1=1, nstypes
        do i_a1_2=1, nstypes                
            a1_out = a1_out + Seg_array(i_a1_1)%xsi * Seg_array(i_a1_2)%xsi * aij(i_a1_1,i_a1_2)     
        end do
        end do

        return
    end function a1
!***************************************************************************************************
    function da1n() result(a1_out)
        implicit none
        
        real(kind=DP)   ::  a1_out
        integer         ::  i_a1_1, i_a1_2
      
        a1_out=0.0e0_DP
      
        do i_a1_1=1, nstypes
        do i_a1_2=1, nstypes                
            a1_out = a1_out                                                                 + &
            &   Seg_array(i_a1_1)%dxsin * Seg_array(i_a1_2)%xsi * aij(i_a1_1,i_a1_2)        + & 
            &   Seg_array(i_a1_1)%xsi * Seg_array(i_a1_2)%dxsin * aij(i_a1_1,i_a1_2)        + &
            &   Seg_array(i_a1_1)%xsi * Seg_array(i_a1_2)%xsi * daijn(i_a1_1,i_a1_2)           
        end do
        end do

        return
    end function da1n
!***************************************************************************************************
    function da1v() result(a1_out)
        implicit none
        
        real(kind=DP)   ::  a1_out
        integer         ::  i_a1_1, i_a1_2
      
        a1_out=0.0e0_DP
      
        do i_a1_1=1, nstypes
        do i_a1_2=1, nstypes                
            a1_out = a1_out                                                                 + &
            &   Seg_array(i_a1_1)%xsi * Seg_array(i_a1_2)%xsi * daijv(i_a1_1,i_a1_2)           
        end do
        end do

        return
    end function da1v
!***************************************************************************************************
    function a2() result(a2_out)
        implicit none

        real(kind=DP)   ::  a2_out
        integer         ::  i_a2_1, i_a2_2
   
        a2_out=0.0e0_DP
      
        do i_a2_1=1, nstypes
        do i_a2_2=1, nstypes   
            a2_out = a2_out + Seg_array(i_a2_1)%xsi * Seg_array(i_a2_2)%xsi * as2ij(i_a2_1, i_a2_2)
        end do
        end do

        return        
    end function 
!***************************************************************************************************
    function da2n() result(a2_out)
        implicit none

        real(kind=DP)   ::  a2_out
        integer         ::  i_a2_1, i_a2_2
   
        a2_out=0.0e0_DP
      
        do i_a2_1=1, nstypes
        do i_a2_2=1, nstypes   
            a2_out = a2_out                                                                   + &
            &   Seg_array(i_a2_1)%dxsin * Seg_array(i_a2_2)%xsi * as2ij(i_a2_1, i_a2_2)       + &
            &   Seg_array(i_a2_1)%xsi * Seg_array(i_a2_2)%dxsin * as2ij(i_a2_1, i_a2_2)       + &
            &   Seg_array(i_a2_1)%xsi * Seg_array(i_a2_2)%xsi * das2ijn(i_a2_1, i_a2_2)       
        end do
        end do

        return        
    end function da2n
!***************************************************************************************************
    function da2v() result(a2_out)
        implicit none

        real(kind=DP)   ::  a2_out
        integer         ::  i_a2_1, i_a2_2
   
        a2_out=0.0e0_DP
      
        do i_a2_1=1, nstypes
        do i_a2_2=1, nstypes   
            a2_out = a2_out                                                                   + &
            &   Seg_array(i_a2_1)%xsi * Seg_array(i_a2_2)%xsi * das2ijv(i_a2_1, i_a2_2)       
        end do
        end do

        return        
    end function da2v
!***************************************************************************************************
    function a3() result(a3_out)    !A25
        implicit none
        
        real(kind=DP)   ::  a3_out
        integer         ::  i_a3_1, i_a3_2
    
        a3_out=0.0e0_DP
        
        do i_a3_1=1, nstypes
        do i_a3_2=1, nstypes
            a3_out = a3_out + Seg_array(i_a3_1)%xsi * Seg_array(i_a3_2)%xsi * as3ij(i_a3_1, i_a3_2)        
        end do
        end do

        return
    end function
!***************************************************************************************************
    function da3n() result(a3_out)    !A25
        implicit none
        
        real(kind=DP)   ::  a3_out
        integer         ::  i_a3_1, i_a3_2
   
        a3_out=0.0e0_DP
        
        do i_a3_1=1, nstypes
        do i_a3_2=1, nstypes
            a3_out = a3_out                                                                   + &
            &   Seg_array(i_a3_1)%dxsin * Seg_array(i_a3_2)%xsi * as3ij(i_a3_1, i_a3_2)       + &
            &   Seg_array(i_a3_1)%xsi * Seg_array(i_a3_2)%dxsin * as3ij(i_a3_1, i_a3_2)       + &
            &   Seg_array(i_a3_1)%xsi * Seg_array(i_a3_2)%xsi * das3ijn(i_a3_1, i_a3_2)       
        end do
        end do

        return
    end function da3n
!***************************************************************************************************
    function da3v() result(a3_out)    !A25
        implicit none
        
        real(kind=DP)   ::  a3_out
        integer         ::  i_a3_1, i_a3_2
   
        a3_out=0.0e0_DP
        
        do i_a3_1=1, nstypes
        do i_a3_2=1, nstypes
            a3_out = a3_out                                                                 + &
            &   Seg_array(i_a3_1)%xsi * Seg_array(i_a3_2)%xsi * das3ijv(i_a3_1, i_a3_2)       
        end do
        end do

        return
    end function da3v
!***************************************************************************************************
    function aij( i_aij1, i_aij2 ) result(aij_out)  !A11
        implicit none
        
        integer,intent(in)  ::  i_aij1, i_aij2
        real(kind=DP)       ::  xodij, lam1, lam2      
        real(kind=DP)       ::  aij_out
     
        lam1 = la(i_aij1, i_aij2)
        lam2 = lr(i_aij1, i_aij2)   
        xodij = sig(i_aij1, i_aij2) / dij(i_aij1, i_aij2)     !A11 +

        aij_out = cij(i_aij1, i_aij2) * ((xodij**la(i_aij1, i_aij2)) *      &
        &       (as1ij( i_aij1, i_aij2, lam1, zeff(lam1) )+                 &
        &       bij( i_aij1, i_aij2, lam1, xodij )) -                       &
        &       (xodij**lam2) *                                             &
        &       (as1ij( i_aij1, i_aij2, lam2, zeff(lam2) ) +                &
        &       bij( i_aij1, i_aij2, lam2, xodij) )        )  

        return
    end function
!***************************************************************************************************
    function daijn( i_aij1, i_aij2 ) result(aij_out)  !A11
        implicit none
        
        integer,intent(in)  ::  i_aij1, i_aij2
        real(kind=DP)       ::  xodij, lam1, lam2      
        real(kind=DP)       ::  aij_out
            
        lam1 = la(i_aij1, i_aij2)
        lam2 = lr(i_aij1, i_aij2)   
        xodij = sig(i_aij1, i_aij2) / dij(i_aij1, i_aij2)     !A11 +

        aij_out = cij(i_aij1, i_aij2) * ((xodij**la(i_aij1, i_aij2))        * &
        &       (das1ijn( i_aij1, i_aij2, lam1, zeff(lam1), dzeffn(lam1) )  + &
        &       dbijn( i_aij1, i_aij2, lam1, xodij ))                       - &
        &       (xodij**lam2)                                               * &
        &       (das1ijn( i_aij1, i_aij2, lam2, zeff(lam2), dzeffn(lam2) )  + &
        &       dbijn( i_aij1, i_aij2, lam2, xodij) )        )  

        return
    end function daijn
!***************************************************************************************************
    function daijv( i_aij1, i_aij2 ) result(aij_out)  !A11
        implicit none
        
        integer,intent(in)  ::  i_aij1, i_aij2
        real(kind=DP)       ::  xodij, lam1, lam2      
        real(kind=DP)       ::  aij_out
            
        lam1 = la(i_aij1, i_aij2)
        lam2 = lr(i_aij1, i_aij2)   
        xodij = sig(i_aij1, i_aij2) / dij(i_aij1, i_aij2)     !A11 +

        aij_out = cij(i_aij1, i_aij2) * ((xodij**la(i_aij1, i_aij2))        * &
        &       (das1ijv( i_aij1, i_aij2, lam1, zeff(lam1), dzeffv(lam1) )  + &
        &       dbijv( i_aij1, i_aij2, lam1, xodij ))                       - &
        &       (xodij**lam2)                                               * &
        &       (das1ijv( i_aij1, i_aij2, lam2, zeff(lam2), dzeffv(lam2) )  + &
        &       dbijv( i_aij1, i_aij2, lam2, xodij) )        )  

        return
    end function daijv
!***************************************************************************************************
    pure function as1ij( asij1, asij2, lam, zxf ) result(as1ij_out)
        implicit none
        
        integer, intent(in)         ::  asij1,asij2
        real(kind=DP), intent(in)   ::  lam,zxf
        real(kind=DP)               ::  as1ij_out
        
        as1ij_out = -2.0e0_DP * rhos * (PI * eps(asij1, asij2) * (dij(asij1, asij2)**3.0e0_DP)  &
        &   / (lam - 3.0e0_DP)) * (( 1.0e0_DP - zxf / 2.0e0_DP) / (1.0e0_DP - zxf)**3.0e0_DP)

        return
    end function as1ij
!***************************************************************************************************
    function das1ijn( asij1, asij2, lam, zxf, dzxfn ) result(as1ij_out)
        implicit none
        
        integer, intent(in)         ::  asij1,asij2
        real(kind=DP), intent(in)   ::  lam,zxf,dzxfn
        real(kind=DP)               ::  as1ij_out
        real(kind=DP)               ::  l_asu,l_asv,l_asup,l_asvp
        
        l_asu  =  1.0e0_DP - zxf / 2.0e0_DP
        l_asv  =  (1.0e0_DP - zxf)**3.0e0_DP
        l_asup =  -dzxfn / 2.0e0_DP
        l_asvp =  -3.0e0_DP * dzxfn * (1.0e0_DP - zxf)**2.0e0_DP    
        
        as1ij_out = -2.0e0_DP * drhosn * (PI * eps(asij1, asij2) * (dij(asij1, asij2)**3.0e0_DP)    /   &
        &     (lam - 3.0e0_DP)) * (( 1.0e0_DP - zxf / 2.0e0_DP) / (1.0e0_DP - zxf)**3.0e0_DP)       -   &
        &            2.0e0_DP * rhos * (PI * eps(asij1, asij2) * (dij(asij1, asij2)**3.0e0_DP)      /   &
        &     (lam - 3.0e0_DP)) * (l_asv * l_asup - l_asu * l_asvp) / l_asv**2.0e0_DP  

        return
    end function das1ijn
!***************************************************************************************************
    function das1ijv( asij1, asij2, lam, zxf, dzxfv ) result(as1ij_out)
        implicit none
        
        integer, intent(in)         ::  asij1, asij2
        real(kind=DP), intent(in)   ::  lam, zxf, dzxfv
        real(kind=DP)               ::  as1ij_out
        real(kind=DP)               ::  l_asu, l_asv, l_asup, l_asvp
        
        l_asu  =  1.0e0_DP - zxf / 2.0e0_DP
        l_asv  =  (1.0e0_DP - zxf)**3.0e0_DP
        l_asup =  -dzxfv / 2.0e0_DP
        l_asvp =  -3.0e0_DP * dzxfv * (1.0e0_DP - zxf)**2.0e0_DP    
        
        as1ij_out = -2.0e0_DP * drhosv * (PI * eps(asij1, asij2) * (dij(asij1, asij2)**3.0e0_DP)    /   &
        &     (lam - 3.0e0_DP)) * (( 1.0e0_DP - zxf / 2.0e0_DP) / (1.0e0_DP - zxf)**3.0e0_DP)       -   &
        &            2.0e0_DP * rhos * (PI * eps(asij1, asij2) * (dij(asij1, asij2)**3.0e0_DP)      /   &
        &     (lam - 3.0e0_DP)) * (l_asv * l_asup - l_asu * l_asvp) / l_asv**2.0e0_DP  

        return
    end function das1ijv
!***************************************************************************************************
    pure function bij( bij1, bij2, lam, xo ) result(bij_out)  
        implicit none
        
        integer,intent(in)          ::  bij1,bij2 
        real(kind=DP),intent(in)    ::  lam,xo
        real(kind=DP)               ::  ilij,jlij
        real(kind=DP)               ::  bij_out
        
        ilij = -(xo**(3.0e0_DP - lam) - 1.0e0_DP) / (lam - 3.0e0_DP)  
        jlij = -(xo**(4.0e0_DP - lam) * (lam - 3.0e0_DP) - xo**(3.0e0_DP - lam) *   &
        &   (lam - 4.0e0_DP) - 1.0e0_DP) / ((lam - 3.0e0_DP) * (lam - 4.0e0_DP)) 

        bij_out = TWOPI * rhos * dij(bij1, bij2)**3.0e0_DP * eps(bij1, bij2) *      &
        &   ( ((1.0e0_DP - zetax / 2.0e0_DP) / (1.0e0_DP - zetax)**3.0e0_DP) *      &
        &   ilij - (9.0e0_DP * zetax * (1.0e0_DP + zetax) /                         &
        &   (2.0e0_DP * (1.0e0_DP - zetax)**3.0e0_DP)) * jlij  ) 

        return
    end function bij
!***************************************************************************************************
    function dbijn( bij1, bij2, lam, xo ) result(bij_out)
        implicit none
        
        integer,intent(in)          ::  bij1,bij2 
        real(kind=DP),intent(in)    ::  lam,xo
        real(kind=DP)               ::  ilij,jlij
        real(kind=DP)               ::  bij_out
        real(kind=DP)               ::  l_bu,l_bv,l_bup,l_bvp
        
        ilij = -(xo**(3.0e0_DP - lam) - 1.0e0_DP) / (lam - 3.0e0_DP)  
        jlij = -(xo**(4.0e0_DP - lam) * (lam - 3.0e0_DP) - xo**(3.0e0_DP - lam) *   &
        &   (lam - 4.0e0_DP) - 1.0e0_DP) / ((lam - 3.0e0_DP) * (lam - 4.0e0_DP)) 

        l_bu  =  TWOPI * rhos * dij(bij1, bij2)**3.0e0_DP * eps(bij1, bij2)
        l_bup =  TWOPI * drhosn * dij(bij1, bij2)**3.0e0_DP * eps(bij1, bij2)

        l_bv  =  ((1.0e0_DP - zetax / 2.0e0_DP) / (1.0e0_DP - zetax)**3.0e0_DP) *   &
        &   ilij - (9.0e0_DP * zetax * (1.0e0_DP + zetax) /                         &
        &   (2.0e0_DP * (1.0e0_DP - zetax)**3.0e0_DP)) * jlij
        
        l_bvp = ( (-(1.0e0_DP - zetax)**3.0e0_DP * dzetaxn/2.0e0_DP                                     + &
        &       (3.0e0_DP * dzetaxn * (1.0e0_DP - zetax)**2.0e0_DP * (1.0e0_DP - zetax / 2.0e0_DP)))    / &
        &       (1.0e0_DP - zetax)**6.0e0_DP ) * ilij                                                   - &
        &       jlij * ( (2.0e0_DP*(1.0e0_DP-zetax)**3.0e0_DP * (18.0e0_DP*zetax*dzetaxn + 9.0e0_DP     * &
        &       dzetaxn )                                                                               + &
        &       54.0e0_DP*zetax * (1.0e0_DP+zetax) * (1.0e0_DP-zetax)**2.0e0_DP * dzetaxn )             / &
        &       (4.0e0_DP * (1.0e0_DP - zetax)**6.0e0_DP) ) 
   
        bij_out = l_bu*l_bvp + l_bv*l_bup
  
        return
    end function dbijn
!***************************************************************************************************
    function dbijv( bij1, bij2, lam, xo ) result(bij_out)
        implicit none
        
        integer,intent(in)          ::  bij1,bij2 
        real(kind=DP),intent(in)    ::  lam,xo
        real(kind=DP)               ::  ilij,jlij
        real(kind=DP)               ::  bij_out
        real(kind=DP)               ::  l_bu,l_bv,l_bup,l_bvp
        
        ilij = -(xo**(3.0e0_DP - lam) - 1.0e0_DP) / (lam - 3.0e0_DP)  
        jlij = -(xo**(4.0e0_DP - lam) * (lam - 3.0e0_DP) - xo**(3.0e0_DP - lam) *   &
        &   (lam - 4.0e0_DP) - 1.0e0_DP) / ((lam - 3.0e0_DP) * (lam - 4.0e0_DP)) 

        l_bu  =  TWOPI * rhos * dij(bij1, bij2)**3.0e0_DP * eps(bij1, bij2)
        l_bup =  TWOPI * drhosv * dij(bij1, bij2)**3.0e0_DP * eps(bij1, bij2)

        l_bv  =  ((1.0e0_DP - zetax / 2.0e0_DP) / (1.0e0_DP - zetax)**3.0e0_DP) *   &
        &   ilij - (9.0e0_DP * zetax * (1.0e0_DP + zetax) /                         &
        &   (2.0e0_DP * (1.0e0_DP - zetax)**3.0e0_DP)) * jlij
        
        l_bvp = ( (-(1.0e0_DP - zetax)**3.0e0_DP * dzetaxv/2.0e0_DP                                     + &
        &       (3.0e0_DP * dzetaxv * (1.0e0_DP - zetax)**2.0e0_DP * (1.0e0_DP - zetax / 2.0e0_DP)))    / &
        &       (1.0e0_DP - zetax)**6.0e0_DP ) * ilij                                                   - &
        &       jlij * ( (2.0e0_DP*(1.0e0_DP-zetax)**3.0e0_DP * (18.0e0_DP*zetax*dzetaxv + 9.0e0_DP     * &
        &       dzetaxv )                                                                               + &
        &       54.0e0_DP*zetax * (1.0e0_DP+zetax) * (1.0e0_DP-zetax)**2.0e0_DP * dzetaxv )             / &
        &       (4.0e0_DP * (1.0e0_DP - zetax)**6.0e0_DP) ) 
   
        bij_out = l_bu*l_bvp + l_bv*l_bup
  
        return
    end function dbijv
!***************************************************************************************************
    pure function as2ij( as2ij1, as2ij2 ) result(as2ij_out)   
        implicit none
        
        integer,intent(in)      ::  as2ij1, as2ij2
        real(kind=DP)           ::  xodij, chiij
        real(kind=DP)           ::  lam1, lam2, lam3, zx1, zx2, zx3
        real(kind=DP)           ::  as2ij_out

        xodij = sig(as2ij1, as2ij2) / dij(as2ij1, as2ij2)     !A11 +
    
        chiij = fkij(1, as2ij1, as2ij2) * zetabar + fkij(2, as2ij1, as2ij2) * zetabar**5.0e0_DP + &
        &       fkij(3, as2ij1, as2ij2) * zetabar**8.0e0_DP

        lam1  = 2.0e0_DP * la(as2ij1, as2ij2)
        lam2  = la(as2ij1, as2ij2) + lr(as2ij1, as2ij2)
        lam3  = 2.0e0_DP * lr(as2ij1, as2ij2)
        zx1   = zeff(lam1)
        zx2   = zeff(lam2)
        zx3   = zeff(lam3)

        as2ij_out=0.5e0_DP * khs * (1.0e0_DP + chiij) * eps(as2ij1, as2ij2) * cij(as2ij1, as2ij2)**2.0e0_DP     &
        &   * ( xodij**lam1 * (as1ij(as2ij1, as2ij2, lam1, zx1) + bij(as2ij1, as2ij2, lam1, xodij))             &
        &   - 2.0e0_DP * xodij**lam2 * (as1ij(as2ij1, as2ij2, lam2, zx2) +                                      &
        &   bij(as2ij1, as2ij2, lam2, xodij))  + xodij**lam3 * (as1ij(as2ij1, as2ij2, lam3, zx3) +              &
        &   bij(as2ij1, as2ij2, lam3, xodij)))  
 
        return
    end function as2ij
!***************************************************************************************************
    function das2ijn( as2ij1, as2ij2 ) result(as2ij_out)   
        implicit none
        
        integer,intent(in)      ::  as2ij1, as2ij2
        real(kind=DP)           ::  xodij, chiij, dchiijn
        real(kind=DP)           ::  lam1, lam2, lam3, zx1, zx2, zx3, dzx1n, dzx2n, dzx3n
        real(kind=DP)           ::  as2ij_out
        real(kind=DP)           ::  las2_u, las2_v, las2_w, las2_up, las2_vp, las2_wp
  
        xodij = sig(as2ij1, as2ij2) / dij(as2ij1, as2ij2)     !A11 +
    
        chiij = fkij(1, as2ij1, as2ij2) * zetabar + fkij(2, as2ij1, as2ij2) * zetabar**5.0e0_DP + &
        &       fkij(3, as2ij1, as2ij2) * zetabar**8.0e0_DP
    
        dchiijn = dzetabarn * (fkij(1, as2ij1, as2ij2) + 5.0e0_DP * fkij(2, as2ij1, as2ij2) * zetabar**4.0e0_DP + &
        &       8.0e0_DP * fkij(3, as2ij1, as2ij2) * zetabar**7.0e0_DP)

        lam1  = 2.0e0_DP * la(as2ij1, as2ij2)
        lam2  = la(as2ij1, as2ij2) + lr(as2ij1, as2ij2)
        lam3  = 2.0e0_DP * lr(as2ij1, as2ij2)
        zx1   = zeff(lam1)
        zx2   = zeff(lam2)
        zx3   = zeff(lam3)
        
        dzx1n   = dzeffn(lam1)
        dzx2n   = dzeffn(lam2)
        dzx3n   = dzeffn(lam3)

        las2_u = 0.5e0_DP * KHS * eps(as2ij1, as2ij2) * cij(as2ij1, as2ij2)**2.0e0_DP
        las2_v = 1.0e0_DP + chiij
        las2_w = xodij**lam1 * (as1ij(as2ij1, as2ij2, lam1, zx1) + bij(as2ij1, as2ij2, lam1, xodij))            &
        &   - 2.0e0_DP * xodij**lam2 * (as1ij(as2ij1, as2ij2, lam2, zx2) +                                      &
        &   bij(as2ij1, as2ij2, lam2, xodij))  + xodij**lam3 * (as1ij(as2ij1, as2ij2, lam3, zx3) +              &
        &   bij(as2ij1, as2ij2, lam3, xodij))

        las2_up = 0.5e0_DP * dKHSn * eps(as2ij1, as2ij2) * cij(as2ij1, as2ij2)**2.0e0_DP
        las2_vp = dchiijn
        las2_wp = xodij**lam1 * (das1ijn(as2ij1, as2ij2, lam1, zx1, dzx1n) + dbijn(as2ij1, as2ij2, lam1, xodij))   - &
        &   2.0e0_DP * xodij**lam2 * (das1ijn(as2ij1, as2ij2, lam2, zx2, dzx2n)                                    + &
        &   dbijn(as2ij1, as2ij2, lam2, xodij))  + xodij**lam3 * (das1ijn(as2ij1, as2ij2, lam3, zx3, dzx3n)        + &
        &   dbijn(as2ij1, as2ij2, lam3, xodij))

        as2ij_out = las2_up*las2_v*las2_w + las2_u*las2_vp*las2_w + las2_u*las2_v*las2_wp  

        return
    end function das2ijn
!***************************************************************************************************
    function das2ijv( as2ij1, as2ij2 ) result(as2ij_out)   
        implicit none
        
        integer,intent(in)      ::  as2ij1, as2ij2
        real(kind=DP)           ::  xodij, chiij, dchiijv
        real(kind=DP)           ::  lam1, lam2, lam3, zx1, zx2, zx3, dzx1v, dzx2v, dzx3v
        real(kind=DP)           ::  as2ij_out
        real(kind=DP)           ::  las2_u, las2_v, las2_w, las2_up, las2_vp, las2_wp
  
        xodij = sig(as2ij1, as2ij2) / dij(as2ij1, as2ij2)     !A11 +
    
        chiij = fkij(1, as2ij1, as2ij2) * zetabar + fkij(2, as2ij1, as2ij2) * zetabar**5.0e0_DP + &
        &       fkij(3, as2ij1, as2ij2) * zetabar**8.0e0_DP
    
        dchiijv = dzetabarv * (fkij(1, as2ij1, as2ij2) + 5.0e0_DP * fkij(2, as2ij1, as2ij2) * zetabar**4.0e0_DP + &
        &       8.0e0_DP * fkij(3, as2ij1, as2ij2) * zetabar**7.0e0_DP)

        lam1  = 2.0e0_DP * la(as2ij1, as2ij2)
        lam2  = la(as2ij1, as2ij2) + lr(as2ij1, as2ij2)
        lam3  = 2.0e0_DP * lr(as2ij1, as2ij2)
        zx1   = zeff(lam1)
        zx2   = zeff(lam2)
        zx3   = zeff(lam3)
        
        dzx1v   = dzeffv(lam1)
        dzx2v   = dzeffv(lam2)
        dzx3v   = dzeffv(lam3)

        las2_u = 0.5e0_DP * KHS * eps(as2ij1, as2ij2) * cij(as2ij1, as2ij2)**2.0e0_DP
        las2_v = 1.0e0_DP + chiij
        las2_w = xodij**lam1 * (as1ij(as2ij1, as2ij2, lam1, zx1) + bij(as2ij1, as2ij2, lam1, xodij))            &
        &   - 2.0e0_DP * xodij**lam2 * (as1ij(as2ij1, as2ij2, lam2, zx2) +                                      &
        &   bij(as2ij1, as2ij2, lam2, xodij))  + xodij**lam3 * (as1ij(as2ij1, as2ij2, lam3, zx3) +              &
        &   bij(as2ij1, as2ij2, lam3, xodij))

        las2_up = 0.5e0_DP * dKHSv * eps(as2ij1, as2ij2) * cij(as2ij1, as2ij2)**2.0e0_DP
        las2_vp = dchiijv
        las2_wp = xodij**lam1 * (das1ijv(as2ij1, as2ij2, lam1, zx1, dzx1v) + dbijv(as2ij1, as2ij2, lam1, xodij))   - &
        &   2.0e0_DP * xodij**lam2 * (das1ijv(as2ij1, as2ij2, lam2, zx2, dzx2v)                                    + &
        &   dbijv(as2ij1, as2ij2, lam2, xodij))  + xodij**lam3 * (das1ijv(as2ij1, as2ij2, lam3, zx3, dzx3v)        + &
        &   dbijv(as2ij1, as2ij2, lam3, xodij))

        as2ij_out = las2_up*las2_v*las2_w + las2_u*las2_vp*las2_w + las2_u*las2_v*las2_wp  

        return
    end function das2ijv
!***************************************************************************************************
    function as3ij( as3ij1, as3ij2 ) result(as3ij_out)
        implicit none
        
        integer,intent(in)      ::  as3ij1, as3ij2
        real(kind=DP)           ::  as3ij_out
        
        as3ij_out = -eps(as3ij1, as3ij2)**3.0e0_DP * fkij(4, as3ij1, as3ij2) * zetabar * &
            & dexp(fkij( 5, as3ij1, as3ij2 ) * zetabar + fkij(6, as3ij1, as3ij2)* zetabar**2.0e0_DP)

        return
    end function as3ij
!***************************************************************************************************
    function das3ijn( as3ij1, as3ij2 ) result(as3ij_out)
        implicit none
        
        integer,intent(in)      ::  as3ij1, as3ij2
        real(kind=DP)           ::  as3ij_out
        real(kind=DP)           ::  las3_u, las3_v, las3_up, las3_vp
      
        las3_u  = -eps(as3ij1, as3ij2)**3.0e0_DP * fkij(4, as3ij1, as3ij2) * zetabar
        las3_up = -eps(as3ij1, as3ij2)**3.0e0_DP * fkij(4, as3ij1, as3ij2) * dzetabarn
        las3_v  = dexp(fkij( 5, as3ij1, as3ij2 ) * zetabar + fkij(6, as3ij1, as3ij2)* zetabar**2.0e0_DP)
        las3_vp =  dzetabarn * (fkij( 5, as3ij1, as3ij2 ) + 2.0e0_DP * fkij(6, as3ij1, as3ij2)* zetabar) * las3_v
    
        as3ij_out = las3_u * las3_vp + las3_v * las3_up

        return
    end function das3ijn
!***************************************************************************************************
    function das3ijv( as3ij1, as3ij2 ) result(as3ij_out)
        implicit none
        
        integer,intent(in)      ::  as3ij1, as3ij2
        real(kind=DP)           ::  as3ij_out
        real(kind=DP)           ::  las3_u, las3_v, las3_up, las3_vp
      
        las3_u  = -eps(as3ij1, as3ij2)**3.0e0_DP * fkij(4, as3ij1, as3ij2) * zetabar
        las3_up = -eps(as3ij1, as3ij2)**3.0e0_DP * fkij(4, as3ij1, as3ij2) * dzetabarv
        las3_v  = dexp(fkij( 5, as3ij1, as3ij2 ) * zetabar + fkij(6, as3ij1, as3ij2)* zetabar**2.0e0_DP)
        las3_vp =  dzetabarv * (fkij( 5, as3ij1, as3ij2 ) + 2.0e0_DP * fkij(6, as3ij1, as3ij2)* zetabar) * las3_v
    
        as3ij_out = las3_u * las3_vp + las3_v * las3_up

        return
    end function das3ijv
!***************************************************************************************************

!***************************************************************************************************
end module Mono_mod
!***************************************************************************************************
!***************************************************************************************************
!   SAFT Module to:
!       1. Calculate A chain
!***************************************************************************************************
!
!***************************************************************************************************
module Chain_mod
!***************************************************************************************************
!Modules
!=======    
    use Types_mod           ! Definitions of types and double precision
    use Global_mod          ! Important global parameters  
    use Zeff_mod            ! Zeff calculations
!***************************************************************************************************
    implicit none
!***************************************************************************************************
    public  ::  A_chain, A_chain_dn, A_chain_dv
!***************************************************************************************************
    contains
!*************************************************************************************************** 
!***************************************************************************************************
!***************************************************************************************************
!***************************************************************************************************
!Note to self, if any Chain error remains, it is in lngii, but not g1 or g2...
    function A_chain() result(a_res)     ![1]  A27     [2] 18      
        implicit none
        
        integer             ::  i_a1, i_a2
        real(kind=DP)       ::  tsum
        real(kind=DP)       ::  a_res  

        a_res = 0.0e0_DP
      
        if(chain_switch) then             
            do i_a1=1, nctypes
                tsum = 0.0e0_DP
            do i_a2=1, nstypes
                if(Comp_array(i_a1)%comp(i_a2) > 0.0e0_DP) THEN                        
                    tsum=tsum+Comp_array(i_a1)%comp(i_a2) * Seg_array(i_a2)%nseg * Seg_array(i_a2)%sf 
                end if
            end do
                a_res = a_res - Comp_array(i_a1)%xi * (tsum - 1.0e0_DP) * lngii(i_a1)                               
            end do
        end if

        return
    end function A_chain    
!***************************************************************************************************   
    function A_chain_dn( i_mu_in ) result(a_res)     ![1]  A27     [2] 18      
        implicit none
        
        integer, intent(in) ::  i_mu_in 
        integer             ::  i_a1, i_a2
        real(kind=DP)       ::  tsum
        real(kind=DP)       ::  a_res  

        a_res = 0.0e0_DP
      
        if(chain_switch) then             
            do i_a1=1, nctypes
                if(Comp_array(i_a1)%xi/=0.0e0_DP) then
                    tsum = 0.0e0_DP
                do i_a2=1, nstypes
                    if(Comp_array(i_a1)%comp(i_a2) > 0.0e0_DP) THEN                        
                        tsum=tsum+Comp_array(i_a1)%comp(i_a2) * Seg_array(i_a2)%nseg * Seg_array(i_a2)%sf 
                    end if
                end do
                    a_res = a_res -  &
                    &       nsum * KB * t * Comp_array(i_a1)%dxin * (tsum - 1.0e0_DP) * lngii(i_a1)     - &
                    &       nsum * KB * t * Comp_array(i_a1)%xi * (tsum - 1.0e0_DP) * dlngiin(i_a1)     - &
                    &       KB * t * Comp_array(i_a1)%xi * (tsum - 1.0e0_DP) * lngii(i_a1)    
                end if
            end do        
        end if

        return
    end function A_chain_dn   
!***************************************************************************************************  
    function A_chain_dv() result(a_res)     ![1]  A27     [2] 18      
        implicit none
        
        integer             ::  i_a1, i_a2
        real(kind=DP)       ::  tsum
        real(kind=DP)       ::  a_res  

        a_res = 0.0e0_DP
      
        if(chain_switch) then             
            do i_a1=1, nctypes
                if(Comp_array(i_a1)%xi/=0.0e0_DP) then
                    tsum = 0.0e0_DP
                do i_a2=1, nstypes
                    if(Comp_array(i_a1)%comp(i_a2) > 0.0e0_DP) THEN                        
                        tsum=tsum+Comp_array(i_a1)%comp(i_a2) * Seg_array(i_a2)%nseg * Seg_array(i_a2)%sf 
                    end if
                end do
                    a_res = a_res - nsum * KB * t * Comp_array(i_a1)%xi * (tsum - 1.0e0_DP) * dlngiiv(i_a1)      
                end if    
            end do
        end if

        return
    end function A_chain_dv   
!*************************************************************************************************** 
    function lngii( gri ) result(lngii_out)
        implicit none
        
        integer,intent(in)      ::  gri
        integer                 ::  i_lngii1, i_lngii2
        
        real(kind=DP)           ::  gii, k0, k1, k2, k3
        real(kind=DP)           ::  l_ghs, l_x0, l_d, l_C, l_la, l_lr
        real(kind=DP)           ::  l_eps, l_a1a, l_a1r, l_Ba, l_Br
        real(kind=DP)           ::  l_G1, l_G2
        real(kind=DP)           ::  g2mca, l_gam, l_alpha, l_theta, l_zetaxs
        real(kind=DP)           ::  l_a12r, l_B2r, l_a1ar, l_Bar, l_a12a, l_B2a
        real(kind=DP)           ::  lngii_out
         
        l_d = d3ij(gri)**(1.0e0_DP / 3.0e0_DP)
        l_x0 = sig3ij(gri,gri)**(1.0e0_DP / 3.0e0_DP)/ l_d

        k0 = -dlog(1.0e0_DP - zetax)   +                                                      &
        &   (42.0e0_DP * zetax - 39.0e0_DP * zetax**2.0e0_DP + 9.0e0_DP * zetax**3.0e0_DP  - &
        &   2.0e0_DP * zetax**4.0e0_DP) / (6.0e0_DP * (1.0e0_DP - zetax)**3.0e0_DP)
 
        k1 = (zetax**4.0e0_DP + 6.0e0_DP * zetax**2.0e0_DP - 12.0e0_DP * zetax) / (2.0e0_DP &
        &   * (1.0e0_DP - zetax)**3.0e0_DP)
        k2 = (-3.0e0_DP * zetax**2.0e0_DP) / (8.0e0_DP * (1.0e0_DP - zetax)**2.0e0_DP)
        k3 = (-zetax**4.0e0_DP + 3.0e0_DP * zetax**2.0e0_DP + 3.0e0_DP * zetax) / (6.0e0_DP &
        &   * (1.0e0_DP - zetax)**3.0e0_DP)

        l_ghs = dexp(k0 + k1 * l_x0 + k2 * l_x0**2.0e0_DP + k3 * l_x0**3.0e0_DP)

        l_la = lachij(gri)
        l_lr = lrchij(gri)
    
        l_C = l_lr / (l_lr - l_la) * (l_lr / l_la) ** (l_la / (l_lr - l_la))       
        l_eps = epschij(gri, gri)

        l_Ba  = bchain(l_la, l_x0, l_d, l_eps)
        l_Br  = bchain(l_lr, l_x0, l_d, l_eps)
        l_B2r = bchain(2.0e0_DP * l_lr, l_x0, l_d, l_eps) 
        l_Bar = bchain(l_la + l_lr, l_x0, l_d, l_eps) 
        l_B2a = bchain(2.0e0_DP * l_la, l_x0, l_d, l_eps) 

        l_a1a  = as1chain(l_la, l_d, l_eps)        
        l_a1r  = as1chain(l_lr, l_d, l_eps)
        l_a12r = as1chain(2.0e0_DP * l_lr, l_d, l_eps)
     
        l_a1ar = as1chain(l_la + l_lr, l_d, l_eps)
        l_a12a = as1chain(2.0e0_DP * l_la, l_d, l_eps)

        l_G1 = 1.0e0_DP / (2.0e0_DP * PI * l_eps * l_d**3.0e0_DP) * (3.0e0_DP *     &
        &   as1grad(l_C, l_x0 ,l_d ,l_la ,l_lr ,l_eps) -                            &
        &   l_C * l_la * l_x0**l_la * (l_a1a + l_Ba) / rhos   +                     &
        &   l_C * l_lr * l_x0**l_lr * (l_a1r + l_Br) / rhos )

        l_alpha = l_C * (1.0e0_DP / (l_la - 3.0e0_DP) - 1.0e0_DP / (l_lr - 3.0e0_DP))
        l_theta = dexp(beta_K * l_eps) - 1.0e0_DP

        l_zetaxs = 0.0e0_DP
      
        do i_lngii1=1, nstypes
        do i_lngii2=1, nstypes
            l_zetaxs = l_zetaxs + Seg_array(i_lngii1)%xsi * Seg_array(i_lngii2)%xsi &
            &   * sig(i_lngii1, i_lngii2)**3.0e0_DP
        end do
        end do
        
        l_zetaxs = l_zetaxs * PI * rhos / 6.0e0_DP * NA

        l_gam = phik(7, 0) * (-1.0e0_DP * dtanh(phik(7, 1) * (phik(7, 2) - l_alpha)) + 1.0e0_DP) * &
        &   l_zetaxs * l_theta * dexp(phik(7, 3) * l_zetaxs + phik(7, 4) * l_zetaxs**2.0e0_DP)

        g2mca = 1.0e0_DP / (2.0e0_DP * PI * l_eps**2.00e0_DP * l_d**3.0e0_DP) * ( 3.0e0_DP *            &
        &   as2grad(rhos, l_C, l_x0, l_d, l_la,l_lr,l_eps) -                                            &
        &   l_eps * khs * l_C**2.0e0_DP * l_lr * l_x0**(2.0e0_DP * l_lr) * (l_a12r + l_B2r ) / rhos +   &
        &   l_eps * khs * l_C**2.0e0_DP * (l_la + l_lr) * l_x0**(l_la + l_lr)                       *   &    
        &   (l_a1ar + l_Bar) / rhos - l_eps * khs * l_C**2.0e0_DP * l_la * l_x0**(2.0e0_DP * l_la)  *   &
        &   (l_a12a + l_B2a) / rhos)

        l_G2 = (1.0e0_DP + l_gam) * g2mca

        gii = l_ghs * dexp(beta_K * l_eps * l_G1 / l_ghs + (beta_K * l_eps)**2.0e0_DP * l_G2 / l_ghs)     
        lngii_out = dlog(gii)
!print*,l_ghs,l_G1,l_G2,gii, lngii_out 
        return
    end function lngii
!***************************************************************************************************  
    function dlngiin( gri ) result(lngii_out)
        implicit none
        
        integer,intent(in)      ::  gri
        integer                 ::  i_lngii1, i_lngii2
        
        real(kind=DP)           ::  gii, k0, k1, k2, k3
        real(kind=DP)           ::  l_ghs, l_x0, l_d, l_C, l_la, l_lr
        real(kind=DP)           ::  l_eps, l_a1a, l_a1r, l_Ba, l_Br
        real(kind=DP)           ::  l_G1, l_G2
        real(kind=DP)           ::  g2mca, l_gam, l_alpha, l_theta, l_zetaxs
        real(kind=DP)           ::  l_a12r, l_B2r, l_a1ar, l_Bar, l_a12a, l_B2a
        real(kind=DP)           ::  lngii_out
        !differential parts
        real(kind=DP)           ::  dk0n, dk1n, dk2n, dk3n, dl_ghsn
        real(kind=DP)           ::  dl_a1an, dl_a1rn, dl_Ban, dl_Brn
        real(kind=DP)           ::  dl_a12rn, dl_B2rn, dl_a1arn, dl_Barn, dl_a12an, dl_B2an
        real(kind=DP)           ::  dl_G1n, dl_gamn, dl_G2n, dgiin, dl_zetaxsn, dg2mcan
             
        l_d = d3ij(gri)**(1.0e0_DP / 3.0e0_DP)

        l_x0 = sig3ij(gri,gri)**(1.0e0_DP / 3.0e0_DP)/ l_d

        k0 = -dlog(1.0e0_DP - zetax)   +                                                      &
        &   (42.0e0_DP * zetax - 39.0e0_DP * zetax**2.0e0_DP + 9.0e0_DP * zetax**3.0e0_DP  - &
        &   2.0e0_DP * zetax**4.0e0_DP) / (6.0e0_DP * (1.0e0_DP - zetax)**3.0e0_DP)
        
        k1 = (zetax**4.0e0_DP + 6.0e0_DP * zetax**2.0e0_DP - 12.0e0_DP * zetax) / (2.0e0_DP &
        &   * (1.0e0_DP - zetax)**3.0e0_DP)
        
        k2 = (-3.0e0_DP * zetax**2.0e0_DP) / (8.0e0_DP * (1.0e0_DP - zetax)**2.0e0_DP)
        
        k3 = (-zetax**4.0e0_DP + 3.0e0_DP * zetax**2.0e0_DP + 3.0e0_DP * zetax) / (6.0e0_DP &
        &   * (1.0e0_DP - zetax)**3.0e0_DP)
        
        dk0n = dzetaxn/(1.0e0_DP - zetax)                                                           + &
        &   (6.0e0_DP * (1.0e0_DP - zetax)**3.0e0_DP * (42.0e0_DP * dzetaxn - 78.0e0_DP * zetax     * &
        &   dzetaxn + 27.0e0_DP * zetax**2.0e0_DP * dzetaxn - 8.0e0_DP * zetax**3.0e0_DP * dzetaxn) + &
        &   18.0e0_DP * dzetaxn * (1.0e0_DP - zetax)**2.0e0_DP                                      * &
        &   (42.0e0_DP * zetax - 39.0e0_DP * zetax**2.0e0_DP + 9.0e0_DP * zetax**3.0e0_DP           - &
        &   2.0e0_DP * zetax**4.0e0_DP))                                                            / &
        &   (36.0e0_DP * (1.0e0_DP - zetax)**6.0e0_DP)  

        dk1n = ( 2.0e0_DP  * (1.0e0_DP - zetax)**3.0e0_DP                                           * &
        &   (4.0e0_DP*dzetaxn*zetax**3.0e0_DP + 12.0e0_DP * dzetaxn*zetax - 12.0e0_DP * dzetaxn)    + &
        &   6.0e0_DP * dzetaxn * (1.0e0_DP - zetax)**2.0e0_DP                                       * &
        &   (zetax**4.0e0_DP + 6.0e0_DP * zetax**2.0e0_DP - 12.0e0_DP * zetax))                     / &
        &   (4.0e0_DP  * (1.0e0_DP - zetax)**6.0e0_DP)

        dk2n = (-6.0e0_DP * dzetaxn * zetax * 8.0e0_DP * (1.0e0_DP - zetax)**2.0e0_DP               - &
        &   3.0e0_DP * zetax**2.0e0_DP * 16.0e0_DP * dzetaxn * (1.0e0_DP - zetax))                  / &
        &   (64.0e0_DP * (1.0e0_DP - zetax)**4.0e0_DP)    

        dk3n = ((6.0e0_DP * (1.0e0_DP - zetax)**3.0e0_DP)                                           * &
        &   (-4.0e0_DP*dzetaxn*zetax**3.0e0_DP + 6.0e0_DP * dzetaxn * zetax + 3.0e0_DP * dzetaxn)   - &
        &   (-zetax**4.0e0_DP + 3.0e0_DP * zetax**2.0e0_DP + 3.0e0_DP * zetax)                      * &
        &   (-18.0e0_DP * dzetaxn * (1.0e0_DP - zetax)**2.0e0_DP))                                  / &
        &   (36.0e0_DP * (1.0e0_DP - zetax)**6.0e0_DP)
       
        l_ghs = dexp(k0 + k1 * l_x0 + k2 * l_x0**2.0e0_DP + k3 * l_x0**3.0e0_DP)

        dl_ghsn = (dk0n + dk1n * l_x0 + dk2n * l_x0**2.0e0_DP + dk3n * l_x0**3.0e0_DP)              * &
        &   dexp(k0 + k1 * l_x0 + k2 * l_x0**2.0e0_DP + k3 * l_x0**3.0e0_DP)

        l_la = lachij(gri)
        l_lr = lrchij(gri)
       
        l_C = l_lr / (l_lr - l_la) * (l_lr / l_la) ** (l_la / (l_lr - l_la))       
        l_eps = epschij(gri, gri)

        l_Ba  = bchain(l_la, l_x0, l_d, l_eps)
        l_Br  = bchain(l_lr, l_x0, l_d, l_eps)
        l_B2r = bchain(2.0e0_DP * l_lr, l_x0, l_d, l_eps) 
        l_Bar = bchain(l_la + l_lr, l_x0, l_d, l_eps) 
        l_B2a = bchain(2.0e0_DP * l_la, l_x0, l_d, l_eps)

        dl_Ban  = dbchainn(l_la, l_x0, l_d, l_eps)
        dl_Brn  = dbchainn(l_lr, l_x0, l_d, l_eps)
        dl_B2rn = dbchainn(2.0e0_DP * l_lr, l_x0, l_d, l_eps) 
        dl_Barn = dbchainn(l_la + l_lr, l_x0, l_d, l_eps) 
        dl_B2an = dbchainn(2.0e0_DP * l_la, l_x0, l_d, l_eps) 

        l_a1a  = as1chain(l_la, l_d, l_eps)        
        l_a1r  = as1chain(l_lr, l_d, l_eps)
        l_a12r = as1chain(2.0e0_DP * l_lr, l_d, l_eps)
        l_a1ar = as1chain(l_la + l_lr, l_d, l_eps)
        l_a12a = as1chain(2.0e0_DP * l_la, l_d, l_eps)
        
        dl_a1an  = das1chainn(l_la, l_d, l_eps)        
        dl_a1rn  = das1chainn(l_lr, l_d, l_eps)
        dl_a12rn = das1chainn(2.0e0_DP * l_lr, l_d, l_eps)
        dl_a1arn = das1chainn(l_la + l_lr, l_d, l_eps)
        dl_a12an = das1chainn(2.0e0_DP * l_la, l_d, l_eps)

        l_G1 = 1.0e0_DP / (2.0e0_DP * PI * l_eps * l_d**3.0e0_DP) * (3.0e0_DP *     &
        &   as1grad(l_C, l_x0 ,l_d ,l_la ,l_lr ,l_eps) -                            &
        &   l_C * l_la * l_x0**l_la * (l_a1a + l_Ba) / rhos   +                     &
        &   l_C * l_lr * l_x0**l_lr * (l_a1r + l_Br) / rhos )

        dl_G1n = 1.0e0_DP / (2.0e0_DP * PI * l_eps * l_d**3.0e0_DP) * (3.0e0_DP     * &
        &   das1gradn(l_C, l_x0 ,l_d ,l_la ,l_lr ,l_eps)                            - &
        &   l_C * l_la * l_x0**l_la                                                 * & 
        &   (rhos*(dl_a1an + dl_Ban) - drhosn*(l_a1a + l_Ba))/rhos**2.0e0_DP        + &
        &   l_C * l_lr * l_x0**l_lr                                                 * &
        &   (rhos*(dl_a1rn + dl_Brn) - drhosn*(l_a1r + l_Br))/rhos**2.0e0_DP) 

        l_alpha = l_C * (1.0e0_DP / (l_la - 3.0e0_DP) - 1.0e0_DP / (l_lr - 3.0e0_DP))
        l_theta = dexp(beta_K * l_eps) - 1.0e0_DP
        
        l_zetaxs    = 0.0e0_DP
        dl_zetaxsn  = 0.0e0_DP

        do i_lngii1=1, nstypes
        do i_lngii2=1, nstypes
            l_zetaxs = l_zetaxs + Seg_array(i_lngii1)%xsi * Seg_array(i_lngii2)%xsi &
            &   * sig(i_lngii1, i_lngii2)**3.0e0_DP
            
            dl_zetaxsn = dl_zetaxsn + (Seg_array(i_lngii1)%dxsin * Seg_array(i_lngii2)%xsi      + &
            &   Seg_array(i_lngii1)%xsi * Seg_array(i_lngii2)%dxsin)                            * &
            &   sig(i_lngii1, i_lngii2)**3.0e0_DP
        end do
        end do
               
        dl_zetaxsn = PI / 6.0e0_DP * NA * (l_zetaxs * drhosn + dl_zetaxsn * rhos)
        l_zetaxs   = l_zetaxs * PI * rhos / 6.0e0_DP * NA
    
        l_gam = phik(7, 0) * (-1.0e0_DP * dtanh(phik(7, 1) * (phik(7, 2) - l_alpha)) + 1.0e0_DP) * &
        &   l_zetaxs * l_theta * dexp(phik(7, 3) * l_zetaxs + phik(7, 4) * l_zetaxs**2.0e0_DP)
        
        dl_gamn = phik(7, 0) * (-1.0e0_DP * dtanh(phik(7, 1) * (phik(7, 2) - l_alpha)) + 1.0e0_DP)   * &
        &   l_theta * dexp(phik(7, 3) * l_zetaxs + phik(7, 4) * l_zetaxs**2.0e0_DP)                  * &
        &   (dl_zetaxsn + l_zetaxs * (phik(7, 3) * dl_zetaxsn + 2.0e0_DP * phik(7, 4) * dl_zetaxsn  * &
        &   l_zetaxs))

        g2mca = 1.0e0_DP / (2.0e0_DP * PI * l_eps**2.00e0_DP * l_d**3.0e0_DP) * ( 3.0e0_DP *            &
        &   as2grad(rhos, l_C, l_x0, l_d, l_la,l_lr,l_eps) -                                            &
        &   l_eps * khs * l_C**2.0e0_DP * l_lr * l_x0**(2.0e0_DP * l_lr) * (l_a12r + l_B2r ) / rhos +   &
        &   l_eps * khs * l_C**2.0e0_DP * (l_la + l_lr) * l_x0**(l_la + l_lr)                       *   &    
        &   (l_a1ar + l_Bar) / rhos - l_eps * khs * l_C**2.0e0_DP * l_la * l_x0**(2.0e0_DP * l_la)  *   &
        &   (l_a12a + l_B2a) / rhos)
        
        dg2mcan = 1.0e0_DP / (2.0e0_DP * PI * l_eps**2.00e0_DP * l_d**3.0e0_DP)                         * &
        &   ( 3.0e0_DP * das2gradn(rhos, l_C, l_x0, l_d, l_la,l_lr,l_eps)                               - &

        &   l_eps * l_C**2.0e0_DP * l_lr * l_x0**(2.0e0_DP * l_lr)                              * ( &
        &   khs * (rhos * (dl_a12rn + dl_B2rn ) - drhosn * (l_a12r + l_B2r )) / rhos**2.0e0_DP  +   &
        &   dkhsn * (l_a12r + l_B2r ) / rhos)                                                   +   &  

        &   l_eps * l_C**2.0e0_DP * (l_la + l_lr) * l_x0**(l_la + l_lr)                         * ( &
        &   khs * (rhos * (dl_a1arn + dl_Barn ) - drhosn * (l_a1ar + l_Bar )) / rhos**2.0e0_DP  +   &
        &   dkhsn * (l_a1ar + l_Bar ) / rhos)                                                   -   &  
        
        &   l_eps * l_C**2.0e0_DP * l_la * l_x0**(2.0e0_DP * l_la)                              * ( &
        &   khs * (rhos * (dl_a12an + dl_B2an ) - drhosn * (l_a12a + l_B2a )) / rhos**2.0e0_DP  +   &
        &   dkhsn * (l_a12a + l_B2a ) / rhos)) 
      
        l_G2 = (1.0e0_DP + l_gam) * g2mca
        
        dl_G2n = (1.0e0_DP + l_gam) * dg2mcan + dl_gamn * g2mca
     
        gii = l_ghs * dexp(beta_K * l_eps * l_G1 / l_ghs + (beta_K * l_eps)**2.0e0_DP * l_G2 / l_ghs)     
        
        dgiin = dexp(beta_K * l_eps * l_G1 / l_ghs + (beta_K * l_eps)**2.0e0_DP * l_G2 / l_ghs)      * &
        &   (dl_ghsn + l_ghs                                                                        * &
        &   (beta_K * l_eps * (dl_G1n*l_ghs - l_G1*dl_ghsn)/l_ghs**2.0e0_DP                         + &
        &   (beta_K * l_eps)**2.0e0_DP * (dl_G2n*l_ghs - l_G2*dl_ghsn)/l_ghs**2.0e0_DP))      
   
        lngii_out = dgiin / gii

        return
    end function dlngiin
!*************************************************************************************************** 
    function dlngiiv( gri ) result(lngii_out)
        implicit none
        
        integer,intent(in)      ::  gri
        integer                 ::  i_lngii1, i_lngii2
        
        real(kind=DP)           ::  gii, k0, k1, k2, k3
        real(kind=DP)           ::  l_ghs, l_x0, l_d, l_C, l_la, l_lr
        real(kind=DP)           ::  l_eps, l_a1a, l_a1r, l_Ba, l_Br
        real(kind=DP)           ::  l_G1, l_G2
        real(kind=DP)           ::  g2mca, l_gam, l_alpha, l_theta, l_zetaxs
        real(kind=DP)           ::  l_a12r, l_B2r, l_a1ar, l_Bar, l_a12a, l_B2a
        real(kind=DP)           ::  lngii_out
        !differential parts
        real(kind=DP)           ::  dk0v, dk1v, dk2v, dk3v, dl_ghsv
        real(kind=DP)           ::  dl_a1av, dl_a1rv, dl_Bav, dl_Brv
        real(kind=DP)           ::  dl_a12rv, dl_B2rv, dl_a1arv, dl_Barv, dl_a12av, dl_B2av
        real(kind=DP)           ::  dl_G1v, dl_gamv, dl_G2v, dgiiv, dl_zetaxsv, dg2mcav
           
        l_d = d3ij(gri)**(1.0e0_DP / 3.0e0_DP)
        l_x0 = sig3ij(gri,gri)**(1.0e0_DP / 3.0e0_DP)/ l_d

        k0 = -dlog(1.0e0_DP - zetax)   +                                                      &
        &   (42.0e0_DP * zetax - 39.0e0_DP * zetax**2.0e0_DP + 9.0e0_DP * zetax**3.0e0_DP  - &
        &   2.0e0_DP * zetax**4.0e0_DP) / (6.0e0_DP * (1.0e0_DP - zetax)**3.0e0_DP)
        
        k1 = (zetax**4.0e0_DP + 6.0e0_DP * zetax**2.0e0_DP - 12.0e0_DP * zetax) / (2.0e0_DP &
        &   * (1.0e0_DP - zetax)**3.0e0_DP)
        
        k2 = (-3.0e0_DP * zetax**2.0e0_DP) / (8.0e0_DP * (1.0e0_DP - zetax)**2.0e0_DP)
        
        k3 = (-zetax**4.0e0_DP + 3.0e0_DP * zetax**2.0e0_DP + 3.0e0_DP * zetax) / (6.0e0_DP &
        &   * (1.0e0_DP - zetax)**3.0e0_DP)
        
        dk0v = dzetaxv/(1.0e0_DP - zetax)                                                           + &
        &   (6.0e0_DP * (1.0e0_DP - zetax)**3.0e0_DP * (42.0e0_DP * dzetaxv - 78.0e0_DP * zetax     * &
        &   dzetaxv + 27.0e0_DP * zetax**2.0e0_DP * dzetaxv - 8.0e0_DP * zetax**3.0e0_DP * dzetaxv) + &
        &   18.0e0_DP * dzetaxv * (1.0e0_DP - zetax)**2.0e0_DP                                      * &
        &   (42.0e0_DP * zetax - 39.0e0_DP * zetax**2.0e0_DP + 9.0e0_DP * zetax**3.0e0_DP           - &
        &   2.0e0_DP * zetax**4.0e0_DP))                                                            / &
        &   (36.0e0_DP * (1.0e0_DP - zetax)**6.0e0_DP)  

        dk1v = ( 2.0e0_DP  * (1.0e0_DP - zetax)**3.0e0_DP                                           * &
        &   (4.0e0_DP*dzetaxv*zetax**3.0e0_DP + 12.0e0_DP * dzetaxv*zetax - 12.0e0_DP * dzetaxv)    + &
        &   6.0e0_DP * dzetaxv * (1.0e0_DP - zetax)**2.0e0_DP                                       * &
        &   (zetax**4.0e0_DP + 6.0e0_DP * zetax**2.0e0_DP - 12.0e0_DP * zetax))                     / &
        &   (4.0e0_DP  * (1.0e0_DP - zetax)**6.0e0_DP)

        dk2v = (-6.0e0_DP * dzetaxv * zetax * 8.0e0_DP * (1.0e0_DP - zetax)**2.0e0_DP               - &
        &   3.0e0_DP * zetax**2.0e0_DP * 16.0e0_DP * dzetaxv * (1.0e0_DP - zetax))                  / &
        &   (64.0e0_DP * (1.0e0_DP - zetax)**4.0e0_DP)    

        dk3v = ((6.0e0_DP * (1.0e0_DP - zetax)**3.0e0_DP)                                           * &
        &   (-4.0e0_DP*dzetaxv*zetax**3.0e0_DP + 6.0e0_DP * dzetaxv * zetax + 3.0e0_DP * dzetaxv)   - &
        &   (-zetax**4.0e0_DP + 3.0e0_DP * zetax**2.0e0_DP + 3.0e0_DP * zetax)                      * &
        &   (-18.0e0_DP * dzetaxv * (1.0e0_DP - zetax)**2.0e0_DP))                                  / &
        &   (36.0e0_DP * (1.0e0_DP - zetax)**6.0e0_DP)
   
        l_ghs = dexp(k0 + k1 * l_x0 + k2 * l_x0**2.0e0_DP + k3 * l_x0**3.0e0_DP)

        dl_ghsv = (dk0v + dk1v * l_x0 + dk2v * l_x0**2.0e0_DP + dk3v * l_x0**3.0e0_DP)              * &
        &   dexp(k0 + k1 * l_x0 + k2 * l_x0**2.0e0_DP + k3 * l_x0**3.0e0_DP)

        l_la = lachij(gri)
        l_lr = lrchij(gri)
       
        l_C = l_lr / (l_lr - l_la) * (l_lr / l_la) ** (l_la / (l_lr - l_la))       
        l_eps = epschij(gri, gri)

        l_Ba  = bchain(l_la, l_x0, l_d, l_eps)
        l_Br  = bchain(l_lr, l_x0, l_d, l_eps)
        l_B2r = bchain(2.0e0_DP * l_lr, l_x0, l_d, l_eps) 
        l_Bar = bchain(l_la + l_lr, l_x0, l_d, l_eps) 
        l_B2a = bchain(2.0e0_DP * l_la, l_x0, l_d, l_eps)

        dl_Bav  = dbchainv(l_la, l_x0, l_d, l_eps)
        dl_Brv  = dbchainv(l_lr, l_x0, l_d, l_eps)
        dl_B2rv = dbchainv(2.0e0_DP * l_lr, l_x0, l_d, l_eps) 
        dl_Barv = dbchainv(l_la + l_lr, l_x0, l_d, l_eps) 
        dl_B2av = dbchainv(2.0e0_DP * l_la, l_x0, l_d, l_eps) 

        l_a1a  = as1chain(l_la, l_d, l_eps)        
        l_a1r  = as1chain(l_lr, l_d, l_eps)
        l_a12r = as1chain(2.0e0_DP * l_lr, l_d, l_eps)
        l_a1ar = as1chain(l_la + l_lr, l_d, l_eps)
        l_a12a = as1chain(2.0e0_DP * l_la, l_d, l_eps)
        
        dl_a1av  = das1chainv(l_la, l_d, l_eps)        
        dl_a1rv  = das1chainv(l_lr, l_d, l_eps)
        dl_a12rv = das1chainv(2.0e0_DP * l_lr, l_d, l_eps)
        dl_a1arv = das1chainv(l_la + l_lr, l_d, l_eps)
        dl_a12av = das1chainv(2.0e0_DP * l_la, l_d, l_eps)

        l_G1 = 1.0e0_DP / (2.0e0_DP * PI * l_eps * l_d**3.0e0_DP) * (3.0e0_DP *     &
        &   as1grad(l_C, l_x0 ,l_d ,l_la ,l_lr ,l_eps) -                            &
        &   l_C * l_la * l_x0**l_la * (l_a1a + l_Ba) / rhos   +                     &
        &   l_C * l_lr * l_x0**l_lr * (l_a1r + l_Br) / rhos )

        dl_G1v = 1.0e0_DP / (2.0e0_DP * PI * l_eps * l_d**3.0e0_DP) * (3.0e0_DP     * &
        &   das1gradv(l_C, l_x0 ,l_d ,l_la ,l_lr ,l_eps)                            - &
        &   l_C * l_la * l_x0**l_la                                                 * & 
        &   (rhos*(dl_a1av + dl_Bav) - drhosv*(l_a1a + l_Ba))/rhos**2.0e0_DP        + &
        &   l_C * l_lr * l_x0**l_lr                                                 * &
        &   (rhos*(dl_a1rv + dl_Brv) - drhosv*(l_a1r + l_Br))/rhos**2.0e0_DP) 

        l_alpha = l_C * (1.0e0_DP / (l_la - 3.0e0_DP) - 1.0e0_DP / (l_lr - 3.0e0_DP))
        l_theta = dexp(beta_K * l_eps) - 1.0e0_DP
        
        l_zetaxs    = 0.0e0_DP
        dl_zetaxsv  = 0.0e0_DP

        do i_lngii1=1, nstypes
        do i_lngii2=1, nstypes
            l_zetaxs = l_zetaxs + Seg_array(i_lngii1)%xsi * Seg_array(i_lngii2)%xsi &
            &   * sig(i_lngii1, i_lngii2)**3.0e0_DP
        end do
        end do
        
        l_zetaxs   = l_zetaxs * PI * rhos / 6.0e0_DP * NA
        dl_zetaxsv = l_zetaxs *  drhosv / rhos

        l_gam = phik(7, 0) * (-1.0e0_DP * dtanh(phik(7, 1) * (phik(7, 2) - l_alpha)) + 1.0e0_DP) * &
        &   l_zetaxs * l_theta * dexp(phik(7, 3) * l_zetaxs + phik(7, 4) * l_zetaxs**2.0e0_DP)
        
        dl_gamv = phik(7, 0) * (-1.0e0_DP * dtanh(phik(7, 1) * (phik(7, 2) - l_alpha)) + 1.0e0_DP)   * &
        &   l_theta * dexp(phik(7, 3) * l_zetaxs + phik(7, 4) * l_zetaxs**2.0e0_DP)                  * &
        &   (dl_zetaxsv + l_zetaxs * (phik(7, 3) * dl_zetaxsv + 2.0e0_DP * phik(7, 4) * dl_zetaxsv  * &
        &   l_zetaxs))

        g2mca = 1.0e0_DP / (2.0e0_DP * PI * l_eps**2.00e0_DP * l_d**3.0e0_DP) * ( 3.0e0_DP *            &
        &   as2grad(rhos, l_C, l_x0, l_d, l_la,l_lr,l_eps) -                                            &
        &   l_eps * khs * l_C**2.0e0_DP * l_lr * l_x0**(2.0e0_DP * l_lr) * (l_a12r + l_B2r ) / rhos +   &
        &   l_eps * khs * l_C**2.0e0_DP * (l_la + l_lr) * l_x0**(l_la + l_lr)                       *   &    
        &   (l_a1ar + l_Bar) / rhos - l_eps * khs * l_C**2.0e0_DP * l_la * l_x0**(2.0e0_DP * l_la)  *   &
        &   (l_a12a + l_B2a) / rhos)
        
        dg2mcav = 1.0e0_DP / (2.0e0_DP * PI * l_eps**2.00e0_DP * l_d**3.0e0_DP)                         * &
        &   ( 3.0e0_DP * das2gradv(rhos, l_C, l_x0, l_d, l_la,l_lr,l_eps)                               - &

        &   l_eps * l_C**2.0e0_DP * l_lr * l_x0**(2.0e0_DP * l_lr)                              * ( &
        &   khs * (rhos * (dl_a12rv + dl_B2rv ) - drhosv * (l_a12r + l_B2r )) / rhos**2.0e0_DP  +   &
        &   dkhsv * (l_a12r + l_B2r ) / rhos)                                                   +   &  

        &   l_eps * l_C**2.0e0_DP * (l_la + l_lr) * l_x0**(l_la + l_lr)                         * ( &
        &   khs * (rhos * (dl_a1arv + dl_Barv ) - drhosv * (l_a1ar + l_Bar )) / rhos**2.0e0_DP  +   &
        &   dkhsv * (l_a1ar + l_Bar ) / rhos)                                                   -   &  
        
        &   l_eps * l_C**2.0e0_DP * l_la * l_x0**(2.0e0_DP * l_la)                              * ( &
        &   khs * (rhos * (dl_a12av + dl_B2av ) - drhosv * (l_a12a + l_B2a )) / rhos**2.0e0_DP  +   &
        &   dkhsv * (l_a12a + l_B2a ) / rhos)) 

        l_G2 = (1.0e0_DP + l_gam) * g2mca
        
        dl_G2v = (1.0e0_DP + l_gam) * dg2mcav + dl_gamv * g2mca

        gii = l_ghs * dexp(beta_K * l_eps * l_G1 / l_ghs + (beta_K * l_eps)**2.0e0_DP * l_G2 / l_ghs)     
        
        dgiiv = dexp(beta_K * l_eps * l_G1 / l_ghs + (beta_K * l_eps)**2.0e0_DP * l_G2 / l_ghs)      * &
        &   (dl_ghsv + l_ghs                                                                        * &
        &   (beta_K * l_eps * (dl_G1v*l_ghs - l_G1*dl_ghsv)/l_ghs**2.0e0_DP                         + &
        &   (beta_K * l_eps)**2.0e0_DP * (dl_G2v*l_ghs - l_G2*dl_ghsv)/l_ghs**2.0e0_DP))           
     
        lngii_out = dgiiv / gii

        return
    end function dlngiiv
!*************************************************************************************************** 
    function bchain( lam, xo, l_dd, l_e ) result(bchain_out)
        implicit none
        
        real(kind=DP), intent(in)   ::  lam, xo, l_dd, l_e
        real(kind=DP)               ::  ilij, jlij
        real(kind=DP)               ::  bchain_out
        
        ilij = -(xo**(3.0e0_DP - lam) - 1.0e0_DP) / (lam - 3.0e0_DP)  
        jlij = -(xo**(4.0e0_DP - lam) * (lam - 3.0e0_DP) - xo**(3.0e0_DP - lam) *  &
        &       (lam - 4.0e0_DP) - 1.0e0_DP) / ((lam - 3.0e0_DP) * (lam - 4.0e0_DP))

        bchain_out = TWOPI * rhos * l_dd**3.0e0_DP * l_e *                          &
        &   (   ((1.0e0_DP - zetax / 2.0e0_DP) / (1.0e0_DP - zetax)**3.0e0_DP) *    &
        &   ilij - (9.0e0_DP * zetax * (1.0e0_DP + zetax) / (2.0e0_DP*(1.0e0_DP -   &
        &   zetax)**3.0e0_DP)) * jlij  ) 

        return
    end function bchain 
!*************************************************************************************************** 
    function dbchainn( lam, xo, l_dd, l_e ) result(bchain_out)
        implicit none
        
        real(kind=DP), intent(in)   ::  lam, xo, l_dd, l_e
        real(kind=DP)               ::  ilij, jlij
        real(kind=DP)               ::  bchain_out
        real(kind=DP)               ::  lb_1,lb_1p,lb_2,lb_2p,lb_3,lb_3p
        
        ilij = -(xo**(3.0e0_DP - lam) - 1.0e0_DP) / (lam - 3.0e0_DP)  
        jlij = -(xo**(4.0e0_DP - lam) * (lam - 3.0e0_DP) - xo**(3.0e0_DP - lam) *  &
        &       (lam - 4.0e0_DP) - 1.0e0_DP) / ((lam - 3.0e0_DP) * (lam - 4.0e0_DP))
        
        lb_1  = TWOPI * rhos * l_dd**3.0e0_DP * l_e
        lb_1p = TWOPI * drhosn * l_dd**3.0e0_DP * l_e
        
        lb_2  = ((1.0e0_DP - zetax / 2.0e0_DP) / (1.0e0_DP - zetax)**3.0e0_DP) * ilij
        lb_2p = ilij * (-1.0e0_DP*dzetaxn / 2.0e0_DP * (1.0e0_DP - zetax)**3.0e0_DP + 3.0e0_DP*dzetaxn   * &
        &   (1.0e0_DP - zetax)**2.0e0_DP *(1.0e0_DP - zetax / 2.0e0_DP))                                / &
        &   (1.0e0_DP - zetax)**6.0e0_DP
        
        lb_3  = (9.0e0_DP * zetax * (1.0e0_DP + zetax) / (2.0e0_DP*(1.0e0_DP - zetax)**3.0e0_DP)) * jlij
        lb_3p = jlij * ((2.0e0_DP*(1.0e0_DP - zetax)**3.0e0_DP) * 9.0e0_DP * (dzetaxn * (1.0e0_DP + zetax) + &
        &   zetax * dzetaxn) + 54.0e0_DP*dzetaxn*(1.0e0_DP - zetax)**2.0e0_DP * zetax * (1.0e0_DP + zetax)) / &
        &   (4.0e0_DP*(1.0e0_DP - zetax)**6.0e0_DP)
        
        bchain_out = lb_1p * (lb_2 - lb_3) + lb_1 * (lb_2p - lb_3p) 

        return
    end function dbchainn
!***************************************************************************************************
    function dbchainv( lam, xo, l_dd, l_e ) result(bchain_out)
        implicit none
        
        real(kind=DP), intent(in)   ::  lam, xo, l_dd, l_e
        real(kind=DP)               ::  ilij, jlij
        real(kind=DP)               ::  bchain_out
        real(kind=DP)               ::  lb_1,lb_1p,lb_2,lb_2p,lb_3,lb_3p
        
        ilij = -(xo**(3.0e0_DP - lam) - 1.0e0_DP) / (lam - 3.0e0_DP)  
        jlij = -(xo**(4.0e0_DP - lam) * (lam - 3.0e0_DP) - xo**(3.0e0_DP - lam) *  &
        &       (lam - 4.0e0_DP) - 1.0e0_DP) / ((lam - 3.0e0_DP) * (lam - 4.0e0_DP))
        
        lb_1  = TWOPI * rhos * l_dd**3.0e0_DP * l_e
        lb_1p = TWOPI * drhosv * l_dd**3.0e0_DP * l_e
        
        lb_2  = ((1.0e0_DP - zetax / 2.0e0_DP) / (1.0e0_DP - zetax)**3.0e0_DP) * ilij
        lb_2p = ilij * (-dzetaxv / 2.0e0_DP * (1.0e0_DP - zetax)**3.0e0_DP + 3.0e0_DP*dzetaxv   * &
        &   (1.0e0_DP - zetax)**2.0e0_DP *(1.0e0_DP - zetax / 2.0e0_DP))                                / &
        &   (1.0e0_DP - zetax)**6.0e0_DP
        
        lb_3  = (9.0e0_DP * zetax * (1.0e0_DP + zetax) / (2.0e0_DP*(1.0e0_DP - zetax)**3.0e0_DP)) * jlij
        lb_3p = jlij * ((2.0e0_DP*(1.0e0_DP - zetax)**3.0e0_DP) * 9.0e0_DP * (dzetaxv * (1.0e0_DP + zetax) + &
        &   zetax * dzetaxv) + 54.0e0_DP*dzetaxv*(1.0e0_DP - zetax)**2.0e0_DP * zetax * (1.0e0_DP + zetax)) / &
        &   (4.0e0_DP*(1.0e0_DP - zetax)**6.0e0_DP)
        
        bchain_out = lb_1p * (lb_2 - lb_3) + lb_1 * (lb_2p - lb_3p) 

        return
    end function dbchainv
!***************************************************************************************************
    function as1chain( lam, l_dd, l_e ) result(as1chain_out)
        implicit none

        real(kind=DP), intent(in)    ::  lam, l_dd, l_e
        real(kind=DP)                ::  zxf
        real(kind=DP)                ::  as1chain_out
        
        zxf = zeff(lam)
      
        as1chain_out = -2.0e0_DP * rhos * (PI * l_e) * (l_dd**3.0e0_DP) / (lam - 3.0e0_DP)   &
        &  * ((1.0e0_DP - zxf/2.0e0_DP) / (1.0e0_DP - zxf)**3.0e0_DP)

        return
    end function as1chain
!***************************************************************************************************
    pure function das1chainn( lam, l_dd, l_e ) result(as1chain_out)
        implicit none

        real(kind=DP), intent(in)    ::  lam, l_dd, l_e
        real(kind=DP)                ::  zxf, dzxfn
        real(kind=DP)                ::  as1chain_out
        real(kind=DP)                ::  las1_1,las1_1p,las1_2,las1_2p
        
        zxf = zeff(lam)
        dzxfn = dzeffn(lam)
        
        las1_1  = -2.0e0_DP * rhos * (PI * l_e) * (l_dd**3.0e0_DP) / (lam - 3.0e0_DP)
        las1_1p = -2.0e0_DP * drhosn * (PI * l_e) * (l_dd**3.0e0_DP) / (lam - 3.0e0_DP)
        las1_2  = (1.0e0_DP - zxf/2.0e0_DP) / (1.0e0_DP - zxf)**3.0e0_DP
        las1_2p =  (-dzxfn/2.0e0_DP * (1.0e0_DP - zxf)**3.0e0_DP + 3.0e0_DP * dzxfn * (1.0e0_DP - zxf)**2.0e0_DP * &
        &   (1.0e0_DP - zxf/2.0e0_DP)) / (1.0e0_DP - zxf)**6.0e0_DP
        
        as1chain_out = las1_1 * las1_2p + las1_1p * las1_2
        
        return
    end function das1chainn
!***************************************************************************************************
    function das1chainv( lam, l_dd, l_e ) result(as1chain_out)
        implicit none

        real(kind=DP), intent(in)    ::  lam, l_dd, l_e
        real(kind=DP)                ::  zxf, dzxfv
        real(kind=DP)                ::  as1chain_out
        real(kind=DP)                ::  las1_1,las1_1p,las1_2,las1_2p
        
        zxf = zeff(lam)
        dzxfv = dzeffv(lam)
      
        las1_1  = -2.0e0_DP * rhos * (PI * l_e) * (l_dd**3.0e0_DP) / (lam - 3.0e0_DP)
        las1_1p = -2.0e0_DP * drhosv * (PI * l_e) * (l_dd**3.0e0_DP) / (lam - 3.0e0_DP)
        las1_2  = (1.0e0_DP - zxf/2.0e0_DP) / (1.0e0_DP - zxf)**3.0e0_DP
        las1_2p =  (-dzxfv/2.0e0_DP * (1.0e0_DP - zxf)**3.0e0_DP + 3.0e0_DP * dzxfv * (1.0e0_DP - zxf)**2.0e0_DP * &
        &   (1.0e0_DP - zxf/2.0e0_DP)) / (1.0e0_DP - zxf)**6.0e0_DP
        
        as1chain_out = las1_1 * las1_2p + las1_1p * las1_2
        
        return
    end function das1chainv
!***************************************************************************************************
    function as1grad( cc, xx, dd, lla, llr, leps) result(as1grad_out) 
        implicit none
        
        real(kind=DP),intent(in)    ::  cc, xx, dd, lla, llr, leps
        real(kind=DP)               ::  as1grad_out
        real(kind=DP)               ::  l_da1, l_dr1, l_const1a, l_const1r, l_za, l_zr, l_zap, l_zrp
        real(kind=DP)               ::  l_dBa, l_dBr, l_Ja, l_Jr, l_Ia, l_Ir, l_const2

        l_const1a = PI * leps*dd**3.0e0_DP / (lla - 3.0e0_DP)
        l_const1r = PI * leps*dd**3.0e0_DP / (llr - 3.0e0_DP)
       
        l_za  = zeff(lla)
        l_zr  = zeff(llr)
        l_zap = zeffprime(lla)
        l_zrp = zeffprime(llr)
               
        l_da1 = -2.0e0_DP * l_const1a * ( (1.0e0_DP - l_za/2.0e0_DP) / (1.0e0_DP - l_za)**3.0e0_DP + rhos *  &
        &   (2.5e0_DP * l_zap - l_zap * l_za) / (1.0e0_DP - l_za)**4.0e0_DP)

        l_dr1 = -2.0e0_DP * l_const1r * ( (1.0e0_DP - l_zr/2.0e0_DP) / (1.0e0_DP - l_zr)**3.0e0_DP + rhos *  &
        &   (2.5e0_DP * l_zrp - l_zrp * l_zr) / (1.0e0_DP - l_zr)**4.0e0_DP)
        
        l_const2 = 2.0e0_DP * PI * dd**3.0e0_DP * leps     
        
        l_Ia = -(xx**(3.0e0_DP - lla) - 1.0e0_DP) / (lla - 3.0e0_DP)  
        l_Ir = -(xx**(3.0e0_DP - llr) - 1.0e0_DP) / (llr - 3.0e0_DP)
   
        l_Ja = -(xx**(4.0e0_DP - lla) * (lla - 3.0e0_DP) - xx**(3.0e0_DP - lla) * (lla - 4.0e0_DP) - 1.0e0_DP) / &
        &   ((lla - 3.0e0_DP) * (lla - 4.0e0_DP))
        l_Jr = -(xx**(4.0e0_DP - llr) * (llr - 3.0e0_DP) - xx**(3.0e0_DP - llr) * (llr - 4.0e0_DP) - 1.0e0_DP) / &
        &   ((llr - 3.0e0_DP) * (llr - 4.0e0_DP))   
        
        l_dBa = l_const2 * (l_Ia * (((1.0e0_DP - zetax / 2.0e0_DP) / (1.0e0_DP - zetax)**3.0e0_DP )                 +   &
                &   zetax * (2.5e0_DP - zetax) / (1.0e0_DP - zetax)**4.0e0_DP)                                      -   &   
                &   l_Ja * (zetax/( 2.0d0 * ( 1.0e0_DP - zetax )**3.0e0_DP ) * (( 9.0e0_DP * (1.0e0_DP + zetax))    +   &
                &   (( 9.0e0_DP  * ( 1.0e0_DP + 2.0e0_DP * zetax ) * 2.0e0_DP * ( 1.0e0_DP - zetax )  )             +   &
                &   ( 9.e0_DP * ( zetax + zetax**2.0e0_DP ) * 6.0e0_DP ) ) / ( 2.0e0_DP * (1.0e0_DP - zetax) ))  ))
        
        l_dBr = l_const2 * (l_Ir * (((1.0e0_DP - zetax / 2.0e0_DP) / (1.0e0_DP - zetax)**3.0e0_DP )                 +   &
                &   zetax * (2.5e0_DP - zetax) / (1.0e0_DP - zetax)**4.0e0_DP)                                      -   &   
                &   l_Jr * (zetax/( 2.0d0 * ( 1.0e0_DP - zetax )**3.0e0_DP ) * (( 9.0e0_DP * (1.0e0_DP + zetax))    +   &
                &   (( 9.0e0_DP  * ( 1.0e0_DP + 2.0e0_DP * zetax ) * 2.0e0_DP * ( 1.0e0_DP - zetax )  )             +   &
                &   ( 9.e0_DP * ( zetax + zetax**2.0e0_DP ) * 6.0e0_DP ) ) / ( 2.0e0_DP * (1.0e0_DP - zetax) ))  ))
    
        as1grad_out = cc * ( xx**lla * (l_da1 + l_dBa) - xx**llr * (l_dr1 + l_dBr) )

        return
    end function as1grad
!***************************************************************************************************
    function das1gradn( cc, xx, dd, lla, llr, leps) result(as1grad_out) 
        implicit none
        
        real(kind=DP),intent(in)    ::  cc, xx, dd, lla, llr, leps
        real(kind=DP)               ::  as1grad_out
        real(kind=DP)               ::  l_da1, l_dr1, l_const1a, l_const1r, l_za, l_zr, l_zap, l_zrp
        real(kind=DP)               ::  l_dBa, l_dBr, l_Ja, l_Jr, l_Ia, l_Ir, l_const2
        !differentials
        real(kind=DP)               ::  dl_zan, dl_zrn, dl_zapn, dl_zrpn
        real(kind=DP)               ::  dl_da1n, dl_dr1n, dl_dBan, dl_dBrn
        real(kind=DP)               ::  temp_1,temp_2,temp_3,temp_4,temp_5,temp_6,temp_7
        real(kind=DP)               ::  temp_1p,temp_2p,temp_3p,temp_4p,temp_5p,temp_6p,temp_7p

        l_const1a = PI * leps*dd**3.0e0_DP / (lla - 3.0e0_DP)
        l_const1r = PI * leps*dd**3.0e0_DP / (llr - 3.0e0_DP)
       
        l_za  = zeff(lla)
        l_zr  = zeff(llr)
        l_zap = zeffprime(lla)
        l_zrp = zeffprime(llr)

        dl_zan  = dzeffn(lla)
        dl_zrn  = dzeffn(llr)
        dl_zapn = dzeffprimen(lla)
        dl_zrpn = dzeffprimen(llr)
             
        l_da1 = -2.0e0_DP * l_const1a * ( (1.0e0_DP - l_za/2.0e0_DP) / (1.0e0_DP - l_za)**3.0e0_DP + rhos *  &
        &   (2.5e0_DP * l_zap - l_zap * l_za) / (1.0e0_DP - l_za)**4.0e0_DP)

        l_dr1 = -2.0e0_DP * l_const1r * ( (1.0e0_DP - l_zr/2.0e0_DP) / (1.0e0_DP - l_zr)**3.0e0_DP + rhos *  &
        &   (2.5e0_DP * l_zrp - l_zrp * l_zr) / (1.0e0_DP - l_zr)**4.0e0_DP)

        dl_da1n = -2.0e0_DP * l_const1a                                                                     * &
        &   ( ((1.0e0_DP - l_za)**3.0e0_DP*(-dl_zan/2.0e0_DP)                                               + &
        &   3.0e0_DP * dl_zan * (1.0e0_DP - l_za)**2.0e0_DP * (1.0e0_DP - l_za/2.0e0_DP) )                  / &
        &   (1.0e0_DP - l_za)**6.0e0_DP                                                                     + &
        &   drhosn * (2.5e0_DP * l_zap - l_zap * l_za) / (1.0e0_DP - l_za)**4.0e0_DP                        + &
        &   rhos * ((1.0e0_DP - l_za)**4.0e0_DP*(2.5e0_DP * dl_zapn - (dl_zapn*l_za + l_zap*dl_zan))        + &
        &   4.0e0_DP * dl_zan                                                                              * &
        &   (1.0e0_DP - l_za)**3.0e0_DP * (2.5e0_DP * l_zap - l_zap * l_za)) / (1.0e0_DP - l_za)**8.0e0_DP)

        dl_dr1n = -2.0e0_DP * l_const1r                                                                     * &
        &   ( ((1.0e0_DP - l_zr)**3.0e0_DP*(-dl_zrn/2.0e0_DP)                                               + &
        &   3.0e0_DP * dl_zrn * (1.0e0_DP - l_zr)**2.0e0_DP * (1.0e0_DP - l_zr/2.0e0_DP) )                  / &
        &   (1.0e0_DP - l_zr)**6.0e0_DP                                                                     + &
        &   drhosn * (2.5e0_DP * l_zrp - l_zrp * l_zr) / (1.0e0_DP - l_zr)**4.0e0_DP                        + &
        &   rhos * ((1.0e0_DP - l_zr)**4.0e0_DP*(2.5e0_DP * dl_zrpn - (dl_zrpn*l_zr + l_zrp*dl_zrn))        + &
        &   4.0e0_DP * dl_zrn                                                                              * &
        &   (1.0e0_DP - l_zr)**3.0e0_DP * (2.5e0_DP * l_zrp - l_zrp * l_zr)) / (1.0e0_DP - l_zr)**8.0e0_DP)

        l_const2 = 2.0e0_DP * PI * dd**3.0e0_DP * leps     
        
        l_Ia = -(xx**(3.0e0_DP - lla) - 1.0e0_DP) / (lla - 3.0e0_DP)  
        l_Ir = -(xx**(3.0e0_DP - llr) - 1.0e0_DP) / (llr - 3.0e0_DP)
   
        l_Ja = -(xx**(4.0e0_DP - lla) * (lla - 3.0e0_DP) - xx**(3.0e0_DP - lla) * (lla - 4.0e0_DP) - 1.0e0_DP) / &
        &   ((lla - 3.0e0_DP) * (lla - 4.0e0_DP))
        l_Jr = -(xx**(4.0e0_DP - llr) * (llr - 3.0e0_DP) - xx**(3.0e0_DP - llr) * (llr - 4.0e0_DP) - 1.0e0_DP) / &
        &   ((llr - 3.0e0_DP) * (llr - 4.0e0_DP))   

        l_dBa = l_const2 * (l_Ia * (((1.0e0_DP - zetax / 2.0e0_DP) / (1.0e0_DP - zetax)**3.0e0_DP )                 +   &
                &   zetax * (2.5e0_DP - zetax) / (1.0e0_DP - zetax)**4.0e0_DP)                                      -   &   
                &   l_Ja * (zetax/( 2.0d0 * ( 1.0e0_DP - zetax )**3.0e0_DP ) * (( 9.0e0_DP * (1.0e0_DP + zetax))    +   &
                &   (( 9.0e0_DP  * ( 1.0e0_DP + 2.0e0_DP * zetax ) * 2.0e0_DP * ( 1.0e0_DP - zetax )  )             +   &
                &   ( 9.e0_DP * ( zetax + zetax**2.0e0_DP ) * 6.0e0_DP ) ) / ( 2.0e0_DP * (1.0e0_DP - zetax) ))  ))
        
        l_dBr = l_const2 * (l_Ir * (((1.0e0_DP - zetax / 2.0e0_DP) / (1.0e0_DP - zetax)**3.0e0_DP )                 +   &
                &   zetax * (2.5e0_DP - zetax) / (1.0e0_DP - zetax)**4.0e0_DP)                                      -   &   
                &   l_Jr * (zetax/( 2.0d0 * ( 1.0e0_DP - zetax )**3.0e0_DP ) * (( 9.0e0_DP * (1.0e0_DP + zetax))    +   &
                &   (( 9.0e0_DP  * ( 1.0e0_DP + 2.0e0_DP * zetax ) * 2.0e0_DP * ( 1.0e0_DP - zetax )  )             +   &
                &   ( 9.e0_DP * ( zetax + zetax**2.0e0_DP ) * 6.0e0_DP ) ) / ( 2.0e0_DP * (1.0e0_DP - zetax) ))  ))
        
        temp_1  = (1.0e0_DP - zetax / 2.0e0_DP) / (1.0e0_DP - zetax)**3.0e0_DP   
        temp_1p = ((-dzetaxn / 2.0e0_DP)*(1.0e0_DP - zetax)**3.0e0_DP + 3.0e0_DP*dzetaxn*(1.0e0_DP              - &
        &   zetax)**2.0e0_DP)*(1.0e0_DP - zetax / 2.0e0_DP) / (1.0e0_DP - zetax)**6.0e0_DP  
        
        temp_2  = zetax * (2.5e0_DP - zetax) / (1.0e0_DP - zetax)**4.0e0_DP
        temp_2p = dzetaxn * (2.5e0_DP - zetax) / (1.0e0_DP - zetax)**4.0e0_DP                                   + &
        &   zetax * (-dzetaxn*(1.0e0_DP - zetax)**4.0e0_DP + 4.0e0_DP*dzetaxn*(1.0e0_DP - zetax)**3.0e0_DP      * &
        &   (2.5e0_DP - zetax)) /  (1.0e0_DP - zetax)**8.0e0_DP
        
        temp_3  = zetax/( 2.0d0 * ( 1.0e0_DP - zetax )**3.0e0_DP )
        temp_3p = (dzetaxn * ( 2.0d0 * ( 1.0e0_DP - zetax )**3.0e0_DP ) + 6.0d0 * dzetaxn * zetax               * &
        &   ( 1.0e0_DP - zetax )**2.0e0_DP) / ( 4.0d0 * ( 1.0e0_DP - zetax )**6.0e0_DP )

        temp_4  = 9.0e0_DP * (1.0e0_DP + zetax) 
        temp_4p = 9.0e0_DP * dzetaxn
 
        temp_5  = 18.0e0_DP * ( 1.0e0_DP + 2.0e0_DP * zetax ) * ( 1.0e0_DP - zetax )
        temp_5p = 18.0e0_DP * (( 2.0e0_DP * dzetaxn ) * ( 1.0e0_DP - zetax )                                    - &
        &    dzetaxn * ( 1.0e0_DP + 2.0e0_DP * zetax ))
        
        temp_6  = 54.0e0_DP * ( zetax + zetax**2.0e0_DP ) 
        temp_6p = 54.0e0_DP * ( dzetaxn + 2.0e0_DP * dzetaxn * zetax ) 
        
        temp_7  = 2.0e0_DP * (1.0e0_DP - zetax)
        temp_7p = -2.0e0_DP * dzetaxn
        
        dl_dBan = l_const2 * (l_Ia * (temp_1p + temp_2p)                                                        - &
        &   l_Ja * (temp_3p*(temp_4 +(temp_5+temp_6)/temp_7)                                                    + &
        &   temp_3*(temp_4p + ((temp_5p+temp_6p)*temp_7 - temp_7p*(temp_5+temp_6))/temp_7**2.0e0_DP)))
        
        dl_dBrn = l_const2 * (l_Ir * (temp_1p + temp_2p)                                                        - &
        &   l_Jr * (temp_3p*(temp_4 +(temp_5+temp_6)/temp_7)                                                    + &
        &   temp_3*(temp_4p + ((temp_5p+temp_6p)*temp_7 - temp_7p*(temp_5+temp_6))/temp_7**2.0e0_DP)))         
      
        as1grad_out = cc * ( xx**lla * (dl_da1n + dl_dBan) - xx**llr * (dl_dr1n + dl_dBrn) )

        return
    end function das1gradn
!***************************************************************************************************
    function das1gradv( cc, xx, dd, lla, llr, leps) result(as1grad_out) 
        implicit none
        
        real(kind=DP),intent(in)    ::  cc, xx, dd, lla, llr, leps
        real(kind=DP)               ::  as1grad_out
        real(kind=DP)               ::  l_da1, l_dr1, l_const1a, l_const1r, l_za, l_zr, l_zap, l_zrp
        real(kind=DP)               ::  l_dBa, l_dBr, l_Ja, l_Jr, l_Ia, l_Ir, l_const2
        !differentials
        real(kind=DP)               ::  dl_zav, dl_zrv, dl_zapv, dl_zrpv
        real(kind=DP)               ::  dl_da1v, dl_dr1v, dl_dBav, dl_dBrv
        real(kind=DP)               ::  temp_1,temp_2,temp_3,temp_4,temp_5,temp_6,temp_7
        real(kind=DP)               ::  temp_1p,temp_2p,temp_3p,temp_4p,temp_5p,temp_6p,temp_7p

        l_const1a = PI * leps*dd**3.0e0_DP / (lla - 3.0e0_DP)
        l_const1r = PI * leps*dd**3.0e0_DP / (llr - 3.0e0_DP)
       
        l_za  = zeff(lla)
        l_zr  = zeff(llr)
        l_zap = zeffprime(lla)
        l_zrp = zeffprime(llr)

        dl_zav  = dzeffv(lla)
        dl_zrv  = dzeffv(llr)
        dl_zapv = dzeffprimev(lla)
        dl_zrpv = dzeffprimev(llr)
     
        l_da1 = -2.0e0_DP * l_const1a * ( (1.0e0_DP - l_za/2.0e0_DP) / (1.0e0_DP - l_za)**3.0e0_DP + rhos *  &
        &   (2.5e0_DP * l_zap - l_zap * l_za) / (1.0e0_DP - l_za)**4.0e0_DP)

        l_dr1 = -2.0e0_DP * l_const1r * ( (1.0e0_DP - l_zr/2.0e0_DP) / (1.0e0_DP - l_zr)**3.0e0_DP + rhos *  &
        &   (2.5e0_DP * l_zrp - l_zrp * l_zr) / (1.0e0_DP - l_zr)**4.0e0_DP)

        dl_da1v = -2.0e0_DP * l_const1a                                                                     * &
        &   ( ((1.0e0_DP - l_za)**3.0e0_DP*(-dl_zav/2.0e0_DP)                                               + &
        &   3.0e0_DP * dl_zav * (1.0e0_DP - l_za)**2.0e0_DP * (1.0e0_DP - l_za/2.0e0_DP) )                  / &
        &   (1.0e0_DP - l_za)**6.0e0_DP                                                                     + &
        &   drhosv * (2.5e0_DP * l_zap - l_zap * l_za) / (1.0e0_DP - l_za)**4.0e0_DP                        + &
        &   rhos * ((1.0e0_DP - l_za)**4.0e0_DP*(2.5e0_DP * dl_zapv - (dl_zapv*l_za + l_zap*dl_zav))        + &
        &   4.0e0_DP * dl_zav  * &!ALTERED THESE dl_zav from dl_zapv                                                                            
        &   (1.0e0_DP - l_za)**3.0e0_DP * (2.5e0_DP * l_zap - l_zap * l_za)) / (1.0e0_DP - l_za)**8.0e0_DP)

        dl_dr1v = -2.0e0_DP * l_const1r                                                                     * &
        &   ( ((1.0e0_DP - l_zr)**3.0e0_DP*(-dl_zrv/2.0e0_DP)                                               + &
        &   3.0e0_DP * dl_zrv * (1.0e0_DP - l_zr)**2.0e0_DP * (1.0e0_DP - l_zr/2.0e0_DP) )                  / &
        &   (1.0e0_DP - l_zr)**6.0e0_DP                                                                     + &
        &   drhosv * (2.5e0_DP * l_zrp - l_zrp * l_zr) / (1.0e0_DP - l_zr)**4.0e0_DP                        + &
        &   rhos * ((1.0e0_DP - l_zr)**4.0e0_DP*(2.5e0_DP * dl_zrpv - (dl_zrpv*l_zr + l_zrp*dl_zrv))        + &
        &   4.0e0_DP * dl_zrv  * &!ALTERED THESE dl_zrv from dl_zrpv  
        &   (1.0e0_DP - l_zr)**3.0e0_DP * (2.5e0_DP * l_zrp - l_zrp * l_zr)) / (1.0e0_DP - l_zr)**8.0e0_DP)

        l_const2 = 2.0e0_DP * PI * dd**3.0e0_DP * leps     
        
        l_Ia = -(xx**(3.0e0_DP - lla) - 1.0e0_DP) / (lla - 3.0e0_DP)  
        l_Ir = -(xx**(3.0e0_DP - llr) - 1.0e0_DP) / (llr - 3.0e0_DP)
   
        l_Ja = -(xx**(4.0e0_DP - lla) * (lla - 3.0e0_DP) - xx**(3.0e0_DP - lla) * (lla - 4.0e0_DP) - 1.0e0_DP) / &
        &   ((lla - 3.0e0_DP) * (lla - 4.0e0_DP))
        l_Jr = -(xx**(4.0e0_DP - llr) * (llr - 3.0e0_DP) - xx**(3.0e0_DP - llr) * (llr - 4.0e0_DP) - 1.0e0_DP) / &
        &   ((llr - 3.0e0_DP) * (llr - 4.0e0_DP))   

        l_dBa = l_const2 * (l_Ia * (((1.0e0_DP - zetax / 2.0e0_DP) / (1.0e0_DP - zetax)**3.0e0_DP )                 +   &
                &   zetax * (2.5e0_DP - zetax) / (1.0e0_DP - zetax)**4.0e0_DP)                                      -   &   
                &   l_Ja * (zetax/( 2.0d0 * ( 1.0e0_DP - zetax )**3.0e0_DP ) * (( 9.0e0_DP * (1.0e0_DP + zetax))    +   &
                &   (( 9.0e0_DP  * ( 1.0e0_DP + 2.0e0_DP * zetax ) * 2.0e0_DP * ( 1.0e0_DP - zetax )  )             +   &
                &   ( 9.e0_DP * ( zetax + zetax**2.0e0_DP ) * 6.0e0_DP ) ) / ( 2.0e0_DP * (1.0e0_DP - zetax) ))  ))
        
        l_dBr = l_const2 * (l_Ir * (((1.0e0_DP - zetax / 2.0e0_DP) / (1.0e0_DP - zetax)**3.0e0_DP )                 +   &
                &   zetax * (2.5e0_DP - zetax) / (1.0e0_DP - zetax)**4.0e0_DP)                                      -   &   
                &   l_Jr * (zetax/( 2.0d0 * ( 1.0e0_DP - zetax )**3.0e0_DP ) * (( 9.0e0_DP * (1.0e0_DP + zetax))    +   &
                &   (( 9.0e0_DP  * ( 1.0e0_DP + 2.0e0_DP * zetax ) * 2.0e0_DP * ( 1.0e0_DP - zetax )  )             +   &
                &   ( 9.e0_DP * ( zetax + zetax**2.0e0_DP ) * 6.0e0_DP ) ) / ( 2.0e0_DP * (1.0e0_DP - zetax) ))  ))
        
        temp_1  = (1.0e0_DP - zetax / 2.0e0_DP) / (1.0e0_DP - zetax)**3.0e0_DP   
        temp_1p = ((-dzetaxv / 2.0e0_DP)*(1.0e0_DP - zetax)**3.0e0_DP + 3.0e0_DP*dzetaxv*(1.0e0_DP              - &
        &   zetax)**2.0e0_DP)*(1.0e0_DP - zetax / 2.0e0_DP) / (1.0e0_DP - zetax)**6.0e0_DP 
        
        temp_2  = zetax * (2.5e0_DP - zetax) / (1.0e0_DP - zetax)**4.0e0_DP
        temp_2p = dzetaxv * (2.5e0_DP - zetax) / (1.0e0_DP - zetax)**4.0e0_DP                                   + &
        &   zetax * (-dzetaxv*(1.0e0_DP - zetax)**4.0e0_DP + 4.0e0_DP*dzetaxv*(1.0e0_DP - zetax)**3.0e0_DP      * &
        &   (2.5e0_DP - zetax)) /  (1.0e0_DP - zetax)**8.0e0_DP
        
        temp_3  = zetax/( 2.0d0 * ( 1.0e0_DP - zetax )**3.0e0_DP )
        temp_3p = (dzetaxv * ( 2.0d0 * ( 1.0e0_DP - zetax )**3.0e0_DP ) + 6.0d0 * dzetaxv * zetax               * &
        &   ( 1.0e0_DP - zetax )**2.0e0_DP) / ( 4.0d0 * ( 1.0e0_DP - zetax )**6.0e0_DP )

        temp_4  = 9.0e0_DP * (1.0e0_DP + zetax) 
        temp_4p = 9.0e0_DP * dzetaxv
 
        temp_5  = 18.0e0_DP * ( 1.0e0_DP + 2.0e0_DP * zetax ) * ( 1.0e0_DP - zetax )
        temp_5p = 18.0e0_DP * (( 2.0e0_DP * dzetaxv ) * ( 1.0e0_DP - zetax )                                    - &
        &    dzetaxv * ( 1.0e0_DP + 2.0e0_DP * zetax ))
        
        temp_6  = 54.0e0_DP * ( zetax + zetax**2.0e0_DP ) 
        temp_6p = 54.0e0_DP * ( dzetaxv + 2.0e0_DP * dzetaxv * zetax ) 
        
        temp_7  = 2.0e0_DP * (1.0e0_DP - zetax)
        temp_7p = -2.0e0_DP * dzetaxv
        
        dl_dBav = l_const2 * (l_Ia * (temp_1p + temp_2p)                                                        - &
        &   l_Ja * (temp_3p*(temp_4 +(temp_5+temp_6)/temp_7)                                                    + &
        &   temp_3*(temp_4p + ((temp_5p+temp_6p)*temp_7 - temp_7p*(temp_5+temp_6))/temp_7**2.0e0_DP)))
        
        dl_dBrv = l_const2 * (l_Ir * (temp_1p + temp_2p)                                                        - &
        &   l_Jr * (temp_3p*(temp_4 +(temp_5+temp_6)/temp_7)                                                    + &
        &   temp_3*(temp_4p + ((temp_5p+temp_6p)*temp_7 - temp_7p*(temp_5+temp_6))/temp_7**2.0e0_DP)))      
      
        as1grad_out = cc * ( xx**lla * (dl_da1v + dl_dBav) - xx**llr * (dl_dr1v + dl_dBrv) )

        return
    end function das1gradv
!***************************************************************************************************
     pure function bchain2( lam, xo, l_dd, l_e, l_r ) result(bchain2_out)
        implicit none
        
        real(kind=DP), intent(in)   ::  lam, xo, l_dd, l_e, l_r
        real(kind=DP)               ::  ilij, jlij
        real(kind=DP)               ::  bchain2_out
        
        ilij = -(xo**(3.0e0_DP - lam) - 1.0e0_DP) / (lam - 3.0e0_DP)  
        jlij = -(xo**(4.0e0_DP - lam) * (lam - 3.0e0_DP) - xo**(3.0e0_DP - lam) * (lam - 4.0e0_DP) - 1.0e0_DP) &
        &   / ((lam - 3.0e0_DP) * (lam - 4.0e0_DP))

        bchain2_out = TWOPI * l_r * l_dd**3.0e0_DP * l_e * &
        &   (   ((1.0e0_DP - zetax / 2.0e0_DP) / (1.0e0_DP - zetax)**3.0e0_DP) * ilij - (9.0e0_DP * zetax *     &
        &   (1.0e0_DP + zetax) / (2.0e0_DP * (1.0e0_DP - zetax)**3.0e0_DP)) * jlij  ) 

        return
    end function bchain2
!***************************************************************************************************   
     pure function dbchain2n( lam, xo, l_dd, l_e, l_r, dl_r ) result(bchain2_out)
        implicit none
        
        real(kind=DP), intent(in)   ::  lam, xo, l_dd, l_e, l_r, dl_r
        real(kind=DP)               ::  ilij, jlij
        real(kind=DP)               ::  bchain2_out
        real(kind=DP)               ::  lb_1,lb_1p,lb_2,lb_2p,lb_3,lb_3p
        
        ilij = -(xo**(3.0e0_DP - lam) - 1.0e0_DP) / (lam - 3.0e0_DP)  
        jlij = -(xo**(4.0e0_DP - lam) * (lam - 3.0e0_DP) - xo**(3.0e0_DP - lam) * (lam - 4.0e0_DP) - 1.0e0_DP) &
        &   / ((lam - 3.0e0_DP) * (lam - 4.0e0_DP))

        lb_1  = TWOPI * l_r * l_dd**3.0e0_DP * l_e 
        lb_1p = TWOPI * dl_r * l_dd**3.0e0_DP * l_e
        
        lb_2  = ((1.0e0_DP - zetax / 2.0e0_DP) / (1.0e0_DP - zetax)**3.0e0_DP) * ilij
        lb_2p = ilij * (-dzetaxn / 2.0e0_DP * (1.0e0_DP - zetax)**3.0e0_DP + 3.0e0_DP * dzetaxn             * &
        &   (1.0e0_DP - zetax)**2.0e0_DP * (1.0e0_DP - zetax / 2.0e0_DP)) / (1.0e0_DP - zetax)**6.0e0_DP
        
        lb_3  = (9.0e0_DP * zetax * (1.0e0_DP + zetax) / (2.0e0_DP * (1.0e0_DP - zetax)**3.0e0_DP)) * jlij
        lb_3p = -jlij * ( 2.0e0_DP * (1.0e0_DP - zetax)**3.0e0_DP * (9.0e0_DP * dzetaxn * (1.0e0_DP + zetax) + &
        &   54.0e0_DP * zetax * (1.0e0_DP + zetax) * dzetaxn * (1.0e0_DP - zetax)**2.0e0_DP     )    )      / &
        &   (4.0e0_DP *  (1.0e0_DP - zetax)**6.0e0_DP )

        bchain2_out = lb_1p * (lb_2 - lb_3) + lb_1 * (lb_2p - lb_3p) 

        return
    end function dbchain2n
!***************************************************************************************************     
     pure function dbchain2v( lam, xo, l_dd, l_e, l_r, dl_r ) result(bchain2_out)
        implicit none
        
        real(kind=DP), intent(in)   ::  lam, xo, l_dd, l_e, l_r, dl_r
        real(kind=DP)               ::  ilij, jlij
        real(kind=DP)               ::  bchain2_out
        real(kind=DP)               ::  lb_1,lb_1p,lb_2,lb_2p,lb_3,lb_3p
        
        ilij = -(xo**(3.0e0_DP - lam) - 1.0e0_DP) / (lam - 3.0e0_DP)  
        jlij = -(xo**(4.0e0_DP - lam) * (lam - 3.0e0_DP) - xo**(3.0e0_DP - lam) * (lam - 4.0e0_DP) - 1.0e0_DP) &
        &   / ((lam - 3.0e0_DP) * (lam - 4.0e0_DP))

        lb_1  = TWOPI * l_r * l_dd**3.0e0_DP * l_e 
        lb_1p = TWOPI * dl_r * l_dd**3.0e0_DP * l_e
        
        lb_2  = ((1.0e0_DP - zetax / 2.0e0_DP) / (1.0e0_DP - zetax)**3.0e0_DP) * ilij
        lb_2p = ilij * (-dzetaxv / 2.0e0_DP * (1.0e0_DP - zetax)**3.0e0_DP + 3.0e0_DP * dzetaxv             * &
        &   (1.0e0_DP - zetax)**2.0e0_DP * (1.0e0_DP - zetax / 2.0e0_DP)) / (1.0e0_DP - zetax)**6.0e0_DP
        
        lb_3  = (9.0e0_DP * zetax * (1.0e0_DP + zetax) / (2.0e0_DP * (1.0e0_DP - zetax)**3.0e0_DP)) * jlij
        lb_3  = -jlij * ( 2.0e0_DP * (1.0e0_DP - zetax)**3.0e0_DP * (9.0e0_DP * dzetaxn * (1.0e0_DP + zetax) + &
        &   54.0e0_DP * zetax * (1.0e0_DP + zetax) * dzetaxv * (1.0e0_DP - zetax)**2.0e0_DP     )    )      / &
        &   (4.0e0_DP *  (1.0e0_DP - zetax)**6.0e0_DP )

        bchain2_out = lb_1p * (lb_2 - lb_3) + lb_1 * (lb_2p - lb_3p) 

        return
    end function dbchain2v
!***************************************************************************************************   
     function as1chain2( lam, l_dd, l_e, l_r ) result(as1chain2_out)
        implicit none

        real(kind=DP), intent(in)   ::  lam, l_dd, l_e, l_r
        real(kind=DP)               ::  zxf
        real(kind=DP)               ::  as1chain2_out
     
        zxf = zeff(lam)

        as1chain2_out = -2.0e0_DP * l_r * (PI * l_e) * (l_dd**3.0e0_DP) / (lam - 3.0e0_DP) *   &
            &  ((1.0e0_DP - zxf / 2.0e0_DP) / (1.0e0_DP - zxf)**3.0e0_DP)

        return
    end function as1chain2
!***************************************************************************************************
     function das1chain2n( lam, l_dd, l_e, l_r, dl_r ) result(as1chain2_out)
        implicit none

        real(kind=DP), intent(in)   ::  lam, l_dd, l_e, l_r, dl_r
        real(kind=DP)               ::  zxf, dzxf
        real(kind=DP)               ::  as1chain2_out
        real(kind=DP)               ::  las1_1,las1_1p,las1_2,las1_2p,las1_3,las1_3p
        
        zxf  = zeff(lam)
        dzxf = dzeffn(lam)
        
        las1_1  = -2.0e0_DP * l_r * (PI * l_e) * (l_dd**3.0e0_DP) / (lam - 3.0e0_DP)
        las1_1p = -2.0e0_DP * dl_r * (PI * l_e) * (l_dd**3.0e0_DP) / (lam - 3.0e0_DP)
        
        las1_2  = 1.0e0_DP - zxf / 2.0e0_DP
        las1_2p = -dzxf / 2.0e0_DP
        
        las1_3  = (1.0e0_DP - zxf)**3.0e0_DP
        las1_3p = -3.0e0_DP * dzxf * (1.0e0_DP - zxf)**2.0e0_DP
        
         
        as1chain2_out = las1_1p * las1_2/las1_3 + las1_1 * (las1_3*las1_2p - las1_2*las1_3p)/las1_3**2.0e0_DP

        return
    end function das1chain2n
!***************************************************************************************************
    function das1chain2v( lam, l_dd, l_e, l_r, dl_r ) result(as1chain2_out)
        implicit none

        real(kind=DP), intent(in)   ::  lam, l_dd, l_e, l_r, dl_r
        real(kind=DP)               ::  zxf, dzxf
        real(kind=DP)               ::  as1chain2_out
        real(kind=DP)               ::  las1_1,las1_1p,las1_2,las1_2p,las1_3,las1_3p
        
        zxf  = zeff(lam)
        dzxf = dzeffv(lam)
       
        las1_1  = -2.0e0_DP * l_r * (PI * l_e) * (l_dd**3.0e0_DP) / (lam - 3.0e0_DP)
        las1_1p = -2.0e0_DP * dl_r * (PI * l_e) * (l_dd**3.0e0_DP) / (lam - 3.0e0_DP)
        
        las1_2  = 1.0e0_DP - zxf / 2.0e0_DP
        las1_2p = -dzxf / 2.0e0_DP
        
        las1_3  = (1.0e0_DP - zxf)**3.0e0_DP
        las1_3p = -3.0e0_DP * dzxf * (1.0e0_DP - zxf)**2.0e0_DP
        
         
        as1chain2_out = las1_1p * las1_2/las1_3 + las1_1 * (las1_3*las1_2p - las1_2*las1_3p)/las1_3**2.0e0_DP
        
        return
    end function das1chain2v
!***************************************************************************************************
    function as2grad( l_rhos, cc, xx, dd, lla, llr, leps ) result(as2grad_out)
        implicit none
        
        real(kind=DP), intent(in) ::  l_rhos,cc,xx,dd,lla,llr,leps
        
        real(kind=DP)       ::  dkhs, d_zetax, l_denom
        real(kind=DP)       ::  l_constaa, l_zaap, l_constar, l_zarp, l_constrr, l_zrrp, l_const2
        real(kind=DP)       ::  l_zaa, l_zar, l_zrr
        real(kind=DP)       ::  l_aa, l_ar, l_rr, l_daa, l_dar, l_drr, l_Baa, l_Bar, l_Brr, l_dBaa, l_dBar, l_dBrr
        real(kind=DP)       ::  diff_v, diff_dv
        real(kind=DP)       ::  l_Iaa, l_Iar, l_Irr, l_Jaa, l_Jar, l_Jrr
        real(kind=DP)       ::  as2grad_out
     
        d_zetax = zetax / l_rhos
    
        l_denom = (1.0e0_DP + 4.0e0_DP * zetax + 4.0e0_DP * zetax**2.0e0_DP - 4.0e0_DP * &
        &   zetax**3.0e0_DP + zetax**4.0e0_DP)
      
        dkhs = -4.0e0_DP * d_zetax * (1.0e0_DP - zetax)**3.0e0_DP   *                    &
        &   (l_denom + (1.0e0_DP + 2.0e0_DP * zetax - 3.0e0_DP * zetax**2.0e0_DP +       &
        &   zetax**3.0e0_DP) * (1.0e0_DP - zetax)) / l_denom**2.0e0_DP

        l_aa    =   as1chain( 2.0e0_DP * lla, dd, leps )
        l_ar    =   as1chain( llr + lla, dd, leps )
        l_rr    =   as1chain( 2.0e0_DP * llr, dd, leps )
       
        l_constaa = PI * leps * dd**3.0e0_DP / (2.0e0_DP * lla - 3.0e0_DP)
        l_zaap    = zeffprime( 2.0e0_DP * lla )
        l_constar = PI * leps * dd**3.0e0_DP / (llr + lla - 3.0e0_DP)
        l_zarp    = zeffprime( llr + lla )
        l_constrr = PI * leps * dd**3.0e0_DP / (2.0e0_DP * llr - 3.0e0_DP)
        l_zrrp    = zeffprime( 2.0e0_DP * llr )
        l_zaa     = zeff( 2.0e0_DP * lla )
        l_zar     = zeff( llr + lla )
        l_zrr     = zeff( 2.0e0_DP * llr )

        l_daa   =   -2.0e0_DP * l_constaa * ( (1.0e0_DP - l_zaa / 2.0e0_DP) / (1.0e0_DP - l_zaa)**3.0e0_DP +  &
        &    l_rhos * (2.5e0_DP * l_zaap - l_zaap * l_zaa) / (1.0e0_DP - l_zaa)**4.0e0_DP)
        l_dar   =   -2.0e0_DP * l_constar * ( (1.0e0_DP - l_zar / 2.0e0_DP) / (1.0e0_DP - l_zar)**3.0e0_DP +  &
        &    l_rhos * (2.5e0_DP * l_zarp - l_zarp * l_zar) / (1.0e0_DP - l_zar)**4.0e0_DP)
        l_drr   =   -2.0e0_DP * l_constrr * ( (1.0e0_DP - l_zrr / 2.0e0_DP) / (1.0e0_DP - l_zrr)**3.0e0_DP +  &
        &    l_rhos * (2.5e0_DP * l_zrrp - l_zrrp * l_zrr) / (1.0e0_DP - l_zrr)**4.0e0_DP)
        
        l_Baa   =   bchain(2.0e0_DP * lla, xx,dd, leps)
        l_Bar   =   bchain(llr + lla, xx, dd, leps)
        l_Brr   =   bchain(2.0e0_DP * llr, xx, dd, leps)

        l_Iaa = -(xx**(3.0e0_DP - 2.0e0_DP * lla) - 1.0e0_DP) / (2.0e0_DP * lla - 3.0e0_DP)  
        l_Iar = -(xx**(3.0e0_DP - llr - lla) - 1.0e0_DP) / (llr + lla - 3.0e0_DP) 
        l_Irr = -(xx**(3.0e0_DP - 2.0e0_DP * llr) - 1.0e0_DP) / (2.0e0_DP * llr - 3.0e0_DP)
        
        l_Jaa = -(xx**(4.0e0_DP - 2.0e0_DP * lla) * (2.0e0_DP * lla - 3.0e0_DP) - xx**(3.0e0_DP - 2.0e0_DP * lla) &
        &   * (2.0e0_DP * lla - 4.0e0_DP) - 1.0e0_DP) / ((2.0e0_DP * lla - 3.0e0_DP) * (2.0e0_DP * lla - 4.0e0_DP))
        l_Jar = -(xx**(4.0e0_DP - lla - llr) * (lla + llr - 3.0e0_DP) - xx**(3.0e0_DP - lla - llr) * (lla + llr - &
        &   4.0e0_DP) - 1.0e0_DP) / ((lla + llr - 3.0e0_DP) * (lla + llr - 4.0e0_DP))
        l_Jrr = -(xx**(4.0e0_DP - 2.0e0_DP * llr) * (2.0e0_DP * llr - 3.0e0_DP) - xx**(3.0e0_DP - 2.0e0_DP * llr) * &
        &   (2.0e0_DP * llr - 4.0e0_DP) - 1.0e0_DP) / ((2.0e0_DP * llr - 3.0e0_DP) * (2.0e0_DP * llr - 4.0e0_DP))   
        
        l_const2 = 2.0e0_DP * PI * dd**3.0e0_DP * leps 
       
        l_dBaa  =   l_const2 * ( (1.0e0_DP - zetax / 2.0e0_DP) / (1.0e0_DP - zetax)**3.0e0_DP * l_Iaa       -   &
        &        9.0e0_DP * zetax * (1.0e0_DP + zetax) / (2.0e0_DP * (1.0e0_DP - zetax)**3.0e0_DP) * l_Jaa  +   &
        &        l_rhos * (((-d_zetax / 2.0e0_DP * (1.0e0_DP - zetax) + 3.0e0_DP * d_zetax * (1.0e0_DP      -   &
        &        zetax /2.0e0_DP)) / (1.0e0_DP - zetax)**4.0e0_DP) * l_Iaa - (((9.0e0_DP * d_zetax          *   &
        &        (1.0e0_DP + zetax) + 9.0e0_DP * zetax * d_zetax) * 2.0e0_DP * (1.0e0_DP - zetax)           +   &
        &        6.0e0_DP * d_zetax * 9.0e0_DP * zetax * (1.0e0_DP + zetax)) / (4.0e0_DP * (1.0e0_DP        -   &
        &        zetax)**4.0e0_DP)) * l_Jaa )) 
        
        l_dBar  =  l_const2 * ( (1.0e0_DP - zetax / 2.0e0_DP) / (1.0e0_DP - zetax)**3.0e0_DP * l_Iar       -   &
        &        9.0e0_DP * zetax * (1.0e0_DP + zetax) / (2.0e0_DP * (1.0e0_DP - zetax)**3.0e0_DP) * l_Jar  +   &
        &        l_rhos * (((-d_zetax / 2.0e0_DP * (1.0e0_DP - zetax) + 3.0e0_DP * d_zetax * (1.0e0_DP      -   &
        &        zetax /2.0e0_DP)) / (1.0e0_DP - zetax)**4.0e0_DP) * l_Iar - (((9.0e0_DP * d_zetax          *   &
        &        (1.0e0_DP + zetax) + 9.0e0_DP * zetax * d_zetax) * 2.0e0_DP * (1.0e0_DP - zetax)           +   &
        &        6.0e0_DP * d_zetax * 9.0e0_DP * zetax * (1.0e0_DP + zetax)) / (4.0e0_DP * (1.0e0_DP        -   &
        &        zetax)**4.0e0_DP)) * l_Jar )) 
     
        l_dBrr  =  l_const2 * ( (1.0e0_DP - zetax / 2.0e0_DP) / (1.0e0_DP - zetax)**3.0e0_DP * l_Irr       -   &
        &        9.0e0_DP * zetax * (1.0e0_DP + zetax) / (2.0e0_DP * (1.0e0_DP - zetax)**3.0e0_DP) * l_Jrr  +   &
        &        l_rhos * (((-d_zetax / 2.0e0_DP * (1.0e0_DP - zetax) + 3.0e0_DP * d_zetax * (1.0e0_DP      -   &
        &        zetax /2.0e0_DP)) / (1.0e0_DP - zetax)**4.0e0_DP) * l_Irr - (((9.0e0_DP * d_zetax          *   &
        &        (1.0e0_DP + zetax) + 9.0e0_DP * zetax * d_zetax) * 2.0e0_DP * (1.0e0_DP - zetax)           +   &
        &        6.0e0_DP * d_zetax * 9.0e0_DP * zetax * (1.0e0_DP + zetax)) / (4.0e0_DP * (1.0e0_DP        -   &
        &        zetax)**4.0e0_DP)) * l_Jrr ))             

        diff_v  = xx**(2.0e0_DP * lla) * (l_aa + l_Baa) - 2.0e0_DP * xx**(lla + llr) * (l_ar + l_Bar) + &
        &       xx**(2.0e0_DP * llr) * (l_rr + l_Brr) 
        
        diff_dv = xx**(2.0e0_DP * lla) * (l_daa + l_dBaa) - 2.0e0_DP * xx**(lla + llr) * (l_dar + l_dBar) + &
        &       xx**(2.0e0_DP * llr) * (l_drr + l_dBrr)
       
        as2grad_out = 0.5e0_DP * leps * cc**2.0e0_DP *(diff_v * dkhs + diff_dv * khs)
       
        return
    end function as2grad
!***************************************************************************************************
    function das2gradn( l_rhos, cc, xx, dd, lla, llr, leps ) result(as2grad_out)
        implicit none
        
        real(kind=DP), intent(in) ::  l_rhos,cc,xx,dd,lla,llr,leps
        
        real(kind=DP)       ::  dkhs, d_zetax, l_denom
        real(kind=DP)       ::  l_constaa, l_zaap, l_constar, l_zarp, l_constrr, l_zrrp, l_const2
        real(kind=DP)       ::  l_zaa, l_zar, l_zrr
        real(kind=DP)       ::  l_aa, l_ar, l_rr, l_daa, l_dar, l_drr, l_Baa, l_Bar, l_Brr, l_dBaa, l_dBar, l_dBrr
        real(kind=DP)       ::  diff_v, diff_dv
        real(kind=DP)       ::  l_Iaa, l_Iar, l_Irr, l_Jaa, l_Jar, l_Jrr
        real(kind=DP)       ::  as2grad_out
        !differentials
        real(kind=DP)       ::  dd_zetaxn, dl_denomn, ddkhsn
        real(kind=DP)       ::  temp_1,temp_2,temp_3,temp_4,temp_5,temp_6,temp_7
        real(kind=DP)       ::  temp_1p,temp_2p,temp_3p,temp_4p,temp_5p,temp_6p,temp_7p
        real(kind=DP)       ::  dl_aan, dl_arn, dl_rrn
        real(kind=DP)       ::  dl_zaapn, dl_zarpn, dl_zrrpn
        real(kind=DP)       ::  dl_zaan, dl_zarn, dl_zrrn
        real(kind=DP)       ::  dl_daan, dl_darn, dl_drrn
        real(kind=DP)       ::  dl_Baan, dl_Barn, dl_Brrn, dl_dBaan, dl_dBarn, dl_dBrrn
        real(kind=DP)       ::  ddiff_dvn, ddiff_vn
        real(kind=DP)       ::  temp_as2u, temp_as2v, temp_as2up, temp_as2vp
        real(kind=DP)       ::  temp_as2_a, temp_as2_ap, temp_as2_b, temp_as2_bp, temp_as2_cp
        
        d_zetax = zetax / l_rhos
        dd_zetaxn = (dzetaxn*l_rhos - drhosn*zetax) / l_rhos**2.0e0_DP
        
        l_denom = 1.0e0_DP + 4.0e0_DP * zetax + 4.0e0_DP * zetax**2.0e0_DP - 4.0e0_DP * &
        &   zetax**3.0e0_DP + zetax**4.0e0_DP
        
        dl_denomn = 4.0e0_DP * dzetaxn + 8.0e0_DP * dzetaxn * zetax - 12.0e0_DP * &
        &   dzetaxn * zetax**2.0e0_DP + 4.0e0_DP * dzetaxn * zetax**3.0e0_DP
        
        dkhs = -4.0e0_DP * d_zetax * (1.0e0_DP - zetax)**3.0e0_DP   *                    &
        &   (l_denom + (1.0e0_DP + 2.0e0_DP * zetax - 3.0e0_DP * zetax**2.0e0_DP +       &
        &   zetax**3.0e0_DP) * (1.0e0_DP - zetax)) / l_denom**2.0e0_DP

        temp_1  = -4.0e0_DP * d_zetax
        temp_1p =  -4.0e0_DP * dd_zetaxn
        
        temp_2  = (1.0e0_DP - zetax)**3.0e0_DP
        temp_2p = -3.0e0_DP * dzetaxn * (1.0e0_DP - zetax)**2.0e0_DP
        
        temp_3  = l_denom + (1.0e0_DP + 2.0e0_DP * zetax - 3.0e0_DP * zetax**2.0e0_DP      + &
        &   zetax**3.0e0_DP) * (1.0e0_DP - zetax)
        temp_3p = dl_denomn + (1.0e0_DP + 2.0e0_DP * zetax - 3.0e0_DP * zetax**2.0e0_DP    + &
        &   zetax**3.0e0_DP) * (-dzetaxn)                                                  + &
        &   (2.0e0_DP * dzetaxn - 6.0e0_DP * dzetaxn * zetax                               + &
        &   3.0e0_DP * dzetaxn * zetax**2.0e0_DP) * (1.0e0_DP - zetax)
        
        temp_4  = l_denom**2.0e0_DP
        temp_4p = 2.0e0_DP * dl_denomn * l_denom
        
        ddkhsn = temp_1p*temp_2*temp_3/temp_4 + temp_1*temp_2p*temp_3/temp_4              + &
        &   temp_1*temp_2 * (temp_4*temp_3p - temp_4p*temp_3)/temp_4**2.0e0_DP

        l_aa    =   as1chain( 2.0e0_DP * lla, dd, leps )
        l_ar    =   as1chain( llr + lla, dd, leps )
        l_rr    =   as1chain( 2.0e0_DP * llr, dd, leps )
        
        l_constaa = PI * leps * dd**3.0e0_DP / (2.0e0_DP * lla - 3.0e0_DP)
        l_zaap    = zeffprime( 2.0e0_DP * lla )
        l_constar = PI * leps * dd**3.0e0_DP / (llr + lla - 3.0e0_DP)
        l_zarp    = zeffprime( llr + lla )
        l_constrr = PI * leps * dd**3.0e0_DP / (2.0e0_DP * llr - 3.0e0_DP)
        l_zrrp    = zeffprime( 2.0e0_DP * llr )
        l_zaa     = zeff( 2.0e0_DP * lla )
        l_zar     = zeff( llr + lla )
        l_zrr     = zeff( 2.0e0_DP * llr )

        dl_aan    =   das1chainn( 2.0e0_DP * lla, dd, leps )
        dl_arn    =   das1chainn( llr + lla, dd, leps )
        dl_rrn    =   das1chainn( 2.0e0_DP * llr, dd, leps )
        
        dl_zaapn    = dzeffprimen( 2.0e0_DP * lla )
        dl_zarpn    = dzeffprimen( llr + lla )
        dl_zrrpn    = dzeffprimen( 2.0e0_DP * llr )
        
        dl_zaan     = dzeffn( 2.0e0_DP * lla )
        dl_zarn     = dzeffn( llr + lla )
        dl_zrrn     = dzeffn( 2.0e0_DP * llr )

        l_daa   =   -2.0e0_DP * l_constaa * ( (1.0e0_DP - l_zaa / 2.0e0_DP) / (1.0e0_DP - l_zaa)**3.0e0_DP +  &
        &    l_rhos * (2.5e0_DP * l_zaap - l_zaap * l_zaa) / (1.0e0_DP - l_zaa)**4.0e0_DP)
        l_dar   =   -2.0e0_DP * l_constar * ( (1.0e0_DP - l_zar / 2.0e0_DP) / (1.0e0_DP - l_zar)**3.0e0_DP +  &
        &    l_rhos * (2.5e0_DP * l_zarp - l_zarp * l_zar) / (1.0e0_DP - l_zar)**4.0e0_DP)
        l_drr   =   -2.0e0_DP * l_constrr * ( (1.0e0_DP - l_zrr / 2.0e0_DP) / (1.0e0_DP - l_zrr)**3.0e0_DP +  &
        &    l_rhos * (2.5e0_DP * l_zrrp - l_zrrp * l_zrr) / (1.0e0_DP - l_zrr)**4.0e0_DP)
     
        temp_1  = 1.0e0_DP - l_zaa / 2.0e0_DP
        temp_1p = -dl_zaan / 2.0e0_DP
        
        temp_2  = (1.0e0_DP - l_zaa)**3.0e0_DP
        temp_2p = -3.0e0_DP * dl_zaan * (1.0e0_DP - l_zaa)**2.0e0_DP
        
        temp_3  = l_rhos * (2.5e0_DP * l_zaap - l_zaap * l_zaa)
        temp_3p = drhosn * (2.5e0_DP * l_zaap - l_zaap * l_zaa) + l_rhos * (2.5e0_DP * dl_zaapn - (dl_zaapn     * &
        &    l_zaa + l_zaap * dl_zaan))
        
        temp_4  = (1.0e0_DP - l_zaa)**4.0e0_DP
        temp_4p = -4.0e0_DP * dl_zaan * (1.0e0_DP - l_zaa)**3.0e0_DP
        
        dl_daan = -2.0e0_DP * l_constaa * ( (temp_2*temp_1p - temp_2p*temp_1)/temp_2**2.0e0_DP                  + &
        &   (temp_4*temp_3p - temp_4p*temp_3)/temp_4**2.0e0_DP)
               
        temp_1  = 1.0e0_DP - l_zar / 2.0e0_DP
        temp_1p = -dl_zarn / 2.0e0_DP
        
        temp_2  = (1.0e0_DP - l_zar)**3.0e0_DP
        temp_2p = -3.0e0_DP * dl_zarn * (1.0e0_DP - l_zar)**2.0e0_DP
        
        temp_3  = l_rhos * (2.5e0_DP * l_zarp - l_zarp * l_zar)
        temp_3p = drhosn * (2.5e0_DP * l_zarp - l_zarp * l_zar) + l_rhos * (2.5e0_DP * dl_zarpn - (dl_zarpn     * &
        &    l_zar + l_zarp * dl_zarn))
        
        temp_4  = (1.0e0_DP - l_zar)**4.0e0_DP
        temp_4p = -4.0e0_DP * dl_zarn * (1.0e0_DP - l_zar)**3.0e0_DP
        
        dl_darn = -2.0e0_DP * l_constar * ( (temp_2*temp_1p - temp_2p*temp_1)/temp_2**2.0e0_DP                  + &
        &   (temp_4*temp_3p - temp_4p*temp_3)/temp_4**2.0e0_DP)                
        
        temp_1  = 1.0e0_DP - l_zrr / 2.0e0_DP
        temp_1p = -dl_zrrn / 2.0e0_DP
        
        temp_2  = (1.0e0_DP - l_zrr)**3.0e0_DP
        temp_2p = -3.0e0_DP * dl_zrrn * (1.0e0_DP - l_zrr)**2.0e0_DP
        
        temp_3  = l_rhos * (2.5e0_DP * l_zrrp - l_zrrp * l_zrr)
        temp_3p = drhosn * (2.5e0_DP * l_zrrp - l_zrrp * l_zrr) + l_rhos * (2.5e0_DP * dl_zrrpn - (dl_zrrpn     * &
        &    l_zrr + l_zrrp * dl_zrrn))
        
        temp_4  = (1.0e0_DP - l_zrr)**4.0e0_DP
        temp_4p = -4.0e0_DP * dl_zrrn * (1.0e0_DP - l_zrr)**3.0e0_DP
        
        dl_drrn = -2.0e0_DP * l_constrr * ( (temp_2*temp_1p - temp_2p*temp_1)/temp_2**2.0e0_DP                  + &
        &   (temp_4*temp_3p - temp_4p*temp_3)/temp_4**2.0e0_DP) 

        l_Baa   =   bchain(2.0e0_DP * lla, xx,dd, leps)
        l_Bar   =   bchain(llr + lla, xx, dd, leps)
        l_Brr   =   bchain(2.0e0_DP * llr, xx, dd, leps)
        
        dl_Baan   =   dbchainn(2.0e0_DP * lla, xx,dd, leps)
        dl_Barn   =   dbchainn(llr + lla, xx, dd, leps)
        dl_Brrn   =   dbchainn(2.0e0_DP * llr, xx, dd, leps)

        l_Iaa = -(xx**(3.0e0_DP - 2.0e0_DP * lla) - 1.0e0_DP) / (2.0e0_DP * lla - 3.0e0_DP)  
        l_Iar = -(xx**(3.0e0_DP - llr - lla) - 1.0e0_DP) / (llr + lla - 3.0e0_DP) 
        l_Irr = -(xx**(3.0e0_DP - 2.0e0_DP * llr) - 1.0e0_DP) / (2.0e0_DP * llr - 3.0e0_DP)
        
        l_Jaa = -(xx**(4.0e0_DP - 2.0e0_DP * lla) * (2.0e0_DP * lla - 3.0e0_DP) - xx**(3.0e0_DP - 2.0e0_DP * lla) &
        &   * (2.0e0_DP * lla - 4.0e0_DP) - 1.0e0_DP) / ((2.0e0_DP * lla - 3.0e0_DP) * (2.0e0_DP * lla - 4.0e0_DP))
        l_Jar = -(xx**(4.0e0_DP - lla - llr) * (lla + llr - 3.0e0_DP) - xx**(3.0e0_DP - lla - llr) * (lla + llr - &
        &   4.0e0_DP) - 1.0e0_DP) / ((lla + llr - 3.0e0_DP) * (lla + llr - 4.0e0_DP))
        l_Jrr = -(xx**(4.0e0_DP - 2.0e0_DP * llr) * (2.0e0_DP * llr - 3.0e0_DP) - xx**(3.0e0_DP - 2.0e0_DP * llr) * &
        &   (2.0e0_DP * llr - 4.0e0_DP) - 1.0e0_DP) / ((2.0e0_DP * llr - 3.0e0_DP) * (2.0e0_DP * llr - 4.0e0_DP))   
        
        l_const2 = 2.0e0_DP * PI * dd**3.0e0_DP * leps 

        l_dBaa  =   l_const2 * ( (1.0e0_DP - zetax / 2.0e0_DP) / (1.0e0_DP - zetax)**3.0e0_DP * l_Iaa       -   &
        &        9.0e0_DP * zetax * (1.0e0_DP + zetax) / (2.0e0_DP * (1.0e0_DP - zetax)**3.0e0_DP) * l_Jaa  +   &
        &        l_rhos * (((-d_zetax / 2.0e0_DP * (1.0e0_DP - zetax) + 3.0e0_DP * d_zetax * (1.0e0_DP      -   &
        &        zetax /2.0e0_DP)) / (1.0e0_DP - zetax)**4.0e0_DP) * l_Iaa - (((9.0e0_DP * d_zetax          *   &
        &        (1.0e0_DP + zetax) + 9.0e0_DP * zetax * d_zetax) * 2.0e0_DP * (1.0e0_DP - zetax)           +   &
        &        6.0e0_DP * d_zetax * 9.0e0_DP * zetax * (1.0e0_DP + zetax)) / (4.0e0_DP * (1.0e0_DP        -   &
        &        zetax)**4.0e0_DP)) * l_Jaa )) 

        l_dBar  =  l_const2 * ( (1.0e0_DP - zetax / 2.0e0_DP) / (1.0e0_DP - zetax)**3.0e0_DP * l_Iar       -   &
        &        9.0e0_DP * zetax * (1.0e0_DP + zetax) / (2.0e0_DP * (1.0e0_DP - zetax)**3.0e0_DP) * l_Jar  +   &
        &        l_rhos * (((-d_zetax / 2.0e0_DP * (1.0e0_DP - zetax) + 3.0e0_DP * d_zetax * (1.0e0_DP      -   &
        &        zetax /2.0e0_DP)) / (1.0e0_DP - zetax)**4.0e0_DP) * l_Iar - (((9.0e0_DP * d_zetax          *   &
        &        (1.0e0_DP + zetax) + 9.0e0_DP * zetax * d_zetax) * 2.0e0_DP * (1.0e0_DP - zetax)           +   &
        &        6.0e0_DP * d_zetax * 9.0e0_DP * zetax * (1.0e0_DP + zetax)) / (4.0e0_DP * (1.0e0_DP        -   &
        &        zetax)**4.0e0_DP)) * l_Jar )) 
     
        l_dBrr  =  l_const2 * ( (1.0e0_DP - zetax / 2.0e0_DP) / (1.0e0_DP - zetax)**3.0e0_DP * l_Irr       -   &
        &        9.0e0_DP * zetax * (1.0e0_DP + zetax) / (2.0e0_DP * (1.0e0_DP - zetax)**3.0e0_DP) * l_Jrr  +   &
        &        l_rhos * (((-d_zetax / 2.0e0_DP * (1.0e0_DP - zetax) + 3.0e0_DP * d_zetax * (1.0e0_DP      -   &
        &        zetax /2.0e0_DP)) / (1.0e0_DP - zetax)**4.0e0_DP) * l_Irr - (((9.0e0_DP * d_zetax          *   &
        &        (1.0e0_DP + zetax) + 9.0e0_DP * zetax * d_zetax) * 2.0e0_DP * (1.0e0_DP - zetax)           +   &
        &        6.0e0_DP * d_zetax * 9.0e0_DP * zetax * (1.0e0_DP + zetax)) / (4.0e0_DP * (1.0e0_DP        -   &
        &        zetax)**4.0e0_DP)) * l_Jrr )) 
      
        temp_1  = (1.0e0_DP - zetax / 2.0e0_DP) / (1.0e0_DP - zetax)**3.0e0_DP 
        temp_1p = ((1.0e0_DP - zetax)**3.0e0_DP * (-dzetaxn / 2.0e0_DP)                         + &
        &   3.0e0_DP * dzetaxn * (1.0e0_DP - zetax / 2.0e0_DP) * (1.0e0_DP - zetax)**2.0e0_DP)  / &
        &   (1.0e0_DP - zetax)**6.0e0_DP 
        
        temp_2  = -9.0e0_DP * zetax * (1.0e0_DP + zetax) / (2.0e0_DP * (1.0e0_DP - zetax)**3.0e0_DP) 
        temp_2p = ( (2.0e0_DP * (1.0e0_DP - zetax)**3.0e0_DP) * (-9.0e0_DP * dzetaxn * (1.0e0_DP        + &
        &   zetax) -9.0e0_DP * zetax * dzetaxn) - (9.0e0_DP * zetax * (1.0e0_DP + zetax))               * &
        &   (6.0e0_DP * dzetaxn * (1.0e0_DP - zetax)**2.0e0_DP)) / (4.0e0_DP * (1.0e0_DP - zetax)**6.0e0_DP)
        
        temp_3 = (-d_zetax / 2.0e0_DP * (1.0e0_DP - zetax) + 3.0e0_DP * d_zetax * (1.0e0_DP      -   &
        &        zetax /2.0e0_DP)) / (1.0e0_DP - zetax)**4.0e0_DP
        temp_as2u  = -d_zetax / 2.0e0_DP * (1.0e0_DP - zetax) + 3.0e0_DP * d_zetax * (1.0e0_DP - zetax /2.0e0_DP)
        temp_as2up = -dd_zetaxn/2.0e0_DP * (1.0e0_DP - zetax) + d_zetax * dzetaxn / 2.0e0_DP  + &
        &   3.0e0_DP*dd_zetaxn * (1.0e0_DP - zetax/2.0e0_DP) - 1.5e0_DP*d_zetax * dzetaxn
        temp_as2v  = (1.0e0_DP - zetax)**4.0e0_DP
        temp_as2vp = -4.0e0_DP * dzetaxn*(1.0e0_DP - zetax)**3.0e0_DP     
        temp_3p = (temp_as2v*temp_as2up - temp_as2u*temp_as2vp) / temp_as2v**2.0e0_DP 
        
        temp_4 = ((9.0e0_DP * d_zetax          *   &
        &        (1.0e0_DP + zetax) + 9.0e0_DP * zetax * d_zetax) * 2.0e0_DP * (1.0e0_DP - zetax)           +   &
        &        6.0e0_DP * d_zetax * 9.0e0_DP * zetax * (1.0e0_DP + zetax)) / (4.0e0_DP * (1.0e0_DP        -   &
        &        zetax)**4.0e0_DP)
        temp_as2u  = (9.0e0_DP * d_zetax * (1.0e0_DP + zetax) + 9.0e0_DP * zetax * d_zetax) * 2.0e0_DP      *   &
        &   (1.0e0_DP - zetax) + 6.0e0_DP * d_zetax * 9.0e0_DP * zetax * (1.0e0_DP + zetax)
        
        temp_as2_a  = 9.0e0_DP * d_zetax * (1.0e0_DP + zetax) + 9.0e0_DP * zetax * d_zetax
        temp_as2_ap = 9.0e0_DP*dd_zetaxn * (1.0e0_DP + zetax) +  9.0e0_DP * d_zetax * dzetaxn + &
        &  9.0e0_DP * dzetaxn * d_zetax + 9.0e0_DP * zetax * dd_zetaxn
        temp_as2_b  = 2.0e0_DP * (1.0e0_DP - zetax)
        temp_as2_bp = -2.0e0_DP*dzetaxn
        temp_as2_cp = 54.0e0_DP*(dd_zetaxn * zetax * (1.0e0_DP + zetax) + d_zetax * dzetaxn *   &
        &   (1.0e0_DP + zetax) + d_zetax * zetax * dzetaxn)
        
        temp_as2up = temp_as2_ap*temp_as2_b + temp_as2_a*temp_as2_bp + temp_as2_cp
        
        temp_as2v  = 4.0e0_DP * (1.0e0_DP - zetax)**4.0e0_DP
        temp_as2vp = -16.0e0_DP * dzetaxn * (1.0e0_DP - zetax)**3.0e0_DP
        temp_4p = (temp_as2v*temp_as2up - temp_as2u*temp_as2vp) / temp_as2v**2.0e0_DP

        dl_dBaan = l_const2 * (temp_1p*l_Iaa + temp_2p*l_Jaa + l_rhos * (temp_3p * l_Iaa - temp_4p * l_Jaa) + &
        &           drhosn * (temp_3 * l_Iaa - temp_4 * l_Jaa))
        dl_dBarn = l_const2 * (temp_1p*l_Iar + temp_2p*l_Jar + l_rhos * (temp_3p * l_Iar - temp_4p * l_Jar) + &
        &           drhosn * (temp_3 * l_Iar - temp_4 * l_Jar))
        dl_dBrrn = l_const2 * (temp_1p*l_Irr + temp_2p*l_Jrr + l_rhos * (temp_3p * l_Irr - temp_4p * l_Jrr) + &
        &           drhosn * (temp_3 * l_Irr - temp_4 * l_Jrr))     

        diff_v  = xx**(2.0e0_DP * lla) * (l_aa + l_Baa) - 2.0e0_DP * xx**(lla + llr) * (l_ar + l_Bar) + &
        &       xx**(2.0e0_DP * llr) * (l_rr + l_Brr) 
        diff_dv = xx**(2.0e0_DP * lla) * (l_daa + l_dBaa) - 2.0e0_DP * xx**(lla + llr) * (l_dar + l_dBar) + &
        &       xx**(2.0e0_DP * llr) * (l_drr + l_dBrr)
        
        ddiff_vn  = xx**(2.0e0_DP * lla) * (dl_aan + dl_Baan) - 2.0e0_DP * xx**(lla + llr) * (dl_arn + dl_Barn) + &
        &       xx**(2.0e0_DP * llr) * (dl_rrn + dl_Brrn) 
        ddiff_dvn = xx**(2.0e0_DP * lla) * (dl_daan + dl_dBaan) - 2.0e0_DP*xx**(lla + llr) * (dl_darn + dl_dBarn) + &
        &       xx**(2.0e0_DP * llr) * (dl_drrn + dl_dBrrn)
   
        as2grad_out = 0.5e0_DP * leps * cc**2.0e0_DP *(ddiff_vn * dkhs + diff_v * ddkhsn  + &
        &   ddiff_dvn * khs + diff_dv * dkhsn)
   
        return
    end function das2gradn
!***************************************************************************************************
    function das2gradv( l_rhos, cc, xx, dd, lla, llr, leps ) result(as2grad_out)
        implicit none
        
        real(kind=DP), intent(in) ::  l_rhos,cc,xx,dd,lla,llr,leps
        
        real(kind=DP)       ::  dkhs, d_zetax, l_denom
        real(kind=DP)       ::  l_constaa, l_zaap, l_constar, l_zarp, l_constrr, l_zrrp, l_const2
        real(kind=DP)       ::  l_zaa, l_zar, l_zrr
        real(kind=DP)       ::  l_aa, l_ar, l_rr, l_daa, l_dar, l_drr, l_Baa, l_Bar, l_Brr, l_dBaa, l_dBar, l_dBrr
        real(kind=DP)       ::  diff_v, diff_dv
        real(kind=DP)       ::  l_Iaa, l_Iar, l_Irr, l_Jaa, l_Jar, l_Jrr
        real(kind=DP)       ::  as2grad_out
        !differentials
        real(kind=DP)       ::  dd_zetaxv, dl_denomv, ddkhsv
        real(kind=DP)       ::  temp_1,temp_2,temp_3,temp_4,temp_5,temp_6,temp_7
        real(kind=DP)       ::  temp_1p,temp_2p,temp_3p,temp_4p,temp_5p,temp_6p,temp_7p
        real(kind=DP)       ::  dl_aav, dl_arv, dl_rrv
        real(kind=DP)       ::  dl_zaapv, dl_zarpv, dl_zrrpv
        real(kind=DP)       ::  dl_zaav, dl_zarv, dl_zrrv
        real(kind=DP)       ::  dl_daav, dl_darv, dl_drrv
        real(kind=DP)       ::  dl_Baav, dl_Barv, dl_Brrv, dl_dBaav, dl_dBarv, dl_dBrrv
        real(kind=DP)       ::  ddiff_dvv, ddiff_vv
        real(kind=DP)       ::  temp_as2u, temp_as2v, temp_as2up, temp_as2vp
        real(kind=DP)       ::  temp_as2_a, temp_as2_ap, temp_as2_b, temp_as2_bp, temp_as2_cp
        
        d_zetax = zetax / l_rhos
        dd_zetaxv = (drhosv*zetax - dzetaxv*l_rhos) / l_rhos**2.0e0_DP
      
        l_denom = 1.0e0_DP + 4.0e0_DP * zetax + 4.0e0_DP * zetax**2.0e0_DP - 4.0e0_DP * &
        &   zetax**3.0e0_DP + zetax**4.0e0_DP
        
        dl_denomv = 4.0e0_DP * dzetaxv + 8.0e0_DP * dzetaxv * zetax - 12.0e0_DP * &
        &   dzetaxv * zetax**2.0e0_DP + 4.0e0_DP * dzetaxv * zetax**3.0e0_DP
      
        dkhs = -4.0e0_DP * d_zetax * (1.0e0_DP - zetax)**3.0e0_DP   *                    &
        &   (l_denom + (1.0e0_DP + 2.0e0_DP * zetax - 3.0e0_DP * zetax**2.0e0_DP +       &
        &   zetax**3.0e0_DP) * (1.0e0_DP - zetax)) / l_denom**2.0e0_DP

        temp_1  = -4.0e0_DP * d_zetax
        temp_1p =  -4.0e0_DP * dd_zetaxv
        
        temp_2  = (1.0e0_DP - zetax)**3.0e0_DP
        temp_2p = -3.0e0_DP * dzetaxv * (1.0e0_DP - zetax)**2.0e0_DP
        
        temp_3  = l_denom + (1.0e0_DP + 2.0e0_DP * zetax - 3.0e0_DP * zetax**2.0e0_DP      + &
        &   zetax**3.0e0_DP) * (1.0e0_DP - zetax)
        temp_3p = dl_denomv + (1.0e0_DP + 2.0e0_DP * zetax - 3.0e0_DP * zetax**2.0e0_DP    + &
        &   zetax**3.0e0_DP) * (-dzetaxv)                                                  + &
        &   (2.0e0_DP * dzetaxv - 6.0e0_DP * dzetaxv * zetax                               + &
        &   3.0e0_DP * dzetaxv * zetax**2.0e0_DP) * (1.0e0_DP - zetax)
        
        temp_4  = l_denom**2.0e0_DP
        temp_4p = 2.0e0_DP * dl_denomv * l_denom
        
        ddkhsv = temp_1p*temp_2*temp_3/temp_4 + temp_1*temp_2p*temp_3/temp_4              + &
        &   temp_1*temp_2 * (temp_4*temp_3p - temp_4p*temp_3)/temp_4**2.0e0_DP

        l_aa    =   as1chain( 2.0e0_DP * lla, dd, leps )
        l_ar    =   as1chain( llr + lla, dd, leps )
        l_rr    =   as1chain( 2.0e0_DP * llr, dd, leps )
        
        l_constaa = PI * leps * dd**3.0e0_DP / (2.0e0_DP * lla - 3.0e0_DP)
        l_zaap    = zeffprime( 2.0e0_DP * lla )
        l_constar = PI * leps * dd**3.0e0_DP / (llr + lla - 3.0e0_DP)
        l_zarp    = zeffprime( llr + lla )
        l_constrr = PI * leps * dd**3.0e0_DP / (2.0e0_DP * llr - 3.0e0_DP)
        l_zrrp    = zeffprime( 2.0e0_DP * llr )
        l_zaa     = zeff( 2.0e0_DP * lla )
        l_zar     = zeff( llr + lla )
        l_zrr     = zeff( 2.0e0_DP * llr )

        dl_aav    =   das1chainv( 2.0e0_DP * lla, dd, leps )
        dl_arv    =   das1chainv( llr + lla, dd, leps )
        dl_rrv    =   das1chainv( 2.0e0_DP * llr, dd, leps )
   
        dl_zaapv    = dzeffprimev( 2.0e0_DP * lla )
        dl_zarpv    = dzeffprimev( llr + lla )
        dl_zrrpv    = dzeffprimev( 2.0e0_DP * llr )
        
        dl_zaav     = dzeffv( 2.0e0_DP * lla )
        dl_zarv     = dzeffv( llr + lla )
        dl_zrrv     = dzeffv( 2.0e0_DP * llr )

        l_daa   =   -2.0e0_DP * l_constaa * ( (1.0e0_DP - l_zaa / 2.0e0_DP) / (1.0e0_DP - l_zaa)**3.0e0_DP +  &
        &    l_rhos * (2.5e0_DP * l_zaap - l_zaap * l_zaa) / (1.0e0_DP - l_zaa)**4.0e0_DP)
        l_dar   =   -2.0e0_DP * l_constar * ( (1.0e0_DP - l_zar / 2.0e0_DP) / (1.0e0_DP - l_zar)**3.0e0_DP +  &
        &    l_rhos * (2.5e0_DP * l_zarp - l_zarp * l_zar) / (1.0e0_DP - l_zar)**4.0e0_DP)
        l_drr   =   -2.0e0_DP * l_constrr * ( (1.0e0_DP - l_zrr / 2.0e0_DP) / (1.0e0_DP - l_zrr)**3.0e0_DP +  &
        &    l_rhos * (2.5e0_DP * l_zrrp - l_zrrp * l_zrr) / (1.0e0_DP - l_zrr)**4.0e0_DP)
        
        temp_1  = 1.0e0_DP - l_zaa / 2.0e0_DP
        temp_1p = -dl_zaav / 2.0e0_DP
        
        temp_2  = (1.0e0_DP - l_zaa)**3.0e0_DP
        temp_2p = -3.0e0_DP * dl_zaav * (1.0e0_DP - l_zaa)**2.0e0_DP
        
        temp_3  = l_rhos * (2.5e0_DP * l_zaap - l_zaap * l_zaa)
        temp_3p = drhosv * (2.5e0_DP * l_zaap - l_zaap * l_zaa) + l_rhos * (2.5e0_DP * dl_zaapv - (dl_zaapv     * &
        &    l_zaa + l_zaap * dl_zaav))
        
        temp_4  = (1.0e0_DP - l_zaa)**4.0e0_DP
        temp_4p = -4.0e0_DP * dl_zaav * (1.0e0_DP - l_zaa)**3.0e0_DP
        
        dl_daav = -2.0e0_DP * l_constaa * ( (temp_2*temp_1p - temp_2p*temp_1)/temp_2**2.0e0_DP                  + &
        &   (temp_4*temp_3p - temp_4p*temp_3)/temp_4**2.0e0_DP)
             
        temp_1  = 1.0e0_DP - l_zar / 2.0e0_DP
        temp_1p = -dl_zarv / 2.0e0_DP
        
        temp_2  = (1.0e0_DP - l_zar)**3.0e0_DP
        temp_2p = -3.0e0_DP * dl_zarv * (1.0e0_DP - l_zar)**2.0e0_DP
        
        temp_3  = l_rhos * (2.5e0_DP * l_zarp - l_zarp * l_zar)
        temp_3p = drhosv * (2.5e0_DP * l_zarp - l_zarp * l_zar) + l_rhos * (2.5e0_DP * dl_zarpv - (dl_zarpv     * &
        &    l_zar + l_zarp * dl_zarv))
        
        temp_4  = (1.0e0_DP - l_zar)**4.0e0_DP
        temp_4p = -4.0e0_DP * dl_zarv * (1.0e0_DP - l_zar)**3.0e0_DP
        
        dl_darv = -2.0e0_DP * l_constar * ( (temp_2*temp_1p - temp_2p*temp_1)/temp_2**2.0e0_DP                  + &
        &   (temp_4*temp_3p - temp_4p*temp_3)/temp_4**2.0e0_DP)        
          
        temp_1  = 1.0e0_DP - l_zrr / 2.0e0_DP
        temp_1p = -dl_zrrv / 2.0e0_DP
        
        temp_2  = (1.0e0_DP - l_zrr)**3.0e0_DP
        temp_2p = -3.0e0_DP * dl_zrrv * (1.0e0_DP - l_zrr)**2.0e0_DP
        
        temp_3  = l_rhos * (2.5e0_DP * l_zrrp - l_zrrp * l_zrr)
        temp_3p = drhosv * (2.5e0_DP * l_zrrp - l_zrrp * l_zrr) + l_rhos * (2.5e0_DP * dl_zrrpv - (dl_zrrpv     * &
        &    l_zrr + l_zrrp * dl_zrrv))
        
        temp_4  = (1.0e0_DP - l_zrr)**4.0e0_DP
        temp_4p = -4.0e0_DP * dl_zrrv * (1.0e0_DP - l_zrr)**3.0e0_DP
        
        dl_drrv = -2.0e0_DP * l_constrr * ( (temp_2*temp_1p - temp_2p*temp_1)/temp_2**2.0e0_DP                  + &
        &   (temp_4*temp_3p - temp_4p*temp_3)/temp_4**2.0e0_DP) 

        l_Baa   =   bchain(2.0e0_DP * lla, xx,dd, leps)
        l_Bar   =   bchain(llr + lla, xx, dd, leps)
        l_Brr   =   bchain(2.0e0_DP * llr, xx, dd, leps)
        
        dl_Baav   =   dbchainv(2.0e0_DP * lla, xx,dd, leps)
        dl_Barv   =   dbchainv(llr + lla, xx, dd, leps)
        dl_Brrv   =   dbchainv(2.0e0_DP * llr, xx, dd, leps)

        l_Iaa = -(xx**(3.0e0_DP - 2.0e0_DP * lla) - 1.0e0_DP) / (2.0e0_DP * lla - 3.0e0_DP)  
        l_Iar = -(xx**(3.0e0_DP - llr - lla) - 1.0e0_DP) / (llr + lla - 3.0e0_DP) 
        l_Irr = -(xx**(3.0e0_DP - 2.0e0_DP * llr) - 1.0e0_DP) / (2.0e0_DP * llr - 3.0e0_DP)
        
        l_Jaa = -(xx**(4.0e0_DP - 2.0e0_DP * lla) * (2.0e0_DP * lla - 3.0e0_DP) - xx**(3.0e0_DP - 2.0e0_DP * lla) &
        &   * (2.0e0_DP * lla - 4.0e0_DP) - 1.0e0_DP) / ((2.0e0_DP * lla - 3.0e0_DP) * (2.0e0_DP * lla - 4.0e0_DP))
        l_Jar = -(xx**(4.0e0_DP - lla - llr) * (lla + llr - 3.0e0_DP) - xx**(3.0e0_DP - lla - llr) * (lla + llr - &
        &   4.0e0_DP) - 1.0e0_DP) / ((lla + llr - 3.0e0_DP) * (lla + llr - 4.0e0_DP))
        l_Jrr = -(xx**(4.0e0_DP - 2.0e0_DP * llr) * (2.0e0_DP * llr - 3.0e0_DP) - xx**(3.0e0_DP - 2.0e0_DP * llr) * &
        &   (2.0e0_DP * llr - 4.0e0_DP) - 1.0e0_DP) / ((2.0e0_DP * llr - 3.0e0_DP) * (2.0e0_DP * llr - 4.0e0_DP))   
        
        l_const2 = 2.0e0_DP * PI * dd**3.0e0_DP * leps 

        l_dBaa  =   l_const2 * ( (1.0e0_DP - zetax / 2.0e0_DP) / (1.0e0_DP - zetax)**3.0e0_DP * l_Iaa       -   &
        &        9.0e0_DP * zetax * (1.0e0_DP + zetax) / (2.0e0_DP * (1.0e0_DP - zetax)**3.0e0_DP) * l_Jaa  +   &
        &        l_rhos * (((-d_zetax / 2.0e0_DP * (1.0e0_DP - zetax) + 3.0e0_DP * d_zetax * (1.0e0_DP      -   &
        &        zetax /2.0e0_DP)) / (1.0e0_DP - zetax)**4.0e0_DP) * l_Iaa - (((9.0e0_DP * d_zetax          *   &
        &        (1.0e0_DP + zetax) + 9.0e0_DP * zetax * d_zetax) * 2.0e0_DP * (1.0e0_DP - zetax)           +   &
        &        6.0e0_DP * d_zetax * 9.0e0_DP * zetax * (1.0e0_DP + zetax)) / (4.0e0_DP * (1.0e0_DP        -   &
        &        zetax)**4.0e0_DP)) * l_Jaa )) 

        l_dBar  =  l_const2 * ( (1.0e0_DP - zetax / 2.0e0_DP) / (1.0e0_DP - zetax)**3.0e0_DP * l_Iar       -   &
        &        9.0e0_DP * zetax * (1.0e0_DP + zetax) / (2.0e0_DP * (1.0e0_DP - zetax)**3.0e0_DP) * l_Jar  +   &
        &        l_rhos * (((-d_zetax / 2.0e0_DP * (1.0e0_DP - zetax) + 3.0e0_DP * d_zetax * (1.0e0_DP      -   &
        &        zetax /2.0e0_DP)) / (1.0e0_DP - zetax)**4.0e0_DP) * l_Iar - (((9.0e0_DP * d_zetax          *   &
        &        (1.0e0_DP + zetax) + 9.0e0_DP * zetax * d_zetax) * 2.0e0_DP * (1.0e0_DP - zetax)           +   &
        &        6.0e0_DP * d_zetax * 9.0e0_DP * zetax * (1.0e0_DP + zetax)) / (4.0e0_DP * (1.0e0_DP        -   &
        &        zetax)**4.0e0_DP)) * l_Jar )) 
     
        l_dBrr  =  l_const2 * ( (1.0e0_DP - zetax / 2.0e0_DP) / (1.0e0_DP - zetax)**3.0e0_DP * l_Irr       -   &
        &        9.0e0_DP * zetax * (1.0e0_DP + zetax) / (2.0e0_DP * (1.0e0_DP - zetax)**3.0e0_DP) * l_Jrr  +   &
        &        l_rhos * (((-d_zetax / 2.0e0_DP * (1.0e0_DP - zetax) + 3.0e0_DP * d_zetax * (1.0e0_DP      -   &
        &        zetax /2.0e0_DP)) / (1.0e0_DP - zetax)**4.0e0_DP) * l_Irr - (((9.0e0_DP * d_zetax          *   &
        &        (1.0e0_DP + zetax) + 9.0e0_DP * zetax * d_zetax) * 2.0e0_DP * (1.0e0_DP - zetax)           +   &
        &        6.0e0_DP * d_zetax * 9.0e0_DP * zetax * (1.0e0_DP + zetax)) / (4.0e0_DP * (1.0e0_DP        -   &
        &        zetax)**4.0e0_DP)) * l_Jrr )) 
        
        temp_1  = (1.0e0_DP - zetax / 2.0e0_DP) / (1.0e0_DP - zetax)**3.0e0_DP 
        temp_1p = ((1.0e0_DP - zetax)**3.0e0_DP * (-dzetaxv / 2.0e0_DP)                         + &
        &   3.0e0_DP * dzetaxv * (1.0e0_DP - zetax / 2.0e0_DP) * (1.0e0_DP - zetax)**2.0e0_DP)  / &
        &   (1.0e0_DP - zetax)**6.0e0_DP 
        
        temp_2  = -9.0e0_DP * zetax * (1.0e0_DP + zetax) / (2.0e0_DP * (1.0e0_DP - zetax)**3.0e0_DP) 
        temp_2p = ( (2.0e0_DP * (1.0e0_DP - zetax)**3.0e0_DP) * (-9.0e0_DP * dzetaxv * (1.0e0_DP        + &
        &   zetax) -9.0e0_DP * zetax * dzetaxv) - (9.0e0_DP * zetax * (1.0e0_DP + zetax))               * &
        &   (6.0e0_DP * dzetaxv * (1.0e0_DP - zetax)**2.0e0_DP)) / (4.0e0_DP * (1.0e0_DP - zetax)**6.0e0_DP)
        
        temp_3 = (-d_zetax / 2.0e0_DP * (1.0e0_DP - zetax) + 3.0e0_DP * d_zetax * (1.0e0_DP      -   &
        &        zetax /2.0e0_DP)) / (1.0e0_DP - zetax)**4.0e0_DP
        temp_as2u  = -d_zetax / 2.0e0_DP * (1.0e0_DP - zetax) + 3.0e0_DP * d_zetax * (1.0e0_DP - zetax /2.0e0_DP)
        temp_as2up = -dd_zetaxv/2.0e0_DP * (1.0e0_DP - zetax) + d_zetax * dzetaxv / 2.0e0_DP  + &
        &   3.0e0_DP*dd_zetaxv * (1.0e0_DP - zetax/2.0e0_DP) - 1.5e0_DP*d_zetax * dzetaxv
        temp_as2v  = (1.0e0_DP - zetax)**4.0e0_DP
        temp_as2vp = -4.0e0_DP * dzetaxv*(1.0e0_DP - zetax)**3.0e0_DP     
        temp_3p = (temp_as2v*temp_as2up - temp_as2u*temp_as2vp) / temp_as2v**2.0e0_DP 

        temp_4 = ((9.0e0_DP * d_zetax          *   &
        &        (1.0e0_DP + zetax) + 9.0e0_DP * zetax * d_zetax) * 2.0e0_DP * (1.0e0_DP - zetax)           +   &
        &        6.0e0_DP * d_zetax * 9.0e0_DP * zetax * (1.0e0_DP + zetax)) / (4.0e0_DP * (1.0e0_DP        -   &
        &        zetax)**4.0e0_DP)
        temp_as2u  = (9.0e0_DP * d_zetax * (1.0e0_DP + zetax) + 9.0e0_DP * zetax * d_zetax) * 2.0e0_DP      *   &
        &   (1.0e0_DP - zetax) + 6.0e0_DP * d_zetax * 9.0e0_DP * zetax * (1.0e0_DP + zetax)

        temp_as2_a  = 9.0e0_DP * d_zetax * (1.0e0_DP + zetax) + 9.0e0_DP * zetax * d_zetax
        temp_as2_ap = 9.0e0_DP*dd_zetaxv * (1.0e0_DP + zetax) +  9.0e0_DP * d_zetax * dzetaxv + &
        &  9.0e0_DP * dzetaxv * d_zetax + 9.0e0_DP * zetax * dd_zetaxv
        temp_as2_b  = 2.0e0_DP * (1.0e0_DP - zetax)
        temp_as2_bp = -2.0e0_DP*dzetaxv
        temp_as2_cp = 54.0e0_DP*(dd_zetaxv * zetax * (1.0e0_DP + zetax) + d_zetax * dzetaxv *   &
        &   (1.0e0_DP + zetax) + d_zetax * zetax * dzetaxv)

        temp_as2up = temp_as2_ap*temp_as2_b + temp_as2_a*temp_as2_bp + temp_as2_cp
        
        temp_as2v  = 4.0e0_DP * (1.0e0_DP - zetax)**4.0e0_DP
        temp_as2vp = -16.0e0_DP * dzetaxv * (1.0e0_DP - zetax)**3.0e0_DP
        temp_4p = (temp_as2v*temp_as2up - temp_as2u*temp_as2vp) / temp_as2v**2.0e0_DP
        
        dl_dBaav = l_const2 * (temp_1p*l_Iaa + temp_2p*l_Jaa + l_rhos * (temp_3p * l_Iaa - temp_4p * l_Jaa) + &
        &           drhosv * (temp_3 * l_Iaa - temp_4 * l_Jaa))
        dl_dBarv = l_const2 * (temp_1p*l_Iar + temp_2p*l_Jar + l_rhos * (temp_3p * l_Iar - temp_4p * l_Jar) + &
        &           drhosv * (temp_3 * l_Iar - temp_4 * l_Jar))
        dl_dBrrv = l_const2 * (temp_1p*l_Irr + temp_2p*l_Jrr + l_rhos * (temp_3p * l_Irr - temp_4p * l_Jrr) + &
        &           drhosv * (temp_3 * l_Irr - temp_4 * l_Jrr))     

        diff_v  = xx**(2.0e0_DP * lla) * (l_aa + l_Baa) - 2.0e0_DP * xx**(lla + llr) * (l_ar + l_Bar) + &
        &       xx**(2.0e0_DP * llr) * (l_rr + l_Brr) 
        diff_dv = xx**(2.0e0_DP * lla) * (l_daa + l_dBaa) - 2.0e0_DP * xx**(lla + llr) * (l_dar + l_dBar) + &
        &       xx**(2.0e0_DP * llr) * (l_drr + l_dBrr)
        
        ddiff_vv  = xx**(2.0e0_DP * lla) * (dl_aav + dl_Baav) - 2.0e0_DP * xx**(lla + llr) * (dl_arv + dl_Barv) + &
        &       xx**(2.0e0_DP * llr) * (dl_rrv + dl_Brrv) 
        ddiff_dvv = xx**(2.0e0_DP * lla) * (dl_daav + dl_dBaav) - 2.0e0_DP*xx**(lla + llr) * (dl_darv + dl_dBarv) + &
        &       xx**(2.0e0_DP * llr) * (dl_drrv + dl_dBrrv)
   
        as2grad_out = 0.5e0_DP * leps * cc**2.0e0_DP *(ddiff_vv * dkhs + diff_v * ddkhsv  + &
        &   ddiff_dvv * khs + diff_dv * dkhsv)

        return
    end function das2gradv
!***************************************************************************************************
end module Chain_mod
!***************************************************************************************************
!***************************************************************************************************
!   SAFT Module to:
!       1. Calculate A assoc
!
!       Analtical derivatives added March 2017
!
!***************************************************************************************************
!
!***************************************************************************************************
!***************************************************************************************************
!   SAFT Module to:
!       1. Calculate A assoc
!
!       Analtical derivatives added March 2017
!
!***************************************************************************************************
!
!***************************************************************************************************
module Assoc_mod
!***************************************************************************************************
!Modules
!=======    
    use Types_mod           ! Definitions of types and double precision
    use Global_mod          ! Important global parameters 
!***************************************************************************************************
    implicit none
!***************************************************************************************************
    contains
!***************************************************************************************************    
    function A_assoc() result(a_res)      ! [2] 21
        implicit none
        
        integer             ::  i_a1,i_a2,i_a3
        real(kind=DP)       ::  a_res
        
        a_res = 0.0e0_DP

        if(assoc_switch) then        
            call Mass_action()
      
            do i_a1=1, nctypes
            do i_a2=1, nstypes
            
            if(Comp_array(i_a1)%comp(i_a2)>0) then      !added Jan 2017   
               
            do i_a3=1, Seg_array(i_a2)%nassoc
                  a_res = a_res + Comp_array(i_a1)%xi * Comp_array(i_a1)%comp(i_a2) *       &
                  & Seg_array(i_a2)%nassoctyp(i_a3) * (dlog(x(i_a1, i_a2,i_a3)) + (1.0e0_DP  &
                  & -x(i_a1,i_a2,i_a3))/2.0e0_DP)            
            end do

            end if
            
            end do
            end do
        end if   

        return
    end function A_assoc
!***************************************************************************************************
    function A_assoc_dn() result(adn_res)      ! [2] 21
        implicit none
        
        integer             ::  i_a1,i_a2,i_a3
        real(kind=DP)       ::  adn_res
  
        adn_res = 0.0e0_DP

        if(assoc_switch) then        
            call Mass_action()
            call Mass_action_dn()
      
            do i_a1=1, nctypes
            do i_a2=1, nstypes
            
            if(Comp_array(i_a1)%comp(i_a2)>0) then      !added Jan 2017   
               
            do i_a3=1, Seg_array(i_a2)%nassoc                  
                  adn_res = adn_res +   &
                  & nsum * Comp_array(i_a1)%xi * Comp_array(i_a1)%comp(i_a2) * Seg_array(i_a2)%nassoctyp(i_a3)   * &
                  & (dxn(i_a1, i_a2,i_a3)/x(i_a1,i_a2,i_a3) - dxn(i_a1, i_a2,i_a3)/2.0e0_DP )                    + &
                  & nsum * Comp_array(i_a1)%dxin * Comp_array(i_a1)%comp(i_a2) * Seg_array(i_a2)%nassoctyp(i_a3) * &
                  & (dlog(x(i_a1, i_a2,i_a3)) + (1.0e0_DP - x(i_a1,i_a2,i_a3))/2.0e0_DP)                          + &
                  & Comp_array(i_a1)%xi * Comp_array(i_a1)%comp(i_a2) * Seg_array(i_a2)%nassoctyp(i_a3)          * &
                  & (dlog(x(i_a1, i_a2,i_a3)) + (1.0e0_DP - x(i_a1,i_a2,i_a3))/2.0e0_DP) 
            end do

            end if
            
            end do
            end do
        end if   
     
        return
    end function A_assoc_dn
!***************************************************************************************************
    function A_assoc_dv() result(adv_res)      ! [2] 21
        implicit none
        
        integer             ::  i_a1,i_a2,i_a3
        real(kind=DP)       ::  adv_res
  
        adv_res = 0.0e0_DP

        if(assoc_switch) then               
            call Mass_action()
            call Mass_action_dv()
      
            do i_a1=1, nctypes
            do i_a2=1, nstypes
          
            if(Comp_array(i_a1)%comp(i_a2)>0) then      !added Jan 2017   
               
            do i_a3=1, Seg_array(i_a2)%nassoc                  
                  adv_res = adv_res + Comp_array(i_a1)%xi * Comp_array(i_a1)%comp(i_a2)             * &
                  & Seg_array(i_a2)%nassoctyp(i_a3) * (dxv(i_a1, i_a2,i_a3) / x(i_a1, i_a2,i_a3)    - &
                  & dxv(i_a1,i_a2,i_a3)/2.0e0_DP)     
            end do

            end if
            
            end do
            end do
        end if   
     
        adv_res = -adv_res * KB * NA * t
     
        return
    end function A_assoc_dv
!***************************************************************************************************
    subroutine mass_action() !BRUTE FORCE ITERATION AT PRESENT
        implicit none
        
        integer                         ::  mac1, mac2, mac3, mbc1, mbc2, mbc3, n_it
        real(kind=DP),parameter         ::  CRITERIA=1.0e-15_DP   
        logical                         ::  converged
        real(kind=DP),allocatable       ::  x_old(:,:,:)
        real(kind=DP)                   ::  xsum
        
        if(allocated(x))        deallocate(x)
        if(allocated(delx))     deallocate(delx)
        if(allocated(x_old))    deallocate(x_old)   
        
        allocate(x(nctypes, nstypes, 3), delx(nctypes, nstypes, 3, nctypes, nstypes, 3), &
        &   x_old(nctypes, nstypes, 3))

        call delta_x( )                  
        
        x = 0.5e0_DP  
  
        n_it=0
        converged=.false.
      
        do while(.not.converged)      
            x_old=x     
            converged=.true.
            
            do mac1=1, nctypes
            do mac2=1, nstypes
            
            if(Comp_array(mac1)%comp(mac2)>0) then    !added Jan 2017   
    
            do mac3=1, Seg_array(mac2)%nassoc
        
                xsum = 1.0e0_DP
                
                do mbc1=1, nctypes
                do mbc2=1, nstypes
                do mbc3=1, Seg_array(mbc2)%nassoc
                
                    xsum = xsum + rho * Comp_array(mbc1)%xi * Comp_array(mbc1)%comp(mbc2) * &
                    &   Seg_array(mbc2)%nassoctyp(mbc3) * x_old(mbc1, mbc2, mbc3) *         &
                    &   delx( mac1, mac2, mac3, mbc1, mbc2, mbc3 )                
                end do
                end do
                end do

                x(mac1, mac2, mac3) = 1.0e0_DP / xsum
                if(dabs(x(mac1, mac2, mac3) - x_old(mac1, mac2, mac3)) > CRITERIA) converged=.false.                    
            end do
            
            end if   
                     
            end do
            end do
   
            n_it = n_it + 1
        end do  

        return
    end subroutine Mass_action
!***************************************************************************************************
    subroutine mass_action_dn() 
        implicit none
        
        integer                         ::  mac1, mac2, mac3, mbc1, mbc2, mbc3
        real(kind=DP),allocatable       ::  dxn_old(:,:,:)
        real(kind=DP)                   ::  xpart2, xpart2p  
        real(kind=DP),parameter         ::  CRITERIA=1.0e-20_DP 
        logical                         ::  CONVERGED

        if(allocated(dxn))          deallocate(dxn)
        if(allocated(dxn_old))      deallocate(dxn_old)
        if(allocated(ddelxn))       deallocate(ddelxn)  
        
        allocate(ddelxn(nctypes, nstypes, 3, nctypes, nstypes, 3), dxn(nctypes, nstypes, 3), dxn_old(nctypes, nstypes, 3))

        call ddelta_xn( )                  

        xpart2  = 0.0e0_DP
        xpart2p = 0.0e0_DP    
           
        dxn(:,:,:) = 1.0e0_DP   
        CONVERGED=.false.   
           
        do while(.not.CONVERGED)
           
            dxn_old=dxn  
            CONVERGED=.true.

        do mac1=1, nctypes
        do mac2=1, nstypes
        
        if(Comp_array(mac1)%comp(mac2)>0) then    !added Jan 2017   

        do mac3=1, Seg_array(mac2)%nassoc
            
            xpart2  = 0.0e0_DP
            xpart2p = 0.0e0_DP
            
            do mbc1=1, nctypes
            do mbc2=1, nstypes
            do mbc3=1, Seg_array(mbc2)%nassoc
                xpart2 = xpart2 + Comp_array(mbc1)%xi * Comp_array(mbc1)%comp(mbc2) * &
                &   Seg_array(mbc2)%nassoctyp(mbc3) * x(mbc1, mbc2, mbc3) * delx( mac1, mac2, mac3, mbc1, mbc2, mbc3 )    
                
                xpart2p = xpart2p + Comp_array(mbc1)%comp(mbc2) * Seg_array(mbc2)%nassoctyp(mbc3) * &
                & (Comp_array(mbc1)%dxin * x(mbc1, mbc2, mbc3) *  delx( mac1, mac2, mac3, mbc1, mbc2, mbc3) + &
                &  Comp_array(mbc1)%xi * dxn(mbc1, mbc2, mbc3) *  delx( mac1, mac2, mac3, mbc1, mbc2, mbc3) + &
                &  Comp_array(mbc1)%xi * x(mbc1, mbc2, mbc3) *  ddelxn( mac1, mac2, mac3, mbc1, mbc2, mbc3))                         
            end do
            end do
            end do

            dxn(mac1, mac2, mac3) = (-1.0e0_DP * x(mac1, mac2, mac3) * drhon * (xpart2 +  Comp_array(mac1)%xi           * &
            &   Comp_array(mac1)%comp(mac2)                                                                             * &
            &   Seg_array(mac2)%nassoctyp(mac3) * x(mac1, mac2, mac3) * delx( mac1, mac2, mac3, mac1, mac2, mac3 ))     - &
            &   x(mac1, mac2, mac3) * rho *  (xpart2p + Comp_array(mac1)%comp(mac2) * Seg_array(mac2)%nassoctyp(mac3)   * &
            &   (Comp_array(mac1)%xi * x(mac1, mac2, mac3) * ddelxn( mac1, mac2, mac3, mac1, mac2, mac3 )               + &
            &   Comp_array(mac1)%dxin * x(mac1, mac2, mac3) * delx( mac1, mac2, mac3, mac1, mac2, mac3 )) ))            / &
            &   (1.0e0_DP + rho * (xpart2 +  Comp_array(mac1)%xi * Comp_array(mac1)%comp(mac2)                          * &
            &   Seg_array(mac2)%nassoctyp(mac3) * x(mac1, mac2, mac3) * delx( mac1, mac2, mac3, mac1, mac2, mac3 ))     + &
            &   x(mac1, mac2, mac3) * rho * Comp_array(mac1)%comp(mac2) * Seg_array(mac2)%nassoctyp(mac3)               * &
            &   Comp_array(mac1)%xi * delx( mac1, mac2, mac3, mac1, mac2, mac3 )        ) 

            if(dabs((dxn(mac1, mac2, mac3) - dxn_old(mac1, mac2, mac3))/dxn(mac1, mac2, mac3)) > CRITERIA) converged=.false.               
        end do
        
        end if   
                 
        end do
        end do
     
        end do

        return
    end subroutine mass_action_dn    
!***************************************************************************************************
    subroutine mass_action_dv() !BRUTE FORCE ITERATION AT PRESENT
        implicit none
        
        integer                         ::  mac1, mac2, mac3, mbc1, mbc2, mbc3, n_it
        real(kind=DP),parameter         ::  CRITERIA=1.0e-15_DP   
        logical                         ::  converged
        real(kind=DP),allocatable       ::  dxv_old(:,:,:)
        real(kind=DP)                   ::  xsum, xsum2
        
        if(allocated(dxv))          deallocate(dxv)
        if(allocated(dxv_old))      deallocate(dxv_old)
        if(allocated(ddelxv))       deallocate(ddelxv)  
        
        allocate(ddelxv(nctypes, nstypes, 3, nctypes, nstypes, 3), dxv(nctypes, nstypes, 3), dxv_old(nctypes, nstypes, 3))

        call ddelta_xv( )                  
        
        dxv = 1.0e0_DP  
  
        n_it=0
        converged=.false.
       
        do while(.not.converged)
           
            dxv_old=dxv     
            converged=.true.
            
            do mac1=1, nctypes
            do mac2=1, nstypes
            
            if(Comp_array(mac1)%comp(mac2)>0) then    !added Jan 2017   
    
            do mac3=1, Seg_array(mac2)%nassoc
            
                xsum = 1.0e0_DP
                xsum2 = 0.0e0_DP
                
                do mbc1=1, nctypes
                do mbc2=1, nstypes
                do mbc3=1, Seg_array(mbc2)%nassoc
                    xsum = xsum + rho * Comp_array(mbc1)%xi * Comp_array(mbc1)%comp(mbc2)   * &
                    &   Seg_array(mbc2)%nassoctyp(mbc3) * x(mbc1, mbc2, mbc3)               * &
                    &   delx( mac1, mac2, mac3, mbc1, mbc2, mbc3 )      
                    
                    xsum2 = xsum2 - Comp_array(mbc1)%xi * Comp_array(mbc1)%comp(mbc2)                   * &
                    &   Seg_array(mbc2)%nassoctyp(mbc3)                                                 * &
                    &   (drhov * x(mbc1, mbc2, mbc3) * delx( mac1, mac2, mac3, mbc1, mbc2, mbc3 )       + &
                    &    rho * dxv(mbc1, mbc2, mbc3) * delx( mac1, mac2, mac3, mbc1, mbc2, mbc3 )       + &
                    &    rho * x(mbc1, mbc2, mbc3) * ddelxv( mac1, mac2, mac3, mbc1, mbc2, mbc3 ))
                end do
                end do
                end do

                dxv(mac1, mac2, mac3) = xsum2 / xsum**2.0e0_DP
                
                if(dabs(dxv(mac1, mac2, mac3) - dxv_old(mac1, mac2, mac3)) > CRITERIA) converged=.false.                    
            end do
           
            end if   
                     
            end do
            end do
  
            n_it = n_it + 1
        end do  

        return
    end subroutine Mass_action_dv
!***************************************************************************************************
    subroutine delta_x( )
        implicit none
        
        integer          ::  i_dx1, i_dx2, i_dx3, i_dx4, i_dx5, i_dx6, dxp, dxq
        real(kind=DP)    ::  f_x(1:nstypes, 1:3, 1:nstypes, 1:3)
        real(kind=DP)    ::  i_x(1:nctypes, 1:nstypes, 1:3, 1:nctypes, 1:nstypes, 1:3)
        
        do i_dx1=1, nstypes
        do i_dx2=1, Seg_array(i_dx1)%nassoc
            do i_dx3=1, nstypes
            do i_dx4=1, Seg_array(i_dx3)%nassoc               
                f_x(i_dx1, i_dx2, i_dx3, i_dx4) = dexp( ehb(i_dx1, i_dx2, i_dx3, i_dx4) / T) - 1.0e0_DP                                
            end do
            end do
        end do
        end do

        i_x = 0.0e0_DP
          
        do i_dx1=1, nctypes
        do i_dx2=1, nstypes
        do i_dx3=1, Seg_array(i_dx2)%nassoc
            do i_dx4=1, nctypes
            do i_dx5=1, nstypes
            do i_dx6=1, Seg_array(i_dx5)%nassoc
           
                do dxp=0, 10
                do dxq=0, 10 - dxp              
                    i_x(i_dx1, i_dx2, i_dx3, i_dx4, i_dx5, i_dx6) =         &
                    &   i_x(i_dx1, i_dx2, i_dx3, i_dx4, i_dx5, i_dx6)   +   &

                    &   (c_x(dxp, dxq) * (rhos * NA * sigx3) ** dxp      *   &
                    &   (t / epschij(i_dx1, i_dx4))**dxq)     
                end do
                end do
                
                delx(i_dx1, i_dx2, i_dx3, i_dx4, i_dx5, i_dx6) =  f_x(i_dx2, i_dx3, i_dx5, i_dx6) * &
                &   i_x(i_dx1, i_dx2, i_dx3, i_dx4, i_dx5, i_dx6) * khb(i_dx2, i_dx3, i_dx5, i_dx6)
                          
                !SHOULDN'T HAPPEN, BUT.... CAUSED BY ISSUE WITH i_x               
                if(delx(i_dx1, i_dx2, i_dx3, i_dx4, i_dx5, i_dx6) < 0.0e0_DP) then             
                    delx(i_dx1, i_dx2, i_dx3, i_dx4, i_dx5, i_dx6)=0.0e0_DP
                end if

            end do
            end do
            end do
        end do
        end do
        end do 
!print*, rhos, drhosv,i_x

        return
    end subroutine          
!***************************************************************************************************
    subroutine ddelta_xn( )
        implicit none
        
        integer          ::  i_dx1, i_dx2, i_dx3, i_dx4, i_dx5, i_dx6, dxp, dxq
        real(kind=DP)    ::  f_x(1:nstypes, 1:3, 1:nstypes, 1:3)
        real(kind=DP)    ::  di_xn(1:nctypes, 1:nstypes, 1:3, 1:nctypes, 1:nstypes, 1:3)
        
        do i_dx1=1, nstypes
        do i_dx2=1, Seg_array(i_dx1)%nassoc
            do i_dx3=1, nstypes
            do i_dx4=1, Seg_array(i_dx3)%nassoc               
                f_x(i_dx1, i_dx2, i_dx3, i_dx4) = dexp( ehb(i_dx1, i_dx2, i_dx3, i_dx4) / T) - 1.0e0_DP                                
            end do
            end do
        end do
        end do

        di_xn = 0.0e0_DP
          
        do i_dx1=1, nctypes
        do i_dx2=1, nstypes
        do i_dx3=1, Seg_array(i_dx2)%nassoc
            do i_dx4=1, nctypes
            do i_dx5=1, nstypes
            do i_dx6=1, Seg_array(i_dx5)%nassoc
           
                do dxp=0, 10
                do dxq=0, 10 - dxp              
                    di_xn(i_dx1, i_dx2, i_dx3, i_dx4, i_dx5, i_dx6) =         &
                    &   di_xn(i_dx1, i_dx2, i_dx3, i_dx4, i_dx5, i_dx6)   +   &

                    &   (c_x(dxp, dxq) * &                
                                       
                    &   (rhos*dsigx3n + drhosn*sigx3) * dxp * NA *  &
                    &   (rhos * NA * sigx3) ** (dxp-1)  *   &                 
                                        
                    &   (t / epschij(i_dx1, i_dx4))**dxq)  
                                      
                                                                   
                end do
                end do

                ddelxn(i_dx1, i_dx2, i_dx3, i_dx4, i_dx5, i_dx6) =  f_x(i_dx2, i_dx3, i_dx5, i_dx6) * &
                &   di_xn(i_dx1, i_dx2, i_dx3, i_dx4, i_dx5, i_dx6) * khb(i_dx2, i_dx3, i_dx5, i_dx6)
            end do
            end do
            end do
        end do
        end do
        end do 

        return
    end subroutine ddelta_xn        
!***************************************************************************************************
    subroutine ddelta_xv( )
        implicit none
        
        integer          ::  i_dx1, i_dx2, i_dx3, i_dx4, i_dx5, i_dx6, dxp, dxq
        real(kind=DP)    ::  f_x(1:nstypes, 1:3, 1:nstypes, 1:3)
        real(kind=DP)    ::  di_xv(1:nctypes, 1:nstypes, 1:3, 1:nctypes, 1:nstypes, 1:3)
        
        do i_dx1=1, nstypes
        do i_dx2=1, Seg_array(i_dx1)%nassoc
            do i_dx3=1, nstypes
            do i_dx4=1, Seg_array(i_dx3)%nassoc               
                f_x(i_dx1, i_dx2, i_dx3, i_dx4) = dexp( ehb(i_dx1, i_dx2, i_dx3, i_dx4) / T) - 1.0e0_DP                                
            end do
            end do
        end do
        end do

        di_xv = 0.0e0_DP
          
        do i_dx1=1, nctypes
        do i_dx2=1, nstypes
        do i_dx3=1, Seg_array(i_dx2)%nassoc
            do i_dx4=1, nctypes
            do i_dx5=1, nstypes
            do i_dx6=1, Seg_array(i_dx5)%nassoc
           
                do dxp=0, 10
                do dxq=0, 10 - dxp              
                    di_xv(i_dx1, i_dx2, i_dx3, i_dx4, i_dx5, i_dx6) =         &
                    &   di_xv(i_dx1, i_dx2, i_dx3, i_dx4, i_dx5, i_dx6)   +   &

                    &   (c_x(dxp, dxq)                                    *   &
                    &   real(dxp) * drhosv * sigx3                        *   &
                    &   (rhos * NA * sigx3) ** (dxp-1)            *   &
                    &   (t / epschij(i_dx1, i_dx4))**dxq)                                                 
                end do
                end do
                
                ddelxv(i_dx1, i_dx2, i_dx3, i_dx4, i_dx5, i_dx6) =  f_x(i_dx2, i_dx3, i_dx5, i_dx6) * &
                &   di_xv(i_dx1, i_dx2, i_dx3, i_dx4, i_dx5, i_dx6) * khb(i_dx2, i_dx3, i_dx5, i_dx6) * NA
            end do
            end do
            end do
        end do
        end do
        end do 

        return
    end subroutine          
!***************************************************************************************************

end module Assoc_mod
!***************************************************************************************************
!***************************************************************************************************
!***************************************************************************************************
!   SAFT Module to:
!       1. Calculate A ion and A born
!
!       Based on:
!           D.K. Eriksen, G. Lazarou, A. Galindo, G. Jackson, C.S. Adjiman, A.J. Haslam,           
!           "Development of intermolecular potential models for electrolyte
!           solutions using an electrolyte SAFT-VR Mie equation of state", 
!           Mol. Phys. 114(18), 2016 2724-2749
!           
!       Analytical derivatives added March 2017
!
!***************************************************************************************************
!
!***************************************************************************************************
module Ion_mod
!***************************************************************************************************
!Modules
!=======    
    use Types_mod           ! Definitions of types and double precision
    use Global_mod          ! Important global parameters     
!***************************************************************************************************
    implicit none
!***************************************************************************************************
    contains
!***************************************************************************************************
    function A_ion( ) result(a_res)     !AJ Haslam et al Mol Phys 2016            
        implicit none
    
        real(kind=DP)       ::  a_res
        real(kind=DP)       ::  dielec

        a_res=0.0e0_DP

        if(ion_switch) then   
            dielec = Diel()                 
            a_res  = Aion( dielec ) + Aborn( dielec )
        end if

        return
    end function A_ion
!***************************************************************************************************
    function A_ion_only( ) result(a_res)     !As part of pressure calculation              
        implicit none
    
        real(kind=DP)       ::  a_res
        real(kind=DP)       ::  dielec

        a_res=0.0e0_DP

        if(ion_switch) then
            dielec = Diel()          
            a_res = Aion( dielec )
        end if

        return
    end function A_ion_only  
!***************************************************************************************************
    function A_born_only( ) result(a_res)     !As part of mu inf calculation              
        implicit none
    
        real(kind=DP)       ::  a_res
        real(kind=DP)       ::  dielec

        a_res=0.0e0_DP

        if(ion_switch) then 
            dielec = Diel()          
            a_res = Aborn( dielec )
        end if

        return
    end function A_born_only  
!*************************************************************************************************** 
    function Diel( ) result( diel_res )
        implicit none

        real(kind=DP)               ::  diel_res
        real(kind=DP)               ::  x_solv, ddiel 
        integer                     ::  i_diel1, i_diel2

        x_solv=0.0e0_DP
        ddiel=0.0e0_DP
        
        do i_diel1=1, nctypes
            if(Comp_array(i_diel1)%solv) x_solv = x_solv + Comp_array(i_diel1)%xi
        end do
        
        do i_diel1=1, nctypes
        do i_diel2=i_diel1, nctypes
            if(( Comp_array(i_diel1)%solv ).and.( Comp_array(i_diel2)%solv )) then           
                ddiel = ddiel + ( Comp_array(i_diel1)%xi / x_solv ) * ( Comp_array(i_diel2)%xi / x_solv ) * &
                & ( Comp_array(i_diel1)%dv * ( Comp_array(i_diel1)%dt / t - 1.0e0_DP ) + &
                & Comp_array(i_diel2)%dv * ( Comp_array(i_diel2)%dt / t - 1.0e0_DP ) ) / 2.0e0_DP             
            end if 
        end do
        end do

        diel_res= 1.0e0_DP + rho * ddiel * x_solv 
     
        return
    end function Diel
!***************************************************************************************************
    function dDielN( ) result( diel_res )
        implicit none

        real(kind=DP)               ::  diel_res
        real(kind=DP)               ::  x_solv, ddiel 
        real(kind=DP)               ::  x_solv_dn, ddiel_dn 
        integer                     ::  i_diel1, i_diel2

        x_solv=0.0e0_DP
        ddiel=0.0e0_DP
        
        x_solv_dn=0.0e0_DP
        ddiel_dn=0.0e0_DP
        
        do i_diel1=1, nctypes
            if(Comp_array(i_diel1)%solv) then
                x_solv = x_solv + Comp_array(i_diel1)%xi
                x_solv_dn = x_solv_dn + Comp_array(i_diel1)%dxin
            end if
        end do
        
        do i_diel1=1, nctypes
        do i_diel2=i_diel1, nctypes
            if(( Comp_array(i_diel1)%solv ).and.( Comp_array(i_diel2)%solv )) then           
                ddiel = ddiel + ( Comp_array(i_diel1)%xi / x_solv ) * ( Comp_array(i_diel2)%xi / x_solv ) * &
                & ( Comp_array(i_diel1)%dv * ( Comp_array(i_diel1)%dt / t - 1.0e0_DP ) + &
                & Comp_array(i_diel2)%dv * ( Comp_array(i_diel2)%dt / t - 1.0e0_DP ) ) / 2.0e0_DP   
                
                ddiel_dn = ddiel_dn + ( Comp_array(i_diel1)%xi / x_solv ) * ( Comp_array(i_diel2)%xi / x_solv ) * &
                
                &   ((x_solv*(Comp_array(i_diel1)%xi*Comp_array(i_diel2)%dxin + Comp_array(i_diel1)%dxin*Comp_array(i_diel2)%xi) - &
                & 2.0e0_DP*Comp_array(i_diel1)%xi*Comp_array(i_diel2)%xi*x_solv_dn ) / x_solv**3.0e0_DP) * &
                
                & ( Comp_array(i_diel1)%dv * ( Comp_array(i_diel1)%dt / t - 1.0e0_DP ) + &
                & Comp_array(i_diel2)%dv * ( Comp_array(i_diel2)%dt / t - 1.0e0_DP ) ) / 2.0e0_DP 
            end if 
        end do
        end do

        diel_res= x_solv*rho*ddiel_dn + x_solv_dn*rho*ddiel_dn + x_solv*drhon*ddiel + x_solv_dn*drhon*ddiel  
      
        return
    end function dDielN
!***************************************************************************************************
    function dDielv( ) result( diel_res )
        implicit none

        real(kind=DP)               ::  diel_res
        real(kind=DP)               ::  x_solv, ddiel 
        integer                     ::  i_diel1, i_diel2

        x_solv=0.0e0_DP
        ddiel=0.0e0_DP
        
        do i_diel1=1, nctypes
            if(Comp_array(i_diel1)%solv) x_solv = x_solv + Comp_array(i_diel1)%xi
        end do
        
        do i_diel1=1, nctypes
        do i_diel2=i_diel1, nctypes
            if(( Comp_array(i_diel1)%solv ).and.( Comp_array(i_diel2)%solv )) then           
                ddiel = ddiel + ( Comp_array(i_diel1)%xi / x_solv ) * ( Comp_array(i_diel2)%xi / x_solv ) * &
                & ( Comp_array(i_diel1)%dv * ( Comp_array(i_diel1)%dt / t - 1.0e0_DP ) + &
                & Comp_array(i_diel2)%dv * ( Comp_array(i_diel2)%dt / t - 1.0e0_DP ) ) / 2.0e0_DP             
            end if 
        end do
        end do

        diel_res= 1.0e0_DP + drhov * ddiel * x_solv 
   
        return
    end function dDielv
!***************************************************************************************************
    function Aborn( l_dielec ) result( born_res )
        implicit none
        
        real(kind=DP), intent(in)   ::  l_dielec
        real(kind=DP)   ::  born_res
        real(kind=DP)   ::  bsum, born_const, born_diel
        integer         ::  i_born1, i_born2, i_born3

        bsum=0.0e0_DP
        
        born_diel = 1.0e0_DP - 1.0e0_DP / l_dielec
        born_const = -1.0e0 * QE**2.0e0_DP / 4.0e0_DP / PI / E0
     
        do i_born1=1, nctypes
            do i_born2=1, nstypes
                if(Seg_array(i_born2)%charged) then
                    do i_born3=1, Comp_array(i_born1)%comp(i_born2)
!This changed - xi * NA to  nm                    
                        bsum = bsum + Comp_array(i_born1)%nm * Seg_array(i_born2)%q**2.0e0_DP /     &                      
                        &       Seg_array(i_born2)%sigb                     
                    end do
                end if
            end do           
        end do
!This changed - NA to  nm sum      
        born_res = born_const * born_diel * bsum
        born_res = born_res / sum(Comp_array(:)%nm) / KB / t         

        return
    end function Aborn
!***************************************************************************************************
    function A_born_dv(  ) result( born_res )      
        implicit none
        
        real(kind=DP)   ::  born_res
        real(kind=DP)   ::  bsum, born_const, born_diel, x_solv, ddiel
        integer         ::  i_born1, i_born2, i_born3

        bsum=0.0e0_DP
        
        if(ion_switch) then 
            do i_born1=1, nctypes
                do i_born2=1, nstypes
                    if(Seg_array(i_born2)%charged) then
                        do i_born3=1, Comp_array(i_born1)%comp(i_born2)
                            bsum = bsum + Comp_array(i_born1)%xi * NA * Seg_array(i_born2)%q**2.0e0_DP /     &                      
                            &       Seg_array(i_born2)%sigb
                        end do
                    end if
                end do
            end do
            
            x_solv=0.0e0_DP
            ddiel=0.0e0_DP
           
            do i_born1=1, nctypes
                if(Comp_array(i_born1)%solv) x_solv = x_solv + Comp_array(i_born1)%xi
            end do
           
            do i_born1=1, nctypes
            do i_born2=i_born1, nctypes
                if(( Comp_array(i_born1)%solv ).and.( Comp_array(i_born2)%solv )) then           
                    ddiel = ddiel + ( Comp_array(i_born1)%xi / x_solv ) * ( Comp_array(i_born2)%xi / x_solv ) * &
                    & ( Comp_array(i_born1)%dv * ( Comp_array(i_born1)%dt / t - 1.0e0_DP ) + &
                    & Comp_array(i_born2)%dv * ( Comp_array(i_born2)%dt / t - 1.0e0_DP ) ) / 2.0e0_DP
                end if 
            end do
            end do
        end if
        
        born_diel = x_solv * ddiel / (v + x_solv * ddiel)**2.0e0_DP
        born_const = QE**2.0e0_DP / 4.0e0_DP / PI / E0
        
        born_res = born_const * born_diel * bsum
        born_res = born_res / NA / KB / t         

        return
    end function A_born_dv
!***************************************************************************************************
    function A_born_dn( i_ion ) result( born_res )      
        implicit none
        
        integer, intent(in) ::  i_ion
        real(kind=DP)       ::  born_res
        real(kind=DP)       ::  diel_l, born_const, diel_l_dn
        integer     ::  i_born1, i_born2, i_born3
        
        born_res = 0.0e0_DP
        
        if(ion_switch) then 
            
            born_const = -1.0e0 * QE**2.0e0_DP / 4.0e0_DP / PI / E0
            
            if((i_ion==properties%cation).or.(i_ion==properties%anion)) then
                diel_l = Diel( )             
                born_res = born_const * (1.0e0_DP - 1.0e0_DP / diel_l)      *   &
                &   Seg_array(i_ion)%q**2.0e0_DP / Seg_array(i_ion)%sigb
            else            
                diel_l = Diel( ) 
                diel_l_dn = dDielN()
                                
                do i_born1=1, nctypes
                do i_born2=1, nstypes
                if(Seg_array(i_born2)%charged) then
                do i_born3=1, Comp_array(i_born1)%comp(i_born2)
                born_res = born_const * (-diel_l_dn/diel_l**2.0e0_DP)      *   &
                      &   Seg_array(i_born2)%q**2.0e0_DP / Seg_array(i_born2)%sigb                          
                end do
                end if
                end do           
                end do            
            end if

        end if
                                                                        
        born_res = born_res / KB / t / sum(Comp_array(:)%nm)  
    
        return
    end function A_born_dn
!***************************************************************************************************
    function Aion( l_dielec ) result( ion_res )        
        implicit none

        real(kind=DP)               ::  ion_res
        real(kind=DP), intent(in)   ::  l_dielec
        real(kind=DP)               ::  c_shield, shield
        
        real(kind=DP), parameter    ::  SHIELD_CONV=1.0e-12_DP
        
        real(kind=DP)               ::  qsum_l, del_l, omega_l, pn_l, shield_l
        real(kind=DP)               ::  umsa, qtemp
        
        integer                     ::  i_ion1, i_ion2, i_ion3
        integer                     ::  it_tot
!Here also changed all to nm        
        c_shield = QE**2.0e0_DP / (4.0e0_DP * l_dielec * E0 *  KB * t * v)        
        qsum_l = 0.0e0_DP
        
        do i_ion1=1, nctypes   
            do i_ion2=1, nstypes
                if( Seg_array(i_ion2)%charged ) then
                    do i_ion3=1, Comp_array(i_ion1)%comp(i_ion2)
                        qsum_l = qsum_l + Comp_array(i_ion1)%nm * Seg_array(i_ion2)%q**2.0e0_DP                        
                    end do
                end if 
            end do
        end do
        
        !initial shield guess
        shield = 0.5e0_DP * dsqrt( c_shield * qsum_l )

        !set delta  - this is constant for iteration       
        del_l = 0.0e0_DP

        do i_ion1=1, nctypes
            do i_ion2=1, nstypes
                if(Seg_array(i_ion2)%charged) then
                    do i_ion3=1, Comp_array(i_ion1)%comp(i_ion2)
                        del_l = del_l + Comp_array(i_ion1)%nm * Seg_array(i_ion2)%sig**3.0e0_DP                   
                    end do
                end if 
            end do
        end do
        
        del_l = 1.0e0_DP - PI / 6.0e0_DP / v * del_l            
     
        !iterate parameters to find shield
        it_tot = 0       
        
        do             
            omega_l = 0.0e0_DP
            
            do i_ion1=1, nctypes
                do i_ion2=1, nstypes
                    if(Seg_array(i_ion2)%charged) then
                        do i_ion3=1, Comp_array(i_ion1)%comp(i_ion2)
                            omega_l = omega_l + Comp_array(i_ion1)%nm * Seg_array(i_ion2)%sig**3.0e0_DP / &
                                &   (1.0e0_DP + shield * Seg_array(i_ion2)%sig)
                        end do
                    end if 
                end do
            end do

            omega_l = 1.0e0_DP + PI / v / 2.e0_DP / del_l * omega_l

            pn_l = 0.0e0_DP
            
            do i_ion1=1, nctypes
                do i_ion2=1, nstypes
                    if(Seg_array(i_ion2)%charged) then
                        do i_ion3=1, Comp_array(i_ion1)%comp(i_ion2)
!added *QE
                            pn_l = pn_l + Comp_array(i_ion1)%nm * Seg_array(i_ion2)%sig *       &
                            &   Seg_array(i_ion2)%q*QE / (1.0e0_DP + shield * Seg_array(i_ion2)%sig )
                           
                        end do
                    end if 
                end do
            end do
           
            pn_l = 1.0e0_DP / omega_l / v * pn_l 

            shield_l = 0.0e0_DP
            
            do i_ion1=1, nctypes
                do i_ion2=1, nstypes
                    if(Seg_array(i_ion2)%charged) then
                        do i_ion3=1, Comp_array(i_ion1)%comp(i_ion2)
                            qtemp = ( Seg_array(i_ion2)%q * QE - (PI / 2.0e0_DP / del_l) *  &
                            &       Seg_array(i_ion2)%sig**2.0e0_DP * pn_l ) / (1.0e0_DP + shield * Seg_array(i_ion2)%sig)                     
                            shield_l = shield_l + Comp_array(i_ion1)%nm * qtemp**2.0e0_DP                                
                        end do
                    end if 
                end do
            end do
    
            shield_l = dsqrt( PI  / 4.0e0_DP / PI / E0 / l_dielec / KB / T / v * shield_l )
        
            it_tot = it_tot + 1
                    
            if(dabs( shield - shield_l ) / dabs( shield_l )  < SHIELD_CONV) then
                shield = shield_l
                exit
            else
                shield = shield_l
            end if    
        end do 

        umsa = 0.0e0_DP

        do i_ion1=1, nctypes
            do i_ion2=1, nstypes
                if(Seg_array(i_ion2)%charged) then
                    do i_ion3=1, Comp_array(i_ion1)%comp(i_ion2)
                        umsa = umsa + (Comp_array(i_ion1)%nm * Seg_array(i_ion2)%q**2.0e0_DP* QE**2.0e0_DP )  /   &
                        &   ( 1.0e0_DP + shield * Seg_array(i_ion2)%sig )  
                    end do
                end if 
            end do
        end do 

        umsa = -1.0e0_DP * v / 4.0e0_DP / PI / E0 / l_dielec * (shield / v * umsa +  &
        &       PI / 2.0e0_DP / del_l * omega_l * pn_l**2.0e0_DP )

        ion_res = umsa + shield**3.0e0_DP * KB * v * t / 3.0e0_DP / PI

        !Make it in equivalent units...
        ion_res = ion_res / sum(Comp_array(:)%nm) / KB / t
                 
        return
    end function Aion
!***************************************************************************************************  
    function A_ion_dn( i_part ) result(a_res)
        implicit none
        
        integer, intent(in)     ::  i_part
        real(kind=DP)           ::  a_res

        real(kind=DP)           ::  l_dielec
        real(kind=DP)               ::  c_shield, shield
        
        real(kind=DP), parameter    ::  SHIELD_CONV=1.0e-12_DP
        
        real(kind=DP)               ::  qsum_l, del_l, omega_l, pn_l, shield_l
        real(kind=DP)               ::  umsa, qtemp(1:nctypes,1:nstypes,1:20)
        
        real(kind=DP)               ::  del_lp, omega_lp, pn_lp, shield_lp, shield_temp
        real(kind=DP)               ::  qtempp(1:nctypes,1:nstypes,1:20), umsa_p
        
        real(kind=DP)               ::  temp_u, temp_v, temp_up, temp_vp
        
        integer                     ::  i_ion1, i_ion2, i_ion3
        integer                     ::  it_tot

        if(.not.ion_switch) then
            a_res = 0.0e0_DP
            return
        end if
        
        l_dielec = Diel()
        
        qtemp = 0.0e0_DP
!Do full A_ion calc first...      
        c_shield = QE**2.0e0_DP / (4.0e0_DP * l_dielec * E0 *  KB * t * v)        

        qsum_l = 0.0e0_DP
             
        do i_ion1=1, nctypes
            do i_ion2=1, nstypes
                if( Seg_array(i_ion2)%charged ) then
                    do i_ion3=1, Comp_array(i_ion1)%comp(i_ion2)
                        qsum_l = qsum_l + Comp_array(i_ion1)%nm * Seg_array(i_ion2)%q**2.0e0_DP                        
                    end do
                end if 
            end do
        end do
        
        !initial shield guess
        shield = 0.5e0_DP * dsqrt( c_shield * qsum_l )

        !set delta  - this is constant for iteration       
        del_l = 0.0e0_DP

        do i_ion1=1, nctypes
            do i_ion2=1, nstypes
                if(Seg_array(i_ion2)%charged) then
                    do i_ion3=1, Comp_array(i_ion1)%comp(i_ion2)
                        del_l = del_l + Comp_array(i_ion1)%nm * Seg_array(i_ion2)%sig**3.0e0_DP                   
                    end do
                end if 
            end do
        end do
        
        del_l = 1.0e0_DP - PI / 6.0e0_DP / v * del_l            

        !iterate parameters to find shield
        it_tot = 0       
        
        do             
            omega_l = 0.0e0_DP
            
            do i_ion1=1, nctypes
                do i_ion2=1, nstypes
                    if(Seg_array(i_ion2)%charged) then
                        do i_ion3=1, Comp_array(i_ion1)%comp(i_ion2)
                            omega_l = omega_l + Comp_array(i_ion1)%nm * Seg_array(i_ion2)%sig**3.0e0_DP / &
                                &   (1.0e0_DP + shield * Seg_array(i_ion2)%sig)
                        end do
                    end if 
                end do
            end do

            omega_l = 1.0e0_DP + PI / v / 2.e0_DP / del_l * omega_l

            pn_l = 0.0e0_DP
            
            do i_ion1=1, nctypes
                do i_ion2=1, nstypes
                    if(Seg_array(i_ion2)%charged) then
                        do i_ion3=1, Comp_array(i_ion1)%comp(i_ion2)
!added *QE
                            pn_l = pn_l + Comp_array(i_ion1)%nm * Seg_array(i_ion2)%sig *       &
                            &   Seg_array(i_ion2)%q*QE / (1.0e0_DP + shield * Seg_array(i_ion2)%sig )
                           
                        end do
                    end if 
                end do
            end do
           
            pn_l = 1.0e0_DP / omega_l / v * pn_l 

            shield_l = 0.0e0_DP
            
            do i_ion1=1, nctypes
                do i_ion2=1, nstypes
                    if(Seg_array(i_ion2)%charged) then
                        do i_ion3=1, Comp_array(i_ion1)%comp(i_ion2)
                            qtemp(i_ion1,i_ion2,i_ion3) =                                   &
                            &       ( Seg_array(i_ion2)%q * QE - (PI / 2.0e0_DP / del_l) *  &
                            &       Seg_array(i_ion2)%sig**2.0e0_DP * pn_l ) / (1.0e0_DP + shield * Seg_array(i_ion2)%sig)                     

                            shield_l = shield_l + Comp_array(i_ion1)%nm * qtemp(i_ion1,i_ion2,i_ion3)**2.0e0_DP                     
                        end do
                    end if 
                end do
            end do
       
            shield_l = dsqrt( PI  / 4.0e0_DP / PI / E0 / l_dielec / KB / T / v * shield_l )
          
            it_tot = it_tot + 1
                    
            if(dabs( shield - shield_l ) / dabs( shield_l )  < SHIELD_CONV) then
                shield = shield_l
                exit
            else
                shield = shield_l
            end if    
        end do 
        
        umsa = 0.0e0_DP

        do i_ion1=1, nctypes
            do i_ion2=1, nstypes
                if(Seg_array(i_ion2)%charged) then
                    do i_ion3=1, Comp_array(i_ion1)%comp(i_ion2)
                        umsa = umsa + (Comp_array(i_ion1)%nm * Seg_array(i_ion2)%q**2.0e0_DP* QE**2.0e0_DP )  /   &
                        &   ( 1.0e0_DP + shield * Seg_array(i_ion2)%sig )  
                    end do
                end if 
            end do
        end do 
  
        umsa = -1.0e0_DP * v / 4.0e0_DP / PI / E0 / l_dielec * (shield / v * umsa +  &
        &       PI / 2.0e0_DP / del_l * omega_l * pn_l**2.0e0_DP )
    
!shield converged, do differential
        del_lp = 0.0e0_DP
        i_ion1=i_part
        
        do i_ion2=1, nstypes
            if(Seg_array(i_ion2)%charged) then
                do i_ion3=1, Comp_array(i_ion1)%comp(i_ion2)
                    del_lp = del_lp - Seg_array(i_ion2)%sig**3.0e0_DP                   
                end do
            end if 
        end do
      
        del_lp = del_lp * PI / 6.0e0_DP / v

!initial shield' guess        
        qsum_l = 0.0e0_DP           

        temp_u  = 0.0e0_DP
        do i_ion1=1, nctypes
            do i_ion2=1, nstypes
                if( Seg_array(i_ion2)%charged ) then
                    do i_ion3=1, Comp_array(i_ion1)%comp(i_ion2)
                        temp_u = temp_u + Comp_array(i_ion1)%nm * Seg_array(i_ion2)%q**2.0e0_DP                        
                    end do
                end if 
            end do
        end do
     
        shield_lp = 0.5e0_DP / (0.25e0_DP * c_shield * temp_u)**0.5e0_DP * 0.25e0_DP * c_shield * &
        &   Seg_array(i_part)%q**2.0e0_DP       
        
        do !shield differential
            temp_v  = 0.0e0_DP 
            temp_vp = 0.0e0_DP
            temp_up = -PI * del_lp / 2.0e0_DP / v / del_l**2.0e0_DP
            temp_u = PI / 2.0e0_DP / del_l / v
    
            do i_ion1=1, nctypes
                do i_ion2=1, nstypes
                    if(Seg_array(i_ion2)%charged) then
                        do i_ion3=1, Comp_array(i_ion1)%comp(i_ion2)
                        
                            temp_v = temp_v + Comp_array(i_ion1)%nm * Seg_array(i_ion2)%sig**3.0e0_DP / &
                                &   (1.0e0_DP + shield * Seg_array(i_ion2)%sig)
                        
                            if(i_ion1==i_part) then
                                temp_vp = temp_vp + (Seg_array(i_ion2)%sig**3.0e0_DP * (1.0e0_DP + shield   * & 
                                & Seg_array(i_ion2)%sig) - Comp_array(i_ion1)%nm * shield_lp                * &
                                & Seg_array(i_ion2)%sig**4.0e0_DP)                                          / &
                                & (1.0e0_DP + shield * Seg_array(i_ion2)%sig)**2.0e0_DP
                            else
                                temp_vp = temp_vp - (Comp_array(i_ion1)%nm * shield_lp                      * &
                                & Seg_array(i_ion2)%sig**4.0e0_DP)                                          / &
                                & (1.0e0_DP + shield * Seg_array(i_ion2)%sig)**2.0e0_DP
                            end if
    
                        end do
                    end if 
                end do
            end do     
    
            omega_lp = temp_u * temp_vp + temp_v * temp_up
     
            temp_up = -omega_lp / v / omega_l**2.0e0_DP
            temp_vp = 0.0e0_DP
            temp_v  = 0.0e0_DP
            temp_u  = 1.0e0_DP /omega_l / v
            
            do i_ion1=1, nctypes
                do i_ion2=1, nstypes
                    if(Seg_array(i_ion2)%charged) then
                        do i_ion3=1, Comp_array(i_ion1)%comp(i_ion2)
    
                            temp_v = temp_v + Comp_array(i_ion1)%nm * Seg_array(i_ion2)%sig *       &
                            &   Seg_array(i_ion2)%q*QE / (1.0e0_DP + shield * Seg_array(i_ion2)%sig )
                            
                            if(i_ion1==i_part) then
                                temp_vp = temp_vp + (Seg_array(i_ion2)%sig * Seg_array(i_ion2)%q*QE * &
                                &   (1.0e0_DP + shield * Seg_array(i_ion2)%sig)                     - &
                                &   shield_lp * Seg_array(i_ion2)%sig * (Comp_array(i_ion1)%nm      * &
                                &   Seg_array(i_ion2)%sig * Seg_array(i_ion2)%q*QE))                / &
                                &   (1.0e0_DP + shield * Seg_array(i_ion2)%sig)**2.0e0_DP
                            else
                                temp_vp = temp_vp - (shield_lp * Seg_array(i_ion2)%sig * (Comp_array(i_ion1)%nm * &
                                &   Seg_array(i_ion2)%sig * Seg_array(i_ion2)%q*QE))                            / &
                                &   (1.0e0_DP + shield * Seg_array(i_ion2)%sig)**2.0e0_DP                        
                            end if
                                
                        end do
                    end if 
                end do
            end do        
            
            pn_lp = temp_u * temp_vp + temp_v * temp_up
      
            do i_ion1=1, nctypes
                do i_ion2=1, nstypes
                    if(Seg_array(i_ion2)%charged) then
                        do i_ion3=1, Comp_array(i_ion1)%comp(i_ion2)
                        
                            temp_u  = Seg_array(i_ion2)%q * QE - (PI / 2.0e0_DP / del_l) *  &
                            &   Seg_array(i_ion2)%sig**2.0e0_DP * pn_l 
                            temp_up = (Seg_array(i_ion2)%sig**2.0e0_DP * PI * (pn_l*del_lp - pn_lp*del_l))  / &
                            &   2.0e0_DP / del_l**2.0e0_DP
                            temp_v  = 1.0e0_DP + shield * Seg_array(i_ion2)%sig
                            temp_vp = shield_lp * Seg_array(i_ion2)%sig 
    
                            qtempp(i_ion1,i_ion2,i_ion3) = (temp_v * temp_up - temp_u * temp_vp) / temp_v**2.0e0_DP                              
                        end do
                    end if 
                end do
            end do        
         
            temp_u = 0.0e0_DP
            temp_v = 0.0e0_DP           
            
            do i_ion1=1, nctypes
                do i_ion2=1, nstypes
                    if(Seg_array(i_ion2)%charged) then
                        do i_ion3=1, Comp_array(i_ion1)%comp(i_ion2)
    
                        
                            temp_u = temp_u + Comp_array(i_ion1)%nm * qtemp(i_ion1,i_ion2,i_ion3)**2.0e0_DP 
    
                            if(i_ion1==i_part) then
                                temp_v = temp_v + qtemp(i_ion1,i_ion2,i_ion3)**2.0e0_DP             + &
                                2.0e0_DP * Comp_array(i_ion1)%nm * qtemp(i_ion1,i_ion2,i_ion3)      * &
                                & qtempp(i_ion1,i_ion2,i_ion3)
                            else
                                temp_v = temp_v + qtemp(i_ion1,i_ion2,i_ion3)**2.0e0_DP             + &
                                2.0e0_DP * Comp_array(i_ion1)%nm * qtemp(i_ion1,i_ion2,i_ion3)      * &
                                & qtempp(i_ion1,i_ion2,i_ion3)
                            end if
                        end do
                    end if 
                end do
            end do     
    
            shield_temp = 0.5e0_DP * (PI  / 4.0e0_DP / PI / E0 / l_dielec / KB / T / v * temp_u)**(-0.5e0_DP)     * &
            &   (PI  / 4.0e0_DP / PI / E0 / l_dielec / KB / T / v * temp_v)
            
            if(dabs( shield_lp - shield_temp ) / dabs( shield_temp )  < SHIELD_CONV) then
                shield_lp = shield_temp
                exit
            else
                shield_lp = shield_temp
            end if      
                        
        end do
        
        umsa_p = 0.0e0_DP
        
        temp_u  = shield / v
        temp_v  = 0.0e0_DP
        temp_up = shield_lp / v
        temp_vp = 0.0e0_DP

        do i_ion1=1, nctypes
            do i_ion2=1, nstypes
                if(Seg_array(i_ion2)%charged) then
                    do i_ion3=1, Comp_array(i_ion1)%comp(i_ion2)
                        temp_v = temp_v + (Comp_array(i_ion1)%nm * Seg_array(i_ion2)%q**2.0e0_DP    * &
                        &   QE**2.0e0_DP ) / ( 1.0e0_DP + shield * Seg_array(i_ion2)%sig )  
                    
                        if(i_ion1==i_part) then                    
                            temp_vp = temp_vp                                                                   + &
                            &  (( 1.0e0_DP + shield * Seg_array(i_ion2)%sig ) * Seg_array(i_ion2)%q**2.0e0_DP   * &
                            &   QE**2.0e0_DP                                                                    - &
                            &   Comp_array(i_ion1)%nm * Seg_array(i_ion2)%q**2.0e0_DP * QE**2.0e0_DP            * &
                            &   Seg_array(i_ion2)%sig * shield_lp )                                             / &
                            &   ( 1.0e0_DP + shield * Seg_array(i_ion2)%sig )**2.0e0_DP
                        else
                            temp_vp = temp_vp                                                                   - &
                            &   (Comp_array(i_ion1)%nm * Seg_array(i_ion2)%q**2.0e0_DP * QE**2.0e0_DP           * &
                            &   Seg_array(i_ion2)%sig * shield_lp )                                             / &
                            &   ( 1.0e0_DP + shield * Seg_array(i_ion2)%sig )**2.0e0_DP
                        end if
                    
                    end do
                end if 
            end do
        end do 
       
        umsa_p = temp_u * temp_vp + temp_v * temp_up
        
        temp_u  = PI * omega_l * pn_l**2.0e0_DP 
        temp_up = PI * (omega_lp*pn_l**2.0e0_DP + 2.0e0_DP * pn_l * pn_lp * omega_l)
        temp_v  = 2.0e0_DP * del_l
        temp_vp = 2.0e0_DP * del_lp

        umsa_p = umsa_p + (temp_v * temp_up - temp_u * temp_vp) / temp_v**2.0e0_DP
        umsa_p = -umsa_p * v / (4.0e0_DP * PI * E0 * l_dielec) 
    
        a_res = umsa_p + 3.0e0_DP * shield**2.0e0_DP * shield_lp * KB * t * v / 3.0e0_DP / PI

        a_res = a_res / KB / t

        return
    end function A_ion_dn
!***************************************************************************************************  
    function A_ion_dv( ) result(a_res)
        implicit none

        real(kind=DP)           ::  a_res

        real(kind=DP)           ::  l_dielec, dl_dielecv
        real(kind=DP)               ::  c_shield, shield
        
        real(kind=DP), parameter    ::  SHIELD_CONV=1.0e-12_DP
        
        real(kind=DP)               ::  qsum_l, del_l, omega_l, pn_l, shield_l
        real(kind=DP)               ::  umsa, qtemp(1:nctypes,1:nstypes,1:20)
        
        real(kind=DP)               ::  del_lp, omega_lp, pn_lp, shield_lp, shield_temp
        real(kind=DP)               ::  qtempp(1:nctypes,1:nstypes,1:20), umsa_p
        
        real(kind=DP)               ::  temp_u, temp_v, temp_up, temp_vp, temp_w, temp_wp, temp_c
        
        integer                     ::  i_ion1, i_ion2, i_ion3
        integer                     ::  it_tot

        if(.not.ion_switch) then
            a_res = 0.0e0_DP
            return
        end if

        l_dielec   = Diel()
        dl_dielecv = dDielv()
        
        qtemp = 0.0e0_DP
!Do full A_ion calc first...      
        c_shield = QE**2.0e0_DP / (4.0e0_DP * l_dielec * E0 *  KB * t * v)        

        qsum_l = 0.0e0_DP
             
        do i_ion1=1, nctypes
            do i_ion2=1, nstypes
                if( Seg_array(i_ion2)%charged ) then
                    do i_ion3=1, Comp_array(i_ion1)%comp(i_ion2)
                        qsum_l = qsum_l + Comp_array(i_ion1)%nm * Seg_array(i_ion2)%q**2.0e0_DP                        
                    end do
                end if 
            end do
        end do
        
        !initial shield guess
        shield = 0.5e0_DP * dsqrt( c_shield * qsum_l )

        !set delta  - this is constant for iteration       
        del_l = 0.0e0_DP

        do i_ion1=1, nctypes
            do i_ion2=1, nstypes
                if(Seg_array(i_ion2)%charged) then
                    do i_ion3=1, Comp_array(i_ion1)%comp(i_ion2)
                        del_l = del_l + Comp_array(i_ion1)%nm * Seg_array(i_ion2)%sig**3.0e0_DP                   
                    end do
                end if 
            end do
        end do
        
        del_l = 1.0e0_DP - PI / 6.0e0_DP / v * del_l            

        !iterate parameters to find shield
        it_tot = 0       
        
        do             
            omega_l = 0.0e0_DP
            
            do i_ion1=1, nctypes
                do i_ion2=1, nstypes
                    if(Seg_array(i_ion2)%charged) then
                        do i_ion3=1, Comp_array(i_ion1)%comp(i_ion2)
                            omega_l = omega_l + Comp_array(i_ion1)%nm * Seg_array(i_ion2)%sig**3.0e0_DP / &
                                &   (1.0e0_DP + shield * Seg_array(i_ion2)%sig)
                        end do
                    end if 
                end do
            end do

            omega_l = 1.0e0_DP + PI / v / 2.e0_DP / del_l * omega_l

            pn_l = 0.0e0_DP
            
            do i_ion1=1, nctypes
                do i_ion2=1, nstypes
                    if(Seg_array(i_ion2)%charged) then
                        do i_ion3=1, Comp_array(i_ion1)%comp(i_ion2)
!added *QE
                            pn_l = pn_l + Comp_array(i_ion1)%nm * Seg_array(i_ion2)%sig *       &
                            &   Seg_array(i_ion2)%q*QE / (1.0e0_DP + shield * Seg_array(i_ion2)%sig )
                           
                        end do
                    end if 
                end do
            end do
           
            pn_l = 1.0e0_DP / omega_l / v * pn_l 

            shield_l = 0.0e0_DP
            
            do i_ion1=1, nctypes
                do i_ion2=1, nstypes
                    if(Seg_array(i_ion2)%charged) then
                        do i_ion3=1, Comp_array(i_ion1)%comp(i_ion2)
                            qtemp(i_ion1,i_ion2,i_ion3) =                                   &
                            &       ( Seg_array(i_ion2)%q * QE - (PI / 2.0e0_DP / del_l) *  &
                            &       Seg_array(i_ion2)%sig**2.0e0_DP * pn_l ) / (1.0e0_DP + shield * Seg_array(i_ion2)%sig)                     

                            shield_l = shield_l + Comp_array(i_ion1)%nm * qtemp(i_ion1,i_ion2,i_ion3)**2.0e0_DP                     
                        end do
                    end if 
                end do
            end do
       
            shield_l = dsqrt( PI  / 4.0e0_DP / PI / E0 / l_dielec / KB / T / v * shield_l )
          
            it_tot = it_tot + 1
                    
            if(dabs( shield - shield_l ) / dabs( shield_l )  < SHIELD_CONV) then
                shield = shield_l
                exit
            else
                shield = shield_l
            end if    
        end do 
        
        umsa = 0.0e0_DP

        do i_ion1=1, nctypes
            do i_ion2=1, nstypes
                if(Seg_array(i_ion2)%charged) then
                    do i_ion3=1, Comp_array(i_ion1)%comp(i_ion2)
                        umsa = umsa + (Comp_array(i_ion1)%nm * Seg_array(i_ion2)%q**2.0e0_DP* QE**2.0e0_DP )  /   &
                        &   ( 1.0e0_DP + shield * Seg_array(i_ion2)%sig )  
                    end do
                end if 
            end do
        end do 
  
        umsa = -1.0e0_DP * v / 4.0e0_DP / PI / E0 / l_dielec * (shield / v * umsa +  &
        &       PI / 2.0e0_DP / del_l * omega_l * pn_l**2.0e0_DP )
    
!shield converged, do differential
        del_lp = 0.0e0_DP
        
        do i_ion1=1, nctypes
            do i_ion2=1, nstypes
                if(Seg_array(i_ion2)%charged) then
                    do i_ion3=1, Comp_array(i_ion1)%comp(i_ion2)
                        del_lp = del_lp + Comp_array(i_ion1)%nm * Seg_array(i_ion2)%sig**3.0e0_DP                   
                    end do
                end if 
            end do
        end do
      
        del_lp = del_lp * PI / 6.0e0_DP / v**2.0e0_DP

!initial shield guess        
        temp_u  = 0.0e0_DP
        
        do i_ion1=1, nctypes
            do i_ion2=1, nstypes
                if( Seg_array(i_ion2)%charged ) then
                    do i_ion3=1, Comp_array(i_ion1)%comp(i_ion2)
                        temp_u = temp_u + Comp_array(i_ion1)%nm * Seg_array(i_ion2)%q**2.0e0_DP                        
                    end do
                end if 
            end do
        end do
        
        temp_v = QE**2.0e0_DP * temp_u / (4.0e0_DP * E0 *  KB * t )
        shield_lp = 0.25e0_DP * (temp_v / (l_dielec*v))**(-0.5e0_DP)             * &
        &   (-temp_v * (dl_dielecv*v + l_dielec) / (l_dielec**2.0e0_DP * v**2.0e0_DP))       

        do !shield differential
            temp_u = PI / 2.0e0_DP / del_l / v
            temp_up = -PI * (del_l + del_lp*v) / (2.0e0_DP * v**2.0e0_DP * del_l**2.0e0_DP)

            temp_v  = 0.0e0_DP 
            temp_vp = 0.0e0_DP

            do i_ion1=1, nctypes
                do i_ion2=1, nstypes
                    if(Seg_array(i_ion2)%charged) then
                        do i_ion3=1, Comp_array(i_ion1)%comp(i_ion2)
                        
                            temp_v = temp_v + Comp_array(i_ion1)%nm * Seg_array(i_ion2)%sig**3.0e0_DP / &
                                &   (1.0e0_DP + shield * Seg_array(i_ion2)%sig)
                        
                            temp_vp = temp_vp - Comp_array(i_ion1)%nm * shield_lp                * &
                            & Seg_array(i_ion2)%sig**4.0e0_DP                                    / &
                            & (1.0e0_DP + shield * Seg_array(i_ion2)%sig)**2.0e0_DP
                        end do
                    end if 
                end do
            end do     
    
            omega_lp = temp_u * temp_vp + temp_v * temp_up
 
            temp_u  = 1.0e0_DP /omega_l / v  
            temp_up = -(omega_lp*v + omega_l) / (v**2.0e0_DP * omega_l**2.0e0_DP)
            temp_vp = 0.0e0_DP
            temp_v  = 0.0e0_DP
            
            
            do i_ion1=1, nctypes
                do i_ion2=1, nstypes
                    if(Seg_array(i_ion2)%charged) then
                        do i_ion3=1, Comp_array(i_ion1)%comp(i_ion2)
    
                            temp_v = temp_v + Comp_array(i_ion1)%nm * Seg_array(i_ion2)%sig *       &
                            &   Seg_array(i_ion2)%q*QE / (1.0e0_DP + shield * Seg_array(i_ion2)%sig )
                            
                            temp_vp = temp_vp - shield_lp * Seg_array(i_ion2)%sig**2.0e0_DP         * &
                            &   Comp_array(i_ion1)%nm * Seg_array(i_ion2)%q*QE                      / &
                            &   (1.0e0_DP + shield * Seg_array(i_ion2)%sig)**2.0e0_DP                        
                        end do
                    end if 
                end do
            end do        
            
            pn_lp = temp_u * temp_vp + temp_v * temp_up
     
            do i_ion1=1, nctypes
                do i_ion2=1, nstypes
                    if(Seg_array(i_ion2)%charged) then
                        do i_ion3=1, Comp_array(i_ion1)%comp(i_ion2)
                        
                        
                            temp_v  = 1.0e0_DP + shield * Seg_array(i_ion2)%sig
                            temp_vp = shield_lp * Seg_array(i_ion2)%sig 

                            temp_u  = Seg_array(i_ion2)%q * QE - PI             * &
                            &   Seg_array(i_ion2)%sig**2.0e0_DP * pn_l          / &
                            &   (2.0e0_DP * del_l)

                            temp_up = PI*Seg_array(i_ion2)%sig**2.0e0_DP        * &
                            &   (pn_l*del_lp - del_l*pn_lp)                     / &
                            &   (2.0e0_DP * del_l**2.0e0_DP)

                            qtempp(i_ion1,i_ion2,i_ion3) = (temp_v * temp_up - temp_u * temp_vp) / temp_v**2.0e0_DP                              
                        end do
                    end if 
                end do
            end do  

            temp_c  = 1.0e0_DP / (4.0e0_DP * E0 * KB * T)
            temp_u  = 1.0e0_DP / (l_dielec * v)
            temp_up = -1.0e0_DP * (dl_dielecv * v + l_dielec) / (l_dielec**2.0e0_DP * v**2.0e0_DP) 
  
            temp_v  = 0.0e0_DP
            temp_vp = 0.0e0_DP           
            
            do i_ion1=1, nctypes
                do i_ion2=1, nstypes
                    if(Seg_array(i_ion2)%charged) then
                        do i_ion3=1, Comp_array(i_ion1)%comp(i_ion2)
                            temp_v = temp_v + Comp_array(i_ion1)%nm * qtemp(i_ion1,i_ion2,i_ion3)**2.0e0_DP 
                            
                            temp_vp = temp_vp + 2.0e0_DP * Comp_array(i_ion1)%nm        * &
                            &   qtemp(i_ion1,i_ion2,i_ion3) * qtempp(i_ion1,i_ion2,i_ion3)
                        end do
                    end if 
                end do
            end do 
            
            temp_w = temp_c*(temp_v*temp_up + temp_u*temp_vp)  
    
            shield_temp = 0.5e0_DP * (temp_c*temp_v/(l_dielec*v))**(-0.5e0_DP) * temp_w

            if(dabs( shield_lp - shield_temp ) / dabs( shield_temp )  < SHIELD_CONV) then
                shield_lp = shield_temp
                exit
            else
                shield_lp = shield_temp
            end if      
                        
        end do

        temp_u  = -1.0e0_DP * v /(4.0e0_DP * PI * E0 * l_dielec)
        temp_up = -1.0e0_DP /(4.0e0_DP * PI * E0)                             * &
        &   (l_dielec - v * dl_dielecv) / l_dielec**2.0e0_DP    

        temp_v  = 0.0e0_DP
        temp_vp = 0.0e0_DP
        
        temp_w  = PI/(2.0e0_DP*del_l) * omega_l*pn_l**2.0e0_DP
        temp_wp =  -PI*del_lp/(2.0e0_DP*del_l**2.0e0_DP) * omega_l*pn_l**2.0e0_DP   + &
        &           PI/(2.0e0_DP*del_l) * omega_lp*pn_l**2.0e0_DP                   + &
        &           PI/del_l * omega_l * pn_l*pn_lp

        do i_ion1=1, nctypes
            do i_ion2=1, nstypes
                if(Seg_array(i_ion2)%charged) then
                    do i_ion3=1, Comp_array(i_ion1)%comp(i_ion2)
                        temp_v = temp_v + (Comp_array(i_ion1)%nm * Seg_array(i_ion2)%q**2.0e0_DP    * &
                        &   QE**2.0e0_DP ) / ( 1.0e0_DP + shield * Seg_array(i_ion2)%sig )  
                    
                        temp_vp = -shield*Seg_array(i_ion2)%sig*Seg_array(i_ion2)%q**2.0e0_DP*QE**2.0e0_DP      / &
                        &   ( 1.0e0_DP + shield * Seg_array(i_ion2)%sig )    
                    end do
                end if 
            end do
        end do 
       
        temp_vp = shield/v * temp_vp  +  temp_v * (v*shield_lp - shield)/v**2.0e0_DP
        temp_v = temp_v * shield / v
        umsa_p = temp_u * (temp_vp + temp_wp) + temp_up * (temp_v + temp_w)
    
        a_res = -(umsa_p + 3.0e0_DP * shield**2.0e0_DP * shield_lp * KB * t * v / 3.0e0_DP / PI       + &
        &       shield**3.0e0_DP * KB * t / 3.0e0_DP / PI)

        return
    end function A_ion_dv
!***************************************************************************************************  

!***************************************************************************************************
end module Ion_mod
!***************************************************************************************************
!***************************************************************************************************
!   SAFT Module to:
!       1. Calculate common differentials
!           a. with respect to Ni -> chemical potentials
!           b. with respect to V  -> pressure
!***************************************************************************************************
!
!***************************************************************************************************
module Diff_mod
!***************************************************************************************************
!Modules
!=======    
    use Types_mod           ! Definitions of types and double precision
    use Global_mod          ! Important global parameters
    use Setup_mod           ! To setup the system
!***************************************************************************************************
    implicit none
!***************************************************************************************************
    contains
!***************************************************************************************************    
    subroutine Ndiffs( i_part ) !For mu calculations
        implicit none
        
        integer, intent(in) ::  i_part
        real(kind=DP)       ::  t_dii
        real(kind=DP)       ::  diffv, diffu, diffvp, diffup
        real(kind=DP)       ::  tdiff1, tdiff2, tdiff3
        integer             ::  i_ndiff1, i_ndiff2
       
        call Sys_setup() !ensure eveything initialised   
            
        !sum of particles
        nsum = sum(Comp_array(:)%nm)
        
        !mole fractions
        do i_ndiff1=1,nctypes
            if(i_ndiff1==i_part) then
                Comp_array(i_ndiff1)%dxin = (1.0e0_DP - Comp_array(i_ndiff1)%xi) / nsum
            else
                Comp_array(i_ndiff1)%dxin = -Comp_array(i_ndiff1)%xi / nsum
            end if
        end do
      
        !rho
        drhon = 1.0e0_DP / v / NA
      
        !mole fraction of spherical segments
        diffv  = 0.0e0_DP
        diffvp = 0.0e0_DP
        
        do i_ndiff1=1, nctypes
        do i_ndiff2=1, nstypes
            diffv = diffv + Comp_array(i_ndiff1)%comp(i_ndiff2) * Seg_array(i_ndiff2)%nseg *Comp_array(i_ndiff1)%xi &
            &   * Seg_array(i_ndiff2)%sf
            
            diffvp = diffvp + Comp_array(i_ndiff1)%comp(i_ndiff2) * Seg_array(i_ndiff2)%nseg *Comp_array(i_ndiff1)%dxin &
            &   * Seg_array(i_ndiff2)%sf
        end do
        end do

        do i_ndiff1=1, nstypes       
            diffu  = 0.0e0_DP
            diffup = 0.0e0_DP
            
            do i_ndiff2=1, nctypes
                diffu = diffu +  (Comp_array(i_ndiff2)%comp(i_ndiff1) * Seg_array(i_ndiff1)%sf *      &  
                &                  Seg_array(i_ndiff1)%nseg * Comp_array(i_ndiff2)%xi) 
                
                diffup = diffup +  (Comp_array(i_ndiff2)%comp(i_ndiff1) * Seg_array(i_ndiff1)%sf *      &  
                &                  Seg_array(i_ndiff1)%nseg * Comp_array(i_ndiff2)%dxin) 
            end do

            Seg_array(i_ndiff1)%dxsin = (diffv*diffup - diffu*diffvp) / diffv**2.0e0_DP                    
        end do 

        dsigx3n = 0.0e0_DP
        
        !associaion sigx3
        do i_ndiff1=1, nstypes
        do i_ndiff2=1, nstypes
            dsigx3n = dsigx3n + ( Seg_array(i_ndiff1)%xsi * Seg_array(i_ndiff2)%dxsin    +   &
            &                     Seg_array(i_ndiff2)%xsi * Seg_array(i_ndiff1)%dxsin )  *   &
            &                     sig(i_ndiff1, i_ndiff2)**3.0e0_DP
        end do
        end do
       
        !Mono pieces
        !rhos
        diffu  = 0.0e0_DP
        diffup = 0.0e0_DP

        do i_ndiff1=1,nctypes
        do i_ndiff2=1,nstypes
            diffu = diffu + Comp_array(i_ndiff1)%xi * Comp_array(i_ndiff1)%comp(i_ndiff2) &  
            &   * Seg_array(i_ndiff2)%sf * Seg_array(i_ndiff2)%nseg
            diffup = diffup + Comp_array(i_ndiff1)%dxin * Comp_array(i_ndiff1)%comp(i_ndiff2) &  
            &   * Seg_array(i_ndiff2)%sf * Seg_array(i_ndiff2)%nseg
        end do
        end do

        drhosn = diffup*rho + diffu*drhon  !moles per m^3

        !zl
        dzln(:) = 0.0e0_DP
        
        do i_ndiff1=1, nstypes
            t_dii = dij(i_ndiff1,i_ndiff1) 

            dzln(0) = dzln(0) + PI * rhos / 6.0e0_DP * Seg_array(i_ndiff1)%dxsin                      + &
            &   PI * drhosn / 6.0e0_DP * Seg_array(i_ndiff1)%xsi

            dzln(1) = dzln(1) + PI * rhos / 6.0e0_DP * Seg_array(i_ndiff1)%dxsin * t_dii              + &
            &   PI * drhosn / 6.0e0_DP * Seg_array(i_ndiff1)%xsi * t_dii

            dzln(2) = dzln(2) + PI * rhos / 6.0e0_DP * Seg_array(i_ndiff1)%dxsin * t_dii**2.0e0_DP    + &
            &   PI * drhosn / 6.0e0_DP * Seg_array(i_ndiff1)%xsi * t_dii**2.0e0_DP

            dzln(3) = dzln(3) + PI * rhos / 6.0e0_DP * Seg_array(i_ndiff1)%dxsin * t_dii**3.0e0_DP    + &
            &   PI * drhosn / 6.0e0_DP * Seg_array(i_ndiff1)%xsi * t_dii**3.0e0_DP              

        end do
     
        dzln = dzln * NA
   
        !zetax
        tdiff1 = 0.0e0_DP
        tdiff2 = 0.0e0_DP
        tdiff3 = 0.0e0_DP
        
        do i_ndiff1=1, nstypes
        do i_ndiff2=1, nstypes
            tdiff1 = tdiff1 + Seg_array(i_ndiff1)%xsi * Seg_array(i_ndiff2)%xsi * dij(i_ndiff1,i_ndiff2)**3.0e0_DP
            tdiff2 = tdiff2 + Seg_array(i_ndiff1)%dxsin * Seg_array(i_ndiff2)%xsi * dij(i_ndiff1,i_ndiff2)**3.0e0_DP
            tdiff3 = tdiff3 + Seg_array(i_ndiff1)%xsi * Seg_array(i_ndiff2)%dxsin * dij(i_ndiff1,i_ndiff2)**3.0e0_DP
        end do
        end do
        
        dzetax_sumn = tdiff2 + tdiff3
        dzetaxn = (PI*drhosn/6.0e0_DP * tdiff1 + PI*rhos/6.0e0_DP * (tdiff2 + tdiff3)) * NA
        
        !zetabar
        diffu  = 0.0e0_DP
        diffup = 0.0e0_DP

        do i_ndiff1=1, nstypes
        do i_ndiff2=1, nstypes
           diffu  = diffu + Seg_array(i_ndiff1)%xsi * Seg_array(i_ndiff2)%xsi * sig(i_ndiff1,i_ndiff2)**3.0e0_DP  
           diffup = diffup + Seg_array(i_ndiff1)%dxsin * Seg_array(i_ndiff2)%xsi * sig(i_ndiff1,i_ndiff2)**3.0e0_DP  + &
           &                 Seg_array(i_ndiff1)%xsi * Seg_array(i_ndiff2)%dxsin * sig(i_ndiff1,i_ndiff2)**3.0e0_DP
        end do
        end do
        
        !Averaged packing fractions
        dzetabarn = NA * PI * rhos / 6.0e0_DP * diffup + NA * PI * drhosn / 6.0e0_DP * diffu

        !KHS
        diffu  = (1.d0 - zetax)**4.0e0_DP
        diffv  = 1.d0 + 4.0e0_DP * zetax + 4.0e0_DP * zetax**2 - 4.0e0_DP * zetax**3.0e0_DP + zetax**4.0e0_DP
        diffup = -4.0e0_DP * dzetaxn * (1.d0 - zetax)**3.0e0_DP
        diffvp = dzetaxn *( 4.0e0_DP + 8.0e0_DP * zetax - 12.0e0_DP * zetax**2.0e0_DP + 4.0e0_DP * zetax**3.0e0_DP)
         
        dKHSn = (diffv*diffup - diffu*diffvp)/diffv**2.0e0_DP

        return
    end subroutine Ndiffs
!*************************************************************************************************** 
    subroutine Vdiffs( ) !For pressure calculations
        implicit none
        
        real(kind=DP)       ::  t_dii
        real(kind=DP)       ::  diffv, diffu, diffvp, diffup
        real(kind=DP)       ::  tdiff1, tdiff2, tdiff3
        integer             ::  i_vdiff1, i_vdiff2
    
        call Sys_setup() !ensure eveything initialised   

        !sum of particles
        nsum = sum(Comp_array(:)%nm)
               
        !rho
        drhov = -rho / v 

        !Mono pieces
        !rhos
        drhosv = 0.e0_DP 

        do i_vdiff1=1,nctypes
        do i_vdiff2=1,nstypes
            drhosv = drhosv + Comp_array(i_vdiff1)%xi * Comp_array(i_vdiff1)%comp(i_vdiff2) &  
            &   * Seg_array(i_vdiff2)%sf * Seg_array(i_vdiff2)%nseg
        end do
        end do
        
        drhosv = drhosv * drhov

        !zl
        dzlv(:) = 0.0e0_DP
        
        do i_vdiff1=1,nstypes
            t_dii = dij(i_vdiff1,i_vdiff1)     
           
            dzlv(0) = dzlv(0) + Seg_array(i_vdiff1)%xsi
            dzlv(1) = dzlv(1) + Seg_array(i_vdiff1)%xsi * t_dii
            dzlv(2) = dzlv(2) + Seg_array(i_vdiff1)%xsi * t_dii**2.0e0_DP
            dzlv(3) = dzlv(3) + Seg_array(i_vdiff1)%xsi * t_dii**3.0e0_DP          
        end do

        dzlv(0) = dzlv(0) * PI * drhosv / 6.0e0_DP * NA              
        dzlv(1) = dzlv(1) * PI * drhosv / 6.0e0_DP * NA              
        dzlv(2) = dzlv(2) * PI * drhosv / 6.0e0_DP * NA            
        dzlv(3) = dzlv(3) * PI * drhosv / 6.0e0_DP * NA  

        !zetax
        tdiff1 = 0.0e0_DP

        do i_vdiff1=1, nstypes
        do i_vdiff2=1, nstypes
            tdiff1 = tdiff1 + Seg_array(i_vdiff1)%xsi * Seg_array(i_vdiff2)%xsi * dij(i_vdiff1,i_vdiff2)**3.0e0_DP
        end do
        end do      

        dzetaxv = tdiff1 * PI * drhosv / 6.0e0_DP * NA

        !zetabar
        tdiff1 = 0.0e0_DP

        do i_vdiff1=1, nstypes
        do i_vdiff2=1, nstypes
            tdiff1 = tdiff1 + Seg_array(i_vdiff1)%xsi * Seg_array(i_vdiff2)%xsi * sig(i_vdiff1,i_vdiff2)**3.0e0_DP 
        end do
        end do
        
        !Averaged packing fractions
        dzetabarv = NA * tdiff1 * PI * drhosv / 6.0e0_DP

        !KHS
        diffu  = (1.d0 - zetax)**4.0e0_DP
        diffv  = 1.d0 + 4.0e0_DP * zetax + 4.0e0_DP * zetax**2 - 4.0e0_DP * zetax**3.0e0_DP + zetax**4.0e0_DP
        diffup = -4.0e0_DP * dzetaxv * (1.d0 - zetax)**3.0e0_DP
        diffvp = dzetaxv *( 4.0e0_DP + 8.0e0_DP * zetax - 12.0e0_DP * zetax**2.0e0_DP + 4.0e0_DP * zetax**3.0e0_DP)
         
        dKHSv = (diffv*diffup - diffu*diffvp)/diffv**2.0e0_DP

        return
    end subroutine Vdiffs
!***************************************************************************************************

!*************************************************************************************************** 
end module Diff_mod
!***************************************************************************************************
!***************************************************************************************************
!   SAFT Module to:
!       1. Calculate pressure
!
!       Updated March 2017 to use analytical derivatives
!
!***************************************************************************************************
!
!***************************************************************************************************
module Press_mod
!***************************************************************************************************
!Modules
!=======    
    use Types_mod           ! Definitions of types and double precision
    use Global_mod          ! Important global parameters
    use Setup_mod           ! To setup the system
    use Ideal_mod           ! Calculate ideal A
    use Mono_mod            ! Calculate mono A
    use Chain_mod           ! Calculate chain A
    use Assoc_mod           ! Calculate assoc A
    use Ion_mod             ! Calculate ion A   
    use Diff_mod            ! Differentials setup
!***************************************************************************************************
    implicit none
!***************************************************************************************************
    public  ::  Press
!***************************************************************************************************
    contains
!***************************************************************************************************    
    function Press( ) result(p_out)
        implicit none

        real(kind=DP)           ::  p_out, p_ideal, p_mono, p_chain, p_assoc, p_born, p_ion

! DETERMINE PRESSURE BY CONTRIBUTION      
        call Vdiffs()

        p_ideal =  A_ideal_dv()* KB * NA * T             
        p_mono  = -A_mono_dv() * KB * NA * T     
        p_chain = -A_chain_dv()    
        p_assoc =  A_assoc_dv()      
        p_born  = -A_born_dv() * KB * NA * T     
        p_ion   =  A_ion_dv()

        !hack to avoid p=NaN
        if(p_ideal/=p_ideal) p_ideal=0.0e0_DP
        if(p_mono/=p_mono)   p_mono=0.0e0_DP
        if(p_chain/=p_chain) p_chain=0.0e0_DP
        if(p_assoc/=p_assoc) p_assoc=0.0e0_DP
        if(p_born/=p_born)   p_born=0.0e0_DP
        if(p_ion/=p_ion)     p_ion=0.0e0_DP

        p_out   =  p_ideal + p_mono + p_chain + p_assoc + p_born + p_ion 
       
 99     return
    end function Press
!***************************************************************************************************
!***************************************************************************************************
!   SAFT Module to:
!       1. Calculate chemical potential
!
!       Updated March 2017 to use analytical integrals
!
!***************************************************************************************************
!
end module Press_mod
!***************************************************************************************************
module Mu_mod
!***************************************************************************************************
!Modules
!=======    
    use Types_mod           ! Definitions of types and double precision
    use Global_mod          ! Important global parameters
    use Setup_mod           ! To setup the system
    use Ideal_mod           ! Calculate ideal A
    use Mono_mod            ! Calculate mono A
    use Chain_mod           ! Calculate chain A
    use Assoc_mod           ! Calculate assoc A
    use Ion_mod             ! Calculate ion A
    use Diff_mod            ! Differentials module  
!***************************************************************************************************
    implicit none
!***************************************************************************************************
    contains
!***************************************************************************************************    
    function Mu( i_comp ) result(mu_out)
        implicit none
        
        integer, intent(in)      ::  i_comp
        
        real(kind=DP)            ::  mu_out
        real(kind=DP)            ::  mu_ideal, mu_mono, mu_chain, &
                                    &   mu_ion, mu_born, mu_assoc   
       
        !ensure nsum is set
        Comp_array(:)%nm = Comp_array(:)%xi
        Comp_array(:)%nm = Comp_array(:)%nm * NA
        nsum = sum(Comp_array(:)%nm)

        !call differential set up
        call Ndiffs( i_comp ) 
        
        mu_ideal = A_ideal_dn( i_comp ) * KB * t * NA
        mu_mono  = A_mono_dn( ) * NA * KB * t       
        mu_assoc = A_assoc_dn( ) * NA * KB * t
        mu_born  = A_born_dn( i_comp ) * KB * t * nsum * NA
        mu_ion   = A_ion_dn( i_comp ) * NA * KB * t
        mu_chain = A_chain_dn( i_comp ) * NA

        mu_out = mu_ideal + mu_mono + mu_assoc + mu_chain + mu_born + mu_ion
!        print*,i_comp,mu_out,mu_ideal,mu_mono,mu_chain,mu_assoc,mu_born,mu_ion
        
        return
    end function Mu
!***************************************************************************************************
    function Mu_neut( i_comp ) result(mu_out)
        implicit none
        
        integer, intent(in)      ::  i_comp
        
        real(kind=DP)            ::  mu_out
        real(kind=DP)            ::  mu_ideal, mu_mono, mu_chain, mu_assoc   
       
        !ensure nsum is set
        Comp_array(:)%nm = Comp_array(:)%xi
        Comp_array(:)%nm = Comp_array(:)%nm * NA
        nsum = sum(Comp_array(:)%nm)

        !call differential set up
        call Ndiffs( i_comp ) 
        
        mu_ideal = A_ideal_dn( i_comp ) * KB * t * NA
        mu_mono  = A_mono_dn( ) * NA * KB * t       
        mu_assoc = A_assoc_dn( ) * NA * KB * t
        mu_chain = A_chain_dn( i_comp ) * NA

        mu_out = mu_ideal + mu_mono + mu_assoc + mu_chain 
!        print*,i_comp,mu_out,mu_ideal,mu_mono,mu_chain,mu_assoc
        
        return
    end function Mu_neut
!***************************************************************************************************
    function Mu_id( i_comp ) result(mu_out)
        implicit none
        
        integer, intent(in)      ::  i_comp
        
        real(kind=DP)            ::  mu_out
        real(kind=DP)            ::  mu_ideal, mu_chain

       
        !ensure nsum is set
        Comp_array(:)%nm = Comp_array(:)%xi
        Comp_array(:)%nm = Comp_array(:)%nm * NA
        nsum = sum(Comp_array(:)%nm)

        !call differential set up
        call Ndiffs( i_comp ) 
        
        mu_ideal = A_ideal_dn( i_comp ) * KB * t * NA
        mu_chain = A_chain_dn( i_comp ) * NA

        mu_out = mu_ideal + mu_chain 
        !print*,i_comp,mu_out,mu_ideal,mu_mono,mu_chain,mu_born,mu_ion,mu_assoc
        
        return
    end function Mu_id
!***************************************************************************************************
    function Mu_res( i_comp ) result(mu_out) !As above, sans ideal
        implicit none
        
        integer, intent(in)      ::  i_comp
        
        real(kind=DP)            ::  mu_out
        real(kind=DP)            ::  mu_mono, mu_chain, &
                                    &   mu_ion, mu_born, mu_assoc   
    
        !ensure nsum is set
        Comp_array(:)%nm = Comp_array(:)%xi
        Comp_array(:)%nm = Comp_array(:)%nm * NA
        nsum = sum(Comp_array(:)%nm)

        !call differential set up
        call Ndiffs( i_comp ) 
        
        mu_mono  = A_mono_dn( ) * NA * KB * t 
        mu_assoc = A_assoc_dn( ) * NA * KB * t
        mu_born  = A_born_dn( i_comp ) * KB * t * nsum * NA
        mu_ion   = A_ion_dn( i_comp ) * NA * KB * t
        mu_chain = A_chain_dn( i_comp ) * NA

        mu_out = mu_mono + mu_assoc + mu_chain + mu_born + mu_ion
!        print*,i_comp,mu_out,mu_mono,mu_chain,mu_born,mu_ion,mu_assoc
        
        return
    end function Mu_res
!***************************************************************************************************  
    function Mu_res_inf( i_comp ) result(mu_out) !As above, sans ideal
        implicit none
        
        integer, intent(in)      ::  i_comp
        
        real(kind=DP)            ::  mu_out
        real(kind=DP)            ::  mu_mono, mu_chain, &
                                    &   mu_born, mu_assoc   
                
        !ensure nsum is set
        Comp_array(:)%nm = Comp_array(:)%xi
        Comp_array(:)%nm = Comp_array(:)%nm * NA
        nsum = sum(Comp_array(:)%nm)

        !call differential set up
        call Ndiffs( i_comp ) 
        
        mu_mono  = A_mono_dn( ) * NA * KB * t 
        mu_assoc = A_assoc_dn( ) * NA * KB * t
        mu_born  = A_born_dn( i_comp ) * KB * t * nsum * NA
        if(mu_born/=mu_born) mu_born=0.0e0_DP        
        mu_chain = A_chain_dn( i_comp ) * NA

        mu_out = mu_mono + mu_assoc + mu_chain + mu_born

        return
    end function Mu_res_inf
!***************************************************************************************************

!***************************************************************************************************
end module Mu_mod
!***************************************************************************************************
!***************************************************************************************************
!   SAFT Module to:
!       1. Calculate volume
!***************************************************************************************************
!
!***************************************************************************************************
module Vol_mod
!***************************************************************************************************
!Modules
!=======    
    use Types_mod           ! Definitions of types and double precision
    use Global_mod          ! Important global parameters
    use Setup_mod           ! To setup the system
    use Press_mod           ! Calculate pressure
    use Ideal_mod           ! Calculate ideal A
    use Mono_mod            ! Calculate mono A
    use Chain_mod           ! Calculate chain A
    use Assoc_mod           ! Calculate assoc A
    use Ion_mod             ! Calculate ion A
    use Mu_mod              ! Calculate chem pot
!***************************************************************************************************
    implicit none
!***************************************************************************************************
    contains
!***************************************************************************************************
    function Vol_dens_g(limited) result(v_out)   
        implicit none
        
        logical, optional   ::  limited
        real(kind=DP)       ::  v_out
        real(kind=DP)       ::  lv_array(1:10), lg_array(1:10)
        real(kind=DP)       ::  v_temp1, v_temp2, p_temp1, p_temp2, d_v
        real(kind=DP)       ::  lv1, lv2, lp1, lp2, lr, g_local, grad_local
        integer             ::  i_rho, nroots, piter, froot, imusum
        integer, parameter          :: plimit1=100, plimit2=500, i_rho_limit=5000
        real(kind=DP), parameter    ::  P_crit=1.0e-6_DP,P_crit2=1.0e-4

        !set limit if required
        if(present(limited)) then
            limited = .true.
        end if
        
        nroots = 0
        
        lr  = 0.0e0_DP
        i_rho = 1
                  
        do      
            if(i_rho<100) then
                lr  = lr+1.0e-4_DP 
            else if(i_rho<200) then  
                lr  = lr+1.0e-3_DP 
            else if(i_rho<300) then  
                lr  = lr+1.0e-2_DP   
            else if(i_rho<400) then  
                lr  = lr+1.0e-1_DP 
            else if(i_rho<500) then  
                lr  = lr+1.0e0_DP  
            else if(i_rho<1000) then  
                lr  = lr+5.0e0_DP 
            else if(i_rho<1500) then  
                lr  = lr+10.0e0_DP 
            else if(i_rho<2000) then  
                lr  = lr+50.0e0_DP 
            else if(i_rho<2500) then  
                lr  = lr+100.0e0_DP 
            else                             
                go to 30
            end if
                             
            lv2 = 1.0e0_DP/lr
            if(lv2<=1.0e-5_DP) go to 30
            v   = lv2   
            lp2 = Press()    

            if(i_rho>0) then

                !sensible p and within desired p
                if(((lp1<p).and.(lp2>p)).or.((lp1>p).and.(lp2<p))) then
                if((dabs(lp1)<1.0e200_DP).and.(dabs(lp2)<1.0e200_DP))then
                !check g is real
                g_local = (A_ideal( ) + A_mono( ) + A_chain( ) +  A_assoc( ) + A_ion( )) * NA * KB * t + lv2 * lp2 
                if(g_local==g_local) then
                    !interpolate                   
                    grad_local = (lp2 - lp1) / (lv2 - lv1)  
                    v = (p - lp1) / grad_local + lv1
                
                    !Converge the volume   
                    piter=0
                                                       
                    do              
                        v_temp2 = v                                  
                        p_temp2 = Press( )
    
                        if(piter<=plimit1) then                  
                            if(dabs((p_temp2 - p) / p)<P_crit) go to 20
                        else
                            if(dabs((p_temp2 - p) / p)<P_crit2) go to 20             
                        end if
                        
                        d_v = 1.0e-6_DP * v
                        
                        v_temp1 = v - d_v 
                        v = v_temp1
                        p_temp1 = Press( )
                          
                        !interpolate v
                        v = (p-p_temp1) * (v_temp2-v_temp1) / (p_temp2-p_temp1) + v_temp1
    
                        !escape if spurious volume
                        if(v<0.0e0_DP) then 
                            go to 10
                        else if(v/=v) then
                            go to 10
                        end if
                        
                        !escape if too many iterations
                        if(piter>plimit2) go to 20
                        piter = piter + 1                   
                    end do
                    
 20                 nroots = nroots+1
                    lv_array(nroots) = v
                    !lg_array(nroots) = (A_ideal( ) + A_mono( ) + A_chain( ) +  A_assoc( ) + A_ion( )) * NA * KB * t + lv2 * lp2     
                    !better g calc:
                    lg_array(nroots) = 0      
                    do imusum=1,nctypes
                        lg_array(nroots) = lg_array(nroots) + Comp_array(imusum)%xi * Mu(imusum)
                    end do               
                    
                end if
                end if
                end if
            end if
         
 10         lv1 = lv2
            lp1 = lp2
            
            i_rho = i_rho + 1
            
            if(present(limited)) then
                if (i_rho>i_rho_limit) then
                    limited=.false.
                    go to 30
                end if
            end if
            
            if((i_rho>1000).and.(lp2>5e8_DP)) then
                go to 30
            end if
        end do

        !determine most stable volume
 30     if(nroots<1) stop "no volume roots?"
        froot=1        

        if(nroots>1) then
            do i_rho = 2, nroots
                if(lg_array(i_rho)<lg_array(froot)) froot = i_rho
            end do
        end if
       
        v_out = lv_array(froot)

 99     return
    end function Vol_dens_g
!***************************************************************************************************
    function Vol_dens_g2(v_in1, v_in2) result(v_out)  
        implicit none
        
        real(kind=DP), intent(in)   ::  v_in1, v_in2
        real(kind=DP)       ::  r_low, r_high
        real(kind=DP)       ::  v_out
        real(kind=DP)       ::  lv_array(1:10), lg_array(1:10)
        real(kind=DP)       ::  v_temp1, v_temp2, p_temp1, p_temp2, d_v
        real(kind=DP)       ::  lv1, lv2, lp1, lp2, lr, g_local, grad_local
        integer             ::  i_rho, nroots, piter, froot, imusum
        integer, parameter          :: plimit1=100, plimit2=500
        real(kind=DP), parameter    ::  P_crit=1.0e-6_DP,P_crit2=1.0e-4

        nroots = 0
        
        if(1.0e0_DP/v_in1>1.0e0_DP/v_in2) then
            r_high = 2.0e0_DP/v_in1
            r_low  = 0.5e0_DP/v_in2   
        else 
            r_high = 2.0e0_DP/v_in2
            r_low  = 0.5e0_DP/v_in1
        end if
            
        lr  = r_low
        i_rho = 1
                  
        do      
            if(i_rho==1) then
                continue
            else if(i_rho<100) then
                lr  = lr+1.0e-4_DP 
            else if(i_rho<200) then  
                lr  = lr+1.0e-3_DP 
            else if(i_rho<300) then  
                lr  = lr+1.0e-2_DP   
            else if(i_rho<400) then  
                lr  = lr+1.0e-1_DP 
            else if(i_rho<500) then  
                lr  = lr+1.0e0_DP  
            else if(i_rho<1000) then  
                lr  = lr+5.0e0_DP 
            else if(i_rho<1500) then  
                lr  = lr+10.0e0_DP 
            else if(i_rho<2000) then  
                lr  = lr+50.0e0_DP 
            else if(i_rho<2500) then  
                lr  = lr+100.0e0_DP 
            else                             
                go to 30
            end if
            
            if(lr>r_high) go to 30
                            
            lv2 = 1.0e0_DP/lr
            if(lv2<=1.0e-5_DP) go to 30
            v   = lv2   
            lp2 = Press()    

            if(i_rho>0) then

                !sensible p and within desired p
                if(((lp1<p).and.(lp2>p)).or.((lp1>p).and.(lp2<p))) then
                if((dabs(lp1)<1.0e200_DP).and.(dabs(lp2)<1.0e200_DP))then
                !check g is real
                g_local = (A_ideal( ) + A_mono( ) + A_chain( ) +  A_assoc( ) + A_ion( )) * NA * KB * t + lv2 * lp2 
                if(g_local==g_local) then
                    !interpolate                   
                    grad_local = (lp2 - lp1) / (lv2 - lv1)  
                    v = (p - lp1) / grad_local + lv1
                
                    !Converge the volume   
                    piter=0
                                                       
                    do              
                        v_temp2 = v                                  
                        p_temp2 = Press( )
    
                        if(piter<=plimit1) then                  
                            if(dabs((p_temp2 - p) / p)<P_crit) go to 20
                        else
                            if(dabs((p_temp2 - p) / p)<P_crit2) go to 20             
                        end if
                        
                        d_v = 1.0e-6_DP * v
                        
                        v_temp1 = v - d_v 
                        v = v_temp1
                        p_temp1 = Press( )
                          
                        !interpolate v
                        v = (p-p_temp1) * (v_temp2-v_temp1) / (p_temp2-p_temp1) + v_temp1
    
                        !escape if spurious volume
                        if(v<0.0e0_DP) then 
                            go to 10
                        else if(v/=v) then
                            go to 10
                        end if
                        
                        !escape if too many iterations
                        if(piter>plimit2) go to 20
                        piter = piter + 1                   
                    end do
                    
 20                 nroots = nroots+1
                    lv_array(nroots) = v
                    !lg_array(nroots) = (A_ideal( ) + A_mono( ) + A_chain( ) +  A_assoc( ) + A_ion( )) * NA * KB * t + lv2 * lp2     
                    !better g calc:
                    lg_array(nroots) = 0      
                    do imusum=1,nctypes
                        lg_array(nroots) = lg_array(nroots) + Comp_array(imusum)%xi * Mu(imusum)
                    end do               
                    
                end if
                end if
                end if
            end if
         
 10         lv1 = lv2
            lp1 = lp2
            
            i_rho = i_rho + 1
        end do

        !determine most stable volume
 30     if(nroots<1) then
            v_out=1.e10
        else
            froot=1
            
            if(nroots>1) then
                do i_rho = 2, nroots
                    if(lg_array(i_rho)<lg_array(froot)) froot = i_rho
                end do
            end if
            
            v_out = lv_array(froot)
        end if

        return
    end function Vol_dens_g2
!***************************************************************************************************
!***************************************************************************************************
end module Vol_mod
!***************************************************************************************************
!***************************************************************************************************
!
!   Module to add a quote - please add more as and when used
!
!***************************************************************************************************
module Quote_mod
!***************************************************************************************************
contains
!***************************************************************************************************
    subroutine Quote( )
!***************************************************************************************************
!Modules
!=======
        use Zig_mod 
!***************************************************************************************************
    !Variables
    !=========
        implicit none
        
        integer     ::  n_quote, i_quote
        integer     ::  time(1:8), seed
!***************************************************************************************************    
        n_quote = 108  
        time(:)=0        
        
        call DATE_AND_TIME(values=time)     ! Get the current time
        seed = time(1) * (360000*time(5) + 6000*time(6) + 100*time(7) + time(8)) 
  
        call zigset(seed)
        i_quote = anint((n_quote-1)*uni())+1 

        if(i_quote==1) then
            print*, "I am Guybrush Threepwood, mighty pirate. (Monkey Island)"
        else if(i_quote==2) then    
            print*, "You fight like a dairy farmer. (Monkey Island)"
        else if(i_quote==3) then    
            print*, "Look behind you, a three-headed monkey!. (Monkey Island)"
        else if(i_quote==4) then    
            print*, "In my experience, there is no such thing as luck.  "
            print*, " (Obi-Wan Kenobi, Star Wars Episode IV: A New Hope)"
        else if(i_quote==5) then    
            print*, "The needs of the many outweigh the needs of the few. "
            print*, " (Spock, Star Trek II: The Wrath of Khan)"
        else if(i_quote==6) then    
            print*, "You are what you do. A man is defined by his actions, not his memory. "
            print*, " (Kuato, Total Recall)"
        else if(i_quote==7) then    
            print*, "Do what I do. Hold tight and pretend it's a plan! (The Doctor, Doctor Who)"
        else if(i_quote==8) then    
            print*, "Don't panic! (The Hitchhiker's Guide to the Galaxy)"
        else if(i_quote==9) then    
            print*, "Frak! (Battlestar Galactica)"
        else if(i_quote==10) then    
            print*, "Yeah, I got a plan B: making plan A work! (Stingray)"
        else if(i_quote==11) then    
            print*, "Why don't you smegging well smeg off you annoying little smeggy smegging smegger! "
            print*, " (Red Dwarf)"
        else if(i_quote==12) then    
            print*, "Make it so. (Star Trek TNG)"
        else if(i_quote==13) then    
            print*, "Hokey religions and ancient weapons are no match for a good blaster at your side, kid. "
            print*, " (Hans Solo, Star Wars Episode IV: A New Hope)"
        else if(i_quote==14) then    
            print*, "Roads? Where we're going, we don't need roads. (Back to the Future)"
        else if(i_quote==15) then    
            print*, "Come with me if you want to live. (The Terminator)"
        else if(i_quote==16) then    
            print*, "We'd better get back, 'cause it'll be dark soon, and they mostly come at night... mostly. "
            print*, " (Aliens)"
        else if(i_quote==17) then    
            print*, "Dead or alive, you're coming with me! (Robocop)"
        else if(i_quote==18) then    
            print*, "I mean, have you ever actually seen your brain? (Ghost in the Shell)"
        else if(i_quote==19) then    
            print*, "We spend $250 billion a year on defence, and here we are! The fate of the planet "
            print*, " is in the hands of a bunch of retards I wouldn't trust with a potato gun! (Armageddon)"
        else if(i_quote==20) then    
            print*, "Those who know do not speak. Those who speak do not know. (Lao Tzu, Tao Teh Ching)" 
        else if(i_quote==21) then    
            print*, "The truth is not always beautiful, nor beautiful words the truth. (Lao Tzu, Tao Teh Ching)" 
        else if(i_quote==22) then    
            print*, "When I let go of what I am, I become what I might be. (Lao Tzu, Tao Teh Ching)" 
        else if(i_quote==23) then    
            print*, "Care about what other people think and you will always be their prisoner. "
            print*, " (Lao Tzu, Tao Teh Ching)" 
        else if(i_quote==24) then    
            print*, "A man with outward courage dares to die; a man with inner courage dares to live. "
            print*, " (Lao Tzu, Tao Teh Ching)" 
        else if(i_quote==25) then    
            print*, "Teenage angst has paid off well, Now I'm bored and old. (Nirvana)"
        else if(i_quote==26) then    
            print*, "Hey! Wait! I got a new complaint. (Nirvana)"
        else if(i_quote==27) then    
            print*, "I don't care what you think unless it's about me. (Nirvana)"
        else if(i_quote==28) then    
            print*, "Just because you're paranoid. Don't mean they're now after you. (Nirvana)"
        else if(i_quote==29) then    
            print*, "I'm so happy cause today I've found my friends. They're in my head. (Nirvana)"
        else if(i_quote==30) then    
            print*, "Please allow me to introduce myself / I'm a man of wealth and taste/ I've been"
            print*, " around for a long, long year / Stole many a man's soul and faith. (The Rolling Stones)"
        else if(i_quote==31) then    
            print*, "There's no time for us. There's no place for us. What is this thing that builds our "
            print*, " dreams yet slips away from us? (Queen)"
        else if(i_quote==32) then    
            print*, "And when the brokenhearted people living in the world agree, there will be an answer,"
            print*, " let it be. (The Beatles)"
        else if(i_quote==33) then    
            print*, "The fewer the facts, the stronger the opinion. (Arnold H. Glasow)"
        else if(i_quote==34) then    
            print*, "The distance between insanity and genius is measured only by success. (Bruce Feirstein)"
        else if(i_quote==35) then    
            print*, "Your theory is crazy, but it's not crazy enough to be true. (Niels Bohr)"
        else if(i_quote==36) then    
            print*, "Science is nothing but perception. (Plato)"
        else if(i_quote==37) then    
            print*, "If your experiment needs statistics, you ought to have done a better experiment. "
            print*, " (Ernest Rutherford)"
        else if(i_quote==38) then    
            print*, "A physicist is an atom's way of knowing about atoms. (George Wald)"
        else if(i_quote==39) then    
            print*, "Valid criticism does you a favour. (Carl Sagan)"
        else if(i_quote==40) then    
            print*, "Chemistry, unlike other sciences, sprang originally from delusions and superstitions, "
            print*, " and was at its commencement exactly on a par with magic and astrology. (Thomas Thomson)"
        else if(i_quote==41) then    
            print*, "To the electron -- may it never be of any use to anybody. (JJ. Thomson)"
        else if(i_quote==42) then    
            print*, "After climbing a great hill, one only finds there are many more hills to climb. "
            print*, " (Nelson Mandela)"
        else if(i_quote==43) then    
            print*, "If you cannot do great things, do small things in a great way. (Napoleon Hill)"
        else if(i_quote==44) then    
            print*, "It is a rough road that leads to the heights of greatness. (Lucius Annaeus Seneca)"
        else if(i_quote==45) then    
            print*, "Only those who dare to fail greatly can ever achieve greatly. (Robert Kennedy)"
        else if(i_quote==46) then    
            print*, "It is hard to be humble when you're as great as I am. (Muhammad Ali)"
        else if(i_quote==47) then    
            print*, "Two roads diverged in a wood, and I - I took the one less traveled by, and that "
            print*, " has made all the difference. (Robert Frost)"
        else if(i_quote==48) then    
            print*, "Life is what happens to you while you're busy making other plans. (John Lennon)"
        else if(i_quote==49) then    
            print*, "The most common way people give up their power is by thinking they don't have any."
            print*, " (Alice Walker)"
        else if(i_quote==50) then    
            print*, "The best time to plant a tree was 20 years ago. The second best time is now. "
            print*, " (Chinese Proverb)"
        else if(i_quote==51) then    
            print*, "Every child is an artist.  The problem is how to remain an artist once he grows up. "
            print*, " (Pablo Picasso)"
        else if(i_quote==52) then    
            print*, "Whether you think you can or you think you can't, you're right. (Henry Ford)"
        else if(i_quote==53) then    
            print*, "People often say that motivation doesn't last. Well, neither does bathing.  "
            print*, " That's why we recommend it daily. (Zig Ziglar)"
        else if(i_quote==54) then    
            print*, "There is only one way to avoid criticism: do nothing, say nothing, and be nothing. "
            print*, " (Aristotle)"
        else if(i_quote==55) then    
            print*, "Everything you've ever wanted is on the other side of fear. (George Addair)"
        else if(i_quote==56) then    
            print*, "Teach thy tongue to say, I do not know, and thou shalt progress. (Maimonides)"
        else if(i_quote==57) then    
            print*, "How wonderful it is that nobody need wait a single moment before starting to "
            print*, " improve the world. (Anne Frank)"
        else if(i_quote==58) then    
            print*, "If the wind will not serve, take to the oars. (Latin Proverb)"
        else if(i_quote==59) then    
            print*, "Challenges are what make life interesting and overcoming them is what makes life"
            print*, " meaningful. (Joshua J. Marine)"
        else if(i_quote==60) then    
            print*, "The person who says it cannot be done should not interrupt the person who is "
            print*, " doing it. (Chinese Proverb)"
        else if(i_quote==61) then    
            print*, "Build your own dreams, or someone else will hire you to build theirs. (Farrah Gray)"
        else if(i_quote==62) then    
            print*, "It does not matter how slowly you go as long as you do not stop. (Confucius)"
        else if(i_quote==63) then    
            print*, "Everything in moderation, including moderation. (Oscar Wilde)"
        else if(i_quote==64) then    
            print*, "Just taught my kids about taxes by eating 38% of their ice cream. (Conan O'Brien)"
        else if(i_quote==65) then    
            print*, "I dream of a better tomorrow, where chickens can cross the road and not be questioned"
            print*, " about their motives. (Anonymous)"
        else if(i_quote==66) then    
            print*, "It is the mark of an educated mind to be able to entertain a thought without "
            print*, " accepting it. (Aristotle)"
        else if(i_quote==67) then    
            print*, "You see things and you say Why? But I dream things that never were and I say "
            print*, " Why not? (George Bernard Shaw)"
        else if(i_quote==68) then    
            print*, "I have never in my life learned anything from a man who agreed with me. "
            print*, " (Dudley Field Malone)"
        else if(i_quote==69) then    
            print*, "We have all heard that a million monkeys banging on a million typewriters will "
            print*, " eventually reproduce the entire works of Shakespeare. Now, thanks to the internet, "
            print*, " we know this is not true. (Robert Wilensky)"
        else if(i_quote==70) then    
            print*, "Every truth passes through three stages before it is recognised. In the first,"
            print*, " it is ridiculed. In the second, it is opposed. In the third it is regarded as "
            print*, " self-evident. (Arthur Schopenhauer)"
        else if(i_quote==71) then    
            print*, "I always arrive late at the office, but I make up for it by leaving early. "
            print*, " (Charles Lamb)"
        else if(i_quote==72) then    
            print*, "The elevator to success is out of order. You'll have to use the stairs... "
            print*, " one step at a time. (Joe Girard)"
        else if(i_quote==73) then    
            print*, "The greatest way to live with honour in this world is to be what we"
            print*, " pretend to be. (Socrates)"
        else if(i_quote==74) then    
            print*, "Laws control the lesser man... Right conduct controls the greater one. (Mark Twain)"
        else if(i_quote==75) then    
            print*, "Great men are seldom over-scrupulous in the arrangement of their attire. "
            print*, " (Charles Dickens)"
        else if(i_quote==76) then    
            print*, "Time's fun when you're having flies. (Kermit the Frog )" 
        else if(i_quote==77) then    
            print*, "You have to be odd to be number one. (Dr. Seuss)"
        else if(i_quote==78) then    
            print*, "You have brains in your head. You have feet in your shoes. You can steer "
            print*, " yourself any direction you choose. (Dr. Seuss)"
        else if(i_quote==79) then    
            print*, "Sometimes the questions are complicated and the answers are simple. (Dr. Seuss)"
        else if(i_quote==80) then    
            print*, "Today you are you. That is truer than true. There is no one alive who is "
            print*, " youer than you. (Dr. Seuss)"
        else if(i_quote==81) then    
            print*, "A little nonsense now and then, is relished by the wisest men. (Roald Dahl)"
        else if(i_quote==82) then    
            print*, "In ancient times cats were worshipped as gods. They have not forgotten this. "
            print*, "(Terry Pratchett)"
        else if(i_quote==83) then    
            print*, "The trouble with having an open mind, of course, is that people will insist on"
            print*, " coming along and trying to put things in it. (Terry Pratchett)"
        else if(i_quote==84) then    
            print*, "The pen is mightier than the sword if the sword is very short, and the pen is "
            print*, " very sharp. (Terry Pratchett)"
        else if(i_quote==85) then    
            print*, "The presence of those seeking the truth is infinitely to be preferred to the "
            print*, " presence of those who think they've found it. (Terry Pratchett)"
        else if(i_quote==86) then      
            print*, " We have created amazing machines... and used them to destroy our fellow men. "
            print*, " (Metropolis)"
        else if(i_quote==87) then      
            print*, "You're mad. Bonkers. Off your head... but I'll tell you a secret... "
            print*, "all of the best people are. (Alice in Wonderland)"
        else if(i_quote==88) then      
            print*, "Don't cry because it's over, smile because it happened. (Dr. Seuss)"
        else if(i_quote==89) then      
            print*, "We are all in the gutter, but some of us are looking at the stars."
            print*, "(Oscar Wilde)"
        else if(i_quote==90) then      
            print*, "I may not have gone where I intended to go, but I think I have ended"
            print*, "up where I needed to be. (Douglas Adams)"
        else if(i_quote==91) then      
            print*, "For every minute you are angry you lose sixty seconds of happiness."
            print*, "(Ralph Walso Emerson)"
        else if(i_quote==92) then      
            print*, "If you can't explain it to a six year old, you don't understand it "
            print*, "yourself. (Albert Einstein)"
        else if(i_quote==93) then      
            print*, "Success is not final, failure is not fatal: it is the courage to "
            print*, "continue that counts. (Winston Churchill)"
        else if(i_quote==94) then      
            print*, "Reality continues to ruin my life. (Bill Watterson)"
        else if(i_quote==95) then      
            print*, "I am free of all prejudice. I hate everyone equally. (W.C. Fields)"
        else if(i_quote==96) then      
            print*, "It's no use going back to yesterday, because I was a different"
            print*, "person then. (Lewis Carroll)"
        else if(i_quote==97) then      
            print*, "What you're supposed to do when you don't like a thing is change"
            print*, "it. If you can't change it, change the way you think about it. "
            print*, "Don't complain. (Maya Angelou)"
        else if(i_quote==98) then 
            print*, "Experience is that marvelous thing that enables you to recognize "
            print*, "a mistake when you make it again. (F. P. Jones)"
        else if(i_quote==99) then 
            print*, "Human beings, who are almost unique in having the ability to learn "
            print*, "from the experience of others, are also remarkable for their apparent "
            print*, "disinclination to do so. (Douglas Adams)"
        else if(i_quote==100) then
            print*, "My opinions may have changed, but not the fact that I am right. (Ashleigh Brilliant)"
        else if(i_quote==101) then
            print*, "Every man is born as many men and dies as a single one. (Martin Heidegger)"
        else if(i_quote==102) then
            print*, "Many are stubborn in pursuit of the path they have chosen,"
            print*, "few in pursuit of the goal. (Friedrich Nietzsche)"
        else if(i_quote==103) then
            print*, "If you are lonely when you're alone, you are in bad company. (Jean-Paul Sartre)"
        else if(i_quote==104) then
            print*, "There are two things a person should never be angry at, what they can help,"
            print*, "and what they cannot. (Plato)"
        else if(i_quote==105) then
            print*, "No tree, it is said, can grow to heaven unless its roots reach down to hell. (Carl Jung)"
        else if(i_quote==106) then
            print*, "Truth often suffers more by the heat of its defenders than "
            print*, "the arguments of its opposers. (William Penn)"
        else if(i_quote==107) then
            print*, "You can't depend on your eyes when your imagination is out of focus. (Mark Twain)"
        else if(i_quote==108) then
            print*, "Sometimes you face difficulties not because you're doing something wrong,"
            print*, "but because you're doing something right. (Joel Osteen)"
        end if    
        
        return
    end subroutine Quote
end module Quote_mod

!***************************************************************************************************

    Module c05qbfe_mod

!     C05QBF Example Program Module:
!            Parameters and User-defined Routines

!     .. Use Statements ..
      Use nag_library, Only: nag_wp
	  use Types_mod           ! Definitions of types and double precision
	  use Global_mod          ! Important global parameters
	  use Press_mod
	  use Mu_mod
    use Vol_mod
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
	  use Types_mod           ! Definitions of types and double precision
	  use Global_mod          ! Important global parameters
	  use Press_mod
	  use Input_mod
	  use Input_opt_mod   ! Read optimisation parameters
	  use Vol_mod
	  use c05qbfe_mod
	  use Mu_mod
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
101 Format (/,(2X,A2),(2X,F6.4),(2X,A5),(2X,F6.4),(2X,A6),(2X,F6.4),(2X,A6),(2X,F6.4),(2X,A5),(2X,F10.4),$)   

        Return

      End Subroutine objfun
      Subroutine monfun(n,nf,x,f,rho,iuser,ruser,inform)

!       .. Scalar Arguments ..
		!use Types_mod       ! Definitions of types and double precision
		!use Global_mod,only: opt_num, min_num, max_num, init_values      ! Important global parameters
		
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
	  !use Types_mod       ! Definitions of types and double precision
	  use Global_mod,only: opt_num, min_num, max_num, init_values, param_key, ANG, NA       ! Important global parameters
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
    use Types_mod       ! Definitions of types and double precision
    use Input_mod       ! Read the input
    use Input_opt_mod   ! Read optimisation parameters
    use Global_mod      ! Important global parameters
    !use Simplex_mod     ! Simplex optimiser
    use Quote_mod       ! Random quote
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
    
    call Quote( )
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
