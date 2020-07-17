!***************************************************************************************************
!   SAFT Module to:
!       1. Hold key global variables
!       2. Define constants
!***************************************************************************************************
!
!***************************************************************************************************
module Global
!***************************************************************************************************
!Modules
!=======
    use Types       ! Definitions of types and double precision
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
!Addition VLE
    real(kind=DP),allocatable   ::  x_init(:)      !initial guesses for non-linear solver
    real(kind=DP),allocatable   ::  comp_x(:)     !fix mole fractions in liquid phase in  non-linear solver
    real(kind=DP),allocatable   :: Mu_L(:), Mu_V(:)
    Real (Kind=DP)              :: P_L, P_V






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
END MODULE Global
!***************************************************************************************************