!***************************************************************************************************
!   SAFT Module to:
!       1. Define double precision
!       2. Define types
!***************************************************************************************************
!
!***************************************************************************************************
module Types
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
                                            &   xiv(:,:), xih(:,:), xic(:,:), xip(:,:),ximu(:,:), &
                                            &    t_liqv(:), p_liqv(:), t_liqh(:), p_liqh(:), t_liqc(:), p_liqc(:)      
        character(len=1),allocatable        ::  phase(:)
        character(len=1)                    ::  opt(1:5)
        logical                             ::  opt_l(1:5)
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
end module Types
!***************************************************************************************************