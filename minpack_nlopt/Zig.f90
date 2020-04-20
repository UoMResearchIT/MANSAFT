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
module Zig
    use Types
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
end module Zig
