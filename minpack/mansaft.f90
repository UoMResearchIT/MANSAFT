!***************************************************************************************************
!***************************************************************************************************

!***************************************************************************************************

    Module c05qbfe_mod

!     C05QBF Example Program Module:
!            Parameters and User-defined Routines

!     .. Use Statements ..
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
      Integer, Parameter, Public       :: n = 3, nout = 6
	  Real (Kind=DP), Public       :: P_L, P_V, Mu_L_1, Mu_L_2, Mu_V_1, Mu_V_2, x_1, x_init_1, x_init_2, x_init_3
    Contains
      Subroutine fcn(n,x,fvec,iflag)

!       .. Scalar Arguments ..
        Integer, Intent (Inout)        :: iflag
        Integer, Intent (In)           :: n
!       .. Array Arguments ..
        Real (Kind=DP), Intent (Out) :: fvec(n)
        Real (Kind=DP), Intent (In) :: x(n)
!       .. Executable Statements ..
        v   = x(1)                      ! Liquid volume is x1
        Comp_array(1)%xi = x_1 !x1
        Comp_array(2)%xi = 1.0_DP - x_1 !x2
        p = Press()
        P_L = p
        Mu_L_1 = Mu(1)
        Mu_L_2 = Mu(2)
        v   = x(2)                      ! Vapor volume is x2
        Comp_array(1)%xi = x(3) ! y1
        Comp_array(2)%xi =1.0_DP - x(3) !y2
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
subroutine c05qbfe(output)

!     C05QBF Example Main Program

!     .. Use Statements ..
      Use c05qbfe_mod, Only: fcn, n, nout, x_init_1, x_init_2, x_init_3
      Use types_mod, Only: DP

      Use minpack, Only: dpmpar, enorm, hybrd1

!     .. Implicit None Statement ..
      Implicit None

!     .. Local Scalars ..
      Real (Kind=DP)               :: fnorm, tol
      Integer                      :: i, info, lwa = (n*(3*n+13))/2
!     .. Local Arrays ..
      Real (Kind=DP), Allocatable  :: fvec(:), x(:)
      Real (Kind=DP), allocatable  :: wa(:)
      !     .. Intrinsic Procedures ..
      Intrinsic                        :: sqrt
      Real(Kind=DP)                :: output(3)
!     .. Executable Statements ..

      Allocate (fvec(n),x(n))
      allocate(wa(lwa))

!     The following starting values provide a rough solution.

      x(1) = x_init_1
      x(2) = x_init_2
      x(3) = x_init_3 
      
      tol = sqrt(dpmpar(1))

      Call hybrd1(fcn,n,x,fvec,tol,info,wa,lwa)

      If (info > 0) Then
        If (info == 1) Then
          fnorm = enorm(n,fvec)
        Else
          Write (nout,*)
          Write (nout,*) 'Approximate solution'
        End If
        Do i = 1, n
           output(i) = x(i)
        End Do
      else if (info == 0) then
          write (nout,*) 'Invalid input arguments to hybrd1!'
          stop
      End If
  return
End
!***************************************************************************************************
program pha1
    use Types_mod           ! Definitions of types and double precision
    use Global_mod          ! Important global parameters
    use Press_mod
    use Input_mod
    use Vol_mod
    use Mu_mod
    use c05qbfe_mod
      implicit none
      integer             ::  i
      real(kind=DP)      :: values(3), mu_1, mu_2
      call Read_input( )
        select case(properties%type)
          case('vle','VLE')

                print*, "V-L Phase equilibrium"
                print*, "====================="
                print*, " "
                print*, "  T(K)     P(kPa)    x1     y1     V_L          V_V        "  
                print*, " "
                
    do i = 1, properties%n                    
          t = properties%t(i)
		  p = properties%p(i)*1000
          x_1 = properties%xi(i, 1)
		  Comp_array(1)%xi = x_1
		  Comp_array(2)%xi = 1.0 - x_1
		  P = 1.5*P
		  v = Vol_dens_g( )
		  mu_1 = Mu_res(1)
		  mu_2 = Mu_res(2)
		  x_init_1 = v
		  x_init_2 = 8.314*t/(properties%p(i)*1000)
		  x_init_3 = exp(mu_1/(8.314*t))*x_1/(exp(mu_1/(8.314*t))*x_1+exp(mu_2/(8.314*t))*(1-x_1))
					! print*, mu_1, mu_2, x_init_1, x_init_2, x_init_3
					
          call c05qbfe(values )
          write(*,100) t, (P_L+P_V)/2000, x_1, values(3),values(1), values(2) 
    end do
          end select           
      stop
100 Format ('  ',F7.2,'  ',F8.2,'  ',F5.3,'  ',F5.3,'  'E11.4,'  ',E11.4)  
end
