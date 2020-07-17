Module solver
        ! Module to solve series of non-linear equations
        ! Replaces c05qbf from NAG library

        use Types           ! Definitions of types and double precision
        use Global          ! Important global parameters
        Use Pressure
        use ChemPot
        use Vol

        Implicit None

        Private
        Public                           :: fcn, solve_nle

        Integer, Parameter, Public       :: nout = 6
        Integer, Public                  :: n
        !Real (Kind=DP), Public           :: P_L, P_V   !move to global 
        Contains

        subroutine solve_nle(output)
                ! Replaces C05QBF from NAG library

                Use types, Only: DP
                Use minpack, Only: dpmpar, hybrd1

                Implicit None

                Real (Kind=DP)               :: tol
                Integer                      :: i, info, lwa

                Real (Kind=DP), Allocatable  :: fvec(:), x(:)
                Real (Kind=DP), allocatable  :: wa(:), output(:)

                Intrinsic                        :: sqrt
                !Real(Kind=DP)                :: output(3)

                
                n = nctypes + 1                          ! yichun
                lwa = (n*(3*n+13))/2                     ! yichun
                                
                Allocate (fvec(n),x(n))
                allocate(wa(lwa))
                


                !     The following starting values provide a rough solution.
                do i=1, n
                  x(i) = x_init(i)                  !set initial guesses
                end do
                tol = sqrt(dpmpar(1))

                Call hybrd1(fcn,n,x,fvec,tol,info,wa,lwa)

                if (info == 0) then
                  write (nout,*) 'Invalid input arguments to hybrd1!'
                  stop
                else if (info > 0) Then
                        Do i = 1, n
                           output(i) = x(i)
                        End Do

                        If (info > 1) Then
                          Write (nout,*)
                          Write (nout,*) 'Approximate solution'
                        End If

                End If

                return
        end subroutine solve_nle

        Subroutine fcn(n,x,fvec,iflag)
                ! User subroutine which returns the values of the functions at a point, x.
                Integer, Intent (Inout)        :: iflag
                Integer, Intent (In)           :: n

                Real (Kind=DP), Intent (Out)   :: fvec(n)
                Real (Kind=DP), Intent (In)    :: x(n)
                Integer                        :: i
                Real (Kind=DP)                 :: x_sum

                

                Comp_array(1:nctypes)%xi = comp_x(1:nctypes)
                v   = x(1)                      ! Liquid volume is x1
                p = Press()
                P_L = p
                do i=1, nctypes
                  Mu_L(i) = Mu(i)
                end do
                v   = x(2)                      ! Vapor volume is x2
                x_sum = 0
                do i=1, nctypes-1
                  Comp_array(i)%xi = x(i+2)  ! y1-y(nctypes-1)
                  x_sum = x_sum + Comp_array(i)%xi
                end do
                Comp_array(nctypes)%xi = 1- x_sum  !the last component              
                p = Press()
                P_V = p
                do i=1, nctypes
                  Mu_V(i) = Mu(i)
                end do

                do i=1, n-1
                  fvec(i) = Mu_L(i) - Mu_V(i)
                end do
                fvec(n) = P_L - P_V
                
                ! Set iflag negative to terminate execution for any reason.
                iflag = 0
                Return
        End Subroutine fcn
End Module solver
