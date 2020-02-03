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

        Integer, Parameter, Public       :: n = 3, nout = 6
        Real (Kind=DP), Public           :: P_L, P_V, Mu_L_1, Mu_L_2, Mu_V_1, Mu_V_2, x_1, x_init_1, x_init_2, x_init_3
        Contains

        subroutine solve_nle(output)
                ! Replaces C05QBF from NAG library

                Use types, Only: DP
                Use minpack, Only: dpmpar, enorm, hybrd1

                Implicit None

                Real (Kind=DP)               :: fnorm, tol
                Integer                      :: i, info, lwa = (n*(3*n+13))/2

                Real (Kind=DP), Allocatable  :: fvec(:), x(:)
                Real (Kind=DP), allocatable  :: wa(:)

                Intrinsic                        :: sqrt
                Real(Kind=DP)                :: output(3)

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
        end subroutine solve_nle

        Subroutine fcn(n,x,fvec,iflag)
                ! User subroutine which returns the values of the functions at a point, x.
                Integer, Intent (Inout)        :: iflag
                Integer, Intent (In)           :: n

                Real (Kind=DP), Intent (Out)   :: fvec(n)
                Real (Kind=DP), Intent (In)    :: x(n)

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
                ! Set iflag negative to terminate execution for any reason.
                iflag = 0
                Return
        End Subroutine fcn
End Module solver
