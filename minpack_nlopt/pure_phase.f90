Module pure_phase
        ! Module to solve series of non-linear equations
        ! Replaces c05qbf from NAG library

        use Types           ! Definitions of types and double precision
        use Global          ! Important global parameters
        Use Pressure
        use ChemPot
        use Vol

        Implicit None

        Private
        Public                           :: fcn_pha1, solve_pha1

        Integer, Parameter, Public       :: n_pha1 = 2, nout_pha1 = 6
        !Real (Kind=DP), Public           :: P_L, P_V   !move to global 
        Contains

        subroutine solve_pha1(output_pha1)
                ! Replaces C05QBF from NAG library

                Use types, Only: DP
                Use minpack, Only: dpmpar, enorm, hybrd1

                Implicit None

                Real (Kind=DP)               :: fnorm, tol, output_pha1(2)
                Integer                      :: i, info, lwa=(n_pha1*(3*n_pha1+13))/2

                Real (Kind=DP), Allocatable  :: fvec(:), x(:)
                Real (Kind=DP), allocatable  :: wa(:)

                Intrinsic                        :: sqrt


                                
                Allocate (fvec(n_pha1),x(n_pha1))
                allocate(wa(lwa))
                


                !     The following starting values provide a rough solution.
                do i=1, n_pha1
                  x(i) = x_init(i)
                end do

                tol = sqrt(dpmpar(1))

                Call hybrd1(fcn_pha1,n_pha1,x,fvec,tol,info,wa,lwa)

                If (info > 0) Then
                        If (info == 1) Then
                          fnorm = enorm(n_pha1,fvec)
                        Else
                          Write (nout_pha1,*)
                          Write (nout_pha1,*) 'Approximate solution'
                        End If

                        Do i = 1, n_pha1
                           output_pha1(i) = x(i)
                        End Do
                else if (info == 0) then
                  write (nout_pha1,*) 'Invalid input arguments to hybrd1!'
                  stop
                End If

                return
        end subroutine solve_pha1

        Subroutine fcn_pha1(n,x,fvec,iflag)
                ! User subroutine which returns the values of the functions at a point, x.
                Integer, Intent (Inout)        :: iflag
                Integer, Intent (In)           :: n

                Real (Kind=DP), Intent (Out)   :: fvec(n)
                Real (Kind=DP), Intent (In)    :: x(n)
                Integer                        :: i
                Real (Kind=DP)                 :: x_sum

                

                Comp_array(1)%xi = 1.0e0_DP !x1
                v   = x(1)				! Liquid volume is x1
                P_L = Press()
                Mu_L(1) = Mu(1)
                v   = x(2)				! Vapor volume is x2
                P_V = Press()
                Mu_V(1) = Mu(1)
                
                fvec(1) = P_L - P_V
                fvec(2) = Mu_L(1) - Mu_V(1)
                
                ! Set iflag negative to terminate execution for any reason.
                iflag = 0
                Return
        End Subroutine fcn_pha1
End Module pure_phase
