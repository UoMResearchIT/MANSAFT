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
