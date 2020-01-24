!***************************************************************************************************
!   SAFT Module to:
!       1. Calculate ionic activity coefficient
!       2. Temporarily calculate other activities...
!***************************************************************************************************
!
!***************************************************************************************************
module Act_ion_mod
!***************************************************************************************************
!Modules
!=======    
    use Types_mod           ! Definitions of types and double precision
    use Global_mod          ! Important global parameters
    use Setup_mod           ! To setup the system
    use Vol_mod             ! Calculate volume   
    use Mu_mod              ! Chemical potential  
!***************************************************************************************************
    implicit none
!***************************************************************************************************
    contains
!***************************************************************************************************              
    subroutine Ion_act(  ) 
        implicit none
        
        real(kind=DP)       ::  xi_old(1:nctypes)
        real(kind=DP)       ::  mu1, mu1inf, mu2, mu2inf
        real(kind=DP)       ::  vsys, vinf, act1, act2, actout
        real(kind=DP)       ::  phi1, phi2, phi1inf, phi2inf
        real(kind=DP)       ::  zcomp, zcompinf, tempxi(1:2)
!-----------------------------------------------------------------------------
        xi_old(1:nctypes) = Comp_array(1:nctypes)%xi
     
        vsys = Vol_dens_g()
        v    = vsys

        mu1  = Mu_res(properties%cation)
        mu2  = Mu_res(properties%anion)
        
        zcomp = vsys*p / (NA*KB*t)
     
        phi1 = 1.0e0_DP / zcomp * exp(mu1 / NA/KB/t )
        phi2 = 1.0e0_DP / zcomp * exp(mu2 / NA/KB/t )
!-----------------------------------------------------------------------------
!CURRENTLY ASSUME 1 SOLVENT COMPONENT ONLY, AND IT IS COMPONENT #1
!        Comp_array(1)%xi = 1.0e0_DP
        Comp_array(2)%xi = 1.0e-20_DP
        Comp_array(3)%xi = 1.0e-20_DP
        Comp_array(1)%xi = 1.0e0_DP - Comp_array(2)%xi - Comp_array(3)%xi
        
        ion_switch = .false.         
        vinf    = Vol_dens_g()
        v       = vinf
        ion_switch = .true.
       
        mu1inf  = Mu_res_inf(properties%cation)
        mu2inf  = Mu_res_inf(properties%anion)
     
        zcompinf = vinf*p / (NA*KB*t)
        
        phi1inf = 1.0e0_DP / zcompinf * exp(mu1inf / NA/KB/t )
        phi2inf = 1.0e0_DP / zcompinf * exp(mu2inf / NA/KB/t )
        
        Comp_array(1:nctypes)%xi = xi_old(1:nctypes)
!-----------------------------------------------------------------------------
!AGAIN ASSUMING SOLVENT IS 1 
        tempxi(1) = Comp_array(1)%xi / (Comp_array(1)%xi+Comp_array(2)%xi)
        tempxi(2) = 1.0e0_DP - tempxi(1)
        
        act1    = phi1 / phi1inf * tempxi(1)
        act2    = phi2 / phi2inf * tempxi(1)
        actout  = (act1**properties%cat_stoich * act2**properties%an_stoich)**&
        &   (1.0e0_DP/(properties%cat_stoich + properties%an_stoich)) 
        
        print*, tempxi(2)*1000.0e0_DP/(tempxi(1)*Comp_array(1)%m/AMU) , actout  

        return
    end subroutine Ion_act
!*************************************************************************************************** 



!*************************************************************************************************** 
    function Eopt_act(  ) result(act3) 
        implicit none
        
        real(kind=DP)       ::  xi_old(1:nctypes)
        real(kind=DP)       ::  mu1, mu1inf, mu2, mu2inf
        real(kind=DP)       ::  vsys, vinf, act1, act2, act3
!-----------------------------------------------------------------------------
        xi_old(1:nctypes) = Comp_array(1:nctypes)%xi
      
        vsys = Vol_dens_g()
        v    = vsys

        mu1  = Mu_res(2)
        mu2  = Mu_res(3)
!-----------------------------------------------------------------------------
!CURRENTLY ASSUME 1 SOLVENT COMPONENT ONLY, AND IT IS COMPONENT #1
        Comp_array(1)%xi = 1.0e0_DP
        Comp_array(2)%xi = 0.0e0_DP
        Comp_array(3)%xi = 0.0e0_DP

        ion_switch = .false.         
        vinf    = Vol_dens_g()
        v       = vinf
        ion_switch = .true.
       
        mu1inf  = Mu_res_inf(2)
        mu2inf  = Mu_res_inf(3)
     
        Comp_array(1:nctypes)%xi = xi_old(1:nctypes)
!-----------------------------------------------------------------------------

!AGAIN ASSUMING SOLVENT IS 1
        act1    = vinf / vsys * dexp( (mu1 - mu1inf) / NA / KB / t)
        act1    = act1 * Comp_array(1)%xi
        act2    = vinf / vsys * dexp( (mu2 - mu2inf) / NA / KB / t)
        act2    = act2 * Comp_array(1)%xi
 
!shortcut activity
        act3 = dsqrt(act1*act2)   

        return
    end function Eopt_act
!***************************************************************************************************    
end module Act_ion_mod
!***************************************************************************************************