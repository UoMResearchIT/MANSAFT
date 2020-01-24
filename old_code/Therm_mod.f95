!***************************************************************************************************
!   SAFT Module to:
!       1. Calculate other thermodynamic properties (i.e. not P, V, A)
!           U, H, G, S
!
!
!       Created March 2017
!
!           Sept 2017 - added 6 individual them calculations
!***************************************************************************************************
!
!***************************************************************************************************
module Therm_mod
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
    use Vol_mod             ! Calculate V
    use Press_mod           ! Calculate P
    use Mu_mod              ! Calculate chemical potentials
!***************************************************************************************************
    implicit none
!***************************************************************************************************    
    real(kind=DP)   ::  a_ideal_l, a_mono_l, a_chain_l, a_assoc_l, a_ion_l, a_born_l
    real(kind=DP)   ::  u_ideal_l, u_mono_l, u_chain_l, u_assoc_l, u_ion_l, u_born_l
    real(kind=DP)   ::  h_ideal_l, h_mono_l, h_chain_l, h_assoc_l, h_ion_l, h_born_l
    real(kind=DP)   ::  g_ideal_l, g_mono_l, g_chain_l, g_assoc_l, g_ion_l, g_born_l
    real(kind=DP)   ::  s_ideal_l, s_mono_l, s_chain_l, s_assoc_l, s_ion_l, s_born_l
    
    real(kind=DP)   ::  a_ideal_l1, a_mono_l1, a_chain_l1, a_assoc_l1, a_ion_l1, a_born_l1
    real(kind=DP)   ::  a_ideal_l2, a_mono_l2, a_chain_l2, a_assoc_l2, a_ion_l2, a_born_l2
    real(kind=DP)   ::  a_ideal_l3, a_mono_l3, a_chain_l3, a_assoc_l3, a_ion_l3, a_born_l3
    real(kind=DP)   ::  a_ideal_l4, a_mono_l4, a_chain_l4, a_assoc_l4, a_ion_l4, a_born_l4
    
    real(kind=DP),parameter   ::  dt = 1.0e-6_DP
!***************************************************************************************************
contains

!***************************************************************************************************
    subroutine Thermall()
        implicit none

        !T and P are preset!
        real(kind=DP)   ::  pset, tset, vset, rhoset
        real(kind=DP)   ::  pval(1:4), tval(1:4), rhoval(1:4), tval2(1:9), aval2(1:9), aval(1:4)
        real(kind=DP)   ::  dpdr, dpdt, d2adt2
        real(kind=DP)   ::  isokt, thermalpha, cv, cp, jtcoeff, spsound
        real(kind=DP), parameter    ::  dval=1.0e-3_DP 


        pset   = p
        tset   = t
        vset   = Vol_dens_g()
        rhoset = 1.0e0_DP / vset
        
        tval(1) = tset * (1.0e0_DP - 2.0e0_DP * dval) 
        tval(2) = tset * (1.0e0_DP - dval)
        tval(3) = tset * (1.0e0_DP + dval)
        tval(4) = tset * (1.0e0_DP + 2.0e0_DP * dval)
        
        rhoval(1) = rhoset * (1.0e0_DP - 2.0e0_DP * dval)
        rhoval(2) = rhoset * (1.0e0_DP - dval)
        rhoval(3) = rhoset * (1.0e0_DP + dval)
        rhoval(4) = rhoset * (1.0e0_DP + 2.0e0_DP * dval)
        
        tval2(1) = tset * (1.0e0_DP - 4.0e0_DP * dval) 
        tval2(2) = tset * (1.0e0_DP - 3.0e0_DP * dval)
        tval2(3) = tset * (1.0e0_DP - 2.0e0_DP * dval)
        tval2(4) = tset * (1.0e0_DP - dval)
        tval2(5) = tset 
        tval2(6) = tset * (1.0e0_DP + dval)
        tval2(7) = tset * (1.0e0_DP + 2.0e0_DP * dval)
        tval2(8) = tset * (1.0e0_DP + 3.0e0_DP * dval)
        tval2(9) = tset * (1.0e0_DP + 4.0e0_DP * dval)

        !1. [dP / dRho]T
        t = tset
        v = 1.0e0_DP / rhoval(1)
        pval(1) = Press()
        v = 1.0e0_DP / rhoval(2)
        pval(2) = Press()
        v = 1.0e0_DP / rhoval(3)
        pval(3) = Press()
        v = 1.0e0_DP / rhoval(4)
        pval(4) = Press() 
                
        dpdr = (pval(1) - 8.0e0_DP*pval(2) + 8.0e0_DP*pval(3) - pval(4)) / (12.0e0_DP * rhoset * dval)
          
        !2. [dP / dT]V
        v = vset
        t = tval(1)
        pval(1) = Press()
        t = tval(2)
        pval(2) = Press()
        t = tval(3)
        pval(3) = Press()
        t = tval(4)
        pval(4) = Press() 
                
        dpdt = (pval(1) - 8.0e0_DP*pval(2) + 8.0e0_DP*pval(3) - pval(4)) / (12.0e0_DP * tset * dval) 
        
        !3. [d2a / dT2]V
        v = vset
        t = tval2(1)
        call Sys_setup()
        aval2(1) = A_ideal() + A_mono() + A_chain() + A_assoc() + A_ion()
        t = tval2(2)
        call Sys_setup()
        aval2(2) = A_ideal() + A_mono() + A_chain() + A_assoc() + A_ion()
        t = tval2(3)
        call Sys_setup()
        aval2(3) = A_ideal() + A_mono() + A_chain() + A_assoc() + A_ion()
        t = tval2(4)
        call Sys_setup()
        aval2(4) = A_ideal() + A_mono() + A_chain() + A_assoc() + A_ion()
        t = tval2(5)
        call Sys_setup()
        aval2(5) = A_ideal() + A_mono() + A_chain() + A_assoc() + A_ion()
        t = tval2(6)
        call Sys_setup()
        aval2(6) = A_ideal() + A_mono() + A_chain() + A_assoc() + A_ion()
        t = tval2(7)
        call Sys_setup()
        aval2(7) = A_ideal() + A_mono() + A_chain() + A_assoc() + A_ion()
        t = tval2(8)
        call Sys_setup()
        aval2(8) = A_ideal() + A_mono() + A_chain() + A_assoc() + A_ion()
        t = tval2(9)
        call Sys_setup()
        aval2(9) = A_ideal() + A_mono() + A_chain() + A_assoc() + A_ion()
        
        aval(1) = (aval2(1) - 8.0e0_DP*aval2(2) + 8.0e0_DP*aval2(4) - aval2(5)) / (12.0e0_DP * tset * dval) 
        aval(2) = (aval2(2) - 8.0e0_DP*aval2(3) + 8.0e0_DP*aval2(5) - aval2(6)) / (12.0e0_DP * tset * dval)
        aval(3) = (aval2(4) - 8.0e0_DP*aval2(5) + 8.0e0_DP*aval2(7) - aval2(8)) / (12.0e0_DP * tset * dval)        
        aval(4) = (aval2(5) - 8.0e0_DP*aval2(6) + 8.0e0_DP*aval2(8) - aval2(9)) / (12.0e0_DP * tset * dval)

        d2adt2 = (aval(1) - 8.0e0_DP*aval(2) + 8.0e0_DP*aval(3) - aval(4)) / (12.0e0_DP * tset * dval)

        !properties
        isokt       =   1.0e0_DP / (rhoset*dpdr)
        thermalpha  =   isokt * dpdt
        cv          =   -tset * d2adt2
        cp          =   cv + tset * thermalpha**2.0e0_DP / (isokt * rhoset)
        jtcoeff     =   tset*dpdt - rhoset*dpdr
        spsound     =   dsqrt(cp/cv * dpdr)

        print*,tset,pset,vset,isokt,thermalpha,cv,cp,jtcoeff,spsound

        return
    end subroutine Thermall
!***************************************************************************************************


































!***************************************************************************************************
    subroutine Therm( fixed_var, t_in, pv_in, state_in ) 
        implicit none
        
        character(len=1), intent(in)    ::  fixed_var, state_in
        real(kind=DP), intent(in)       ::  t_in, pv_in
		real(kind=DP)                   ::  t_dt			!make the temperture be right   yichun
        
        if((fixed_var=="p").or.(fixed_var=="P")) then
            !determine V
			t_dt=t_in
            t = t_in
            p = pv_in
            v = Vol_dens_g()
        else if((fixed_var=="v").or.(fixed_var=="V")) then
            !determine V
			t_dt=t_in
            t = t_in
            v = pv_in
            p = Press( )
        else
            stop "unknown fixed variable chosen in Therm subroutine.... should only be v or p"
        end if

        !zero all terms
        a_ideal_l   =   0.0e0_DP
        a_mono_l    =   0.0e0_DP
        a_chain_l   =   0.0e0_DP
        a_assoc_l   =   0.0e0_DP
        a_ion_l     =   0.0e0_DP
        a_born_l    =   0.0e0_DP
        
        u_ideal_l   =   0.0e0_DP
        u_mono_l    =   0.0e0_DP
        u_chain_l   =   0.0e0_DP
        u_assoc_l   =   0.0e0_DP
        u_ion_l     =   0.0e0_DP
        u_born_l    =   0.0e0_DP
        
        h_ideal_l   =   0.0e0_DP
        h_mono_l    =   0.0e0_DP
        h_chain_l   =   0.0e0_DP
        h_assoc_l   =   0.0e0_DP
        h_ion_l     =   0.0e0_DP
        h_born_l    =   0.0e0_DP
        
        g_ideal_l   =   0.0e0_DP
        g_mono_l    =   0.0e0_DP
        g_chain_l   =   0.0e0_DP
        g_assoc_l   =   0.0e0_DP
        g_ion_l     =   0.0e0_DP
        g_born_l    =   0.0e0_DP
        
        s_ideal_l   =   0.0e0_DP
        s_mono_l    =   0.0e0_DP
        s_chain_l   =   0.0e0_DP
        s_assoc_l   =   0.0e0_DP
        s_ion_l     =   0.0e0_DP
        s_born_l    =   0.0e0_DP

!Calculate a values
        call Sys_setup()
        a_ideal_l   =   A_ideal()       * NA * KB * t
        a_mono_l    =   A_mono()        * NA * KB * t
        a_chain_l   =   A_chain()       * NA * KB * t
        a_assoc_l   =   A_assoc()       * NA * KB * t
        a_ion_l     =   A_ion_only()    * NA * KB * t 
        a_born_l    =   A_born_only()   * NA * KB * t        
       
!Calculate u values
        t = t_dt - 2.0e0_DP * dt * t_dt
        call Sys_setup()
        
        a_ideal_l1  =   A_ideal()       * NA * KB 
        a_mono_l1   =   A_mono()        * NA * KB  
        a_chain_l1  =   A_chain()       * NA * KB 
        a_assoc_l1  =   A_assoc()       * NA * KB 
        a_ion_l1    =   A_ion_only()    * NA * KB 
        a_born_l1   =   A_born_only()   * NA * KB 
        
        t = t_dt - dt * t_dt
        call Sys_setup()

        a_ideal_l2  =   A_ideal()       * NA * KB 
        a_mono_l2   =   A_mono()        * NA * KB 
        a_chain_l2  =   A_chain()       * NA * KB 
        a_assoc_l2  =   A_assoc()       * NA * KB 
        a_ion_l2    =   A_ion_only()    * NA * KB 
        a_born_l2   =   A_born_only()   * NA * KB 
        
        t = t_dt+ dt * t_dt
        call Sys_setup()
        
        a_ideal_l3  =   A_ideal()       * NA * KB 
        a_mono_l3   =   A_mono()        * NA * KB 
        a_chain_l3  =   A_chain()       * NA * KB 
        a_assoc_l3  =   A_assoc()       * NA * KB 
        a_ion_l3    =   A_ion_only()    * NA * KB 
        a_born_l3   =   A_born_only()   * NA * KB 
        
        t = t_dt + 2.0e0_DP * dt * t_dt
        call Sys_setup()

        a_ideal_l4  =   A_ideal()       * NA * KB 
        a_mono_l4   =   A_mono()        * NA * KB 
        a_chain_l4  =   A_chain()       * NA * KB 
        a_assoc_l4  =   A_assoc()       * NA * KB 
        a_ion_l4    =   A_ion_only()    * NA * KB 
        a_born_l4   =   A_born_only()   * NA * KB 

        t = t_dt
        
        u_ideal_l   =   (a_ideal_l1 - 8.0e0_DP * a_ideal_l2 + 8.0e0_DP * a_ideal_l3 - a_ideal_l4) &
        &               / (12.0e0_DP * dt * t_dt) * (-1.0e0_DP) * t**2.0e0_DP
        u_mono_l    =   (a_mono_l1 - 8.0e0_DP * a_mono_l2 + 8.0e0_DP * a_mono_l3 - a_mono_l4) &
        &               / (12.0e0_DP * dt * t_dt) * (-1.0e0_DP) * t**2.0e0_DP 
        u_chain_l   =   (a_chain_l1 - 8.0e0_DP * a_chain_l2 + 8.0e0_DP * a_chain_l3 - a_chain_l4) &
        &               / (12.0e0_DP * dt * t_dt) * (-1.0e0_DP) * t**2.0e0_DP
        u_assoc_l   =   (a_assoc_l1 - 8.0e0_DP * a_assoc_l2 + 8.0e0_DP * a_assoc_l3 - a_assoc_l4) &
        &               / (12.0e0_DP * dt * t_dt) * (-1.0e0_DP) * t**2.0e0_DP 
        u_ion_l     =   (a_ion_l1 - 8.0e0_DP * a_ion_l2 + 8.0e0_DP * a_ion_l3 - a_ion_l4) &
        &               / (12.0e0_DP * dt * t_dt) * (-1.0e0_DP) * t**2.0e0_DP 
        u_born_l    =   (a_born_l1 - 8.0e0_DP * a_born_l2 + 8.0e0_DP * a_born_l3 - a_born_l4) &
        &               / (12.0e0_DP * dt * t_dt) * (-1.0e0_DP) * t**2.0e0_DP         
        
!Calculate h values
        h_ideal_l   =   u_ideal_l   
        h_mono_l    =   u_mono_l    
        h_chain_l   =   u_chain_l  
        h_assoc_l   =   u_assoc_l   
        h_ion_l     =   u_ion_l    
        h_born_l    =   u_born_l          
        
!Calculate g values
        g_ideal_l   =   a_ideal_l  
        g_mono_l    =   a_mono_l   
        g_chain_l   =   a_chain_l   
        g_assoc_l   =   a_assoc_l   
        g_ion_l     =   a_ion_l     
        g_born_l    =   a_born_l           
        
!Calculate s values
        s_ideal_l   =   (u_ideal_l - a_ideal_l) / t
        s_mono_l    =   (u_mono_l - a_mono_l) / t
        s_chain_l   =   (u_chain_l - a_chain_l) / t
        s_assoc_l   =   (u_assoc_l - a_assoc_l) / t
        s_ion_l     =   (u_ion_l - a_ion_l) / t
        s_born_l    =   (u_born_l - a_born_l) / t          

        return
    end subroutine Therm
!***************************************************************************************************
   function Cp( ) result(cp_out)
        implicit none
        
        real(kind=DP)                   ::  t_in_cp
        real(kind=DP), parameter        ::  dt_cp = 1.0e-3_DP
        integer                         ::  i_cp
        real(kind=DP)                   ::  t_cp(1:9), g_cp(1:9), h_cp(1:4), cp_out

        t_in_cp = t
      
        t_cp(1) = t_in_cp * (1.0e0_DP - 4.0e0_DP * dt_cp)
        t_cp(2) = t_in_cp * (1.0e0_DP - 3.0e0_DP * dt_cp)
        t_cp(3) = t_in_cp * (1.0e0_DP - 2.0e0_DP * dt_cp)
        t_cp(4) = t_in_cp * (1.0e0_DP - 1.0e0_DP * dt_cp)
        t_cp(5) = t_in_cp 
        t_cp(6) = t_in_cp * (1.0e0_DP + 1.0e0_DP * dt_cp)
        t_cp(7) = t_in_cp * (1.0e0_DP + 2.0e0_DP * dt_cp)
        t_cp(8) = t_in_cp * (1.0e0_DP + 3.0e0_DP * dt_cp)
        t_cp(9) = t_in_cp * (1.0e0_DP + 4.0e0_DP * dt_cp)

        do i_cp = 1,9
            t = t_cp(i_cp)
            v = Vol_dens_g()
            
            call Sys_setup( )
            g_cp(i_cp) = A_ideal() + A_mono() + A_chain() + A_assoc() + A_ion()
            g_cp(i_cp) = g_cp(i_cp) * NA*KB + p * v / t                   
        end do

        h_cp(1) = (g_cp(1) - 8.0e0_DP*g_cp(2) + 8.0e0_DP*g_cp(4) - g_cp(5)) / (12.0e0_DP * dt_cp * t_in_cp )
        h_cp(2) = (g_cp(2) - 8.0e0_DP*g_cp(3) + 8.0e0_DP*g_cp(5) - g_cp(6)) / (12.0e0_DP * dt_cp * t_in_cp )
        h_cp(3) = (g_cp(4) - 8.0e0_DP*g_cp(5) + 8.0e0_DP*g_cp(7) - g_cp(8)) / (12.0e0_DP * dt_cp * t_in_cp )
        h_cp(4) = (g_cp(5) - 8.0e0_DP*g_cp(6) + 8.0e0_DP*g_cp(8) - g_cp(9)) / (12.0e0_DP * dt_cp * t_in_cp )

        h_cp(1) = -1.0e0_DP * h_cp(1) * t_cp(3)**2.0e0_DP
        h_cp(2) = -1.0e0_DP * h_cp(2) * t_cp(4)**2.0e0_DP
        h_cp(3) = -1.0e0_DP * h_cp(3) * t_cp(6)**2.0e0_DP
        h_cp(4) = -1.0e0_DP * h_cp(4) * t_cp(7)**2.0e0_DP
              
        t = t_in_cp
         
        cp_out = (h_cp(1) - 8.0e0_DP*h_cp(2) + 8.0e0_DP*h_cp(3) - h_cp(4)) / (12.0e0_DP * dt_CP * t )

        cp_out = cp_out + properties%A + properties%B * t + properties%C * t**2e0_DP + properties%D    * &
        &   t**3e0_DP + properties%E * t**4e0_DP - 2.5e0_DP * NA * KB

        return
    end function Cp
!***************************************************************************************************
   function Cp_ex( ) result(cp_out)
        implicit none
        
        real(kind=DP)                   ::  t_in_cp
        real(kind=DP), parameter        ::  dt_cp = 1.0e-3_DP
        integer                         ::  i_cp, i_cp2
        real(kind=DP)                   ::  t_cp(1:9), g_cp(1:9), h_cp(1:4), cp_out

        t_in_cp = t
      
        t_cp(1) = t_in_cp * (1.0e0_DP - 4.0e0_DP * dt_cp)
        t_cp(2) = t_in_cp * (1.0e0_DP - 3.0e0_DP * dt_cp)
        t_cp(3) = t_in_cp * (1.0e0_DP - 2.0e0_DP * dt_cp)
        t_cp(4) = t_in_cp * (1.0e0_DP - 1.0e0_DP * dt_cp)
        t_cp(5) = t_in_cp 
        t_cp(6) = t_in_cp * (1.0e0_DP + 1.0e0_DP * dt_cp)
        t_cp(7) = t_in_cp * (1.0e0_DP + 2.0e0_DP * dt_cp)
        t_cp(8) = t_in_cp * (1.0e0_DP + 3.0e0_DP * dt_cp)
        t_cp(9) = t_in_cp * (1.0e0_DP + 4.0e0_DP * dt_cp)

        do i_cp = 1,9
            t = t_cp(i_cp)
            v = Vol_dens_g()
            
            g_cp(i_cp) = 0.0e0_DP
            
            do i_cp2 = 1,nctypes
                g_cp(i_cp) = g_cp(i_cp) + Mu(i_cp2) * Comp_array(i_cp2)%xi
            end do      
                                    
        end do

        h_cp(1) = (g_cp(1) - 8.0e0_DP*g_cp(2) + 8.0e0_DP*g_cp(4) - g_cp(5)) / (12.0e0_DP * dt_cp * t_in_cp )
        h_cp(2) = (g_cp(2) - 8.0e0_DP*g_cp(3) + 8.0e0_DP*g_cp(5) - g_cp(6)) / (12.0e0_DP * dt_cp * t_in_cp )
        h_cp(3) = (g_cp(4) - 8.0e0_DP*g_cp(5) + 8.0e0_DP*g_cp(7) - g_cp(8)) / (12.0e0_DP * dt_cp * t_in_cp )
        h_cp(4) = (g_cp(5) - 8.0e0_DP*g_cp(6) + 8.0e0_DP*g_cp(8) - g_cp(9)) / (12.0e0_DP * dt_cp * t_in_cp )

        h_cp(1) = -1.0e0_DP * h_cp(1) * t_cp(3)**2.0e0_DP
        h_cp(2) = -1.0e0_DP * h_cp(2) * t_cp(4)**2.0e0_DP
        h_cp(3) = -1.0e0_DP * h_cp(3) * t_cp(6)**2.0e0_DP
        h_cp(4) = -1.0e0_DP * h_cp(4) * t_cp(7)**2.0e0_DP
              
        t = t_in_cp
         
        cp_out = (h_cp(1) - 8.0e0_DP*h_cp(2) + 8.0e0_DP*h_cp(3) - h_cp(4)) / (12.0e0_DP * dt_CP * t ) &
        & - 2.5e0_DP * NA * KB

        return
    end function Cp_ex
!***************************************************************************************************
    subroutine Therm_all( i_state_in ) 
        implicit none
        
        integer, intent(in)     ::  i_state_in
        real(kind=DP)           ::  a_sum, u_sum, h_sum, s_sum, g_sum, cv_temp, cp_temp
       
        call Therm( "p", properties%t(i_state_in), properties%p(i_state_in), properties%phase(i_state_in) ) 
    
        a_sum = a_ideal_l + a_mono_l + a_chain_l + a_assoc_l + a_ion_l + a_born_l
        u_sum = u_ideal_l + u_mono_l + u_chain_l + u_assoc_l + u_ion_l + u_born_l
        h_sum = h_ideal_l + h_mono_l + h_chain_l + h_assoc_l + h_ion_l + h_born_l + p * v
        s_sum = s_ideal_l + s_mono_l + s_chain_l + s_assoc_l + s_ion_l + s_born_l
        g_sum = g_ideal_l + g_mono_l + g_chain_l + g_assoc_l + g_ion_l + g_born_l + p * v
               
        cp_temp = Cp( )

        write(*, '(20e16.6)') Comp_array(1:nctypes)%xi, t, p / 1.e6_DP, v,                        &
        &   a_sum / 1.0e3_DP, u_sum / 1.0e3_DP, h_sum / 1.0e3_DP, s_sum , g_sum / 1.0e3_DP,   &
        &   cp_temp 
        
        return
    end subroutine Therm_all
!***************************************************************************************************
    function Enth(  ) result(enth_out) 
        implicit none
        
        real(kind=DP)           ::  spec_mu, enth_out

        call Therm( "p", t, p, "l" )
        enth_out = h_ideal_l + h_mono_l + h_chain_l + h_assoc_l + h_ion_l + h_born_l  + p * v
        
        return
    end function Enth
!***************************************************************************************************
    function ThermG(  ) result(gibbs_out) 
        implicit none
        
        real(kind=DP)           ::  gibbs_out

        call Therm( "p", t, p, "l" )
        gibbs_out = g_ideal_l + g_mono_l + g_chain_l + g_assoc_l + g_ion_l + g_born_l  + p * v
        
        return
    end function ThermG
!***************************************************************************************************
    function EnthV(  ) result(enth_out) 
        implicit none
        
        real(kind=DP)           ::  spec_mu, enth_out

        call Therm( "v", t, v, "l" )
        enth_out = h_ideal_l + h_mono_l + h_chain_l + h_assoc_l + h_ion_l + h_born_l  + p * v
       
        return
    end function EnthV
!***************************************************************************************************
    function ThermGV(  ) result(gibbs_out) 
        implicit none
        
        real(kind=DP)           ::  gibbs_out

        call Therm( "v", t, v, "l" )      
        gibbs_out = g_ideal_l + g_mono_l + g_chain_l + g_assoc_l + g_ion_l + g_born_l  + p * v
        
        return
    end function ThermGV
!***************************************************************************************************
    function Mu_ex( i_mu_ex ) result(gibbs_out) 
        implicit none
        
        real(kind=DP)           ::  gibbs_out
        integer,intent(in)      ::  i_mu_ex

        gibbs_out = Mu( i_mu_ex )
        
        return
    end function Mu_ex
!***************************************************************************************************
end module Therm_mod
!***************************************************************************************************