!***************************************************************************************************
!   SAFT Module to:
!       1. Calculate A ion and A born
!
!       Based on:
!           D.K. Eriksen, G. Lazarou, A. Galindo, G. Jackson, C.S. Adjiman, A.J. Haslam,           
!           "Development of intermolecular potential models for electrolyte
!           solutions using an electrolyte SAFT-VR Mie equation of state", 
!           Mol. Phys. 114(18), 2016 2724-2749
!           
!       Analytical derivatives added March 2017
!
!***************************************************************************************************
!
!***************************************************************************************************
module Ion
!***************************************************************************************************
!Modules
!=======    
    use Types           ! Definitions of types and double precision
    use Global          ! Important global parameters     
!***************************************************************************************************
    implicit none
!***************************************************************************************************
    contains
!***************************************************************************************************
    function A_ion( ) result(a_res)     !AJ Haslam et al Mol Phys 2016            
        implicit none
    
        real(kind=DP)       ::  a_res
        real(kind=DP)       ::  dielec

        a_res=0.0e0_DP

        if(ion_switch) then   
            dielec = Diel()                 
            a_res  = Aion( dielec ) + Aborn( dielec )
        end if

        return
    end function A_ion
!***************************************************************************************************
    function A_ion_only( ) result(a_res)     !As part of pressure calculation              
        implicit none
    
        real(kind=DP)       ::  a_res
        real(kind=DP)       ::  dielec

        a_res=0.0e0_DP

        if(ion_switch) then
            dielec = Diel()          
            a_res = Aion( dielec )
        end if

        return
    end function A_ion_only  
!***************************************************************************************************
    function A_born_only( ) result(a_res)     !As part of mu inf calculation              
        implicit none
    
        real(kind=DP)       ::  a_res
        real(kind=DP)       ::  dielec

        a_res=0.0e0_DP

        if(ion_switch) then 
            dielec = Diel()          
            a_res = Aborn( dielec )
        end if

        return
    end function A_born_only  
!*************************************************************************************************** 
    function Diel( ) result( diel_res )
        implicit none

        real(kind=DP)               ::  diel_res
        real(kind=DP)               ::  x_solv, ddiel 
        integer                     ::  i_diel1, i_diel2

        x_solv=0.0e0_DP
        ddiel=0.0e0_DP
        
        do i_diel1=1, nctypes
            if(Comp_array(i_diel1)%solv) x_solv = x_solv + Comp_array(i_diel1)%xi
        end do
        
        do i_diel1=1, nctypes
        do i_diel2=i_diel1, nctypes
            if(( Comp_array(i_diel1)%solv ).and.( Comp_array(i_diel2)%solv )) then           
                ddiel = ddiel + ( Comp_array(i_diel1)%xi / x_solv ) * ( Comp_array(i_diel2)%xi / x_solv ) * &
                & ( Comp_array(i_diel1)%dv * ( Comp_array(i_diel1)%dt / t - 1.0e0_DP ) + &
                & Comp_array(i_diel2)%dv * ( Comp_array(i_diel2)%dt / t - 1.0e0_DP ) ) / 2.0e0_DP             
            end if 
        end do
        end do

        diel_res= 1.0e0_DP + rho * ddiel * x_solv 
     
        return
    end function Diel
!***************************************************************************************************
    function dDielN( ) result( diel_res )
        implicit none

        real(kind=DP)               ::  diel_res
        real(kind=DP)               ::  x_solv, ddiel 
        real(kind=DP)               ::  x_solv_dn, ddiel_dn 
        integer                     ::  i_diel1, i_diel2

        x_solv=0.0e0_DP
        ddiel=0.0e0_DP
        
        x_solv_dn=0.0e0_DP
        ddiel_dn=0.0e0_DP
        
        do i_diel1=1, nctypes
            if(Comp_array(i_diel1)%solv) then
                x_solv = x_solv + Comp_array(i_diel1)%xi
                x_solv_dn = x_solv_dn + Comp_array(i_diel1)%dxin
            end if
        end do
        
        do i_diel1=1, nctypes
        do i_diel2=i_diel1, nctypes
            if(( Comp_array(i_diel1)%solv ).and.( Comp_array(i_diel2)%solv )) then           
                ddiel = ddiel + ( Comp_array(i_diel1)%xi / x_solv ) * ( Comp_array(i_diel2)%xi / x_solv ) * &
                & ( Comp_array(i_diel1)%dv * ( Comp_array(i_diel1)%dt / t - 1.0e0_DP ) + &
                & Comp_array(i_diel2)%dv * ( Comp_array(i_diel2)%dt / t - 1.0e0_DP ) ) / 2.0e0_DP   
                
                ddiel_dn = ddiel_dn + ( Comp_array(i_diel1)%xi / x_solv ) * ( Comp_array(i_diel2)%xi / x_solv ) * &
                
                &   ((x_solv*(Comp_array(i_diel1)%xi*Comp_array(i_diel2)%dxin + Comp_array(i_diel1)%dxin*Comp_array(i_diel2)%xi) - &
                & 2.0e0_DP*Comp_array(i_diel1)%xi*Comp_array(i_diel2)%xi*x_solv_dn ) / x_solv**3.0e0_DP) * &
                
                & ( Comp_array(i_diel1)%dv * ( Comp_array(i_diel1)%dt / t - 1.0e0_DP ) + &
                & Comp_array(i_diel2)%dv * ( Comp_array(i_diel2)%dt / t - 1.0e0_DP ) ) / 2.0e0_DP 
            end if 
        end do
        end do

        diel_res= x_solv*rho*ddiel_dn + x_solv_dn*rho*ddiel_dn + x_solv*drhon*ddiel + x_solv_dn*drhon*ddiel  
      
        return
    end function dDielN
!***************************************************************************************************
    function dDielv( ) result( diel_res )
        implicit none

        real(kind=DP)               ::  diel_res
        real(kind=DP)               ::  x_solv, ddiel 
        integer                     ::  i_diel1, i_diel2

        x_solv=0.0e0_DP
        ddiel=0.0e0_DP
        
        do i_diel1=1, nctypes
            if(Comp_array(i_diel1)%solv) x_solv = x_solv + Comp_array(i_diel1)%xi
        end do
        
        do i_diel1=1, nctypes
        do i_diel2=i_diel1, nctypes
            if(( Comp_array(i_diel1)%solv ).and.( Comp_array(i_diel2)%solv )) then           
                ddiel = ddiel + ( Comp_array(i_diel1)%xi / x_solv ) * ( Comp_array(i_diel2)%xi / x_solv ) * &
                & ( Comp_array(i_diel1)%dv * ( Comp_array(i_diel1)%dt / t - 1.0e0_DP ) + &
                & Comp_array(i_diel2)%dv * ( Comp_array(i_diel2)%dt / t - 1.0e0_DP ) ) / 2.0e0_DP             
            end if 
        end do
        end do

        diel_res= 1.0e0_DP + drhov * ddiel * x_solv 
   
        return
    end function dDielv
!***************************************************************************************************
    function Aborn( l_dielec ) result( born_res )
        implicit none
        
        real(kind=DP), intent(in)   ::  l_dielec
        real(kind=DP)   ::  born_res
        real(kind=DP)   ::  bsum, born_const, born_diel
        integer         ::  i_born1, i_born2, i_born3

        bsum=0.0e0_DP
        
        born_diel = 1.0e0_DP - 1.0e0_DP / l_dielec
        born_const = -1.0e0 * QE**2.0e0_DP / 4.0e0_DP / PI / E0
     
        do i_born1=1, nctypes
            do i_born2=1, nstypes
                if(Seg_array(i_born2)%charged) then
                    do i_born3=1, Comp_array(i_born1)%comp(i_born2)
!This changed - xi * NA to  nm                    
                        bsum = bsum + Comp_array(i_born1)%nm * Seg_array(i_born2)%q**2.0e0_DP /     &                      
                        &       Seg_array(i_born2)%sigb                     
                    end do
                end if
            end do           
        end do
!This changed - NA to  nm sum      
        born_res = born_const * born_diel * bsum
        born_res = born_res / sum(Comp_array(:)%nm) / KB / t         

        return
    end function Aborn
!***************************************************************************************************
    function A_born_dv(  ) result( born_res )      
        implicit none
        
        real(kind=DP)   ::  born_res
        real(kind=DP)   ::  bsum, born_const, born_diel, x_solv, ddiel
        integer         ::  i_born1, i_born2, i_born3

        bsum=0.0e0_DP
        x_solv=0.0e0_DP
        ddiel=0.0e0_DP
        
        if(ion_switch) then 
            do i_born1=1, nctypes
                do i_born2=1, nstypes
                    if(Seg_array(i_born2)%charged) then
                        do i_born3=1, Comp_array(i_born1)%comp(i_born2)
                            bsum = bsum + Comp_array(i_born1)%xi * NA * Seg_array(i_born2)%q**2.0e0_DP /     &                      
                            &       Seg_array(i_born2)%sigb
                        end do
                    end if
                end do
            end do
            
            do i_born1=1, nctypes
                if(Comp_array(i_born1)%solv) x_solv = x_solv + Comp_array(i_born1)%xi
            end do
           
            do i_born1=1, nctypes
            do i_born2=i_born1, nctypes
                if(( Comp_array(i_born1)%solv ).and.( Comp_array(i_born2)%solv )) then           
                    ddiel = ddiel + ( Comp_array(i_born1)%xi / x_solv ) * ( Comp_array(i_born2)%xi / x_solv ) * &
                    & ( Comp_array(i_born1)%dv * ( Comp_array(i_born1)%dt / t - 1.0e0_DP ) + &
                    & Comp_array(i_born2)%dv * ( Comp_array(i_born2)%dt / t - 1.0e0_DP ) ) / 2.0e0_DP
                end if 
            end do
            end do
        end if
        
        born_diel = x_solv * ddiel / (v + x_solv * ddiel)**2.0e0_DP
        born_const = QE**2.0e0_DP / 4.0e0_DP / PI / E0
        
        born_res = born_const * born_diel * bsum
        born_res = born_res / NA / KB / t         

        return
    end function A_born_dv
!***************************************************************************************************
    function A_born_dn( i_ion ) result( born_res )      
        implicit none
        
        integer, intent(in) ::  i_ion
        real(kind=DP)       ::  born_res
        real(kind=DP)       ::  diel_l, born_const, diel_l_dn
        integer     ::  i_born1, i_born2, i_born3
        
        born_res = 0.0e0_DP
        
        if(ion_switch) then 
            
            born_const = -1.0e0 * QE**2.0e0_DP / 4.0e0_DP / PI / E0
            
            if((i_ion==properties%cation).or.(i_ion==properties%anion)) then
                diel_l = Diel( )             
                born_res = born_const * (1.0e0_DP - 1.0e0_DP / diel_l)      *   &
                &   Seg_array(i_ion)%q**2.0e0_DP / Seg_array(i_ion)%sigb
            else            
                diel_l = Diel( ) 
                diel_l_dn = dDielN()
                                
                do i_born1=1, nctypes
                do i_born2=1, nstypes
                if(Seg_array(i_born2)%charged) then
                do i_born3=1, Comp_array(i_born1)%comp(i_born2)
                born_res = born_const * (-diel_l_dn/diel_l**2.0e0_DP)      *   &
                      &   Seg_array(i_born2)%q**2.0e0_DP / Seg_array(i_born2)%sigb                          
                end do
                end if
                end do           
                end do            
            end if

        end if
                                                                        
        born_res = born_res / KB / t / sum(Comp_array(:)%nm)  
    
        return
    end function A_born_dn
!***************************************************************************************************
    function Aion( l_dielec ) result( ion_res )        
        implicit none

        real(kind=DP)               ::  ion_res
        real(kind=DP), intent(in)   ::  l_dielec
        real(kind=DP)               ::  c_shield, shield
        
        real(kind=DP), parameter    ::  SHIELD_CONV=1.0e-12_DP
        
        real(kind=DP)               ::  qsum_l, del_l, omega_l, pn_l, shield_l
        real(kind=DP)               ::  umsa, qtemp
        
        integer                     ::  i_ion1, i_ion2, i_ion3
        integer                     ::  it_tot
!Here also changed all to nm        
        c_shield = QE**2.0e0_DP / (4.0e0_DP * l_dielec * E0 *  KB * t * v)        
        qsum_l = 0.0e0_DP
        
        do i_ion1=1, nctypes   
            do i_ion2=1, nstypes
                if( Seg_array(i_ion2)%charged ) then
                    do i_ion3=1, Comp_array(i_ion1)%comp(i_ion2)
                        qsum_l = qsum_l + Comp_array(i_ion1)%nm * Seg_array(i_ion2)%q**2.0e0_DP                        
                    end do
                end if 
            end do
        end do
        
        !initial shield guess
        shield = 0.5e0_DP * dsqrt( c_shield * qsum_l )

        !set delta  - this is constant for iteration       
        del_l = 0.0e0_DP

        do i_ion1=1, nctypes
            do i_ion2=1, nstypes
                if(Seg_array(i_ion2)%charged) then
                    do i_ion3=1, Comp_array(i_ion1)%comp(i_ion2)
                        del_l = del_l + Comp_array(i_ion1)%nm * Seg_array(i_ion2)%sig**3.0e0_DP                   
                    end do
                end if 
            end do
        end do
        
        del_l = 1.0e0_DP - PI / 6.0e0_DP / v * del_l            
     
        !iterate parameters to find shield
        it_tot = 0       
        
        do             
            omega_l = 0.0e0_DP
            
            do i_ion1=1, nctypes
                do i_ion2=1, nstypes
                    if(Seg_array(i_ion2)%charged) then
                        do i_ion3=1, Comp_array(i_ion1)%comp(i_ion2)
                            omega_l = omega_l + Comp_array(i_ion1)%nm * Seg_array(i_ion2)%sig**3.0e0_DP / &
                                &   (1.0e0_DP + shield * Seg_array(i_ion2)%sig)
                        end do
                    end if 
                end do
            end do

            omega_l = 1.0e0_DP + PI / v / 2.e0_DP / del_l * omega_l

            pn_l = 0.0e0_DP
            
            do i_ion1=1, nctypes
                do i_ion2=1, nstypes
                    if(Seg_array(i_ion2)%charged) then
                        do i_ion3=1, Comp_array(i_ion1)%comp(i_ion2)
!added *QE
                            pn_l = pn_l + Comp_array(i_ion1)%nm * Seg_array(i_ion2)%sig *       &
                            &   Seg_array(i_ion2)%q*QE / (1.0e0_DP + shield * Seg_array(i_ion2)%sig )
                           
                        end do
                    end if 
                end do
            end do
           
            pn_l = 1.0e0_DP / omega_l / v * pn_l 

            shield_l = 0.0e0_DP
            
            do i_ion1=1, nctypes
                do i_ion2=1, nstypes
                    if(Seg_array(i_ion2)%charged) then
                        do i_ion3=1, Comp_array(i_ion1)%comp(i_ion2)
                            qtemp = ( &
                                     Seg_array(i_ion2)%q * QE - (PI / 2.0e0_DP / del_l)  &
                                     * Seg_array(i_ion2)%sig**2.0e0_DP * pn_l &
                                    ) &
                                    / (1.0e0_DP + shield * Seg_array(i_ion2)%sig)                     
                            shield_l = shield_l + Comp_array(i_ion1)%nm * qtemp**2.0e0_DP
                        end do
                    end if 
                end do
            end do
    
            shield_l = dsqrt( PI  / 4.0e0_DP / PI / E0 / l_dielec / KB / T / v * shield_l )
        
            it_tot = it_tot + 1
                    
            if(dabs( shield - shield_l ) / dabs( shield_l )  < SHIELD_CONV) then
                shield = shield_l
                exit
            else
                shield = shield_l
            end if    
        end do 

        umsa = 0.0e0_DP

        do i_ion1=1, nctypes
            do i_ion2=1, nstypes
                if(Seg_array(i_ion2)%charged) then
                    do i_ion3=1, Comp_array(i_ion1)%comp(i_ion2)
                        umsa = umsa + (Comp_array(i_ion1)%nm * Seg_array(i_ion2)%q**2.0e0_DP* QE**2.0e0_DP )  /   &
                        &   ( 1.0e0_DP + shield * Seg_array(i_ion2)%sig )  
                    end do
                end if 
            end do
        end do 

        umsa = -1.0e0_DP * v / 4.0e0_DP / PI / E0 / l_dielec * (shield / v * umsa +  &
        &       PI / 2.0e0_DP / del_l * omega_l * pn_l**2.0e0_DP )

        ion_res = umsa + shield**3.0e0_DP * KB * v * t / 3.0e0_DP / PI

        !Make it in equivalent units...
        ion_res = ion_res / sum(Comp_array(:)%nm) / KB / t
                 
        return
    end function Aion
!***************************************************************************************************  
    function A_ion_dn( i_part ) result(a_res)
        implicit none
        
        integer, intent(in)     ::  i_part
        real(kind=DP)           ::  a_res

        real(kind=DP)           ::  l_dielec
        real(kind=DP)               ::  c_shield, shield
        
        real(kind=DP), parameter    ::  SHIELD_CONV=1.0e-12_DP
        
        real(kind=DP)               ::  qsum_l, del_l, omega_l, pn_l, shield_l
        real(kind=DP)               ::  umsa, qtemp(1:nctypes,1:nstypes,1:20)
        
        real(kind=DP)               ::  del_lp, omega_lp, pn_lp, shield_lp, shield_temp
        real(kind=DP)               ::  qtempp(1:nctypes,1:nstypes,1:20), umsa_p
        
        real(kind=DP)               ::  temp_u, temp_v, temp_up, temp_vp
        
        integer                     ::  i_ion1, i_ion2, i_ion3
        integer                     ::  it_tot

        if(.not.ion_switch) then
            a_res = 0.0e0_DP
            return
        end if
        
        l_dielec = Diel()
        
        qtemp = 0.0e0_DP
!Do full A_ion calc first...      
        c_shield = QE**2.0e0_DP / (4.0e0_DP * l_dielec * E0 *  KB * t * v)        

        qsum_l = 0.0e0_DP
             
        do i_ion1=1, nctypes
            do i_ion2=1, nstypes
                if( Seg_array(i_ion2)%charged ) then
                    do i_ion3=1, Comp_array(i_ion1)%comp(i_ion2)
                        qsum_l = qsum_l + Comp_array(i_ion1)%nm * Seg_array(i_ion2)%q**2.0e0_DP                        
                    end do
                end if 
            end do
        end do
        
        !initial shield guess
        shield = 0.5e0_DP * dsqrt( c_shield * qsum_l )

        !set delta  - this is constant for iteration       
        del_l = 0.0e0_DP

        do i_ion1=1, nctypes
            do i_ion2=1, nstypes
                if(Seg_array(i_ion2)%charged) then
                    do i_ion3=1, Comp_array(i_ion1)%comp(i_ion2)
                        del_l = del_l + Comp_array(i_ion1)%nm * Seg_array(i_ion2)%sig**3.0e0_DP                   
                    end do
                end if 
            end do
        end do
        
        del_l = 1.0e0_DP - PI / 6.0e0_DP / v * del_l            

        !iterate parameters to find shield
        it_tot = 0       
        
        do             
            omega_l = 0.0e0_DP
            
            do i_ion1=1, nctypes
                do i_ion2=1, nstypes
                    if(Seg_array(i_ion2)%charged) then
                        do i_ion3=1, Comp_array(i_ion1)%comp(i_ion2)
                            omega_l = omega_l + Comp_array(i_ion1)%nm * Seg_array(i_ion2)%sig**3.0e0_DP / &
                                &   (1.0e0_DP + shield * Seg_array(i_ion2)%sig)
                        end do
                    end if 
                end do
            end do

            omega_l = 1.0e0_DP + PI / v / 2.e0_DP / del_l * omega_l

            pn_l = 0.0e0_DP
            
            do i_ion1=1, nctypes
                do i_ion2=1, nstypes
                    if(Seg_array(i_ion2)%charged) then
                        do i_ion3=1, Comp_array(i_ion1)%comp(i_ion2)
!added *QE
                            pn_l = pn_l + Comp_array(i_ion1)%nm * Seg_array(i_ion2)%sig *       &
                            &   Seg_array(i_ion2)%q*QE / (1.0e0_DP + shield * Seg_array(i_ion2)%sig )
                           
                        end do
                    end if 
                end do
            end do
           
            pn_l = 1.0e0_DP / omega_l / v * pn_l 

            shield_l = 0.0e0_DP
            
            do i_ion1=1, nctypes
                do i_ion2=1, nstypes
                    if(Seg_array(i_ion2)%charged) then
                        do i_ion3=1, Comp_array(i_ion1)%comp(i_ion2)
                            qtemp(i_ion1,i_ion2,i_ion3) = ( Seg_array(i_ion2)%q * QE - (PI / 2.0e0_DP / del_l) &
                                    * Seg_array(i_ion2)%sig**2.0e0_DP * pn_l ) &
                                    / (1.0e0_DP + shield * Seg_array(i_ion2)%sig)                     

                            shield_l = shield_l + Comp_array(i_ion1)%nm * qtemp(i_ion1,i_ion2,i_ion3)**2.0e0_DP                     
                        end do
                    end if 
                end do
            end do
       
            shield_l = dsqrt( PI  / 4.0e0_DP / PI / E0 / l_dielec / KB / T / v * shield_l )
          
            it_tot = it_tot + 1
                    
            if(dabs( shield - shield_l ) / dabs( shield_l )  < SHIELD_CONV) then
                shield = shield_l
                exit
            else
                shield = shield_l
            end if    
        end do 
        
        umsa = 0.0e0_DP

        do i_ion1=1, nctypes
            do i_ion2=1, nstypes
                if(Seg_array(i_ion2)%charged) then
                    do i_ion3=1, Comp_array(i_ion1)%comp(i_ion2)
                        umsa = umsa + (Comp_array(i_ion1)%nm * Seg_array(i_ion2)%q**2.0e0_DP* QE**2.0e0_DP )  /   &
                        &   ( 1.0e0_DP + shield * Seg_array(i_ion2)%sig )  
                    end do
                end if 
            end do
        end do 
  
        umsa = -1.0e0_DP * v / 4.0e0_DP / PI / E0 / l_dielec * (shield / v * umsa +  &
        &       PI / 2.0e0_DP / del_l * omega_l * pn_l**2.0e0_DP )
    
!shield converged, do differential
        del_lp = 0.0e0_DP
        i_ion1=i_part
        
        do i_ion2=1, nstypes
            if(Seg_array(i_ion2)%charged) then
                do i_ion3=1, Comp_array(i_ion1)%comp(i_ion2)
                    del_lp = del_lp - Seg_array(i_ion2)%sig**3.0e0_DP                   
                end do
            end if 
        end do
      
        del_lp = del_lp * PI / 6.0e0_DP / v

!initial shield' guess        
        qsum_l = 0.0e0_DP           

        temp_u  = 0.0e0_DP
        do i_ion1=1, nctypes
            do i_ion2=1, nstypes
                if( Seg_array(i_ion2)%charged ) then
                    do i_ion3=1, Comp_array(i_ion1)%comp(i_ion2)
                        temp_u = temp_u + Comp_array(i_ion1)%nm * Seg_array(i_ion2)%q**2.0e0_DP                        
                    end do
                end if 
            end do
        end do
     
        shield_lp = 0.5e0_DP / (0.25e0_DP * c_shield * temp_u)**0.5e0_DP * 0.25e0_DP * c_shield * &
        &   Seg_array(i_part)%q**2.0e0_DP       
        
        do !shield differential
            temp_v  = 0.0e0_DP 
            temp_vp = 0.0e0_DP
            temp_up = -PI * del_lp / 2.0e0_DP / v / del_l**2.0e0_DP
            temp_u = PI / 2.0e0_DP / del_l / v
    
            do i_ion1=1, nctypes
                do i_ion2=1, nstypes
                    if(Seg_array(i_ion2)%charged) then
                        do i_ion3=1, Comp_array(i_ion1)%comp(i_ion2)
                        
                            temp_v = temp_v + Comp_array(i_ion1)%nm * Seg_array(i_ion2)%sig**3.0e0_DP / &
                                &   (1.0e0_DP + shield * Seg_array(i_ion2)%sig)
                        
                            if(i_ion1==i_part) then
                                temp_vp = temp_vp + (Seg_array(i_ion2)%sig**3.0e0_DP * (1.0e0_DP + shield   * & 
                                & Seg_array(i_ion2)%sig) - Comp_array(i_ion1)%nm * shield_lp                * &
                                & Seg_array(i_ion2)%sig**4.0e0_DP)                                          / &
                                & (1.0e0_DP + shield * Seg_array(i_ion2)%sig)**2.0e0_DP
                            else
                                temp_vp = temp_vp - (Comp_array(i_ion1)%nm * shield_lp                      * &
                                & Seg_array(i_ion2)%sig**4.0e0_DP)                                          / &
                                & (1.0e0_DP + shield * Seg_array(i_ion2)%sig)**2.0e0_DP
                            end if
    
                        end do
                    end if 
                end do
            end do     
    
            omega_lp = temp_u * temp_vp + temp_v * temp_up
     
            temp_up = -omega_lp / v / omega_l**2.0e0_DP
            temp_vp = 0.0e0_DP
            temp_v  = 0.0e0_DP
            temp_u  = 1.0e0_DP /omega_l / v
            
            do i_ion1=1, nctypes
                do i_ion2=1, nstypes
                    if(Seg_array(i_ion2)%charged) then
                        do i_ion3=1, Comp_array(i_ion1)%comp(i_ion2)
    
                            temp_v = temp_v + Comp_array(i_ion1)%nm * Seg_array(i_ion2)%sig *       &
                            &   Seg_array(i_ion2)%q*QE / (1.0e0_DP + shield * Seg_array(i_ion2)%sig )
                            
                            if(i_ion1==i_part) then
                                temp_vp = temp_vp + (Seg_array(i_ion2)%sig * Seg_array(i_ion2)%q*QE * &
                                &   (1.0e0_DP + shield * Seg_array(i_ion2)%sig)                     - &
                                &   shield_lp * Seg_array(i_ion2)%sig * (Comp_array(i_ion1)%nm      * &
                                &   Seg_array(i_ion2)%sig * Seg_array(i_ion2)%q*QE))                / &
                                &   (1.0e0_DP + shield * Seg_array(i_ion2)%sig)**2.0e0_DP
                            else
                                temp_vp = temp_vp - (shield_lp * Seg_array(i_ion2)%sig * (Comp_array(i_ion1)%nm * &
                                &   Seg_array(i_ion2)%sig * Seg_array(i_ion2)%q*QE))                            / &
                                &   (1.0e0_DP + shield * Seg_array(i_ion2)%sig)**2.0e0_DP                        
                            end if
                                
                        end do
                    end if 
                end do
            end do        
            
            pn_lp = temp_u * temp_vp + temp_v * temp_up
      
            do i_ion1=1, nctypes
                do i_ion2=1, nstypes
                    if(Seg_array(i_ion2)%charged) then
                        do i_ion3=1, Comp_array(i_ion1)%comp(i_ion2)
                        
                            temp_u  = Seg_array(i_ion2)%q * QE - (PI / 2.0e0_DP / del_l) *  &
                            &   Seg_array(i_ion2)%sig**2.0e0_DP * pn_l 
                            temp_up = (Seg_array(i_ion2)%sig**2.0e0_DP * PI * (pn_l*del_lp - pn_lp*del_l))  / &
                            &   2.0e0_DP / del_l**2.0e0_DP
                            temp_v  = 1.0e0_DP + shield * Seg_array(i_ion2)%sig
                            temp_vp = shield_lp * Seg_array(i_ion2)%sig 
    
                            qtempp(i_ion1,i_ion2,i_ion3) = (temp_v * temp_up - temp_u * temp_vp) / temp_v**2.0e0_DP
                        end do
                    end if 
                end do
            end do        
         
            temp_u = 0.0e0_DP
            temp_v = 0.0e0_DP           
            
            do i_ion1=1, nctypes
                do i_ion2=1, nstypes
                    if(Seg_array(i_ion2)%charged) then
                        do i_ion3=1, Comp_array(i_ion1)%comp(i_ion2)
    
                        
                            temp_u = temp_u + Comp_array(i_ion1)%nm * qtemp(i_ion1,i_ion2,i_ion3)**2.0e0_DP 
    
                            if(i_ion1==i_part) then
                                temp_v = temp_v + qtemp(i_ion1,i_ion2,i_ion3)**2.0e0_DP             + &
                                2.0e0_DP * Comp_array(i_ion1)%nm * qtemp(i_ion1,i_ion2,i_ion3)      * &
                                & qtempp(i_ion1,i_ion2,i_ion3)
                            else
                                temp_v = temp_v + qtemp(i_ion1,i_ion2,i_ion3)**2.0e0_DP             + &
                                2.0e0_DP * Comp_array(i_ion1)%nm * qtemp(i_ion1,i_ion2,i_ion3)      * &
                                & qtempp(i_ion1,i_ion2,i_ion3)
                            end if
                        end do
                    end if 
                end do
            end do     
    
            shield_temp = 0.5e0_DP * (PI  / 4.0e0_DP / PI / E0 / l_dielec / KB / T / v * temp_u)**(-0.5e0_DP)     * &
            &   (PI  / 4.0e0_DP / PI / E0 / l_dielec / KB / T / v * temp_v)
            
            if(dabs( shield_lp - shield_temp ) / dabs( shield_temp )  < SHIELD_CONV) then
                shield_lp = shield_temp
                exit
            else
                shield_lp = shield_temp
            end if      
                        
        end do
        
        umsa_p = 0.0e0_DP
        
        temp_u  = shield / v
        temp_v  = 0.0e0_DP
        temp_up = shield_lp / v
        temp_vp = 0.0e0_DP

        do i_ion1=1, nctypes
            do i_ion2=1, nstypes
                if(Seg_array(i_ion2)%charged) then
                    do i_ion3=1, Comp_array(i_ion1)%comp(i_ion2)
                        temp_v = temp_v + (Comp_array(i_ion1)%nm * Seg_array(i_ion2)%q**2.0e0_DP    * &
                        &   QE**2.0e0_DP ) / ( 1.0e0_DP + shield * Seg_array(i_ion2)%sig )  
                    
                        if(i_ion1==i_part) then                    
                            temp_vp = temp_vp                                                                   + &
                            &  (( 1.0e0_DP + shield * Seg_array(i_ion2)%sig ) * Seg_array(i_ion2)%q**2.0e0_DP   * &
                            &   QE**2.0e0_DP                                                                    - &
                            &   Comp_array(i_ion1)%nm * Seg_array(i_ion2)%q**2.0e0_DP * QE**2.0e0_DP            * &
                            &   Seg_array(i_ion2)%sig * shield_lp )                                             / &
                            &   ( 1.0e0_DP + shield * Seg_array(i_ion2)%sig )**2.0e0_DP
                        else
                            temp_vp = temp_vp                                                                   - &
                            &   (Comp_array(i_ion1)%nm * Seg_array(i_ion2)%q**2.0e0_DP * QE**2.0e0_DP           * &
                            &   Seg_array(i_ion2)%sig * shield_lp )                                             / &
                            &   ( 1.0e0_DP + shield * Seg_array(i_ion2)%sig )**2.0e0_DP
                        end if
                    
                    end do
                end if 
            end do
        end do 
       
        umsa_p = temp_u * temp_vp + temp_v * temp_up
        
        temp_u  = PI * omega_l * pn_l**2.0e0_DP 
        temp_up = PI * (omega_lp*pn_l**2.0e0_DP + 2.0e0_DP * pn_l * pn_lp * omega_l)
        temp_v  = 2.0e0_DP * del_l
        temp_vp = 2.0e0_DP * del_lp

        umsa_p = umsa_p + (temp_v * temp_up - temp_u * temp_vp) / temp_v**2.0e0_DP
        umsa_p = -umsa_p * v / (4.0e0_DP * PI * E0 * l_dielec) 
    
        a_res = umsa_p + 3.0e0_DP * shield**2.0e0_DP * shield_lp * KB * t * v / 3.0e0_DP / PI

        a_res = a_res / KB / t

        return
    end function A_ion_dn
!***************************************************************************************************  
    function A_ion_dv( ) result(a_res)
        implicit none

        real(kind=DP)           ::  a_res

        real(kind=DP)           ::  l_dielec, dl_dielecv
        real(kind=DP)               ::  c_shield, shield
        
        real(kind=DP), parameter    ::  SHIELD_CONV=1.0e-12_DP
        
        real(kind=DP)               ::  qsum_l, del_l, omega_l, pn_l, shield_l
        real(kind=DP)               ::  umsa, qtemp(1:nctypes,1:nstypes,1:20)
        
        real(kind=DP)               ::  del_lp, omega_lp, pn_lp, shield_lp, shield_temp
        real(kind=DP)               ::  qtempp(1:nctypes,1:nstypes,1:20), umsa_p
        
        real(kind=DP)               ::  temp_u, temp_v, temp_up, temp_vp, temp_w, temp_wp, temp_c
        
        integer                     ::  i_ion1, i_ion2, i_ion3
        integer                     ::  it_tot

        if(.not.ion_switch) then
            a_res = 0.0e0_DP
            return
        end if

        l_dielec   = Diel()
        dl_dielecv = dDielv()
        
        qtemp = 0.0e0_DP
!Do full A_ion calc first...      
        c_shield = QE**2.0e0_DP / (4.0e0_DP * l_dielec * E0 *  KB * t * v)        

        qsum_l = 0.0e0_DP
             
        do i_ion1=1, nctypes
            do i_ion2=1, nstypes
                if( Seg_array(i_ion2)%charged ) then
                    do i_ion3=1, Comp_array(i_ion1)%comp(i_ion2)
                        qsum_l = qsum_l + Comp_array(i_ion1)%nm * Seg_array(i_ion2)%q**2.0e0_DP                        
                    end do
                end if 
            end do
        end do
        
        !initial shield guess
        shield = 0.5e0_DP * dsqrt( c_shield * qsum_l )

        !set delta  - this is constant for iteration       
        del_l = 0.0e0_DP

        do i_ion1=1, nctypes
            do i_ion2=1, nstypes
                if(Seg_array(i_ion2)%charged) then
                    do i_ion3=1, Comp_array(i_ion1)%comp(i_ion2)
                        del_l = del_l + Comp_array(i_ion1)%nm * Seg_array(i_ion2)%sig**3.0e0_DP                   
                    end do
                end if 
            end do
        end do
        
        del_l = 1.0e0_DP - PI / 6.0e0_DP / v * del_l            

        !iterate parameters to find shield
        it_tot = 0       
        
        do             
            omega_l = 0.0e0_DP
            
            do i_ion1=1, nctypes
                do i_ion2=1, nstypes
                    if(Seg_array(i_ion2)%charged) then
                        do i_ion3=1, Comp_array(i_ion1)%comp(i_ion2)
                            omega_l = omega_l + Comp_array(i_ion1)%nm * Seg_array(i_ion2)%sig**3.0e0_DP / &
                                &   (1.0e0_DP + shield * Seg_array(i_ion2)%sig)
                        end do
                    end if 
                end do
            end do

            omega_l = 1.0e0_DP + PI / v / 2.e0_DP / del_l * omega_l

            pn_l = 0.0e0_DP
            
            do i_ion1=1, nctypes
                do i_ion2=1, nstypes
                    if(Seg_array(i_ion2)%charged) then
                        do i_ion3=1, Comp_array(i_ion1)%comp(i_ion2)
!added *QE
                            pn_l = pn_l + Comp_array(i_ion1)%nm * Seg_array(i_ion2)%sig *       &
                            &   Seg_array(i_ion2)%q*QE / (1.0e0_DP + shield * Seg_array(i_ion2)%sig )
                           
                        end do
                    end if 
                end do
            end do
           
            pn_l = 1.0e0_DP / omega_l / v * pn_l 

            shield_l = 0.0e0_DP
            
            do i_ion1=1, nctypes
                do i_ion2=1, nstypes
                    if(Seg_array(i_ion2)%charged) then
                        do i_ion3=1, Comp_array(i_ion1)%comp(i_ion2)
                            qtemp(i_ion1,i_ion2,i_ion3) =                                   &
                            &       ( Seg_array(i_ion2)%q * QE - (PI / 2.0e0_DP / del_l) *  &
                            &       Seg_array(i_ion2)%sig**2.0e0_DP * pn_l ) / (1.0e0_DP + shield * Seg_array(i_ion2)%sig) 

                            shield_l = shield_l + Comp_array(i_ion1)%nm * qtemp(i_ion1,i_ion2,i_ion3)**2.0e0_DP                     
                        end do
                    end if 
                end do
            end do
       
            shield_l = dsqrt( PI  / 4.0e0_DP / PI / E0 / l_dielec / KB / T / v * shield_l )
          
            it_tot = it_tot + 1
                    
            if(dabs( shield - shield_l ) / dabs( shield_l )  < SHIELD_CONV) then
                shield = shield_l
                exit
            else
                shield = shield_l
            end if    
        end do 
        
        umsa = 0.0e0_DP

        do i_ion1=1, nctypes
            do i_ion2=1, nstypes
                if(Seg_array(i_ion2)%charged) then
                    do i_ion3=1, Comp_array(i_ion1)%comp(i_ion2)
                        umsa = umsa + (Comp_array(i_ion1)%nm * Seg_array(i_ion2)%q**2.0e0_DP* QE**2.0e0_DP )  /   &
                        &   ( 1.0e0_DP + shield * Seg_array(i_ion2)%sig )  
                    end do
                end if 
            end do
        end do 
  
        umsa = -1.0e0_DP * v / 4.0e0_DP / PI / E0 / l_dielec * (shield / v * umsa +  &
        &       PI / 2.0e0_DP / del_l * omega_l * pn_l**2.0e0_DP )
    
!shield converged, do differential
        del_lp = 0.0e0_DP
        
        do i_ion1=1, nctypes
            do i_ion2=1, nstypes
                if(Seg_array(i_ion2)%charged) then
                    do i_ion3=1, Comp_array(i_ion1)%comp(i_ion2)
                        del_lp = del_lp + Comp_array(i_ion1)%nm * Seg_array(i_ion2)%sig**3.0e0_DP                   
                    end do
                end if 
            end do
        end do
      
        del_lp = del_lp * PI / 6.0e0_DP / v**2.0e0_DP

!initial shield guess        
        temp_u  = 0.0e0_DP
        
        do i_ion1=1, nctypes
            do i_ion2=1, nstypes
                if( Seg_array(i_ion2)%charged ) then
                    do i_ion3=1, Comp_array(i_ion1)%comp(i_ion2)
                        temp_u = temp_u + Comp_array(i_ion1)%nm * Seg_array(i_ion2)%q**2.0e0_DP                        
                    end do
                end if 
            end do
        end do
        
        temp_v = QE**2.0e0_DP * temp_u / (4.0e0_DP * E0 *  KB * t )
        shield_lp = 0.25e0_DP * (temp_v / (l_dielec*v))**(-0.5e0_DP)             * &
        &   (-temp_v * (dl_dielecv*v + l_dielec) / (l_dielec**2.0e0_DP * v**2.0e0_DP))       

        do !shield differential
            temp_u = PI / 2.0e0_DP / del_l / v
            temp_up = -PI * (del_l + del_lp*v) / (2.0e0_DP * v**2.0e0_DP * del_l**2.0e0_DP)

            temp_v  = 0.0e0_DP 
            temp_vp = 0.0e0_DP

            do i_ion1=1, nctypes
                do i_ion2=1, nstypes
                    if(Seg_array(i_ion2)%charged) then
                        do i_ion3=1, Comp_array(i_ion1)%comp(i_ion2)
                        
                            temp_v = temp_v + Comp_array(i_ion1)%nm * Seg_array(i_ion2)%sig**3.0e0_DP / &
                                &   (1.0e0_DP + shield * Seg_array(i_ion2)%sig)
                        
                            temp_vp = temp_vp - Comp_array(i_ion1)%nm * shield_lp                * &
                            & Seg_array(i_ion2)%sig**4.0e0_DP                                    / &
                            & (1.0e0_DP + shield * Seg_array(i_ion2)%sig)**2.0e0_DP
                        end do
                    end if 
                end do
            end do     
    
            omega_lp = temp_u * temp_vp + temp_v * temp_up
 
            temp_u  = 1.0e0_DP /omega_l / v  
            temp_up = -(omega_lp*v + omega_l) / (v**2.0e0_DP * omega_l**2.0e0_DP)
            temp_vp = 0.0e0_DP
            temp_v  = 0.0e0_DP
            
            
            do i_ion1=1, nctypes
                do i_ion2=1, nstypes
                    if(Seg_array(i_ion2)%charged) then
                        do i_ion3=1, Comp_array(i_ion1)%comp(i_ion2)
    
                            temp_v = temp_v + Comp_array(i_ion1)%nm * Seg_array(i_ion2)%sig *       &
                            &   Seg_array(i_ion2)%q*QE / (1.0e0_DP + shield * Seg_array(i_ion2)%sig )
                            
                            temp_vp = temp_vp - shield_lp * Seg_array(i_ion2)%sig**2.0e0_DP         * &
                            &   Comp_array(i_ion1)%nm * Seg_array(i_ion2)%q*QE                      / &
                            &   (1.0e0_DP + shield * Seg_array(i_ion2)%sig)**2.0e0_DP                        
                        end do
                    end if 
                end do
            end do        
            
            pn_lp = temp_u * temp_vp + temp_v * temp_up
     
            do i_ion1=1, nctypes
                do i_ion2=1, nstypes
                    if(Seg_array(i_ion2)%charged) then
                        do i_ion3=1, Comp_array(i_ion1)%comp(i_ion2)
                        
                        
                            temp_v  = 1.0e0_DP + shield * Seg_array(i_ion2)%sig
                            temp_vp = shield_lp * Seg_array(i_ion2)%sig 

                            temp_u  = Seg_array(i_ion2)%q * QE - PI             * &
                            &   Seg_array(i_ion2)%sig**2.0e0_DP * pn_l          / &
                            &   (2.0e0_DP * del_l)

                            temp_up = PI*Seg_array(i_ion2)%sig**2.0e0_DP        * &
                            &   (pn_l*del_lp - del_l*pn_lp)                     / &
                            &   (2.0e0_DP * del_l**2.0e0_DP)

                            qtempp(i_ion1,i_ion2,i_ion3) = (temp_v * temp_up - temp_u * temp_vp) / temp_v**2.0e0_DP
                        end do
                    end if 
                end do
            end do  

            temp_c  = 1.0e0_DP / (4.0e0_DP * E0 * KB * T)
            temp_u  = 1.0e0_DP / (l_dielec * v)
            temp_up = -1.0e0_DP * (dl_dielecv * v + l_dielec) / (l_dielec**2.0e0_DP * v**2.0e0_DP) 
  
            temp_v  = 0.0e0_DP
            temp_vp = 0.0e0_DP           
            
            do i_ion1=1, nctypes
                do i_ion2=1, nstypes
                    if(Seg_array(i_ion2)%charged) then
                        do i_ion3=1, Comp_array(i_ion1)%comp(i_ion2)
                            temp_v = temp_v + Comp_array(i_ion1)%nm * qtemp(i_ion1,i_ion2,i_ion3)**2.0e0_DP 
                            
                            temp_vp = temp_vp + 2.0e0_DP * Comp_array(i_ion1)%nm        * &
                            &   qtemp(i_ion1,i_ion2,i_ion3) * qtempp(i_ion1,i_ion2,i_ion3)
                        end do
                    end if 
                end do
            end do 
            
            temp_w = temp_c*(temp_v*temp_up + temp_u*temp_vp)  
    
            shield_temp = 0.5e0_DP * (temp_c*temp_v/(l_dielec*v))**(-0.5e0_DP) * temp_w

            if(dabs( shield_lp - shield_temp ) / dabs( shield_temp )  < SHIELD_CONV) then
                shield_lp = shield_temp
                exit
            else
                shield_lp = shield_temp
            end if      
                        
        end do

        temp_u  = -1.0e0_DP * v /(4.0e0_DP * PI * E0 * l_dielec)
        temp_up = -1.0e0_DP /(4.0e0_DP * PI * E0)                             * &
        &   (l_dielec - v * dl_dielecv) / l_dielec**2.0e0_DP    

        temp_v  = 0.0e0_DP
        temp_vp = 0.0e0_DP
        
        temp_w  = PI/(2.0e0_DP*del_l) * omega_l*pn_l**2.0e0_DP
        temp_wp =  -PI*del_lp/(2.0e0_DP*del_l**2.0e0_DP) * omega_l*pn_l**2.0e0_DP   + &
        &           PI/(2.0e0_DP*del_l) * omega_lp*pn_l**2.0e0_DP                   + &
        &           PI/del_l * omega_l * pn_l*pn_lp

        do i_ion1=1, nctypes
            do i_ion2=1, nstypes
                if(Seg_array(i_ion2)%charged) then
                    do i_ion3=1, Comp_array(i_ion1)%comp(i_ion2)
                        temp_v = temp_v + (Comp_array(i_ion1)%nm * Seg_array(i_ion2)%q**2.0e0_DP    * &
                        &   QE**2.0e0_DP ) / ( 1.0e0_DP + shield * Seg_array(i_ion2)%sig )  
                    
                        temp_vp = -shield*Seg_array(i_ion2)%sig*Seg_array(i_ion2)%q**2.0e0_DP*QE**2.0e0_DP      / &
                        &   ( 1.0e0_DP + shield * Seg_array(i_ion2)%sig )    
                    end do
                end if 
            end do
        end do 
       
        temp_vp = shield/v * temp_vp  +  temp_v * (v*shield_lp - shield)/v**2.0e0_DP
        temp_v = temp_v * shield / v
        umsa_p = temp_u * (temp_vp + temp_wp) + temp_up * (temp_v + temp_w)
    
        a_res = -(umsa_p + 3.0e0_DP * shield**2.0e0_DP * shield_lp * KB * t * v / 3.0e0_DP / PI       + &
        &       shield**3.0e0_DP * KB * t / 3.0e0_DP / PI)

        return
    end function A_ion_dv
!***************************************************************************************************  

!***************************************************************************************************
end module Ion
!***************************************************************************************************
