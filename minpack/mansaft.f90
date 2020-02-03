!***************************************************************************************************
!***************************************************************************************************
!   SAFT Module to:
!       1. Calculate A assoc
!
!       Analtical derivatives added March 2017
!
!***************************************************************************************************
!
!***************************************************************************************************
!***************************************************************************************************
!   SAFT Module to:
!       1. Calculate A assoc
!
!       Analtical derivatives added March 2017
!
!***************************************************************************************************
!
!***************************************************************************************************
module Assoc_mod
!***************************************************************************************************
!Modules
!=======    
    use Types_mod           ! Definitions of types and double precision
    use Global_mod          ! Important global parameters 
!***************************************************************************************************
    implicit none
!***************************************************************************************************
    contains
!***************************************************************************************************    
    function A_assoc() result(a_res)      ! [2] 21
        implicit none
        
        integer             ::  i_a1,i_a2,i_a3
        real(kind=DP)       ::  a_res
        
        a_res = 0.0e0_DP

        if(assoc_switch) then        
            call Mass_action()
      
            do i_a1=1, nctypes
            do i_a2=1, nstypes
            
            if(Comp_array(i_a1)%comp(i_a2)>0) then      !added Jan 2017   
               
            do i_a3=1, Seg_array(i_a2)%nassoc
                  a_res = a_res + Comp_array(i_a1)%xi * Comp_array(i_a1)%comp(i_a2) *       &
                  & Seg_array(i_a2)%nassoctyp(i_a3) * (dlog(x(i_a1, i_a2,i_a3)) + (1.0e0_DP  &
                  & -x(i_a1,i_a2,i_a3))/2.0e0_DP)            
            end do

            end if
            
            end do
            end do
        end if   

        return
    end function A_assoc
!***************************************************************************************************
    function A_assoc_dn() result(adn_res)      ! [2] 21
        implicit none
        
        integer             ::  i_a1,i_a2,i_a3
        real(kind=DP)       ::  adn_res
  
        adn_res = 0.0e0_DP

        if(assoc_switch) then        
            call Mass_action()
            call Mass_action_dn()
      
            do i_a1=1, nctypes
            do i_a2=1, nstypes
            
            if(Comp_array(i_a1)%comp(i_a2)>0) then      !added Jan 2017   
               
            do i_a3=1, Seg_array(i_a2)%nassoc                  
                  adn_res = adn_res +   &
                  & nsum * Comp_array(i_a1)%xi * Comp_array(i_a1)%comp(i_a2) * Seg_array(i_a2)%nassoctyp(i_a3)   * &
                  & (dxn(i_a1, i_a2,i_a3)/x(i_a1,i_a2,i_a3) - dxn(i_a1, i_a2,i_a3)/2.0e0_DP )                    + &
                  & nsum * Comp_array(i_a1)%dxin * Comp_array(i_a1)%comp(i_a2) * Seg_array(i_a2)%nassoctyp(i_a3) * &
                  & (dlog(x(i_a1, i_a2,i_a3)) + (1.0e0_DP - x(i_a1,i_a2,i_a3))/2.0e0_DP)                          + &
                  & Comp_array(i_a1)%xi * Comp_array(i_a1)%comp(i_a2) * Seg_array(i_a2)%nassoctyp(i_a3)          * &
                  & (dlog(x(i_a1, i_a2,i_a3)) + (1.0e0_DP - x(i_a1,i_a2,i_a3))/2.0e0_DP) 
            end do

            end if
            
            end do
            end do
        end if   
     
        return
    end function A_assoc_dn
!***************************************************************************************************
    function A_assoc_dv() result(adv_res)      ! [2] 21
        implicit none
        
        integer             ::  i_a1,i_a2,i_a3
        real(kind=DP)       ::  adv_res
  
        adv_res = 0.0e0_DP

        if(assoc_switch) then               
            call Mass_action()
            call Mass_action_dv()
      
            do i_a1=1, nctypes
            do i_a2=1, nstypes
          
            if(Comp_array(i_a1)%comp(i_a2)>0) then      !added Jan 2017   
               
            do i_a3=1, Seg_array(i_a2)%nassoc                  
                  adv_res = adv_res + Comp_array(i_a1)%xi * Comp_array(i_a1)%comp(i_a2)             * &
                  & Seg_array(i_a2)%nassoctyp(i_a3) * (dxv(i_a1, i_a2,i_a3) / x(i_a1, i_a2,i_a3)    - &
                  & dxv(i_a1,i_a2,i_a3)/2.0e0_DP)     
            end do

            end if
            
            end do
            end do
        end if   
     
        adv_res = -adv_res * KB * NA * t
     
        return
    end function A_assoc_dv
!***************************************************************************************************
    subroutine mass_action() !BRUTE FORCE ITERATION AT PRESENT
        implicit none
        
        integer                         ::  mac1, mac2, mac3, mbc1, mbc2, mbc3, n_it
        real(kind=DP),parameter         ::  CRITERIA=1.0e-15_DP   
        logical                         ::  converged
        real(kind=DP),allocatable       ::  x_old(:,:,:)
        real(kind=DP)                   ::  xsum
        
        if(allocated(x))        deallocate(x)
        if(allocated(delx))     deallocate(delx)
        if(allocated(x_old))    deallocate(x_old)   
        
        allocate(x(nctypes, nstypes, 3), delx(nctypes, nstypes, 3, nctypes, nstypes, 3), &
        &   x_old(nctypes, nstypes, 3))

        call delta_x( )                  
        
        x = 0.5e0_DP  
  
        n_it=0
        converged=.false.
      
        do while(.not.converged)      
            x_old=x     
            converged=.true.
            
            do mac1=1, nctypes
            do mac2=1, nstypes
            
            if(Comp_array(mac1)%comp(mac2)>0) then    !added Jan 2017   
    
            do mac3=1, Seg_array(mac2)%nassoc
        
                xsum = 1.0e0_DP
                
                do mbc1=1, nctypes
                do mbc2=1, nstypes
                do mbc3=1, Seg_array(mbc2)%nassoc
                
                    xsum = xsum + rho * Comp_array(mbc1)%xi * Comp_array(mbc1)%comp(mbc2) * &
                    &   Seg_array(mbc2)%nassoctyp(mbc3) * x_old(mbc1, mbc2, mbc3) *         &
                    &   delx( mac1, mac2, mac3, mbc1, mbc2, mbc3 )                
                end do
                end do
                end do

                x(mac1, mac2, mac3) = 1.0e0_DP / xsum
                if(dabs(x(mac1, mac2, mac3) - x_old(mac1, mac2, mac3)) > CRITERIA) converged=.false.                    
            end do
            
            end if   
                     
            end do
            end do
   
            n_it = n_it + 1
        end do  

        return
    end subroutine Mass_action
!***************************************************************************************************
    subroutine mass_action_dn() 
        implicit none
        
        integer                         ::  mac1, mac2, mac3, mbc1, mbc2, mbc3
        real(kind=DP),allocatable       ::  dxn_old(:,:,:)
        real(kind=DP)                   ::  xpart2, xpart2p  
        real(kind=DP),parameter         ::  CRITERIA=1.0e-20_DP 
        logical                         ::  CONVERGED

        if(allocated(dxn))          deallocate(dxn)
        if(allocated(dxn_old))      deallocate(dxn_old)
        if(allocated(ddelxn))       deallocate(ddelxn)  
        
        allocate(ddelxn(nctypes, nstypes, 3, nctypes, nstypes, 3), dxn(nctypes, nstypes, 3), dxn_old(nctypes, nstypes, 3))

        call ddelta_xn( )                  

        xpart2  = 0.0e0_DP
        xpart2p = 0.0e0_DP    
           
        dxn(:,:,:) = 1.0e0_DP   
        CONVERGED=.false.   
           
        do while(.not.CONVERGED)
           
            dxn_old=dxn  
            CONVERGED=.true.

        do mac1=1, nctypes
        do mac2=1, nstypes
        
        if(Comp_array(mac1)%comp(mac2)>0) then    !added Jan 2017   

        do mac3=1, Seg_array(mac2)%nassoc
            
            xpart2  = 0.0e0_DP
            xpart2p = 0.0e0_DP
            
            do mbc1=1, nctypes
            do mbc2=1, nstypes
            do mbc3=1, Seg_array(mbc2)%nassoc
                xpart2 = xpart2 + Comp_array(mbc1)%xi * Comp_array(mbc1)%comp(mbc2) * &
                &   Seg_array(mbc2)%nassoctyp(mbc3) * x(mbc1, mbc2, mbc3) * delx( mac1, mac2, mac3, mbc1, mbc2, mbc3 )    
                
                xpart2p = xpart2p + Comp_array(mbc1)%comp(mbc2) * Seg_array(mbc2)%nassoctyp(mbc3) * &
                & (Comp_array(mbc1)%dxin * x(mbc1, mbc2, mbc3) *  delx( mac1, mac2, mac3, mbc1, mbc2, mbc3) + &
                &  Comp_array(mbc1)%xi * dxn(mbc1, mbc2, mbc3) *  delx( mac1, mac2, mac3, mbc1, mbc2, mbc3) + &
                &  Comp_array(mbc1)%xi * x(mbc1, mbc2, mbc3) *  ddelxn( mac1, mac2, mac3, mbc1, mbc2, mbc3))                         
            end do
            end do
            end do

            dxn(mac1, mac2, mac3) = (-1.0e0_DP * x(mac1, mac2, mac3) * drhon * (xpart2 +  Comp_array(mac1)%xi           * &
            &   Comp_array(mac1)%comp(mac2)                                                                             * &
            &   Seg_array(mac2)%nassoctyp(mac3) * x(mac1, mac2, mac3) * delx( mac1, mac2, mac3, mac1, mac2, mac3 ))     - &
            &   x(mac1, mac2, mac3) * rho *  (xpart2p + Comp_array(mac1)%comp(mac2) * Seg_array(mac2)%nassoctyp(mac3)   * &
            &   (Comp_array(mac1)%xi * x(mac1, mac2, mac3) * ddelxn( mac1, mac2, mac3, mac1, mac2, mac3 )               + &
            &   Comp_array(mac1)%dxin * x(mac1, mac2, mac3) * delx( mac1, mac2, mac3, mac1, mac2, mac3 )) ))            / &
            &   (1.0e0_DP + rho * (xpart2 +  Comp_array(mac1)%xi * Comp_array(mac1)%comp(mac2)                          * &
            &   Seg_array(mac2)%nassoctyp(mac3) * x(mac1, mac2, mac3) * delx( mac1, mac2, mac3, mac1, mac2, mac3 ))     + &
            &   x(mac1, mac2, mac3) * rho * Comp_array(mac1)%comp(mac2) * Seg_array(mac2)%nassoctyp(mac3)               * &
            &   Comp_array(mac1)%xi * delx( mac1, mac2, mac3, mac1, mac2, mac3 )        ) 

            if(dabs((dxn(mac1, mac2, mac3) - dxn_old(mac1, mac2, mac3))/dxn(mac1, mac2, mac3)) > CRITERIA) converged=.false.               
        end do
        
        end if   
                 
        end do
        end do
     
        end do

        return
    end subroutine mass_action_dn    
!***************************************************************************************************
    subroutine mass_action_dv() !BRUTE FORCE ITERATION AT PRESENT
        implicit none
        
        integer                         ::  mac1, mac2, mac3, mbc1, mbc2, mbc3, n_it
        real(kind=DP),parameter         ::  CRITERIA=1.0e-15_DP   
        logical                         ::  converged
        real(kind=DP),allocatable       ::  dxv_old(:,:,:)
        real(kind=DP)                   ::  xsum, xsum2
        
        if(allocated(dxv))          deallocate(dxv)
        if(allocated(dxv_old))      deallocate(dxv_old)
        if(allocated(ddelxv))       deallocate(ddelxv)  
        
        allocate(ddelxv(nctypes, nstypes, 3, nctypes, nstypes, 3), dxv(nctypes, nstypes, 3), dxv_old(nctypes, nstypes, 3))

        call ddelta_xv( )                  
        
        dxv = 1.0e0_DP  
  
        n_it=0
        converged=.false.
       
        do while(.not.converged)
           
            dxv_old=dxv     
            converged=.true.
            
            do mac1=1, nctypes
            do mac2=1, nstypes
            
            if(Comp_array(mac1)%comp(mac2)>0) then    !added Jan 2017   
    
            do mac3=1, Seg_array(mac2)%nassoc
            
                xsum = 1.0e0_DP
                xsum2 = 0.0e0_DP
                
                do mbc1=1, nctypes
                do mbc2=1, nstypes
                do mbc3=1, Seg_array(mbc2)%nassoc
                    xsum = xsum + rho * Comp_array(mbc1)%xi * Comp_array(mbc1)%comp(mbc2)   * &
                    &   Seg_array(mbc2)%nassoctyp(mbc3) * x(mbc1, mbc2, mbc3)               * &
                    &   delx( mac1, mac2, mac3, mbc1, mbc2, mbc3 )      
                    
                    xsum2 = xsum2 - Comp_array(mbc1)%xi * Comp_array(mbc1)%comp(mbc2)                   * &
                    &   Seg_array(mbc2)%nassoctyp(mbc3)                                                 * &
                    &   (drhov * x(mbc1, mbc2, mbc3) * delx( mac1, mac2, mac3, mbc1, mbc2, mbc3 )       + &
                    &    rho * dxv(mbc1, mbc2, mbc3) * delx( mac1, mac2, mac3, mbc1, mbc2, mbc3 )       + &
                    &    rho * x(mbc1, mbc2, mbc3) * ddelxv( mac1, mac2, mac3, mbc1, mbc2, mbc3 ))
                end do
                end do
                end do

                dxv(mac1, mac2, mac3) = xsum2 / xsum**2.0e0_DP
                
                if(dabs(dxv(mac1, mac2, mac3) - dxv_old(mac1, mac2, mac3)) > CRITERIA) converged=.false.                    
            end do
           
            end if   
                     
            end do
            end do
  
            n_it = n_it + 1
        end do  

        return
    end subroutine Mass_action_dv
!***************************************************************************************************
    subroutine delta_x( )
        implicit none
        
        integer          ::  i_dx1, i_dx2, i_dx3, i_dx4, i_dx5, i_dx6, dxp, dxq
        real(kind=DP)    ::  f_x(1:nstypes, 1:3, 1:nstypes, 1:3)
        real(kind=DP)    ::  i_x(1:nctypes, 1:nstypes, 1:3, 1:nctypes, 1:nstypes, 1:3)
        
        do i_dx1=1, nstypes
        do i_dx2=1, Seg_array(i_dx1)%nassoc
            do i_dx3=1, nstypes
            do i_dx4=1, Seg_array(i_dx3)%nassoc               
                f_x(i_dx1, i_dx2, i_dx3, i_dx4) = dexp( ehb(i_dx1, i_dx2, i_dx3, i_dx4) / T) - 1.0e0_DP                                
            end do
            end do
        end do
        end do

        i_x = 0.0e0_DP
          
        do i_dx1=1, nctypes
        do i_dx2=1, nstypes
        do i_dx3=1, Seg_array(i_dx2)%nassoc
            do i_dx4=1, nctypes
            do i_dx5=1, nstypes
            do i_dx6=1, Seg_array(i_dx5)%nassoc
           
                do dxp=0, 10
                do dxq=0, 10 - dxp              
                    i_x(i_dx1, i_dx2, i_dx3, i_dx4, i_dx5, i_dx6) =         &
                    &   i_x(i_dx1, i_dx2, i_dx3, i_dx4, i_dx5, i_dx6)   +   &

                    &   (c_x(dxp, dxq) * (rhos * NA * sigx3) ** dxp      *   &
                    &   (t / epschij(i_dx1, i_dx4))**dxq)     
                end do
                end do
                
                delx(i_dx1, i_dx2, i_dx3, i_dx4, i_dx5, i_dx6) =  f_x(i_dx2, i_dx3, i_dx5, i_dx6) * &
                &   i_x(i_dx1, i_dx2, i_dx3, i_dx4, i_dx5, i_dx6) * khb(i_dx2, i_dx3, i_dx5, i_dx6)
                          
                !SHOULDN'T HAPPEN, BUT.... CAUSED BY ISSUE WITH i_x               
                if(delx(i_dx1, i_dx2, i_dx3, i_dx4, i_dx5, i_dx6) < 0.0e0_DP) then             
                    delx(i_dx1, i_dx2, i_dx3, i_dx4, i_dx5, i_dx6)=0.0e0_DP
                end if

            end do
            end do
            end do
        end do
        end do
        end do 
!print*, rhos, drhosv,i_x

        return
    end subroutine          
!***************************************************************************************************
    subroutine ddelta_xn( )
        implicit none
        
        integer          ::  i_dx1, i_dx2, i_dx3, i_dx4, i_dx5, i_dx6, dxp, dxq
        real(kind=DP)    ::  f_x(1:nstypes, 1:3, 1:nstypes, 1:3)
        real(kind=DP)    ::  di_xn(1:nctypes, 1:nstypes, 1:3, 1:nctypes, 1:nstypes, 1:3)
        
        do i_dx1=1, nstypes
        do i_dx2=1, Seg_array(i_dx1)%nassoc
            do i_dx3=1, nstypes
            do i_dx4=1, Seg_array(i_dx3)%nassoc               
                f_x(i_dx1, i_dx2, i_dx3, i_dx4) = dexp( ehb(i_dx1, i_dx2, i_dx3, i_dx4) / T) - 1.0e0_DP                                
            end do
            end do
        end do
        end do

        di_xn = 0.0e0_DP
          
        do i_dx1=1, nctypes
        do i_dx2=1, nstypes
        do i_dx3=1, Seg_array(i_dx2)%nassoc
            do i_dx4=1, nctypes
            do i_dx5=1, nstypes
            do i_dx6=1, Seg_array(i_dx5)%nassoc
           
                do dxp=0, 10
                do dxq=0, 10 - dxp              
                    di_xn(i_dx1, i_dx2, i_dx3, i_dx4, i_dx5, i_dx6) =         &
                    &   di_xn(i_dx1, i_dx2, i_dx3, i_dx4, i_dx5, i_dx6)   +   &

                    &   (c_x(dxp, dxq) * &                
                                       
                    &   (rhos*dsigx3n + drhosn*sigx3) * dxp * NA *  &
                    &   (rhos * NA * sigx3) ** (dxp-1)  *   &                 
                                        
                    &   (t / epschij(i_dx1, i_dx4))**dxq)  
                                      
                                                                   
                end do
                end do

                ddelxn(i_dx1, i_dx2, i_dx3, i_dx4, i_dx5, i_dx6) =  f_x(i_dx2, i_dx3, i_dx5, i_dx6) * &
                &   di_xn(i_dx1, i_dx2, i_dx3, i_dx4, i_dx5, i_dx6) * khb(i_dx2, i_dx3, i_dx5, i_dx6)
            end do
            end do
            end do
        end do
        end do
        end do 

        return
    end subroutine ddelta_xn        
!***************************************************************************************************
    subroutine ddelta_xv( )
        implicit none
        
        integer          ::  i_dx1, i_dx2, i_dx3, i_dx4, i_dx5, i_dx6, dxp, dxq
        real(kind=DP)    ::  f_x(1:nstypes, 1:3, 1:nstypes, 1:3)
        real(kind=DP)    ::  di_xv(1:nctypes, 1:nstypes, 1:3, 1:nctypes, 1:nstypes, 1:3)
        
        do i_dx1=1, nstypes
        do i_dx2=1, Seg_array(i_dx1)%nassoc
            do i_dx3=1, nstypes
            do i_dx4=1, Seg_array(i_dx3)%nassoc               
                f_x(i_dx1, i_dx2, i_dx3, i_dx4) = dexp( ehb(i_dx1, i_dx2, i_dx3, i_dx4) / T) - 1.0e0_DP                                
            end do
            end do
        end do
        end do

        di_xv = 0.0e0_DP
          
        do i_dx1=1, nctypes
        do i_dx2=1, nstypes
        do i_dx3=1, Seg_array(i_dx2)%nassoc
            do i_dx4=1, nctypes
            do i_dx5=1, nstypes
            do i_dx6=1, Seg_array(i_dx5)%nassoc
           
                do dxp=0, 10
                do dxq=0, 10 - dxp              
                    di_xv(i_dx1, i_dx2, i_dx3, i_dx4, i_dx5, i_dx6) =         &
                    &   di_xv(i_dx1, i_dx2, i_dx3, i_dx4, i_dx5, i_dx6)   +   &

                    &   (c_x(dxp, dxq)                                    *   &
                    &   real(dxp) * drhosv * sigx3                        *   &
                    &   (rhos * NA * sigx3) ** (dxp-1)            *   &
                    &   (t / epschij(i_dx1, i_dx4))**dxq)                                                 
                end do
                end do
                
                ddelxv(i_dx1, i_dx2, i_dx3, i_dx4, i_dx5, i_dx6) =  f_x(i_dx2, i_dx3, i_dx5, i_dx6) * &
                &   di_xv(i_dx1, i_dx2, i_dx3, i_dx4, i_dx5, i_dx6) * khb(i_dx2, i_dx3, i_dx5, i_dx6) * NA
            end do
            end do
            end do
        end do
        end do
        end do 

        return
    end subroutine          
!***************************************************************************************************

end module Assoc_mod
!***************************************************************************************************
!***************************************************************************************************
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
module Ion_mod
!***************************************************************************************************
!Modules
!=======    
    use Types_mod           ! Definitions of types and double precision
    use Global_mod          ! Important global parameters     
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
            
            x_solv=0.0e0_DP
            ddiel=0.0e0_DP
           
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
                            qtemp = ( Seg_array(i_ion2)%q * QE - (PI / 2.0e0_DP / del_l) *  &
                            &       Seg_array(i_ion2)%sig**2.0e0_DP * pn_l ) / (1.0e0_DP + shield * Seg_array(i_ion2)%sig)                     
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
end module Ion_mod
!***************************************************************************************************
!***************************************************************************************************
!   SAFT Module to:
!       1. Calculate common differentials
!           a. with respect to Ni -> chemical potentials
!           b. with respect to V  -> pressure
!***************************************************************************************************
!
!***************************************************************************************************
module Diff_mod
!***************************************************************************************************
!Modules
!=======    
    use Types_mod           ! Definitions of types and double precision
    use Global_mod          ! Important global parameters
    use Setup_mod           ! To setup the system
!***************************************************************************************************
    implicit none
!***************************************************************************************************
    contains
!***************************************************************************************************    
    subroutine Ndiffs( i_part ) !For mu calculations
        implicit none
        
        integer, intent(in) ::  i_part
        real(kind=DP)       ::  t_dii
        real(kind=DP)       ::  diffv, diffu, diffvp, diffup
        real(kind=DP)       ::  tdiff1, tdiff2, tdiff3
        integer             ::  i_ndiff1, i_ndiff2
       
        call Sys_setup() !ensure eveything initialised   
            
        !sum of particles
        nsum = sum(Comp_array(:)%nm)
        
        !mole fractions
        do i_ndiff1=1,nctypes
            if(i_ndiff1==i_part) then
                Comp_array(i_ndiff1)%dxin = (1.0e0_DP - Comp_array(i_ndiff1)%xi) / nsum
            else
                Comp_array(i_ndiff1)%dxin = -Comp_array(i_ndiff1)%xi / nsum
            end if
        end do
      
        !rho
        drhon = 1.0e0_DP / v / NA
      
        !mole fraction of spherical segments
        diffv  = 0.0e0_DP
        diffvp = 0.0e0_DP
        
        do i_ndiff1=1, nctypes
        do i_ndiff2=1, nstypes
            diffv = diffv + Comp_array(i_ndiff1)%comp(i_ndiff2) * Seg_array(i_ndiff2)%nseg *Comp_array(i_ndiff1)%xi &
            &   * Seg_array(i_ndiff2)%sf
            
            diffvp = diffvp + Comp_array(i_ndiff1)%comp(i_ndiff2) * Seg_array(i_ndiff2)%nseg *Comp_array(i_ndiff1)%dxin &
            &   * Seg_array(i_ndiff2)%sf
        end do
        end do

        do i_ndiff1=1, nstypes       
            diffu  = 0.0e0_DP
            diffup = 0.0e0_DP
            
            do i_ndiff2=1, nctypes
                diffu = diffu +  (Comp_array(i_ndiff2)%comp(i_ndiff1) * Seg_array(i_ndiff1)%sf *      &  
                &                  Seg_array(i_ndiff1)%nseg * Comp_array(i_ndiff2)%xi) 
                
                diffup = diffup +  (Comp_array(i_ndiff2)%comp(i_ndiff1) * Seg_array(i_ndiff1)%sf *      &  
                &                  Seg_array(i_ndiff1)%nseg * Comp_array(i_ndiff2)%dxin) 
            end do

            Seg_array(i_ndiff1)%dxsin = (diffv*diffup - diffu*diffvp) / diffv**2.0e0_DP                    
        end do 

        dsigx3n = 0.0e0_DP
        
        !associaion sigx3
        do i_ndiff1=1, nstypes
        do i_ndiff2=1, nstypes
            dsigx3n = dsigx3n + ( Seg_array(i_ndiff1)%xsi * Seg_array(i_ndiff2)%dxsin    +   &
            &                     Seg_array(i_ndiff2)%xsi * Seg_array(i_ndiff1)%dxsin )  *   &
            &                     sig(i_ndiff1, i_ndiff2)**3.0e0_DP
        end do
        end do
       
        !Mono pieces
        !rhos
        diffu  = 0.0e0_DP
        diffup = 0.0e0_DP

        do i_ndiff1=1,nctypes
        do i_ndiff2=1,nstypes
            diffu = diffu + Comp_array(i_ndiff1)%xi * Comp_array(i_ndiff1)%comp(i_ndiff2) &  
            &   * Seg_array(i_ndiff2)%sf * Seg_array(i_ndiff2)%nseg
            diffup = diffup + Comp_array(i_ndiff1)%dxin * Comp_array(i_ndiff1)%comp(i_ndiff2) &  
            &   * Seg_array(i_ndiff2)%sf * Seg_array(i_ndiff2)%nseg
        end do
        end do

        drhosn = diffup*rho + diffu*drhon  !moles per m^3

        !zl
        dzln(:) = 0.0e0_DP
        
        do i_ndiff1=1, nstypes
            t_dii = dij(i_ndiff1,i_ndiff1) 

            dzln(0) = dzln(0) + PI * rhos / 6.0e0_DP * Seg_array(i_ndiff1)%dxsin                      + &
            &   PI * drhosn / 6.0e0_DP * Seg_array(i_ndiff1)%xsi

            dzln(1) = dzln(1) + PI * rhos / 6.0e0_DP * Seg_array(i_ndiff1)%dxsin * t_dii              + &
            &   PI * drhosn / 6.0e0_DP * Seg_array(i_ndiff1)%xsi * t_dii

            dzln(2) = dzln(2) + PI * rhos / 6.0e0_DP * Seg_array(i_ndiff1)%dxsin * t_dii**2.0e0_DP    + &
            &   PI * drhosn / 6.0e0_DP * Seg_array(i_ndiff1)%xsi * t_dii**2.0e0_DP

            dzln(3) = dzln(3) + PI * rhos / 6.0e0_DP * Seg_array(i_ndiff1)%dxsin * t_dii**3.0e0_DP    + &
            &   PI * drhosn / 6.0e0_DP * Seg_array(i_ndiff1)%xsi * t_dii**3.0e0_DP              

        end do
     
        dzln = dzln * NA
   
        !zetax
        tdiff1 = 0.0e0_DP
        tdiff2 = 0.0e0_DP
        tdiff3 = 0.0e0_DP
        
        do i_ndiff1=1, nstypes
        do i_ndiff2=1, nstypes
            tdiff1 = tdiff1 + Seg_array(i_ndiff1)%xsi * Seg_array(i_ndiff2)%xsi * dij(i_ndiff1,i_ndiff2)**3.0e0_DP
            tdiff2 = tdiff2 + Seg_array(i_ndiff1)%dxsin * Seg_array(i_ndiff2)%xsi * dij(i_ndiff1,i_ndiff2)**3.0e0_DP
            tdiff3 = tdiff3 + Seg_array(i_ndiff1)%xsi * Seg_array(i_ndiff2)%dxsin * dij(i_ndiff1,i_ndiff2)**3.0e0_DP
        end do
        end do
        
        dzetax_sumn = tdiff2 + tdiff3
        dzetaxn = (PI*drhosn/6.0e0_DP * tdiff1 + PI*rhos/6.0e0_DP * (tdiff2 + tdiff3)) * NA
        
        !zetabar
        diffu  = 0.0e0_DP
        diffup = 0.0e0_DP

        do i_ndiff1=1, nstypes
        do i_ndiff2=1, nstypes
           diffu  = diffu + Seg_array(i_ndiff1)%xsi * Seg_array(i_ndiff2)%xsi * sig(i_ndiff1,i_ndiff2)**3.0e0_DP  
           diffup = diffup + Seg_array(i_ndiff1)%dxsin * Seg_array(i_ndiff2)%xsi * sig(i_ndiff1,i_ndiff2)**3.0e0_DP  + &
           &                 Seg_array(i_ndiff1)%xsi * Seg_array(i_ndiff2)%dxsin * sig(i_ndiff1,i_ndiff2)**3.0e0_DP
        end do
        end do
        
        !Averaged packing fractions
        dzetabarn = NA * PI * rhos / 6.0e0_DP * diffup + NA * PI * drhosn / 6.0e0_DP * diffu

        !KHS
        diffu  = (1.d0 - zetax)**4.0e0_DP
        diffv  = 1.d0 + 4.0e0_DP * zetax + 4.0e0_DP * zetax**2 - 4.0e0_DP * zetax**3.0e0_DP + zetax**4.0e0_DP
        diffup = -4.0e0_DP * dzetaxn * (1.d0 - zetax)**3.0e0_DP
        diffvp = dzetaxn *( 4.0e0_DP + 8.0e0_DP * zetax - 12.0e0_DP * zetax**2.0e0_DP + 4.0e0_DP * zetax**3.0e0_DP)
         
        dKHSn = (diffv*diffup - diffu*diffvp)/diffv**2.0e0_DP

        return
    end subroutine Ndiffs
!*************************************************************************************************** 
    subroutine Vdiffs( ) !For pressure calculations
        implicit none
        
        real(kind=DP)       ::  t_dii
        real(kind=DP)       ::  diffv, diffu, diffvp, diffup
        real(kind=DP)       ::  tdiff1, tdiff2, tdiff3
        integer             ::  i_vdiff1, i_vdiff2
    
        call Sys_setup() !ensure eveything initialised   

        !sum of particles
        nsum = sum(Comp_array(:)%nm)
               
        !rho
        drhov = -rho / v 

        !Mono pieces
        !rhos
        drhosv = 0.e0_DP 

        do i_vdiff1=1,nctypes
        do i_vdiff2=1,nstypes
            drhosv = drhosv + Comp_array(i_vdiff1)%xi * Comp_array(i_vdiff1)%comp(i_vdiff2) &  
            &   * Seg_array(i_vdiff2)%sf * Seg_array(i_vdiff2)%nseg
        end do
        end do
        
        drhosv = drhosv * drhov

        !zl
        dzlv(:) = 0.0e0_DP
        
        do i_vdiff1=1,nstypes
            t_dii = dij(i_vdiff1,i_vdiff1)     
           
            dzlv(0) = dzlv(0) + Seg_array(i_vdiff1)%xsi
            dzlv(1) = dzlv(1) + Seg_array(i_vdiff1)%xsi * t_dii
            dzlv(2) = dzlv(2) + Seg_array(i_vdiff1)%xsi * t_dii**2.0e0_DP
            dzlv(3) = dzlv(3) + Seg_array(i_vdiff1)%xsi * t_dii**3.0e0_DP          
        end do

        dzlv(0) = dzlv(0) * PI * drhosv / 6.0e0_DP * NA              
        dzlv(1) = dzlv(1) * PI * drhosv / 6.0e0_DP * NA              
        dzlv(2) = dzlv(2) * PI * drhosv / 6.0e0_DP * NA            
        dzlv(3) = dzlv(3) * PI * drhosv / 6.0e0_DP * NA  

        !zetax
        tdiff1 = 0.0e0_DP

        do i_vdiff1=1, nstypes
        do i_vdiff2=1, nstypes
            tdiff1 = tdiff1 + Seg_array(i_vdiff1)%xsi * Seg_array(i_vdiff2)%xsi * dij(i_vdiff1,i_vdiff2)**3.0e0_DP
        end do
        end do      

        dzetaxv = tdiff1 * PI * drhosv / 6.0e0_DP * NA

        !zetabar
        tdiff1 = 0.0e0_DP

        do i_vdiff1=1, nstypes
        do i_vdiff2=1, nstypes
            tdiff1 = tdiff1 + Seg_array(i_vdiff1)%xsi * Seg_array(i_vdiff2)%xsi * sig(i_vdiff1,i_vdiff2)**3.0e0_DP 
        end do
        end do
        
        !Averaged packing fractions
        dzetabarv = NA * tdiff1 * PI * drhosv / 6.0e0_DP

        !KHS
        diffu  = (1.d0 - zetax)**4.0e0_DP
        diffv  = 1.d0 + 4.0e0_DP * zetax + 4.0e0_DP * zetax**2 - 4.0e0_DP * zetax**3.0e0_DP + zetax**4.0e0_DP
        diffup = -4.0e0_DP * dzetaxv * (1.d0 - zetax)**3.0e0_DP
        diffvp = dzetaxv *( 4.0e0_DP + 8.0e0_DP * zetax - 12.0e0_DP * zetax**2.0e0_DP + 4.0e0_DP * zetax**3.0e0_DP)
         
        dKHSv = (diffv*diffup - diffu*diffvp)/diffv**2.0e0_DP

        return
    end subroutine Vdiffs
!***************************************************************************************************

!*************************************************************************************************** 
end module Diff_mod
!***************************************************************************************************
!***************************************************************************************************
!   SAFT Module to:
!       1. Calculate pressure
!
!       Updated March 2017 to use analytical derivatives
!
!***************************************************************************************************
!
!***************************************************************************************************
module Press_mod
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
    use Diff_mod            ! Differentials setup
!***************************************************************************************************
    implicit none
!***************************************************************************************************
    public  ::  Press
!***************************************************************************************************
    contains
!***************************************************************************************************    
    function Press( ) result(p_out)
        implicit none

        real(kind=DP)           ::  p_out, p_ideal, p_mono, p_chain, p_assoc, p_born, p_ion

! DETERMINE PRESSURE BY CONTRIBUTION      
        call Vdiffs()

        p_ideal =  A_ideal_dv()* KB * NA * T             
        p_mono  = -A_mono_dv() * KB * NA * T     
        p_chain = -A_chain_dv()    
        p_assoc =  A_assoc_dv()      
        p_born  = -A_born_dv() * KB * NA * T     
        p_ion   =  A_ion_dv()

        !hack to avoid p=NaN
        if(p_ideal/=p_ideal) p_ideal=0.0e0_DP
        if(p_mono/=p_mono)   p_mono=0.0e0_DP
        if(p_chain/=p_chain) p_chain=0.0e0_DP
        if(p_assoc/=p_assoc) p_assoc=0.0e0_DP
        if(p_born/=p_born)   p_born=0.0e0_DP
        if(p_ion/=p_ion)     p_ion=0.0e0_DP

        p_out   =  p_ideal + p_mono + p_chain + p_assoc + p_born + p_ion 
       
 99     return
    end function Press
!***************************************************************************************************
!***************************************************************************************************
!   SAFT Module to:
!       1. Calculate chemical potential
!
!       Updated March 2017 to use analytical integrals
!
!***************************************************************************************************
!
end module Press_mod
!***************************************************************************************************
module Mu_mod
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
    use Diff_mod            ! Differentials module  
!***************************************************************************************************
    implicit none
!***************************************************************************************************
    contains
!***************************************************************************************************    
    function Mu( i_comp ) result(mu_out)
        implicit none
        
        integer, intent(in)      ::  i_comp
        
        real(kind=DP)            ::  mu_out
        real(kind=DP)            ::  mu_ideal, mu_mono, mu_chain, &
                                    &   mu_ion, mu_born, mu_assoc   
       
        !ensure nsum is set
        Comp_array(:)%nm = Comp_array(:)%xi
        Comp_array(:)%nm = Comp_array(:)%nm * NA
        nsum = sum(Comp_array(:)%nm)

        !call differential set up
        call Ndiffs( i_comp ) 
        
        mu_ideal = A_ideal_dn( i_comp ) * KB * t * NA
        mu_mono  = A_mono_dn( ) * NA * KB * t       
        mu_assoc = A_assoc_dn( ) * NA * KB * t
        mu_born  = A_born_dn( i_comp ) * KB * t * nsum * NA
        mu_ion   = A_ion_dn( i_comp ) * NA * KB * t
        mu_chain = A_chain_dn( i_comp ) * NA

        mu_out = mu_ideal + mu_mono + mu_assoc + mu_chain + mu_born + mu_ion
!        print*,i_comp,mu_out,mu_ideal,mu_mono,mu_chain,mu_assoc,mu_born,mu_ion
        
        return
    end function Mu
!***************************************************************************************************
    function Mu_neut( i_comp ) result(mu_out)
        implicit none
        
        integer, intent(in)      ::  i_comp
        
        real(kind=DP)            ::  mu_out
        real(kind=DP)            ::  mu_ideal, mu_mono, mu_chain, mu_assoc   
       
        !ensure nsum is set
        Comp_array(:)%nm = Comp_array(:)%xi
        Comp_array(:)%nm = Comp_array(:)%nm * NA
        nsum = sum(Comp_array(:)%nm)

        !call differential set up
        call Ndiffs( i_comp ) 
        
        mu_ideal = A_ideal_dn( i_comp ) * KB * t * NA
        mu_mono  = A_mono_dn( ) * NA * KB * t       
        mu_assoc = A_assoc_dn( ) * NA * KB * t
        mu_chain = A_chain_dn( i_comp ) * NA

        mu_out = mu_ideal + mu_mono + mu_assoc + mu_chain 
!        print*,i_comp,mu_out,mu_ideal,mu_mono,mu_chain,mu_assoc
        
        return
    end function Mu_neut
!***************************************************************************************************
    function Mu_id( i_comp ) result(mu_out)
        implicit none
        
        integer, intent(in)      ::  i_comp
        
        real(kind=DP)            ::  mu_out
        real(kind=DP)            ::  mu_ideal, mu_chain

       
        !ensure nsum is set
        Comp_array(:)%nm = Comp_array(:)%xi
        Comp_array(:)%nm = Comp_array(:)%nm * NA
        nsum = sum(Comp_array(:)%nm)

        !call differential set up
        call Ndiffs( i_comp ) 
        
        mu_ideal = A_ideal_dn( i_comp ) * KB * t * NA
        mu_chain = A_chain_dn( i_comp ) * NA

        mu_out = mu_ideal + mu_chain 
        !print*,i_comp,mu_out,mu_ideal,mu_mono,mu_chain,mu_born,mu_ion,mu_assoc
        
        return
    end function Mu_id
!***************************************************************************************************
    function Mu_res( i_comp ) result(mu_out) !As above, sans ideal
        implicit none
        
        integer, intent(in)      ::  i_comp
        
        real(kind=DP)            ::  mu_out
        real(kind=DP)            ::  mu_mono, mu_chain, &
                                    &   mu_ion, mu_born, mu_assoc   
    
        !ensure nsum is set
        Comp_array(:)%nm = Comp_array(:)%xi
        Comp_array(:)%nm = Comp_array(:)%nm * NA
        nsum = sum(Comp_array(:)%nm)

        !call differential set up
        call Ndiffs( i_comp ) 
        
        mu_mono  = A_mono_dn( ) * NA * KB * t 
        mu_assoc = A_assoc_dn( ) * NA * KB * t
        mu_born  = A_born_dn( i_comp ) * KB * t * nsum * NA
        mu_ion   = A_ion_dn( i_comp ) * NA * KB * t
        mu_chain = A_chain_dn( i_comp ) * NA

        mu_out = mu_mono + mu_assoc + mu_chain + mu_born + mu_ion
!        print*,i_comp,mu_out,mu_mono,mu_chain,mu_born,mu_ion,mu_assoc
        
        return
    end function Mu_res
!***************************************************************************************************  
    function Mu_res_inf( i_comp ) result(mu_out) !As above, sans ideal
        implicit none
        
        integer, intent(in)      ::  i_comp
        
        real(kind=DP)            ::  mu_out
        real(kind=DP)            ::  mu_mono, mu_chain, &
                                    &   mu_born, mu_assoc   
                
        !ensure nsum is set
        Comp_array(:)%nm = Comp_array(:)%xi
        Comp_array(:)%nm = Comp_array(:)%nm * NA
        nsum = sum(Comp_array(:)%nm)

        !call differential set up
        call Ndiffs( i_comp ) 
        
        mu_mono  = A_mono_dn( ) * NA * KB * t 
        mu_assoc = A_assoc_dn( ) * NA * KB * t
        mu_born  = A_born_dn( i_comp ) * KB * t * nsum * NA
        if(mu_born/=mu_born) mu_born=0.0e0_DP        
        mu_chain = A_chain_dn( i_comp ) * NA

        mu_out = mu_mono + mu_assoc + mu_chain + mu_born

        return
    end function Mu_res_inf
!***************************************************************************************************

!***************************************************************************************************
end module Mu_mod
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
