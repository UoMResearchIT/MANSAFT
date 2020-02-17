
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
    use Types           ! Definitions of types and double precision
    use Global          ! Important global parameters 
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
    use Types           ! Definitions of types and double precision
    use Global          ! Important global parameters
    use Setup           ! To setup the system
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
    use Types           ! Definitions of types and double precision
    use Global          ! Important global parameters
    use Setup           ! To setup the system
    use Ideal           ! Calculate ideal A
    use Mono            ! Calculate mono A
    use Chain           ! Calculate chain A
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
    use Types           ! Definitions of types and double precision
    use Global          ! Important global parameters
    use Setup           ! To setup the system
    use Ideal           ! Calculate ideal A
    use Mono            ! Calculate mono A
    use Chain           ! Calculate chain A
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
    use Types           ! Definitions of types and double precision
    use Global          ! Important global parameters
    use Setup           ! To setup the system
    use Press_mod           ! Calculate pressure
    use Ideal           ! Calculate ideal A
    use Mono            ! Calculate mono A
    use Chain           ! Calculate chain A
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
!
!   Module to add a quote - please add more as and when used
!
!***************************************************************************************************
module Quote_mod
!***************************************************************************************************
contains
!***************************************************************************************************
    subroutine Quote( )
!***************************************************************************************************
!Modules
!=======
        use Zig 
!***************************************************************************************************
    !Variables
    !=========
        implicit none
        
        integer     ::  n_quote, i_quote
        integer     ::  time(1:8), seed
!***************************************************************************************************    
        n_quote = 108  
        time(:)=0        
        
        call DATE_AND_TIME(values=time)     ! Get the current time
        seed = time(1) * (360000*time(5) + 6000*time(6) + 100*time(7) + time(8)) 
  
        call zigset(seed)
        i_quote = anint((n_quote-1)*uni())+1 

        if(i_quote==1) then
            print*, "I am Guybrush Threepwood, mighty pirate. (Monkey Island)"
        else if(i_quote==2) then    
            print*, "You fight like a dairy farmer. (Monkey Island)"
        else if(i_quote==3) then    
            print*, "Look behind you, a three-headed monkey!. (Monkey Island)"
        else if(i_quote==4) then    
            print*, "In my experience, there is no such thing as luck.  "
            print*, " (Obi-Wan Kenobi, Star Wars Episode IV: A New Hope)"
        else if(i_quote==5) then    
            print*, "The needs of the many outweigh the needs of the few. "
            print*, " (Spock, Star Trek II: The Wrath of Khan)"
        else if(i_quote==6) then    
            print*, "You are what you do. A man is defined by his actions, not his memory. "
            print*, " (Kuato, Total Recall)"
        else if(i_quote==7) then    
            print*, "Do what I do. Hold tight and pretend it's a plan! (The Doctor, Doctor Who)"
        else if(i_quote==8) then    
            print*, "Don't panic! (The Hitchhiker's Guide to the Galaxy)"
        else if(i_quote==9) then    
            print*, "Frak! (Battlestar Galactica)"
        else if(i_quote==10) then    
            print*, "Yeah, I got a plan B: making plan A work! (Stingray)"
        else if(i_quote==11) then    
            print*, "Why don't you smegging well smeg off you annoying little smeggy smegging smegger! "
            print*, " (Red Dwarf)"
        else if(i_quote==12) then    
            print*, "Make it so. (Star Trek TNG)"
        else if(i_quote==13) then    
            print*, "Hokey religions and ancient weapons are no match for a good blaster at your side, kid. "
            print*, " (Hans Solo, Star Wars Episode IV: A New Hope)"
        else if(i_quote==14) then    
            print*, "Roads? Where we're going, we don't need roads. (Back to the Future)"
        else if(i_quote==15) then    
            print*, "Come with me if you want to live. (The Terminator)"
        else if(i_quote==16) then    
            print*, "We'd better get back, 'cause it'll be dark soon, and they mostly come at night... mostly. "
            print*, " (Aliens)"
        else if(i_quote==17) then    
            print*, "Dead or alive, you're coming with me! (Robocop)"
        else if(i_quote==18) then    
            print*, "I mean, have you ever actually seen your brain? (Ghost in the Shell)"
        else if(i_quote==19) then    
            print*, "We spend $250 billion a year on defence, and here we are! The fate of the planet "
            print*, " is in the hands of a bunch of retards I wouldn't trust with a potato gun! (Armageddon)"
        else if(i_quote==20) then    
            print*, "Those who know do not speak. Those who speak do not know. (Lao Tzu, Tao Teh Ching)" 
        else if(i_quote==21) then    
            print*, "The truth is not always beautiful, nor beautiful words the truth. (Lao Tzu, Tao Teh Ching)" 
        else if(i_quote==22) then    
            print*, "When I let go of what I am, I become what I might be. (Lao Tzu, Tao Teh Ching)" 
        else if(i_quote==23) then    
            print*, "Care about what other people think and you will always be their prisoner. "
            print*, " (Lao Tzu, Tao Teh Ching)" 
        else if(i_quote==24) then    
            print*, "A man with outward courage dares to die; a man with inner courage dares to live. "
            print*, " (Lao Tzu, Tao Teh Ching)" 
        else if(i_quote==25) then    
            print*, "Teenage angst has paid off well, Now I'm bored and old. (Nirvana)"
        else if(i_quote==26) then    
            print*, "Hey! Wait! I got a new complaint. (Nirvana)"
        else if(i_quote==27) then    
            print*, "I don't care what you think unless it's about me. (Nirvana)"
        else if(i_quote==28) then    
            print*, "Just because you're paranoid. Don't mean they're now after you. (Nirvana)"
        else if(i_quote==29) then    
            print*, "I'm so happy cause today I've found my friends. They're in my head. (Nirvana)"
        else if(i_quote==30) then    
            print*, "Please allow me to introduce myself / I'm a man of wealth and taste/ I've been"
            print*, " around for a long, long year / Stole many a man's soul and faith. (The Rolling Stones)"
        else if(i_quote==31) then    
            print*, "There's no time for us. There's no place for us. What is this thing that builds our "
            print*, " dreams yet slips away from us? (Queen)"
        else if(i_quote==32) then    
            print*, "And when the brokenhearted people living in the world agree, there will be an answer,"
            print*, " let it be. (The Beatles)"
        else if(i_quote==33) then    
            print*, "The fewer the facts, the stronger the opinion. (Arnold H. Glasow)"
        else if(i_quote==34) then    
            print*, "The distance between insanity and genius is measured only by success. (Bruce Feirstein)"
        else if(i_quote==35) then    
            print*, "Your theory is crazy, but it's not crazy enough to be true. (Niels Bohr)"
        else if(i_quote==36) then    
            print*, "Science is nothing but perception. (Plato)"
        else if(i_quote==37) then    
            print*, "If your experiment needs statistics, you ought to have done a better experiment. "
            print*, " (Ernest Rutherford)"
        else if(i_quote==38) then    
            print*, "A physicist is an atom's way of knowing about atoms. (George Wald)"
        else if(i_quote==39) then    
            print*, "Valid criticism does you a favour. (Carl Sagan)"
        else if(i_quote==40) then    
            print*, "Chemistry, unlike other sciences, sprang originally from delusions and superstitions, "
            print*, " and was at its commencement exactly on a par with magic and astrology. (Thomas Thomson)"
        else if(i_quote==41) then    
            print*, "To the electron -- may it never be of any use to anybody. (JJ. Thomson)"
        else if(i_quote==42) then    
            print*, "After climbing a great hill, one only finds there are many more hills to climb. "
            print*, " (Nelson Mandela)"
        else if(i_quote==43) then    
            print*, "If you cannot do great things, do small things in a great way. (Napoleon Hill)"
        else if(i_quote==44) then    
            print*, "It is a rough road that leads to the heights of greatness. (Lucius Annaeus Seneca)"
        else if(i_quote==45) then    
            print*, "Only those who dare to fail greatly can ever achieve greatly. (Robert Kennedy)"
        else if(i_quote==46) then    
            print*, "It is hard to be humble when you're as great as I am. (Muhammad Ali)"
        else if(i_quote==47) then    
            print*, "Two roads diverged in a wood, and I - I took the one less traveled by, and that "
            print*, " has made all the difference. (Robert Frost)"
        else if(i_quote==48) then    
            print*, "Life is what happens to you while you're busy making other plans. (John Lennon)"
        else if(i_quote==49) then    
            print*, "The most common way people give up their power is by thinking they don't have any."
            print*, " (Alice Walker)"
        else if(i_quote==50) then    
            print*, "The best time to plant a tree was 20 years ago. The second best time is now. "
            print*, " (Chinese Proverb)"
        else if(i_quote==51) then    
            print*, "Every child is an artist.  The problem is how to remain an artist once he grows up. "
            print*, " (Pablo Picasso)"
        else if(i_quote==52) then    
            print*, "Whether you think you can or you think you can't, you're right. (Henry Ford)"
        else if(i_quote==53) then    
            print*, "People often say that motivation doesn't last. Well, neither does bathing.  "
            print*, " That's why we recommend it daily. (Zig Ziglar)"
        else if(i_quote==54) then    
            print*, "There is only one way to avoid criticism: do nothing, say nothing, and be nothing. "
            print*, " (Aristotle)"
        else if(i_quote==55) then    
            print*, "Everything you've ever wanted is on the other side of fear. (George Addair)"
        else if(i_quote==56) then    
            print*, "Teach thy tongue to say, I do not know, and thou shalt progress. (Maimonides)"
        else if(i_quote==57) then    
            print*, "How wonderful it is that nobody need wait a single moment before starting to "
            print*, " improve the world. (Anne Frank)"
        else if(i_quote==58) then    
            print*, "If the wind will not serve, take to the oars. (Latin Proverb)"
        else if(i_quote==59) then    
            print*, "Challenges are what make life interesting and overcoming them is what makes life"
            print*, " meaningful. (Joshua J. Marine)"
        else if(i_quote==60) then    
            print*, "The person who says it cannot be done should not interrupt the person who is "
            print*, " doing it. (Chinese Proverb)"
        else if(i_quote==61) then    
            print*, "Build your own dreams, or someone else will hire you to build theirs. (Farrah Gray)"
        else if(i_quote==62) then    
            print*, "It does not matter how slowly you go as long as you do not stop. (Confucius)"
        else if(i_quote==63) then    
            print*, "Everything in moderation, including moderation. (Oscar Wilde)"
        else if(i_quote==64) then    
            print*, "Just taught my kids about taxes by eating 38% of their ice cream. (Conan O'Brien)"
        else if(i_quote==65) then    
            print*, "I dream of a better tomorrow, where chickens can cross the road and not be questioned"
            print*, " about their motives. (Anonymous)"
        else if(i_quote==66) then    
            print*, "It is the mark of an educated mind to be able to entertain a thought without "
            print*, " accepting it. (Aristotle)"
        else if(i_quote==67) then    
            print*, "You see things and you say Why? But I dream things that never were and I say "
            print*, " Why not? (George Bernard Shaw)"
        else if(i_quote==68) then    
            print*, "I have never in my life learned anything from a man who agreed with me. "
            print*, " (Dudley Field Malone)"
        else if(i_quote==69) then    
            print*, "We have all heard that a million monkeys banging on a million typewriters will "
            print*, " eventually reproduce the entire works of Shakespeare. Now, thanks to the internet, "
            print*, " we know this is not true. (Robert Wilensky)"
        else if(i_quote==70) then    
            print*, "Every truth passes through three stages before it is recognised. In the first,"
            print*, " it is ridiculed. In the second, it is opposed. In the third it is regarded as "
            print*, " self-evident. (Arthur Schopenhauer)"
        else if(i_quote==71) then    
            print*, "I always arrive late at the office, but I make up for it by leaving early. "
            print*, " (Charles Lamb)"
        else if(i_quote==72) then    
            print*, "The elevator to success is out of order. You'll have to use the stairs... "
            print*, " one step at a time. (Joe Girard)"
        else if(i_quote==73) then    
            print*, "The greatest way to live with honour in this world is to be what we"
            print*, " pretend to be. (Socrates)"
        else if(i_quote==74) then    
            print*, "Laws control the lesser man... Right conduct controls the greater one. (Mark Twain)"
        else if(i_quote==75) then    
            print*, "Great men are seldom over-scrupulous in the arrangement of their attire. "
            print*, " (Charles Dickens)"
        else if(i_quote==76) then    
            print*, "Time's fun when you're having flies. (Kermit the Frog )" 
        else if(i_quote==77) then    
            print*, "You have to be odd to be number one. (Dr. Seuss)"
        else if(i_quote==78) then    
            print*, "You have brains in your head. You have feet in your shoes. You can steer "
            print*, " yourself any direction you choose. (Dr. Seuss)"
        else if(i_quote==79) then    
            print*, "Sometimes the questions are complicated and the answers are simple. (Dr. Seuss)"
        else if(i_quote==80) then    
            print*, "Today you are you. That is truer than true. There is no one alive who is "
            print*, " youer than you. (Dr. Seuss)"
        else if(i_quote==81) then    
            print*, "A little nonsense now and then, is relished by the wisest men. (Roald Dahl)"
        else if(i_quote==82) then    
            print*, "In ancient times cats were worshipped as gods. They have not forgotten this. "
            print*, "(Terry Pratchett)"
        else if(i_quote==83) then    
            print*, "The trouble with having an open mind, of course, is that people will insist on"
            print*, " coming along and trying to put things in it. (Terry Pratchett)"
        else if(i_quote==84) then    
            print*, "The pen is mightier than the sword if the sword is very short, and the pen is "
            print*, " very sharp. (Terry Pratchett)"
        else if(i_quote==85) then    
            print*, "The presence of those seeking the truth is infinitely to be preferred to the "
            print*, " presence of those who think they've found it. (Terry Pratchett)"
        else if(i_quote==86) then      
            print*, " We have created amazing machines... and used them to destroy our fellow men. "
            print*, " (Metropolis)"
        else if(i_quote==87) then      
            print*, "You're mad. Bonkers. Off your head... but I'll tell you a secret... "
            print*, "all of the best people are. (Alice in Wonderland)"
        else if(i_quote==88) then      
            print*, "Don't cry because it's over, smile because it happened. (Dr. Seuss)"
        else if(i_quote==89) then      
            print*, "We are all in the gutter, but some of us are looking at the stars."
            print*, "(Oscar Wilde)"
        else if(i_quote==90) then      
            print*, "I may not have gone where I intended to go, but I think I have ended"
            print*, "up where I needed to be. (Douglas Adams)"
        else if(i_quote==91) then      
            print*, "For every minute you are angry you lose sixty seconds of happiness."
            print*, "(Ralph Walso Emerson)"
        else if(i_quote==92) then      
            print*, "If you can't explain it to a six year old, you don't understand it "
            print*, "yourself. (Albert Einstein)"
        else if(i_quote==93) then      
            print*, "Success is not final, failure is not fatal: it is the courage to "
            print*, "continue that counts. (Winston Churchill)"
        else if(i_quote==94) then      
            print*, "Reality continues to ruin my life. (Bill Watterson)"
        else if(i_quote==95) then      
            print*, "I am free of all prejudice. I hate everyone equally. (W.C. Fields)"
        else if(i_quote==96) then      
            print*, "It's no use going back to yesterday, because I was a different"
            print*, "person then. (Lewis Carroll)"
        else if(i_quote==97) then      
            print*, "What you're supposed to do when you don't like a thing is change"
            print*, "it. If you can't change it, change the way you think about it. "
            print*, "Don't complain. (Maya Angelou)"
        else if(i_quote==98) then 
            print*, "Experience is that marvelous thing that enables you to recognize "
            print*, "a mistake when you make it again. (F. P. Jones)"
        else if(i_quote==99) then 
            print*, "Human beings, who are almost unique in having the ability to learn "
            print*, "from the experience of others, are also remarkable for their apparent "
            print*, "disinclination to do so. (Douglas Adams)"
        else if(i_quote==100) then
            print*, "My opinions may have changed, but not the fact that I am right. (Ashleigh Brilliant)"
        else if(i_quote==101) then
            print*, "Every man is born as many men and dies as a single one. (Martin Heidegger)"
        else if(i_quote==102) then
            print*, "Many are stubborn in pursuit of the path they have chosen,"
            print*, "few in pursuit of the goal. (Friedrich Nietzsche)"
        else if(i_quote==103) then
            print*, "If you are lonely when you're alone, you are in bad company. (Jean-Paul Sartre)"
        else if(i_quote==104) then
            print*, "There are two things a person should never be angry at, what they can help,"
            print*, "and what they cannot. (Plato)"
        else if(i_quote==105) then
            print*, "No tree, it is said, can grow to heaven unless its roots reach down to hell. (Carl Jung)"
        else if(i_quote==106) then
            print*, "Truth often suffers more by the heat of its defenders than "
            print*, "the arguments of its opposers. (William Penn)"
        else if(i_quote==107) then
            print*, "You can't depend on your eyes when your imagination is out of focus. (Mark Twain)"
        else if(i_quote==108) then
            print*, "Sometimes you face difficulties not because you're doing something wrong,"
            print*, "but because you're doing something right. (Joel Osteen)"
        end if    
        
        return
    end subroutine Quote
end module Quote_mod

!***************************************************************************************************

    Module c05qbfe_mod

!     C05QBF Example Program Module:
!            Parameters and User-defined Routines

!     .. Use Statements ..
      Use nag_library, Only: nag_wp
	  use Types           ! Definitions of types and double precision
	  use Global          ! Important global parameters
	  use Press_mod
	  use Mu_mod
    use Vol_mod
        !     .. Implicit None Statement ..
      Implicit None
!     .. Accessibility Statements ..
      Private
      Public                           :: fcn
!     .. Parameters ..
      Integer, Parameter, Public       :: n_vle = 3, nout_vle = 6
	  Real (Kind=nag_wp), Public			   :: P_L, P_V, Mu_L_1, Mu_L_2, Mu_V_1, Mu_V_2, x_1, x_init_1, x_init_2, x_init_3, x_init_4
    Contains
      Subroutine fcn(n_vle,x,fvec,iuser,ruser,iflag)

!       .. Scalar Arguments ..
        Integer, Intent (Inout)        :: iflag
        Integer, Intent (In)           :: n_vle
!       .. Array Arguments ..
        Real (Kind=nag_wp), Intent (Out) :: fvec(n_vle)
        Real (Kind=nag_wp), Intent (Inout) :: ruser(*)
        Real (Kind=nag_wp), Intent (In) :: x(n_vle)
        Integer, Intent (Inout)        :: iuser(*)
!       .. Executable Statements ..
		v   = x(1)				! Liquid volume is x1
		Comp_array(1)%xi = x_1 !x1
		Comp_array(2)%xi = 1.0_nag_wp - x_1 !x2
		p = Press()
		P_L = p
		Mu_L_1 = Mu(1)
		Mu_L_2 = Mu(2)
		v   = x(2)				! Vapor volume is x2
		Comp_array(1)%xi = x(3) ! y1
		Comp_array(2)%xi =1.0_nag_wp - x(3) !y2
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
subroutine c05qbfe(output_vle)

!     C05QBF Example Main Program

!     .. Use Statements ..
      Use c05qbfe_mod, Only: fcn, n_vle, nout_vle, x_init_1, x_init_2, x_init_3, x_init_4

      Use nag_library, Only: c05qbf, dnrm2, nag_wp, x02ajf
!     .. Implicit None Statement ..
      Implicit None
!     .. Local Scalars ..
      Real (Kind=nag_wp)               :: fnorm, xtol
      Integer                          :: i, ifail
!     .. Local Arrays ..
      Real (Kind=nag_wp), Allocatable  :: fvec(:), x(:)
      Real (Kind=nag_wp)               :: ruser(1)
      Integer                          :: iuser(1)
      !     .. Intrinsic Procedures ..
      Intrinsic                        :: sqrt
	  Real  (Kind=nag_wp)              :: output_vle(3)
!     .. Executable Statements ..

      Allocate (fvec(n_vle),x(n_vle))

!     The following starting values provide a rough solution.

      x(1) = x_init_1
      x(2) = x_init_2
      x(3) = x_init_3
      
      xtol = sqrt(x02ajf())

      ifail = -1
      Call c05qbf(fcn,n_vle,x,fvec,xtol,iuser,ruser,ifail)

      If (ifail==0 .Or. ifail==2 .Or. ifail==3 .Or. ifail==4) Then
        If (ifail==0) Then
!         The NAG name equivalent of dnrm2 is f06ejf
          fnorm = dnrm2(n_vle,fvec,1)
        Else
          Write (nout_vle,*)
          Write (nout_vle,*) 'Approximate solution'
		  x(3)=10
        End If
        !Write (nout,99998)(x(i),i=1,n)
		Do i = 1, n_vle
           output_vle(i) = x(i)
        End Do
      End If

99999 Format (1X,A,E12.4)
99998 Format (1X,3E12.4)
  return
End

!   E04JCF Example Program Text
!   Mark 26.1 Release. NAG Copyright 2017.
    Module e04jcfe_mod

!     E04JCF Example Program Module:
!            Parameters and User-defined Routines

!     .. Use Statements ..
      Use nag_library, Only: nag_wp
	  use Types           ! Definitions of types and double precision
	  use Global          ! Important global parameters
	  use Press_mod
	  use Input
	  use Input_opt   ! Read optimisation parameters
	  use Vol_mod
	  use c05qbfe_mod
	  use Mu_mod
!     .. Implicit None Statement ..
      Implicit None
!     .. Accessibility Statements ..
      Private
      Public                           :: monfun, objfun
!     .. Parameters ..
      Integer, Parameter, Public       :: nout = 6 
      ! integer, public                         ::  n_opt, opt_max
      !real(kind=DP),allocatable, public       ::  min_bound(:), max_bound(:)
      
    Contains
      Subroutine objfun(n,x,f,iuser,ruser,inform)

!       .. Parameters ..
        ! Real (Kind=nag_wp), Parameter  :: five = 5.0_nag_wp
        ! Real (Kind=nag_wp), Parameter  :: ten = 1.0E1_nag_wp
        ! Real (Kind=nag_wp), Parameter  :: two = 2.0_nag_wp
!       .. Scalar Arguments ..
        Real (Kind=nag_wp), Intent (Out) :: f
        Integer, Intent (Out)          :: inform
        Integer, Intent (In)           :: n
!       .. Array Arguments ..
        Real (Kind=nag_wp), Intent (Inout) :: ruser(*)
        Real (Kind=nag_wp), Intent (In) :: x(n)
        Integer, Intent (Inout)        :: iuser(*)
		integer 						:: i
		real(Kind=nag_wp)				::  p_oj, y_oj(2),values(3), mu_1, mu_2, &
										& p_squ, y_squ(2)
!       .. Executable Statements ..
        inform = 0
		
		!allocate(min_bound(1:opt_num),max_bound(1:opt_num))
		do i = 1, opt_num
			if(param_key(i)==1) then               
				!min_bound(i)=min_num(i) * ANG
				!max_bound(i)=max_num(i) * ANG
				sig(param_index(i,1),param_index(i,1)) = x(i)* ANG/1000  !1000 is to make nag get converged
        !print*, Seg_array(param_index(i,1))%sig, x(i)
            else if(param_key(i)==2) then
				!min_bound(i)=min_num(i)
				!max_bound(i)=max_num(i)
				eps(param_index(i,1),param_index(i,1)) = x(i)/10     !10 is to make nag get converged 
            else if(param_key(i)==3) then
				!min_bound(i)=min_num(i)
				!max_bound(i)=max_num(i)
				lr(param_index(i,1),param_index(i,1)) = x(i)/100   !100 is to make nag get converged
            else if(param_key(i)==4) then
				!min_bound(i)=min_num(i)
				!max_bound(i)=max_num(i)
				la(param_index(i,1),param_index(i,1)) = x(i)
            else if(param_key(i)==5) then
				!min_bound(i)=min_num(i)
				!max_bound(i)=max_num(i)
				Seg_array(param_index(i,1))%sf = x(i)/10000 !10000 is to make nag get converged
            else if(param_key(i)==6) then
				!min_bound(i)=min_num(i)
				!max_bound(i)=max_num(i)
				Seg_array(param_index(i,1))%nseg = x(i)
			else if(param_key(i)==7)then
				!min_bound(i)=min_num(i)
				!max_bound(i)=max_num(i)
				eps(param_index(i,1), param_index(i,2))=x(i)/10  !10 is to make nag get converged
				eps(param_index(i,2), param_index(i,1))=x(i)/10 !10 is to make nag get converged				
            else if(param_key(i)==8)then
				!min_bound(i)=min_num(i)
				!max_bound(i)=max_num(i)
				lr(param_index(i,1), param_index(i,2))=x(i)/100 !100 is to make nag get converged
				lr(param_index(i,2), param_index(i,1))=x(i)/100 !100 is to make nag get converged
            else if(param_key(i)==9)then
				!min_bound(i)=min_num(i)
				!max_bound(i)=max_num(i)
				ehb(param_index(i,1),param_index(i,2), &
				& param_index2(i,1),param_index2(i,2)) = x(i)
				ehb(param_index2(i,1),param_index2(i,2), &
				& param_index(i,1),param_index(i,2)) = x(i)
            else if(param_key(i)==10)then
				!min_bound(i)=min_num(i) * ANG**3.0_DP * NA
				!max_bound(i)=max_num(i) * ANG**3.0_DP * NA
				khb(param_index(i,1),param_index(i,2), &
				& param_index2(i,1),param_index2(i,2)) = x(i)* ANG**3.0_DP * NA/10 !10 is to make nag get converged
				khb(param_index2(i,1),param_index2(i,2), &
				& param_index(i,1),param_index(i,2)) = x(i)* ANG**3.0_DP * NA/10 !10 is to make nag get converged
			end if
		end do
		
		p_oj = 0
		y_oj(1)=0
		y_oj(2)=0
		p_squ=0
		y_squ(1)=0
		y_squ(2)=0
		


	if (properties%opt_l(6)) then
		do i = 1, properties%nmu 
          Comp_array(1:nctypes)%xi = properties%ximu(i,1:nctypes)
		  t = properties%t_opt(i)
		  p = properties%p_opt(i)
          x_1 = properties%ximu(i,1)
		  P = 2*P
		  v = Vol_dens_g( )
		  mu_1 = Mu_res(1)
		  mu_2 = Mu_res(2)
		  x_init_1 = v
		  x_init_2 = 8.314*t/(properties%p_opt(i))
		  x_init_3 = exp(mu_1/(8.314*t))*x_1/(exp(mu_1/(8.314*t))*x_1+exp(mu_2/(8.314*t))*(1-x_1))
		  p = properties%p_opt(i)

          call c05qbfe(values )
			p_oj=p_oj+ABS(((P_L+P_V)/2-properties%p_opt(i))/properties%p_opt(i))
			p_squ=p_squ+(((P_L+P_V)/2-properties%p_opt(i))/properties%p_opt(i))**2
			y_oj(1)=y_oj(1)+ABS((values(3)-properties%yimu(i,1))/properties%yimu(i,1))
			y_squ(1)=y_squ(1)+((values(3)-properties%yimu(i,1))/properties%yimu(i,1))**2
			y_oj(2)=y_oj(2)+ABS((1.0_nag_wp-values(3)-properties%yimu(i,2))/properties%yimu(i,2))
			y_squ(2)=y_squ(2)+((1.0_nag_wp-values(3)-properties%yimu(i,2))/properties%yimu(i,2))**2
		end do
			p_oj=p_oj/properties%nmu
			p_squ=p_squ/properties%nmu
			y_oj(1)=y_oj(1)/properties%nmu
			y_squ(1)=y_squ(1)/properties%nmu
			y_oj(2)=y_oj(2)/properties%nmu
			y_squ(2)=y_squ(2)/properties%nmu
		end if		
        
		
		f = p_squ + y_squ(1) + y_squ(2)
   !print*, "f          x(i)"
   write(*, 101)"f=",f, "p_oj=", p_oj, "y_oj1=", y_oj(1), "y_oj2=", y_oj(2),"x(i)=",x(1:opt_num)
101 Format ((2X,A2),(2X,F6.4),(2X,A5),(2X,F6.4),(2X,A6),(2X,F6.4),(2X,A6),(2X,F6.4),(2X,A5),*(2X,F10.4))

        Return

      End Subroutine objfun
      Subroutine monfun(n,nf,x,f,rho,iuser,ruser,inform)

!       .. Scalar Arguments ..
		!use Types       ! Definitions of types and double precision
		!use Global,only: opt_num, min_num, max_num, init_values      ! Important global parameters
		
        Real (Kind=nag_wp), Intent (In) :: f, rho
        Integer, Intent (Out)          :: inform
        Integer, Intent (In)           :: n, nf
!       .. Array Arguments ..
        Real (Kind=nag_wp), Intent (Inout) :: ruser(*)
        Real (Kind=nag_wp), Intent (In) :: x(n)
        Integer, Intent (Inout)        :: iuser(*)
!       .. Local Scalars ..
        Logical                        :: verbose_output
!       .. Executable Statements ..
        inform = 0

        Write (nout,Fmt=99999) 'Monitoring: new trust region radius =', rho

!       Set this to .True. to get more detailed output
        verbose_output = .false.

        If (verbose_output) Then
          Write (nout,Fmt=99998) 'Number of function calls =', nf
          Write (nout,Fmt=99997) 'Current function value =', f
          Write (nout,Fmt=99996) 'The corresponding X is:', x(1:n)
        End If

        Return
99999   Format (/,4X,A,1P,E13.3)
99998   Format (4X,A,I16)
99997   Format (4X,A,1P,E12.4)
99996   Format (4X,A,/,(4X,5E12.4))
      End Subroutine monfun
    End Module e04jcfe_mod

module nag_e04jcfe_mod
contains    
	subroutine e04jcfe(output)

!     Example problem for E04JCF.

!     .. Use Statements ..
      Use e04jcfe_mod, Only: monfun, nout, objfun!, min_bound, max_bound
      Use nag_library, Only: e04jcf, nag_wp, x02alf
	  !use Types       ! Definitions of types and double precision
	  use Global,only: opt_num, min_num, max_num, init_values, param_key, ANG, NA       ! Important global parameters
!     .. Implicit None Statement ..
      Implicit None
!     .. Local Scalars ..
      Real (Kind=nag_wp)               :: f, infbnd, rhobeg, rhoend
      Integer                          :: ifail, maxcal, n, nf, npt
!     .. Local Arrays ..
      Real (Kind=nag_wp), Allocatable  :: bl(:), bu(:), x(:), output(:)
      Real (Kind=nag_wp)               :: ruser(1)!, output(opt_num)
      Integer                          :: iuser(1), i
!     .. Executable Statements ..
      Write (nout,*) 'E04JCF Example Program Results'

      maxcal = 5000
      rhobeg = 100
      rhoend = 1.0E-6_nag_wp
      n = opt_num !4
      npt = ((n+1)*(n+2))/2   !2*n + 1

!     x(3) is unconstrained, so we're going to set bl(3) to a large
!     negative number and bu(3) to a large positive number.

      infbnd = x02alf()**0.25_nag_wp
	!print *, infbnd
      Allocate (bl(n),bu(n),x(n))
		do i=1, opt_num
			bl(i) = min_num(i) !(/1.0_nag_wp,-2.0_nag_wp,-infbnd,1.0_nag_wp/)
			bu(i) = max_num(i) !(/3.0_nag_wp,0.0_nag_wp,infbnd,3.0_nag_wp/)
			x(i) = init_values(i) !(/3.0_nag_wp,-1.0_nag_wp,0.0_nag_wp,1.0_nag_wp/)
			!print *, bl(i),bu(i),x(i)
		end do
		!print*, bl(1:n),bu(1:n),x(1:n),output(1:n)
   print*, "initial values"
   write(*,102) x(1:n)
102 Format (*(2X,F9.4))
      ifail = -1
      Call e04jcf(objfun,n,npt,x,bl,bu,rhobeg,rhoend,monfun,maxcal,f,nf,iuser, &
        ruser,ifail)
	!print*, "test"
      Select Case (ifail)
      Case (0,2:5)
	  !print*, "test"
        If (ifail==0) Then
          Write (nout,Fmt=99999) 'Successful exit from E04JCF.',               &
            'Function value at lowest point found =', f
        Else
          Write (nout,Fmt=99998)                                               &
            'On exit from E04JCF, function value at lowest point found =', f
        End If
		!Write (nout,Fmt=99997) 'The corresponding X is:', x(1:n)
		!print*, "test"
		!Allocate (output(n))	
		do i = 1, opt_num
			
			!print*, "test", x(i)
				output(i) = x(i)
			!print*, "test2", x(i)	
		end do
	!print*, "test"
      End Select

99999 Format (2(/,1X,A),1P,E13.3)
99998 Format (/,1X,A,1P,E13.3)
99997 Format (1X,A,/,(2X,5E13.3))
	return
 End
 end 
 
!***************************************************************************************************
!
!       MANSAFT OPTIMISER
!       =================
!
!       Version 2.7
!***************************************************************************************************
!   SAFT gamma mie program to:
!       1. Read general parameter list form input and
!       2. Calculate system parameters
!       3. Optimise parameters according to simplex method
!***************************************************************************************************
!   Created by SJ Halstead Oct 2015
!       Update 2  Jan 2016
!           - tidied and standardised units
!       Update 3  Feb 2016
!           - added simple P / V calculator
!       Update 4  Mar 2016
!           - standardised MONO units
!       Update 5  Mar 2016
!           - reorganised to focus on propety calculation
!       Update 6  May 2016
!           - property calculations enabled for:
!               - P/V (able to select correct V using lowest G)
!               - Chem pot
!               - Phase equilibrium
!       Update 7  Sept 2016
!           - restructured with FORTRAN best practice
!       Update 7
!           - added electrolyte optimiser
!***************************************************************************************************
!       REFERENCES
!       ==========
!
!       [1] Accurate statistical associating fluid theory for chain molecules formed from Mie segments
!               J Chem Phys 193, 154504 (2013)
!               T. Lafitte, A. Apostolakou, C. Avendano, A. Galindo, C. Adjiman, E. Muller, G. Jackson
!
!               (Main reference)
!
!
!       [2] Prediction of thermodynamic propertied and phase behavior of fluids and mixtures with the
!           SAFT-gamma Mie group-contribution equation of state
!               J. Chem. Eng. Data, 59, 3272-3288 (2014)
!               S. Dufal, V. Papaioannou, M. Sadeqzadeh, T> Pogiatzis, A. Chremos, C.S. Adjiman, 
!               G. Jackson, A. Galindo
!
!               (Associations)
!
!       
!       [3] The A in SAFT: developing the contribution of association to the Helmholtz free energy
!           within a Wertheim TPT1 treatment of generic Mie fluids
!               Mol. Phys. 113, 948-984 (2015)
!               S. Dufal, T. Lafitte, A.J. Haslam, A. Galindo, G.N.I. Clark, C.Vega, G. Jackson            
!
!               (Association coefficients)            
!
!
!***************************************************************************************************

program optimiser
!***************************************************************************************************
!Modules
!=======
    use Types       ! Definitions of types and double precision
    use Input       ! Read the input
    use Input_opt   ! Read optimisation parameters
    use Global      ! Important global parameters
    !use Simplex_mod     ! Simplex optimiser
    use Quote_mod       ! Random quote
    use nag_e04jcfe_mod 
	
!***************************************************************************************************
!Variables
!=========
    implicit none
    
    !Clock
    integer             ::  start, finish, rate
    !Property calculation
    character(len=1)    ::  prop
    !Counter
    integer             ::  i
    !Optimiser
    integer                         ::  n_opt !, opt_max
    real(kind=DP),allocatable       ::  min_opt(:), max_opt(:), paramin(:,:), paramout(:), values(:)!, &
									!&	min_init(:), max_init(:)
    real(kind=DP)                   ::  errout!, values(opt_num)
	! real(kind=DP),allocatable,public::  init_values(:)
!***************************************************************************************************
!Program header
!==============
    print*, " "
    print*, "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    print*, "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    print*, "++                                                                             ++"
    print*, "++                                                                             ++"
    print*, "++            M     M    A    N   N   SSSS    A    FFFFF  TTTTT                ++"
    print*, "++            MM   MM   A A   NN  N  SS      A A   F        T                  ++"
    print*, "++            M M M M  A   A  N N N   SSS   A   A  FFF      T                  ++"
    print*, "++            M  M  M  A A A  N  NN     SS  A A A  F        T                  ++"
    print*, "++            M     M  A   A  N   N  SSSS   A   A  F        T                  ++"
    print*, "++                                                                             ++"
    print*, "++                                                                             ++"
    print*, "++                              MANSAFT                                        ++"
    print*, "++                     Manchester SAFT Gamma Mie                               ++"
    print*, "++                                                                             ++"
    print*, "++                                                                             ++"
    print*, "++                      ****     *****     ********                            ++"
    print*, "++                     ******    ******    ********                            ++"
    print*, "++                    **    **   **  ***      **                               ++"
    print*, "++                    **    **   **  ***      **                               ++"
    print*, "++                    **    **   ******       **                               ++"
    print*, "++                    **    **   ***          **                               ++"
    print*, "++                     ******    ***          **                               ++"
    print*, "++                      ****     ***          **                               ++"
    print*, "++                                                                             ++"
    print*, "++                VERSION 2.7 - The Ant on the Elephant                        ++"
    print*, "++                                2018                                         ++"
    print*, "++                                                                             ++"
    print*, "++                                                                             ++"
    print*, "++                       S.J. Halstead & A.J.Masters                           ++"
    print*, "++                                                                             ++"
    print*, "++         School of Chemical Engineering and Analytical Science               ++"
    print*, "++                     The University of Manchester                            ++"
    print*, "++                                                                             ++"
    print*, "++               contact:  simon.halstead@manchester.ac.uk                     ++"
    print*, "++                         simon.halstead@gmail.com                            ++"
    print*, "++                                                                             ++"
    print*, "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    print*, "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    print*, " "
!***************************************************************************************************
!Start clock
!===========
    call System_clock( start, rate )
    print*,"optimiser"    
!***************************************************************************************************
!Read input
!==========
    call Read_input( )    
     
    if((properties%type/='opt').and.(properties%type/='OPT').and.(properties%type/='Opt')   &
    & .and.(properties%type/='eopt').and.(properties%type/='EOPT')                          &
    & .and.(properties%type/='ilopt').and.(properties%type/='ILOPT')) then
        stop "For optimisation, input run type must be OPT, EOPT or ILOPT"
    end if

    call Read_opt( n_opt, max_opt, min_opt, paramin ) 
		!print*, opt_num, min_num(1:opt_num), max_num(1:opt_num)
!***************************************************************************************************    
		!allocate( init_values(1:opt_num))!, min_init(1:opt_num), max_init(1:opt_num))
		
		!do i = 1, opt_num
			!init_values(i) = (min_num(i) + max_num(i))/2
			!print *, init_values(i)
		!end do
	 !print *, init_values(1:opt_num)
	
	allocate(values(1:opt_num))
	call e04jcfe(values)
	do i = 1, opt_num
			if(param_key(i)==1) then
				values(i) = values(i)/1000  !1000 is to make nag get converged
        !print*, Seg_array(param_index(i,1))%sig, x(i)
            else if(param_key(i)==2) then
				values(i) = values(i)/10     !10 is to make nag get converged 
            else if(param_key(i)==3) then
				values(i) = values(i)/100   !100 is to make nag get converged
            else if(param_key(i)==5) then
				values(i) = values(i)/10000 !10000 is to make nag get converged
			else if(param_key(i)==7)then
				values(i) = values(i)/10  !10 is to make nag get converged				
            else if(param_key(i)==8)then
				values(i) = values(i)/100 !100 is to make nag get converged
            else if(param_key(i)==10)then
				values(i) = values(i)/10 !10 is to make nag get converged
			end if
		end do
	write(*,100)  values(1:opt_num)

100 Format (2X,F9.4,$)
!***************************************************************************************************
!End time
!========
    call System_clock(finish)
   
    print*,  " "
    print*,  "##################################################################"
    print*,  "                    ALL DONE  "
    print*,  "    Time taken: ",REAL(finish-start)/real(rate)," seconds"
    print*,  "##################################################################"
    print*,  " "
    print*,  " "
    
    call Quote( )
    print*, " "
    print*, "#############################################################################################"
    print*, " "
    print*, " "
    print*, " "
    print*, " "
    print*, " "

!***************************************************************************************************            
    stop
!***************************************************************************************************

!***************************************************************************************************
end program optimiser
!
