!***************************************************************************************************
!   SAFT Module to:
!       1. Calculate A assoc
!
!       Analtical derivatives added March 2017
!
!***************************************************************************************************
!
!***************************************************************************************************
module Assoc
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

end module Assoc
!***************************************************************************************************
