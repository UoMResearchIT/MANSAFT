!***************************************************************************************************
!   SAFT Module to:
!       1. Calculate common differentials
!           a. with respect to Ni -> chemical potentials
!           b. with respect to V  -> pressure
!***************************************************************************************************
!
!***************************************************************************************************
module Diff
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
end module Diff
!***************************************************************************************************