!***************************************************************************************************
!   SAFT Module to:
!       1. Calculate A mono
!
!       Analytical derivatives added March 2017
!
!***************************************************************************************************
!
!***************************************************************************************************
module Mono
!***************************************************************************************************
!Modules
!=======    
    use Types           ! Definitions of types and double precision
    use Global          ! Important global parameters 
    use Z_eff            ! Zeff calculations
!***************************************************************************************************
    implicit none
!***************************************************************************************************
    public  ::  A_mono, A_mono_dn, A_mono_dv
!***************************************************************************************************    
    contains
!***************************************************************************************************    
    function A_mono() result(a_res)     ![1]  A4   [2] 6 ish
        implicit none
        
        integer         ::  i_a1, i_a2
        real(kind=DP)   ::  a_res
   
        a_res = 0.0e0_DP
        
        do i_a1=1,nctypes
        do i_a2=1,nstypes    
            a_res = a_res + Comp_array(i_a1)%xi * Comp_array(i_a1)%comp(i_a2) * &
            &   Seg_array(i_a2)%nseg * Seg_array(i_a2)%sf
        end do
        end do

        a_res = a_res * am()

        return
    end function A_mono 
!***************************************************************************************************
    !hacked version to separate hs and dispersion contributions
    subroutine A_mono2(a_res, amono_hs, amono_disp)   ![1]  A4   [2] 6 ish
        implicit none
        
        integer                     ::  i_a1, i_a2
        real(kind=DP),intent(out)   ::  a_res, amono_hs, amono_disp
        real(kind=DP)               ::  am_hs, am_disp, am_full
        
   
        a_res = 0.0e0_DP
        
        do i_a1=1,nctypes
        do i_a2=1,nstypes    
            a_res = a_res + Comp_array(i_a1)%xi * Comp_array(i_a1)%comp(i_a2) * &
            &   Seg_array(i_a2)%nseg * Seg_array(i_a2)%sf
        end do
        end do

        am_full     = am(am_hs, am_disp)
        
        amono_hs    = a_res * am_hs
        amono_disp  = a_res * am_disp
        
        a_res = a_res * am_full

        return
    end subroutine A_mono2 
!***************************************************************************************************
    function A_mono_dn() result(a_res)     ![1]  A4   [2] 6 ish
        implicit none
        
        integer         ::  i_a1, i_a2
        real(kind=DP)   ::  a_res, mono_u, mono_up
  
        mono_u   = 0.0e0_DP
        mono_up  = 0.0e0_DP
        
        do i_a1=1,nctypes       
            do i_a2=1,nstypes    
                mono_u = mono_u + Comp_array(i_a1)%xi * Comp_array(i_a1)%comp(i_a2) * &
                &   Seg_array(i_a2)%nseg * Seg_array(i_a2)%sf
                
                mono_up = mono_up + Comp_array(i_a1)%dxin * Comp_array(i_a1)%comp(i_a2) * &
                &   Seg_array(i_a2)%nseg * Seg_array(i_a2)%sf
            end do       
        end do
        
        a_res = mono_u * damn() * nsum + mono_up * am() * nsum + mono_u * am()

        return
    end function A_mono_dn 
!***************************************************************************************************
    !hacked version to separate hs and dispersion contributions
    subroutine A_mono_dn2(a_res, amono_hs, amono_disp) ![1]  A4   [2] 6 ish
        implicit none
        
        integer         ::  i_a1, i_a2
        real(kind=DP)   ::  a_res, amono_hs, amono_disp, mono_u, mono_up
        real(kind=DP)   ::  am_hs, am_disp, am_full
        real(kind=DP)   ::  damn_hs, damn_disp, damn_full
  
        mono_u   = 0.0e0_DP
        mono_up  = 0.0e0_DP
        
        do i_a1=1,nctypes              
            do i_a2=1,nstypes    
                mono_u = mono_u + Comp_array(i_a1)%xi * Comp_array(i_a1)%comp(i_a2) * &
                &   Seg_array(i_a2)%nseg * Seg_array(i_a2)%sf
                
                mono_up = mono_up + Comp_array(i_a1)%dxin * Comp_array(i_a1)%comp(i_a2) * &
                &   Seg_array(i_a2)%nseg * Seg_array(i_a2)%sf
              
            end do       
        end do
      
        am_full   = am(am_hs, am_disp)
        damn_full = damn(damn_hs, damn_disp)
        
        amono_hs    = mono_u * damn_hs * nsum + mono_up * am_hs * nsum + mono_u * am_hs  
        amono_disp  = mono_u * damn_disp * nsum + mono_up * am_disp * nsum + mono_u * am_disp  
  
        a_res = mono_u * damn_full * nsum + mono_up * am_full * nsum + mono_u * am_full

        return
    end subroutine A_mono_dn2 
!***************************************************************************************************
    function A_mono_dv() result(a_res)     ![1]  A4   [2] 6 ish
        implicit none
        
        integer         ::  i_a1, i_a2
        real(kind=DP)   ::  a_res, mono_u, mono_up
  
        mono_u   = 0.0e0_DP
        mono_up  = 0.0e0_DP
        
        do i_a1=1,nctypes
        if(Comp_array(i_a1)%xi/=0.0e0_DP) then        
            do i_a2=1,nstypes    
                mono_u = mono_u + Comp_array(i_a1)%xi * Comp_array(i_a1)%comp(i_a2) * &
                &   Seg_array(i_a2)%nseg * Seg_array(i_a2)%sf
            end do
        end if        
        end do

        a_res = mono_u * damv() 

        return
    end function A_mono_dv
!***************************************************************************************************
    function am(hs_out, disp_out) result(am_out)      !1 A5
        implicit none
 
        real(kind=DP)           ::  am_out
        real(kind=DP),optional  ::  hs_out, disp_out  
        real(kind=DP)           ::  local_hs, local_disp     
      
        local_hs   = ahs()
        local_disp = beta_K * a1() * NA + beta_K**2.0e0_DP * a2()*NA + beta_K**3.0e0_DP * a3() 
        am_out = local_hs + local_disp

        if(present(hs_out)) hs_out=local_hs
        if(present(disp_out)) disp_out=local_disp

        return
    end function am
!***************************************************************************************************
    function damn(hs_out, disp_out) result(am_out)      !1 A5
        implicit none
 
        real(kind=DP)   ::  am_out
        real(kind=DP),optional  ::  hs_out, disp_out  
        real(kind=DP)           ::  local_hs, local_disp  
        
        local_hs   = dahsn()
        local_disp = beta_K * da1n() * NA + beta_K**2.0e0_DP * da2n()*NA + beta_K**3.0e0_DP * da3n()

        am_out = local_hs + local_disp
        
        if(present(hs_out)) hs_out=local_hs
        if(present(disp_out)) disp_out=local_disp

        return
    end function damn
!***************************************************************************************************
    function damv() result(am_out)      !1 A5
        implicit none
 
        real(kind=DP)   ::  am_out

        am_out = dahsv() + beta_K * da1v() * NA + beta_K**2.0e0_DP * da2v()*NA + beta_K**3.0e0_DP * da3v()

        return
    end function damv
!***************************************************************************************************
    function ahs() result(ahs_out)     !1 A6   2 Eq. 7
        implicit none
        
        real(kind=DP)   ::  ahs_out
        
        ahs_out = 6.0e0_DP / PI / (rhos * NA) * (                                           &
        &       (( zl(2)**3.0e0_DP / zl(3)**2.0e0_DP) - zl(0)) * dlog(1.0e0_DP - zl(3)) +   &
        &       3.0e0_DP * zl(1) * zl(2) / (1.0e0_DP - zl(3)) + zl(2)**3.0e0_DP /          &
        &       (zl(3) * (1.0e0_DP - zl(3))**2.0e0_DP)        )  

        return
    end function ahs
!***************************************************************************************************
    function dahsn() result(dahsn_out)     !1 A6   2 Eq. 7
        implicit none
        
        real(kind=DP)       ::  dahsn_out
        real(kind=DP)       ::  ltemp1,ltemp2,ltemp3,ltemp4,ltemp5,ltemp6
        
        ltemp1 = -6.0e0_DP*drhosn/PI/rhos**2.0e0_DP / NA
        ltemp2 = (3.0e0_DP*zl(3)*zl(2)**2.0e0_DP*dzln(2) - 2.0e0_DP*zl(2)**3.0e0_DP*dzln(3))     / &
        &       zl(3)**3.0e0_DP - dzln(0)   
        ltemp3 = -dzln(3) / (1.0e0_DP-zl(3))
        ltemp4 = 3.0e0_DP * ( (1.0e0_DP-zl(3)) * (zl(1)*dzln(2) + dzln(1)*zl(2)) + zl(1)*zl(2)*dzln(3))     / &
        &       (1.0e0_DP - zl(3))**2.0e0_DP
        ltemp6 = zl(3)*(-2.0e0_DP * (1.0e0_DP-zl(3)) * dzln(3)) + dzln(3)*(1.0e0_DP-zl(3))**2.0e0_DP
        ltemp5 = (3.0e0_DP * zl(3) *(1.0e0_DP-zl(3))**2.0e0_DP * zl(2)**2.0e0_DP * dzln(2)   - &
        &       zl(2)**3.0e0_DP * ltemp6) / (zl(3) * (1.0e0_DP - zl(3))**2.0e0_DP )**2.0e0_DP 
        
        dahsn_out = ltemp1 * (                                                                  &
        &       (( zl(2)**3.0e0_DP / zl(3)**2.0e0_DP) - zl(0)) * dlog(1.0e0_DP - zl(3))      +  &
        &       3.0e0_DP * zl(1) * zl(2) / (1.0e0_DP - zl(3)) + zl(2)**3.0e0_DP             /   &
        &       (zl(3) * (1.0e0_DP - zl(3))**2.0e0_DP)        )                             +   &  
        &       6.0e0_DP / PI / (rhos * NA) * (                                                 &
        &       ltemp3 * (( zl(2)**3.0e0_DP / zl(3)**2.0e0_DP) - zl(0))                     +   &
        &       dlog(1.0e0_DP - zl(3)) * ltemp2 + ltemp4 + ltemp5)
      
        return
    end function dahsn
!***************************************************************************************************
    function dahsv() result(dahsv_out)     !1 A6   2 Eq. 7
        implicit none
        
        real(kind=DP)       ::  dahsv_out
        real(kind=DP)       ::  ltemp1,ltemp2,ltemp3,ltemp4,ltemp5,ltemp6
        
        ltemp1 = -6.0e0_DP*drhosv/PI/rhos**2.0e0_DP / NA
        ltemp2 = (3.0e0_DP*zl(3)*zl(2)**2.0e0_DP*dzlv(2) - 2.0e0_DP*zl(2)**3.0e0_DP*dzlv(3))     / &
        &       zl(3)**3.0e0_DP - dzlv(0)   
        ltemp3 = -dzlv(3) / (1.0e0_DP-zl(3))
        ltemp4 = 3.0e0_DP * ( (1.0e0_DP-zl(3)) * (zl(1)*dzlv(2) + dzlv(1)*zl(2)) + zl(1)*zl(2)*dzlv(3))     / &
        &       (1.0e0_DP - zl(3))**2.0e0_DP
        ltemp6 = zl(3)*(-2.0e0_DP * (1.0e0_DP-zl(3)) * dzlv(3)) + dzlv(3)*(1.0e0_DP-zl(3))**2.0e0_DP
        ltemp5 = (3.0e0_DP * zl(3) *(1.0e0_DP-zl(3))**2.0e0_DP * zl(2)**2.0e0_DP * dzlv(2)   - &
        &       zl(2)**3.0e0_DP * ltemp6) / (zl(3) * (1.0e0_DP - zl(3))**2.0e0_DP )**2.0e0_DP 
        
        dahsv_out = ltemp1 * (                                                                  &
        &       (( zl(2)**3.0e0_DP / zl(3)**2.0e0_DP) - zl(0)) * dlog(1.0e0_DP - zl(3))      +  &
        &       3.0e0_DP * zl(1) * zl(2) / (1.0e0_DP - zl(3)) + zl(2)**3.0e0_DP             /   &
        &       (zl(3) * (1.0e0_DP - zl(3))**2.0e0_DP)        )                             +   &  
        &       6.0e0_DP / PI / (rhos * NA) * (                                                 &
        &       ltemp3 * (( zl(2)**3.0e0_DP / zl(3)**2.0e0_DP) - zl(0))                     +   &
        &       dlog(1.0e0_DP - zl(3)) * ltemp2 + ltemp4 + ltemp5)
      
        return
    end function dahsv
!***************************************************************************************************
    function a1() result(a1_out)
        implicit none
        
        real(kind=DP)   ::  a1_out
        integer         ::  i_a1_1, i_a1_2
        
        a1_out=0.0e0_DP
      
        do i_a1_1=1, nstypes
        do i_a1_2=1, nstypes                
            a1_out = a1_out + Seg_array(i_a1_1)%xsi * Seg_array(i_a1_2)%xsi * aij(i_a1_1,i_a1_2)     
        end do
        end do

        return
    end function a1
!***************************************************************************************************
    function da1n() result(a1_out)
        implicit none
        
        real(kind=DP)   ::  a1_out
        integer         ::  i_a1_1, i_a1_2
      
        a1_out=0.0e0_DP
      
        do i_a1_1=1, nstypes
        do i_a1_2=1, nstypes                
            a1_out = a1_out                                                                 + &
            &   Seg_array(i_a1_1)%dxsin * Seg_array(i_a1_2)%xsi * aij(i_a1_1,i_a1_2)        + & 
            &   Seg_array(i_a1_1)%xsi * Seg_array(i_a1_2)%dxsin * aij(i_a1_1,i_a1_2)        + &
            &   Seg_array(i_a1_1)%xsi * Seg_array(i_a1_2)%xsi * daijn(i_a1_1,i_a1_2)           
        end do
        end do

        return
    end function da1n
!***************************************************************************************************
    function da1v() result(a1_out)
        implicit none
        
        real(kind=DP)   ::  a1_out
        integer         ::  i_a1_1, i_a1_2
      
        a1_out=0.0e0_DP
      
        do i_a1_1=1, nstypes
        do i_a1_2=1, nstypes                
            a1_out = a1_out                                                                 + &
            &   Seg_array(i_a1_1)%xsi * Seg_array(i_a1_2)%xsi * daijv(i_a1_1,i_a1_2)           
        end do
        end do

        return
    end function da1v
!***************************************************************************************************
    function a2() result(a2_out)
        implicit none

        real(kind=DP)   ::  a2_out
        integer         ::  i_a2_1, i_a2_2
   
        a2_out=0.0e0_DP
      
        do i_a2_1=1, nstypes
        do i_a2_2=1, nstypes   
            a2_out = a2_out + Seg_array(i_a2_1)%xsi * Seg_array(i_a2_2)%xsi * as2ij(i_a2_1, i_a2_2)
        end do
        end do

        return        
    end function 
!***************************************************************************************************
    function da2n() result(a2_out)
        implicit none

        real(kind=DP)   ::  a2_out
        integer         ::  i_a2_1, i_a2_2
   
        a2_out=0.0e0_DP
      
        do i_a2_1=1, nstypes
        do i_a2_2=1, nstypes   
            a2_out = a2_out                                                                   + &
            &   Seg_array(i_a2_1)%dxsin * Seg_array(i_a2_2)%xsi * as2ij(i_a2_1, i_a2_2)       + &
            &   Seg_array(i_a2_1)%xsi * Seg_array(i_a2_2)%dxsin * as2ij(i_a2_1, i_a2_2)       + &
            &   Seg_array(i_a2_1)%xsi * Seg_array(i_a2_2)%xsi * das2ijn(i_a2_1, i_a2_2)       
        end do
        end do

        return        
    end function da2n
!***************************************************************************************************
    function da2v() result(a2_out)
        implicit none

        real(kind=DP)   ::  a2_out
        integer         ::  i_a2_1, i_a2_2
   
        a2_out=0.0e0_DP
      
        do i_a2_1=1, nstypes
        do i_a2_2=1, nstypes   
            a2_out = a2_out                                                                   + &
            &   Seg_array(i_a2_1)%xsi * Seg_array(i_a2_2)%xsi * das2ijv(i_a2_1, i_a2_2)       
        end do
        end do

        return        
    end function da2v
!***************************************************************************************************
    function a3() result(a3_out)    !A25
        implicit none
        
        real(kind=DP)   ::  a3_out
        integer         ::  i_a3_1, i_a3_2
    
        a3_out=0.0e0_DP
        
        do i_a3_1=1, nstypes
        do i_a3_2=1, nstypes
            a3_out = a3_out + Seg_array(i_a3_1)%xsi * Seg_array(i_a3_2)%xsi * as3ij(i_a3_1, i_a3_2)        
        end do
        end do

        return
    end function
!***************************************************************************************************
    function da3n() result(a3_out)    !A25
        implicit none
        
        real(kind=DP)   ::  a3_out
        integer         ::  i_a3_1, i_a3_2
   
        a3_out=0.0e0_DP
        
        do i_a3_1=1, nstypes
        do i_a3_2=1, nstypes
            a3_out = a3_out                                                                   + &
            &   Seg_array(i_a3_1)%dxsin * Seg_array(i_a3_2)%xsi * as3ij(i_a3_1, i_a3_2)       + &
            &   Seg_array(i_a3_1)%xsi * Seg_array(i_a3_2)%dxsin * as3ij(i_a3_1, i_a3_2)       + &
            &   Seg_array(i_a3_1)%xsi * Seg_array(i_a3_2)%xsi * das3ijn(i_a3_1, i_a3_2)       
        end do
        end do

        return
    end function da3n
!***************************************************************************************************
    function da3v() result(a3_out)    !A25
        implicit none
        
        real(kind=DP)   ::  a3_out
        integer         ::  i_a3_1, i_a3_2
   
        a3_out=0.0e0_DP
        
        do i_a3_1=1, nstypes
        do i_a3_2=1, nstypes
            a3_out = a3_out                                                                 + &
            &   Seg_array(i_a3_1)%xsi * Seg_array(i_a3_2)%xsi * das3ijv(i_a3_1, i_a3_2)       
        end do
        end do

        return
    end function da3v
!***************************************************************************************************
    function aij( i_aij1, i_aij2 ) result(aij_out)  !A11
        implicit none
        
        integer,intent(in)  ::  i_aij1, i_aij2
        real(kind=DP)       ::  xodij, lam1, lam2      
        real(kind=DP)       ::  aij_out
     
        lam1 = la(i_aij1, i_aij2)
        lam2 = lr(i_aij1, i_aij2)   
        xodij = sig(i_aij1, i_aij2) / dij(i_aij1, i_aij2)     !A11 +

        aij_out = cij(i_aij1, i_aij2) * ((xodij**la(i_aij1, i_aij2)) *      &
        &       (as1ij( i_aij1, i_aij2, lam1, zeff(lam1) )+                 &
        &       bij( i_aij1, i_aij2, lam1, xodij )) -                       &
        &       (xodij**lam2) *                                             &
        &       (as1ij( i_aij1, i_aij2, lam2, zeff(lam2) ) +                &
        &       bij( i_aij1, i_aij2, lam2, xodij) )        )  

        return
    end function
!***************************************************************************************************
    function daijn( i_aij1, i_aij2 ) result(aij_out)  !A11
        implicit none
        
        integer,intent(in)  ::  i_aij1, i_aij2
        real(kind=DP)       ::  xodij, lam1, lam2      
        real(kind=DP)       ::  aij_out
            
        lam1 = la(i_aij1, i_aij2)
        lam2 = lr(i_aij1, i_aij2)   
        xodij = sig(i_aij1, i_aij2) / dij(i_aij1, i_aij2)     !A11 +

        aij_out = cij(i_aij1, i_aij2) * ((xodij**la(i_aij1, i_aij2))        * &
        &       (das1ijn( i_aij1, i_aij2, lam1, zeff(lam1), dzeffn(lam1) )  + &
        &       dbijn( i_aij1, i_aij2, lam1, xodij ))                       - &
        &       (xodij**lam2)                                               * &
        &       (das1ijn( i_aij1, i_aij2, lam2, zeff(lam2), dzeffn(lam2) )  + &
        &       dbijn( i_aij1, i_aij2, lam2, xodij) )        )  

        return
    end function daijn
!***************************************************************************************************
    function daijv( i_aij1, i_aij2 ) result(aij_out)  !A11
        implicit none
        
        integer,intent(in)  ::  i_aij1, i_aij2
        real(kind=DP)       ::  xodij, lam1, lam2      
        real(kind=DP)       ::  aij_out
            
        lam1 = la(i_aij1, i_aij2)
        lam2 = lr(i_aij1, i_aij2)   
        xodij = sig(i_aij1, i_aij2) / dij(i_aij1, i_aij2)     !A11 +

        aij_out = cij(i_aij1, i_aij2) * ((xodij**la(i_aij1, i_aij2))        * &
        &       (das1ijv( i_aij1, i_aij2, lam1, zeff(lam1), dzeffv(lam1) )  + &
        &       dbijv( i_aij1, i_aij2, lam1, xodij ))                       - &
        &       (xodij**lam2)                                               * &
        &       (das1ijv( i_aij1, i_aij2, lam2, zeff(lam2), dzeffv(lam2) )  + &
        &       dbijv( i_aij1, i_aij2, lam2, xodij) )        )  

        return
    end function daijv
!***************************************************************************************************
    pure function as1ij( asij1, asij2, lam, zxf ) result(as1ij_out)
        implicit none
        
        integer, intent(in)         ::  asij1,asij2
        real(kind=DP), intent(in)   ::  lam,zxf
        real(kind=DP)               ::  as1ij_out
        
        as1ij_out = -2.0e0_DP * rhos * (PI * eps(asij1, asij2) * (dij(asij1, asij2)**3.0e0_DP)  &
        &   / (lam - 3.0e0_DP)) * (( 1.0e0_DP - zxf / 2.0e0_DP) / (1.0e0_DP - zxf)**3.0e0_DP)

        return
    end function as1ij
!***************************************************************************************************
    function das1ijn( asij1, asij2, lam, zxf, dzxfn ) result(as1ij_out)
        implicit none
        
        integer, intent(in)         ::  asij1,asij2
        real(kind=DP), intent(in)   ::  lam,zxf,dzxfn
        real(kind=DP)               ::  as1ij_out
        real(kind=DP)               ::  l_asu,l_asv,l_asup,l_asvp
        
        l_asu  =  1.0e0_DP - zxf / 2.0e0_DP
        l_asv  =  (1.0e0_DP - zxf)**3.0e0_DP
        l_asup =  -dzxfn / 2.0e0_DP
        l_asvp =  -3.0e0_DP * dzxfn * (1.0e0_DP - zxf)**2.0e0_DP    
        
        as1ij_out = -2.0e0_DP * drhosn * (PI * eps(asij1, asij2) * (dij(asij1, asij2)**3.0e0_DP)    /   &
        &     (lam - 3.0e0_DP)) * (( 1.0e0_DP - zxf / 2.0e0_DP) / (1.0e0_DP - zxf)**3.0e0_DP)       -   &
        &            2.0e0_DP * rhos * (PI * eps(asij1, asij2) * (dij(asij1, asij2)**3.0e0_DP)      /   &
        &     (lam - 3.0e0_DP)) * (l_asv * l_asup - l_asu * l_asvp) / l_asv**2.0e0_DP  

        return
    end function das1ijn
!***************************************************************************************************
    function das1ijv( asij1, asij2, lam, zxf, dzxfv ) result(as1ij_out)
        implicit none
        
        integer, intent(in)         ::  asij1, asij2
        real(kind=DP), intent(in)   ::  lam, zxf, dzxfv
        real(kind=DP)               ::  as1ij_out
        real(kind=DP)               ::  l_asu, l_asv, l_asup, l_asvp
        
        l_asu  =  1.0e0_DP - zxf / 2.0e0_DP
        l_asv  =  (1.0e0_DP - zxf)**3.0e0_DP
        l_asup =  -dzxfv / 2.0e0_DP
        l_asvp =  -3.0e0_DP * dzxfv * (1.0e0_DP - zxf)**2.0e0_DP    
        
        as1ij_out = -2.0e0_DP * drhosv * (PI * eps(asij1, asij2) * (dij(asij1, asij2)**3.0e0_DP)    /   &
        &     (lam - 3.0e0_DP)) * (( 1.0e0_DP - zxf / 2.0e0_DP) / (1.0e0_DP - zxf)**3.0e0_DP)       -   &
        &            2.0e0_DP * rhos * (PI * eps(asij1, asij2) * (dij(asij1, asij2)**3.0e0_DP)      /   &
        &     (lam - 3.0e0_DP)) * (l_asv * l_asup - l_asu * l_asvp) / l_asv**2.0e0_DP  

        return
    end function das1ijv
!***************************************************************************************************
    pure function bij( bij1, bij2, lam, xo ) result(bij_out)  
        implicit none
        
        integer,intent(in)          ::  bij1,bij2 
        real(kind=DP),intent(in)    ::  lam,xo
        real(kind=DP)               ::  ilij,jlij
        real(kind=DP)               ::  bij_out
        
        ilij = -(xo**(3.0e0_DP - lam) - 1.0e0_DP) / (lam - 3.0e0_DP)  
        jlij = -(xo**(4.0e0_DP - lam) * (lam - 3.0e0_DP) - xo**(3.0e0_DP - lam) *   &
        &   (lam - 4.0e0_DP) - 1.0e0_DP) / ((lam - 3.0e0_DP) * (lam - 4.0e0_DP)) 

        bij_out = TWOPI * rhos * dij(bij1, bij2)**3.0e0_DP * eps(bij1, bij2) *      &
        &   ( ((1.0e0_DP - zetax / 2.0e0_DP) / (1.0e0_DP - zetax)**3.0e0_DP) *      &
        &   ilij - (9.0e0_DP * zetax * (1.0e0_DP + zetax) /                         &
        &   (2.0e0_DP * (1.0e0_DP - zetax)**3.0e0_DP)) * jlij  ) 

        return
    end function bij
!***************************************************************************************************
    function dbijn( bij1, bij2, lam, xo ) result(bij_out)
        implicit none
        
        integer,intent(in)          ::  bij1,bij2 
        real(kind=DP),intent(in)    ::  lam,xo
        real(kind=DP)               ::  ilij,jlij
        real(kind=DP)               ::  bij_out
        real(kind=DP)               ::  l_bu,l_bv,l_bup,l_bvp
        
        ilij = -(xo**(3.0e0_DP - lam) - 1.0e0_DP) / (lam - 3.0e0_DP)  
        jlij = -(xo**(4.0e0_DP - lam) * (lam - 3.0e0_DP) - xo**(3.0e0_DP - lam) *   &
        &   (lam - 4.0e0_DP) - 1.0e0_DP) / ((lam - 3.0e0_DP) * (lam - 4.0e0_DP)) 

        l_bu  =  TWOPI * rhos * dij(bij1, bij2)**3.0e0_DP * eps(bij1, bij2)
        l_bup =  TWOPI * drhosn * dij(bij1, bij2)**3.0e0_DP * eps(bij1, bij2)

        l_bv  =  ((1.0e0_DP - zetax / 2.0e0_DP) / (1.0e0_DP - zetax)**3.0e0_DP) *   &
        &   ilij - (9.0e0_DP * zetax * (1.0e0_DP + zetax) /                         &
        &   (2.0e0_DP * (1.0e0_DP - zetax)**3.0e0_DP)) * jlij
        
        l_bvp = ( (-(1.0e0_DP - zetax)**3.0e0_DP * dzetaxn/2.0e0_DP                                     + &
        &       (3.0e0_DP * dzetaxn * (1.0e0_DP - zetax)**2.0e0_DP * (1.0e0_DP - zetax / 2.0e0_DP)))    / &
        &       (1.0e0_DP - zetax)**6.0e0_DP ) * ilij                                                   - &
        &       jlij * ( (2.0e0_DP*(1.0e0_DP-zetax)**3.0e0_DP * (18.0e0_DP*zetax*dzetaxn + 9.0e0_DP     * &
        &       dzetaxn )                                                                               + &
        &       54.0e0_DP*zetax * (1.0e0_DP+zetax) * (1.0e0_DP-zetax)**2.0e0_DP * dzetaxn )             / &
        &       (4.0e0_DP * (1.0e0_DP - zetax)**6.0e0_DP) ) 
   
        bij_out = l_bu*l_bvp + l_bv*l_bup
  
        return
    end function dbijn
!***************************************************************************************************
    function dbijv( bij1, bij2, lam, xo ) result(bij_out)
        implicit none
        
        integer,intent(in)          ::  bij1,bij2 
        real(kind=DP),intent(in)    ::  lam,xo
        real(kind=DP)               ::  ilij,jlij
        real(kind=DP)               ::  bij_out
        real(kind=DP)               ::  l_bu,l_bv,l_bup,l_bvp
        
        ilij = -(xo**(3.0e0_DP - lam) - 1.0e0_DP) / (lam - 3.0e0_DP)  
        jlij = -(xo**(4.0e0_DP - lam) * (lam - 3.0e0_DP) - xo**(3.0e0_DP - lam) *   &
        &   (lam - 4.0e0_DP) - 1.0e0_DP) / ((lam - 3.0e0_DP) * (lam - 4.0e0_DP)) 

        l_bu  =  TWOPI * rhos * dij(bij1, bij2)**3.0e0_DP * eps(bij1, bij2)
        l_bup =  TWOPI * drhosv * dij(bij1, bij2)**3.0e0_DP * eps(bij1, bij2)

        l_bv  =  ((1.0e0_DP - zetax / 2.0e0_DP) / (1.0e0_DP - zetax)**3.0e0_DP) *   &
        &   ilij - (9.0e0_DP * zetax * (1.0e0_DP + zetax) /                         &
        &   (2.0e0_DP * (1.0e0_DP - zetax)**3.0e0_DP)) * jlij
        
        l_bvp = ( (-(1.0e0_DP - zetax)**3.0e0_DP * dzetaxv/2.0e0_DP                                     + &
        &       (3.0e0_DP * dzetaxv * (1.0e0_DP - zetax)**2.0e0_DP * (1.0e0_DP - zetax / 2.0e0_DP)))    / &
        &       (1.0e0_DP - zetax)**6.0e0_DP ) * ilij                                                   - &
        &       jlij * ( (2.0e0_DP*(1.0e0_DP-zetax)**3.0e0_DP * (18.0e0_DP*zetax*dzetaxv + 9.0e0_DP     * &
        &       dzetaxv )                                                                               + &
        &       54.0e0_DP*zetax * (1.0e0_DP+zetax) * (1.0e0_DP-zetax)**2.0e0_DP * dzetaxv )             / &
        &       (4.0e0_DP * (1.0e0_DP - zetax)**6.0e0_DP) ) 
   
        bij_out = l_bu*l_bvp + l_bv*l_bup
  
        return
    end function dbijv
!***************************************************************************************************
    pure function as2ij( as2ij1, as2ij2 ) result(as2ij_out)   
        implicit none
        
        integer,intent(in)      ::  as2ij1, as2ij2
        real(kind=DP)           ::  xodij, chiij
        real(kind=DP)           ::  lam1, lam2, lam3, zx1, zx2, zx3
        real(kind=DP)           ::  as2ij_out

        xodij = sig(as2ij1, as2ij2) / dij(as2ij1, as2ij2)     !A11 +
    
        chiij = fkij(1, as2ij1, as2ij2) * zetabar + fkij(2, as2ij1, as2ij2) * zetabar**5.0e0_DP + &
        &       fkij(3, as2ij1, as2ij2) * zetabar**8.0e0_DP

        lam1  = 2.0e0_DP * la(as2ij1, as2ij2)
        lam2  = la(as2ij1, as2ij2) + lr(as2ij1, as2ij2)
        lam3  = 2.0e0_DP * lr(as2ij1, as2ij2)
        zx1   = zeff(lam1)
        zx2   = zeff(lam2)
        zx3   = zeff(lam3)

        as2ij_out=0.5e0_DP * khs * (1.0e0_DP + chiij) * eps(as2ij1, as2ij2) * cij(as2ij1, as2ij2)**2.0e0_DP     &
        &   * ( xodij**lam1 * (as1ij(as2ij1, as2ij2, lam1, zx1) + bij(as2ij1, as2ij2, lam1, xodij))             &
        &   - 2.0e0_DP * xodij**lam2 * (as1ij(as2ij1, as2ij2, lam2, zx2) +                                      &
        &   bij(as2ij1, as2ij2, lam2, xodij))  + xodij**lam3 * (as1ij(as2ij1, as2ij2, lam3, zx3) +              &
        &   bij(as2ij1, as2ij2, lam3, xodij)))  
 
        return
    end function as2ij
!***************************************************************************************************
    function das2ijn( as2ij1, as2ij2 ) result(as2ij_out)   
        implicit none
        
        integer,intent(in)      ::  as2ij1, as2ij2
        real(kind=DP)           ::  xodij, chiij, dchiijn
        real(kind=DP)           ::  lam1, lam2, lam3, zx1, zx2, zx3, dzx1n, dzx2n, dzx3n
        real(kind=DP)           ::  as2ij_out
        real(kind=DP)           ::  las2_u, las2_v, las2_w, las2_up, las2_vp, las2_wp
  
        xodij = sig(as2ij1, as2ij2) / dij(as2ij1, as2ij2)     !A11 +
    
        chiij = fkij(1, as2ij1, as2ij2) * zetabar + fkij(2, as2ij1, as2ij2) * zetabar**5.0e0_DP + &
        &       fkij(3, as2ij1, as2ij2) * zetabar**8.0e0_DP
    
        dchiijn = dzetabarn * (fkij(1, as2ij1, as2ij2) + 5.0e0_DP * fkij(2, as2ij1, as2ij2) * zetabar**4.0e0_DP + &
        &       8.0e0_DP * fkij(3, as2ij1, as2ij2) * zetabar**7.0e0_DP)

        lam1  = 2.0e0_DP * la(as2ij1, as2ij2)
        lam2  = la(as2ij1, as2ij2) + lr(as2ij1, as2ij2)
        lam3  = 2.0e0_DP * lr(as2ij1, as2ij2)
        zx1   = zeff(lam1)
        zx2   = zeff(lam2)
        zx3   = zeff(lam3)
        
        dzx1n   = dzeffn(lam1)
        dzx2n   = dzeffn(lam2)
        dzx3n   = dzeffn(lam3)

        las2_u = 0.5e0_DP * KHS * eps(as2ij1, as2ij2) * cij(as2ij1, as2ij2)**2.0e0_DP
        las2_v = 1.0e0_DP + chiij
        las2_w = xodij**lam1 * (as1ij(as2ij1, as2ij2, lam1, zx1) + bij(as2ij1, as2ij2, lam1, xodij))            &
        &   - 2.0e0_DP * xodij**lam2 * (as1ij(as2ij1, as2ij2, lam2, zx2) +                                      &
        &   bij(as2ij1, as2ij2, lam2, xodij))  + xodij**lam3 * (as1ij(as2ij1, as2ij2, lam3, zx3) +              &
        &   bij(as2ij1, as2ij2, lam3, xodij))

        las2_up = 0.5e0_DP * dKHSn * eps(as2ij1, as2ij2) * cij(as2ij1, as2ij2)**2.0e0_DP
        las2_vp = dchiijn
        las2_wp = xodij**lam1 * (das1ijn(as2ij1, as2ij2, lam1, zx1, dzx1n) + dbijn(as2ij1, as2ij2, lam1, xodij))   - &
        &   2.0e0_DP * xodij**lam2 * (das1ijn(as2ij1, as2ij2, lam2, zx2, dzx2n)                                    + &
        &   dbijn(as2ij1, as2ij2, lam2, xodij))  + xodij**lam3 * (das1ijn(as2ij1, as2ij2, lam3, zx3, dzx3n)        + &
        &   dbijn(as2ij1, as2ij2, lam3, xodij))

        as2ij_out = las2_up*las2_v*las2_w + las2_u*las2_vp*las2_w + las2_u*las2_v*las2_wp  

        return
    end function das2ijn
!***************************************************************************************************
    function das2ijv( as2ij1, as2ij2 ) result(as2ij_out)   
        implicit none
        
        integer,intent(in)      ::  as2ij1, as2ij2
        real(kind=DP)           ::  xodij, chiij, dchiijv
        real(kind=DP)           ::  lam1, lam2, lam3, zx1, zx2, zx3, dzx1v, dzx2v, dzx3v
        real(kind=DP)           ::  as2ij_out
        real(kind=DP)           ::  las2_u, las2_v, las2_w, las2_up, las2_vp, las2_wp
  
        xodij = sig(as2ij1, as2ij2) / dij(as2ij1, as2ij2)     !A11 +
    
        chiij = fkij(1, as2ij1, as2ij2) * zetabar + fkij(2, as2ij1, as2ij2) * zetabar**5.0e0_DP + &
        &       fkij(3, as2ij1, as2ij2) * zetabar**8.0e0_DP
    
        dchiijv = dzetabarv * (fkij(1, as2ij1, as2ij2) + 5.0e0_DP * fkij(2, as2ij1, as2ij2) * zetabar**4.0e0_DP + &
        &       8.0e0_DP * fkij(3, as2ij1, as2ij2) * zetabar**7.0e0_DP)

        lam1  = 2.0e0_DP * la(as2ij1, as2ij2)
        lam2  = la(as2ij1, as2ij2) + lr(as2ij1, as2ij2)
        lam3  = 2.0e0_DP * lr(as2ij1, as2ij2)
        zx1   = zeff(lam1)
        zx2   = zeff(lam2)
        zx3   = zeff(lam3)
        
        dzx1v   = dzeffv(lam1)
        dzx2v   = dzeffv(lam2)
        dzx3v   = dzeffv(lam3)

        las2_u = 0.5e0_DP * KHS * eps(as2ij1, as2ij2) * cij(as2ij1, as2ij2)**2.0e0_DP
        las2_v = 1.0e0_DP + chiij
        las2_w = xodij**lam1 * (as1ij(as2ij1, as2ij2, lam1, zx1) + bij(as2ij1, as2ij2, lam1, xodij))            &
        &   - 2.0e0_DP * xodij**lam2 * (as1ij(as2ij1, as2ij2, lam2, zx2) +                                      &
        &   bij(as2ij1, as2ij2, lam2, xodij))  + xodij**lam3 * (as1ij(as2ij1, as2ij2, lam3, zx3) +              &
        &   bij(as2ij1, as2ij2, lam3, xodij))

        las2_up = 0.5e0_DP * dKHSv * eps(as2ij1, as2ij2) * cij(as2ij1, as2ij2)**2.0e0_DP
        las2_vp = dchiijv
        las2_wp = xodij**lam1 * (das1ijv(as2ij1, as2ij2, lam1, zx1, dzx1v) + dbijv(as2ij1, as2ij2, lam1, xodij))   - &
        &   2.0e0_DP * xodij**lam2 * (das1ijv(as2ij1, as2ij2, lam2, zx2, dzx2v)                                    + &
        &   dbijv(as2ij1, as2ij2, lam2, xodij))  + xodij**lam3 * (das1ijv(as2ij1, as2ij2, lam3, zx3, dzx3v)        + &
        &   dbijv(as2ij1, as2ij2, lam3, xodij))

        as2ij_out = las2_up*las2_v*las2_w + las2_u*las2_vp*las2_w + las2_u*las2_v*las2_wp  

        return
    end function das2ijv
!***************************************************************************************************
    function as3ij( as3ij1, as3ij2 ) result(as3ij_out)
        implicit none
        
        integer,intent(in)      ::  as3ij1, as3ij2
        real(kind=DP)           ::  as3ij_out
        
        as3ij_out = -eps(as3ij1, as3ij2)**3.0e0_DP * fkij(4, as3ij1, as3ij2) * zetabar * &
            & dexp(fkij( 5, as3ij1, as3ij2 ) * zetabar + fkij(6, as3ij1, as3ij2)* zetabar**2.0e0_DP)

        return
    end function as3ij
!***************************************************************************************************
    function das3ijn( as3ij1, as3ij2 ) result(as3ij_out)
        implicit none
        
        integer,intent(in)      ::  as3ij1, as3ij2
        real(kind=DP)           ::  as3ij_out
        real(kind=DP)           ::  las3_u, las3_v, las3_up, las3_vp
      
        las3_u  = -eps(as3ij1, as3ij2)**3.0e0_DP * fkij(4, as3ij1, as3ij2) * zetabar
        las3_up = -eps(as3ij1, as3ij2)**3.0e0_DP * fkij(4, as3ij1, as3ij2) * dzetabarn
        las3_v  = dexp(fkij( 5, as3ij1, as3ij2 ) * zetabar + fkij(6, as3ij1, as3ij2)* zetabar**2.0e0_DP)
        las3_vp =  dzetabarn * (fkij( 5, as3ij1, as3ij2 ) + 2.0e0_DP * fkij(6, as3ij1, as3ij2)* zetabar) * las3_v
    
        as3ij_out = las3_u * las3_vp + las3_v * las3_up

        return
    end function das3ijn
!***************************************************************************************************
    function das3ijv( as3ij1, as3ij2 ) result(as3ij_out)
        implicit none
        
        integer,intent(in)      ::  as3ij1, as3ij2
        real(kind=DP)           ::  as3ij_out
        real(kind=DP)           ::  las3_u, las3_v, las3_up, las3_vp
      
        las3_u  = -eps(as3ij1, as3ij2)**3.0e0_DP * fkij(4, as3ij1, as3ij2) * zetabar
        las3_up = -eps(as3ij1, as3ij2)**3.0e0_DP * fkij(4, as3ij1, as3ij2) * dzetabarv
        las3_v  = dexp(fkij( 5, as3ij1, as3ij2 ) * zetabar + fkij(6, as3ij1, as3ij2)* zetabar**2.0e0_DP)
        las3_vp =  dzetabarv * (fkij( 5, as3ij1, as3ij2 ) + 2.0e0_DP * fkij(6, as3ij1, as3ij2)* zetabar) * las3_v
    
        as3ij_out = las3_u * las3_vp + las3_v * las3_up

        return
    end function das3ijv
!***************************************************************************************************

!***************************************************************************************************
end module Mono
!***************************************************************************************************
