!***************************************************************************************************
!***************************************************************************************************
!   SAFT Module to:
!       1. Calculate A chain
!***************************************************************************************************
!
!***************************************************************************************************
module Chain
!***************************************************************************************************
!Modules
!=======    
    use Types           ! Definitions of types and double precision
    use Global          ! Important global parameters  
    use Z_eff            ! Zeff calculations
!***************************************************************************************************
    implicit none
!***************************************************************************************************
    public  ::  A_chain, A_chain_dn, A_chain_dv
!***************************************************************************************************
    contains
!*************************************************************************************************** 
!***************************************************************************************************
!***************************************************************************************************
!***************************************************************************************************
!Note to self, if any Chain error remains, it is in lngii, but not g1 or g2...
    function A_chain() result(a_res)     ![1]  A27     [2] 18      
        implicit none
        
        integer             ::  i_a1, i_a2
        real(kind=DP)       ::  tsum
        real(kind=DP)       ::  a_res  

        a_res = 0.0e0_DP
      
        if(chain_switch) then             
            do i_a1=1, nctypes
                tsum = 0.0e0_DP
            do i_a2=1, nstypes
                if(Comp_array(i_a1)%comp(i_a2) > 0.0e0_DP) THEN                        
                    tsum=tsum+Comp_array(i_a1)%comp(i_a2) * Seg_array(i_a2)%nseg * Seg_array(i_a2)%sf 
                end if
            end do
                a_res = a_res - Comp_array(i_a1)%xi * (tsum - 1.0e0_DP) * lngii(i_a1)                               
            end do
        end if

        return
    end function A_chain    
!***************************************************************************************************   
    function A_chain_dn( i_mu_in ) result(a_res)     ![1]  A27     [2] 18      
        implicit none
        
        integer, intent(in) ::  i_mu_in 
        integer             ::  i_a1, i_a2
        real(kind=DP)       ::  tsum
        real(kind=DP)       ::  a_res  

        a_res = 0.0e0_DP
      
        if(chain_switch) then             
            do i_a1=1, nctypes
                if(Comp_array(i_a1)%xi/=0.0e0_DP) then
                    tsum = 0.0e0_DP
                do i_a2=1, nstypes
                    if(Comp_array(i_a1)%comp(i_a2) > 0.0e0_DP) THEN                        
                        tsum=tsum+Comp_array(i_a1)%comp(i_a2) * Seg_array(i_a2)%nseg * Seg_array(i_a2)%sf 
                    end if
                end do
                    a_res = a_res -  &
                    &       nsum * KB * t * Comp_array(i_a1)%dxin * (tsum - 1.0e0_DP) * lngii(i_a1)     - &
                    &       nsum * KB * t * Comp_array(i_a1)%xi * (tsum - 1.0e0_DP) * dlngiin(i_a1)     - &
                    &       KB * t * Comp_array(i_a1)%xi * (tsum - 1.0e0_DP) * lngii(i_a1)    
                end if
            end do        
        end if

        return
    end function A_chain_dn   
!***************************************************************************************************  
    function A_chain_dv() result(a_res)     ![1]  A27     [2] 18      
        implicit none
        
        integer             ::  i_a1, i_a2
        real(kind=DP)       ::  tsum
        real(kind=DP)       ::  a_res  

        a_res = 0.0e0_DP
      
        if(chain_switch) then             
            do i_a1=1, nctypes
                if(Comp_array(i_a1)%xi/=0.0e0_DP) then
                    tsum = 0.0e0_DP
                do i_a2=1, nstypes
                    if(Comp_array(i_a1)%comp(i_a2) > 0.0e0_DP) THEN                        
                        tsum=tsum+Comp_array(i_a1)%comp(i_a2) * Seg_array(i_a2)%nseg * Seg_array(i_a2)%sf 
                    end if
                end do
                    a_res = a_res - nsum * KB * t * Comp_array(i_a1)%xi * (tsum - 1.0e0_DP) * dlngiiv(i_a1)      
                end if    
            end do
        end if

        return
    end function A_chain_dv   
!*************************************************************************************************** 
    function lngii( gri ) result(lngii_out)
        implicit none
        
        integer,intent(in)      ::  gri
        integer                 ::  i_lngii1, i_lngii2
        
        real(kind=DP)           ::  gii, k0, k1, k2, k3
        real(kind=DP)           ::  l_ghs, l_x0, l_d, l_C, l_la, l_lr
        real(kind=DP)           ::  l_eps, l_a1a, l_a1r, l_Ba, l_Br
        real(kind=DP)           ::  l_G1, l_G2
        real(kind=DP)           ::  g2mca, l_gam, l_alpha, l_theta, l_zetaxs
        real(kind=DP)           ::  l_a12r, l_B2r, l_a1ar, l_Bar, l_a12a, l_B2a
        real(kind=DP)           ::  lngii_out
         
        l_d = d3ij(gri)**(1.0e0_DP / 3.0e0_DP)
        l_x0 = sig3ij(gri,gri)**(1.0e0_DP / 3.0e0_DP)/ l_d

        k0 = -dlog(1.0e0_DP - zetax)   +                                                      &
        &   (42.0e0_DP * zetax - 39.0e0_DP * zetax**2.0e0_DP + 9.0e0_DP * zetax**3.0e0_DP  - &
        &   2.0e0_DP * zetax**4.0e0_DP) / (6.0e0_DP * (1.0e0_DP - zetax)**3.0e0_DP)
 
        k1 = (zetax**4.0e0_DP + 6.0e0_DP * zetax**2.0e0_DP - 12.0e0_DP * zetax) / (2.0e0_DP &
        &   * (1.0e0_DP - zetax)**3.0e0_DP)
        k2 = (-3.0e0_DP * zetax**2.0e0_DP) / (8.0e0_DP * (1.0e0_DP - zetax)**2.0e0_DP)
        k3 = (-zetax**4.0e0_DP + 3.0e0_DP * zetax**2.0e0_DP + 3.0e0_DP * zetax) / (6.0e0_DP &
        &   * (1.0e0_DP - zetax)**3.0e0_DP)

        l_ghs = dexp(k0 + k1 * l_x0 + k2 * l_x0**2.0e0_DP + k3 * l_x0**3.0e0_DP)

        l_la = lachij(gri)
        l_lr = lrchij(gri)
    
        l_C = l_lr / (l_lr - l_la) * (l_lr / l_la) ** (l_la / (l_lr - l_la))       
        l_eps = epschij(gri, gri)

        l_Ba  = bchain(l_la, l_x0, l_d, l_eps)
        l_Br  = bchain(l_lr, l_x0, l_d, l_eps)
        l_B2r = bchain(2.0e0_DP * l_lr, l_x0, l_d, l_eps) 
        l_Bar = bchain(l_la + l_lr, l_x0, l_d, l_eps) 
        l_B2a = bchain(2.0e0_DP * l_la, l_x0, l_d, l_eps) 

        l_a1a  = as1chain(l_la, l_d, l_eps)        
        l_a1r  = as1chain(l_lr, l_d, l_eps)
        l_a12r = as1chain(2.0e0_DP * l_lr, l_d, l_eps)
     
        l_a1ar = as1chain(l_la + l_lr, l_d, l_eps)
        l_a12a = as1chain(2.0e0_DP * l_la, l_d, l_eps)

        l_G1 = 1.0e0_DP / (2.0e0_DP * PI * l_eps * l_d**3.0e0_DP) * (3.0e0_DP *     &
        &   as1grad(l_C, l_x0 ,l_d ,l_la ,l_lr ,l_eps) -                            &
        &   l_C * l_la * l_x0**l_la * (l_a1a + l_Ba) / rhos   +                     &
        &   l_C * l_lr * l_x0**l_lr * (l_a1r + l_Br) / rhos )

        l_alpha = l_C * (1.0e0_DP / (l_la - 3.0e0_DP) - 1.0e0_DP / (l_lr - 3.0e0_DP))
        l_theta = dexp(beta_K * l_eps) - 1.0e0_DP

        l_zetaxs = 0.0e0_DP
      
        do i_lngii1=1, nstypes
        do i_lngii2=1, nstypes
            l_zetaxs = l_zetaxs + Seg_array(i_lngii1)%xsi * Seg_array(i_lngii2)%xsi &
            &   * sig(i_lngii1, i_lngii2)**3.0e0_DP
        end do
        end do
        
        l_zetaxs = l_zetaxs * PI * rhos / 6.0e0_DP * NA

        l_gam = phik(7, 0) * (-1.0e0_DP * dtanh(phik(7, 1) * (phik(7, 2) - l_alpha)) + 1.0e0_DP) * &
        &   l_zetaxs * l_theta * dexp(phik(7, 3) * l_zetaxs + phik(7, 4) * l_zetaxs**2.0e0_DP)

        g2mca = 1.0e0_DP / (2.0e0_DP * PI * l_eps**2.00e0_DP * l_d**3.0e0_DP) * ( 3.0e0_DP *            &
        &   as2grad(rhos, l_C, l_x0, l_d, l_la,l_lr,l_eps) -                                            &
        &   l_eps * khs * l_C**2.0e0_DP * l_lr * l_x0**(2.0e0_DP * l_lr) * (l_a12r + l_B2r ) / rhos +   &
        &   l_eps * khs * l_C**2.0e0_DP * (l_la + l_lr) * l_x0**(l_la + l_lr)                       *   &    
        &   (l_a1ar + l_Bar) / rhos - l_eps * khs * l_C**2.0e0_DP * l_la * l_x0**(2.0e0_DP * l_la)  *   &
        &   (l_a12a + l_B2a) / rhos)

        l_G2 = (1.0e0_DP + l_gam) * g2mca

        gii = l_ghs * dexp(beta_K * l_eps * l_G1 / l_ghs + (beta_K * l_eps)**2.0e0_DP * l_G2 / l_ghs)     
        lngii_out = dlog(gii)
!print*,l_ghs,l_G1,l_G2,gii, lngii_out 
        return
    end function lngii
!***************************************************************************************************  
    function dlngiin( gri ) result(lngii_out)
        implicit none
        
        integer,intent(in)      ::  gri
        integer                 ::  i_lngii1, i_lngii2
        
        real(kind=DP)           ::  gii, k0, k1, k2, k3
        real(kind=DP)           ::  l_ghs, l_x0, l_d, l_C, l_la, l_lr
        real(kind=DP)           ::  l_eps, l_a1a, l_a1r, l_Ba, l_Br
        real(kind=DP)           ::  l_G1, l_G2
        real(kind=DP)           ::  g2mca, l_gam, l_alpha, l_theta, l_zetaxs
        real(kind=DP)           ::  l_a12r, l_B2r, l_a1ar, l_Bar, l_a12a, l_B2a
        real(kind=DP)           ::  lngii_out
        !differential parts
        real(kind=DP)           ::  dk0n, dk1n, dk2n, dk3n, dl_ghsn
        real(kind=DP)           ::  dl_a1an, dl_a1rn, dl_Ban, dl_Brn
        real(kind=DP)           ::  dl_a12rn, dl_B2rn, dl_a1arn, dl_Barn, dl_a12an, dl_B2an
        real(kind=DP)           ::  dl_G1n, dl_gamn, dl_G2n, dgiin, dl_zetaxsn, dg2mcan
             
        l_d = d3ij(gri)**(1.0e0_DP / 3.0e0_DP)

        l_x0 = sig3ij(gri,gri)**(1.0e0_DP / 3.0e0_DP)/ l_d

        k0 = -dlog(1.0e0_DP - zetax)   +                                                      &
        &   (42.0e0_DP * zetax - 39.0e0_DP * zetax**2.0e0_DP + 9.0e0_DP * zetax**3.0e0_DP  - &
        &   2.0e0_DP * zetax**4.0e0_DP) / (6.0e0_DP * (1.0e0_DP - zetax)**3.0e0_DP)
        
        k1 = (zetax**4.0e0_DP + 6.0e0_DP * zetax**2.0e0_DP - 12.0e0_DP * zetax) / (2.0e0_DP &
        &   * (1.0e0_DP - zetax)**3.0e0_DP)
        
        k2 = (-3.0e0_DP * zetax**2.0e0_DP) / (8.0e0_DP * (1.0e0_DP - zetax)**2.0e0_DP)
        
        k3 = (-zetax**4.0e0_DP + 3.0e0_DP * zetax**2.0e0_DP + 3.0e0_DP * zetax) / (6.0e0_DP &
        &   * (1.0e0_DP - zetax)**3.0e0_DP)
        
        dk0n = dzetaxn/(1.0e0_DP - zetax)                                                           + &
        &   (6.0e0_DP * (1.0e0_DP - zetax)**3.0e0_DP * (42.0e0_DP * dzetaxn - 78.0e0_DP * zetax     * &
        &   dzetaxn + 27.0e0_DP * zetax**2.0e0_DP * dzetaxn - 8.0e0_DP * zetax**3.0e0_DP * dzetaxn) + &
        &   18.0e0_DP * dzetaxn * (1.0e0_DP - zetax)**2.0e0_DP                                      * &
        &   (42.0e0_DP * zetax - 39.0e0_DP * zetax**2.0e0_DP + 9.0e0_DP * zetax**3.0e0_DP           - &
        &   2.0e0_DP * zetax**4.0e0_DP))                                                            / &
        &   (36.0e0_DP * (1.0e0_DP - zetax)**6.0e0_DP)  

        dk1n = ( 2.0e0_DP  * (1.0e0_DP - zetax)**3.0e0_DP                                           * &
        &   (4.0e0_DP*dzetaxn*zetax**3.0e0_DP + 12.0e0_DP * dzetaxn*zetax - 12.0e0_DP * dzetaxn)    + &
        &   6.0e0_DP * dzetaxn * (1.0e0_DP - zetax)**2.0e0_DP                                       * &
        &   (zetax**4.0e0_DP + 6.0e0_DP * zetax**2.0e0_DP - 12.0e0_DP * zetax))                     / &
        &   (4.0e0_DP  * (1.0e0_DP - zetax)**6.0e0_DP)

        dk2n = (-6.0e0_DP * dzetaxn * zetax * 8.0e0_DP * (1.0e0_DP - zetax)**2.0e0_DP               - &
        &   3.0e0_DP * zetax**2.0e0_DP * 16.0e0_DP * dzetaxn * (1.0e0_DP - zetax))                  / &
        &   (64.0e0_DP * (1.0e0_DP - zetax)**4.0e0_DP)    

        dk3n = ((6.0e0_DP * (1.0e0_DP - zetax)**3.0e0_DP)                                           * &
        &   (-4.0e0_DP*dzetaxn*zetax**3.0e0_DP + 6.0e0_DP * dzetaxn * zetax + 3.0e0_DP * dzetaxn)   - &
        &   (-zetax**4.0e0_DP + 3.0e0_DP * zetax**2.0e0_DP + 3.0e0_DP * zetax)                      * &
        &   (-18.0e0_DP * dzetaxn * (1.0e0_DP - zetax)**2.0e0_DP))                                  / &
        &   (36.0e0_DP * (1.0e0_DP - zetax)**6.0e0_DP)
       
        l_ghs = dexp(k0 + k1 * l_x0 + k2 * l_x0**2.0e0_DP + k3 * l_x0**3.0e0_DP)

        dl_ghsn = (dk0n + dk1n * l_x0 + dk2n * l_x0**2.0e0_DP + dk3n * l_x0**3.0e0_DP)              * &
        &   dexp(k0 + k1 * l_x0 + k2 * l_x0**2.0e0_DP + k3 * l_x0**3.0e0_DP)

        l_la = lachij(gri)
        l_lr = lrchij(gri)
       
        l_C = l_lr / (l_lr - l_la) * (l_lr / l_la) ** (l_la / (l_lr - l_la))       
        l_eps = epschij(gri, gri)

        l_Ba  = bchain(l_la, l_x0, l_d, l_eps)
        l_Br  = bchain(l_lr, l_x0, l_d, l_eps)
        l_B2r = bchain(2.0e0_DP * l_lr, l_x0, l_d, l_eps) 
        l_Bar = bchain(l_la + l_lr, l_x0, l_d, l_eps) 
        l_B2a = bchain(2.0e0_DP * l_la, l_x0, l_d, l_eps)

        dl_Ban  = dbchainn(l_la, l_x0, l_d, l_eps)
        dl_Brn  = dbchainn(l_lr, l_x0, l_d, l_eps)
        dl_B2rn = dbchainn(2.0e0_DP * l_lr, l_x0, l_d, l_eps) 
        dl_Barn = dbchainn(l_la + l_lr, l_x0, l_d, l_eps) 
        dl_B2an = dbchainn(2.0e0_DP * l_la, l_x0, l_d, l_eps) 

        l_a1a  = as1chain(l_la, l_d, l_eps)        
        l_a1r  = as1chain(l_lr, l_d, l_eps)
        l_a12r = as1chain(2.0e0_DP * l_lr, l_d, l_eps)
        l_a1ar = as1chain(l_la + l_lr, l_d, l_eps)
        l_a12a = as1chain(2.0e0_DP * l_la, l_d, l_eps)
        
        dl_a1an  = das1chainn(l_la, l_d, l_eps)        
        dl_a1rn  = das1chainn(l_lr, l_d, l_eps)
        dl_a12rn = das1chainn(2.0e0_DP * l_lr, l_d, l_eps)
        dl_a1arn = das1chainn(l_la + l_lr, l_d, l_eps)
        dl_a12an = das1chainn(2.0e0_DP * l_la, l_d, l_eps)

        l_G1 = 1.0e0_DP / (2.0e0_DP * PI * l_eps * l_d**3.0e0_DP) * (3.0e0_DP *     &
        &   as1grad(l_C, l_x0 ,l_d ,l_la ,l_lr ,l_eps) -                            &
        &   l_C * l_la * l_x0**l_la * (l_a1a + l_Ba) / rhos   +                     &
        &   l_C * l_lr * l_x0**l_lr * (l_a1r + l_Br) / rhos )

        dl_G1n = 1.0e0_DP / (2.0e0_DP * PI * l_eps * l_d**3.0e0_DP) * (3.0e0_DP     * &
        &   das1gradn(l_C, l_x0 ,l_d ,l_la ,l_lr ,l_eps)                            - &
        &   l_C * l_la * l_x0**l_la                                                 * & 
        &   (rhos*(dl_a1an + dl_Ban) - drhosn*(l_a1a + l_Ba))/rhos**2.0e0_DP        + &
        &   l_C * l_lr * l_x0**l_lr                                                 * &
        &   (rhos*(dl_a1rn + dl_Brn) - drhosn*(l_a1r + l_Br))/rhos**2.0e0_DP) 

        l_alpha = l_C * (1.0e0_DP / (l_la - 3.0e0_DP) - 1.0e0_DP / (l_lr - 3.0e0_DP))
        l_theta = dexp(beta_K * l_eps) - 1.0e0_DP
        
        l_zetaxs    = 0.0e0_DP
        dl_zetaxsn  = 0.0e0_DP

        do i_lngii1=1, nstypes
        do i_lngii2=1, nstypes
            l_zetaxs = l_zetaxs + Seg_array(i_lngii1)%xsi * Seg_array(i_lngii2)%xsi &
            &   * sig(i_lngii1, i_lngii2)**3.0e0_DP
            
            dl_zetaxsn = dl_zetaxsn + (Seg_array(i_lngii1)%dxsin * Seg_array(i_lngii2)%xsi      + &
            &   Seg_array(i_lngii1)%xsi * Seg_array(i_lngii2)%dxsin)                            * &
            &   sig(i_lngii1, i_lngii2)**3.0e0_DP
        end do
        end do
               
        dl_zetaxsn = PI / 6.0e0_DP * NA * (l_zetaxs * drhosn + dl_zetaxsn * rhos)
        l_zetaxs   = l_zetaxs * PI * rhos / 6.0e0_DP * NA
    
        l_gam = phik(7, 0) * (-1.0e0_DP * dtanh(phik(7, 1) * (phik(7, 2) - l_alpha)) + 1.0e0_DP) * &
        &   l_zetaxs * l_theta * dexp(phik(7, 3) * l_zetaxs + phik(7, 4) * l_zetaxs**2.0e0_DP)
        
        dl_gamn = phik(7, 0) * (-1.0e0_DP * dtanh(phik(7, 1) * (phik(7, 2) - l_alpha)) + 1.0e0_DP)   * &
        &   l_theta * dexp(phik(7, 3) * l_zetaxs + phik(7, 4) * l_zetaxs**2.0e0_DP)                  * &
        &   (dl_zetaxsn + l_zetaxs * (phik(7, 3) * dl_zetaxsn + 2.0e0_DP * phik(7, 4) * dl_zetaxsn  * &
        &   l_zetaxs))

        g2mca = 1.0e0_DP / (2.0e0_DP * PI * l_eps**2.00e0_DP * l_d**3.0e0_DP) * ( 3.0e0_DP *            &
        &   as2grad(rhos, l_C, l_x0, l_d, l_la,l_lr,l_eps) -                                            &
        &   l_eps * khs * l_C**2.0e0_DP * l_lr * l_x0**(2.0e0_DP * l_lr) * (l_a12r + l_B2r ) / rhos +   &
        &   l_eps * khs * l_C**2.0e0_DP * (l_la + l_lr) * l_x0**(l_la + l_lr)                       *   &    
        &   (l_a1ar + l_Bar) / rhos - l_eps * khs * l_C**2.0e0_DP * l_la * l_x0**(2.0e0_DP * l_la)  *   &
        &   (l_a12a + l_B2a) / rhos)
        
        dg2mcan = 1.0e0_DP / (2.0e0_DP * PI * l_eps**2.00e0_DP * l_d**3.0e0_DP)                         * &
        &   ( 3.0e0_DP * das2gradn(rhos, l_C, l_x0, l_d, l_la,l_lr,l_eps)                               - &

        &   l_eps * l_C**2.0e0_DP * l_lr * l_x0**(2.0e0_DP * l_lr)                              * ( &
        &   khs * (rhos * (dl_a12rn + dl_B2rn ) - drhosn * (l_a12r + l_B2r )) / rhos**2.0e0_DP  +   &
        &   dkhsn * (l_a12r + l_B2r ) / rhos)                                                   +   &  

        &   l_eps * l_C**2.0e0_DP * (l_la + l_lr) * l_x0**(l_la + l_lr)                         * ( &
        &   khs * (rhos * (dl_a1arn + dl_Barn ) - drhosn * (l_a1ar + l_Bar )) / rhos**2.0e0_DP  +   &
        &   dkhsn * (l_a1ar + l_Bar ) / rhos)                                                   -   &  
        
        &   l_eps * l_C**2.0e0_DP * l_la * l_x0**(2.0e0_DP * l_la)                              * ( &
        &   khs * (rhos * (dl_a12an + dl_B2an ) - drhosn * (l_a12a + l_B2a )) / rhos**2.0e0_DP  +   &
        &   dkhsn * (l_a12a + l_B2a ) / rhos)) 
      
        l_G2 = (1.0e0_DP + l_gam) * g2mca
        
        dl_G2n = (1.0e0_DP + l_gam) * dg2mcan + dl_gamn * g2mca
     
        gii = l_ghs * dexp(beta_K * l_eps * l_G1 / l_ghs + (beta_K * l_eps)**2.0e0_DP * l_G2 / l_ghs)     
        
        dgiin = dexp(beta_K * l_eps * l_G1 / l_ghs + (beta_K * l_eps)**2.0e0_DP * l_G2 / l_ghs)      * &
        &   (dl_ghsn + l_ghs                                                                        * &
        &   (beta_K * l_eps * (dl_G1n*l_ghs - l_G1*dl_ghsn)/l_ghs**2.0e0_DP                         + &
        &   (beta_K * l_eps)**2.0e0_DP * (dl_G2n*l_ghs - l_G2*dl_ghsn)/l_ghs**2.0e0_DP))      
   
        lngii_out = dgiin / gii

        return
    end function dlngiin
!*************************************************************************************************** 
    function dlngiiv( gri ) result(lngii_out)
        implicit none
        
        integer,intent(in)      ::  gri
        integer                 ::  i_lngii1, i_lngii2
        
        real(kind=DP)           ::  gii, k0, k1, k2, k3
        real(kind=DP)           ::  l_ghs, l_x0, l_d, l_C, l_la, l_lr
        real(kind=DP)           ::  l_eps, l_a1a, l_a1r, l_Ba, l_Br
        real(kind=DP)           ::  l_G1, l_G2
        real(kind=DP)           ::  g2mca, l_gam, l_alpha, l_theta, l_zetaxs
        real(kind=DP)           ::  l_a12r, l_B2r, l_a1ar, l_Bar, l_a12a, l_B2a
        real(kind=DP)           ::  lngii_out
        !differential parts
        real(kind=DP)           ::  dk0v, dk1v, dk2v, dk3v, dl_ghsv
        real(kind=DP)           ::  dl_a1av, dl_a1rv, dl_Bav, dl_Brv
        real(kind=DP)           ::  dl_a12rv, dl_B2rv, dl_a1arv, dl_Barv, dl_a12av, dl_B2av
        real(kind=DP)           ::  dl_G1v, dl_gamv, dl_G2v, dgiiv, dl_zetaxsv, dg2mcav
           
        l_d = d3ij(gri)**(1.0e0_DP / 3.0e0_DP)
        l_x0 = sig3ij(gri,gri)**(1.0e0_DP / 3.0e0_DP)/ l_d

        k0 = -dlog(1.0e0_DP - zetax)   +                                                      &
        &   (42.0e0_DP * zetax - 39.0e0_DP * zetax**2.0e0_DP + 9.0e0_DP * zetax**3.0e0_DP  - &
        &   2.0e0_DP * zetax**4.0e0_DP) / (6.0e0_DP * (1.0e0_DP - zetax)**3.0e0_DP)
        
        k1 = (zetax**4.0e0_DP + 6.0e0_DP * zetax**2.0e0_DP - 12.0e0_DP * zetax) / (2.0e0_DP &
        &   * (1.0e0_DP - zetax)**3.0e0_DP)
        
        k2 = (-3.0e0_DP * zetax**2.0e0_DP) / (8.0e0_DP * (1.0e0_DP - zetax)**2.0e0_DP)
        
        k3 = (-zetax**4.0e0_DP + 3.0e0_DP * zetax**2.0e0_DP + 3.0e0_DP * zetax) / (6.0e0_DP &
        &   * (1.0e0_DP - zetax)**3.0e0_DP)
        
        dk0v = dzetaxv/(1.0e0_DP - zetax)                                                           + &
        &   (6.0e0_DP * (1.0e0_DP - zetax)**3.0e0_DP * (42.0e0_DP * dzetaxv - 78.0e0_DP * zetax     * &
        &   dzetaxv + 27.0e0_DP * zetax**2.0e0_DP * dzetaxv - 8.0e0_DP * zetax**3.0e0_DP * dzetaxv) + &
        &   18.0e0_DP * dzetaxv * (1.0e0_DP - zetax)**2.0e0_DP                                      * &
        &   (42.0e0_DP * zetax - 39.0e0_DP * zetax**2.0e0_DP + 9.0e0_DP * zetax**3.0e0_DP           - &
        &   2.0e0_DP * zetax**4.0e0_DP))                                                            / &
        &   (36.0e0_DP * (1.0e0_DP - zetax)**6.0e0_DP)  

        dk1v = ( 2.0e0_DP  * (1.0e0_DP - zetax)**3.0e0_DP                                           * &
        &   (4.0e0_DP*dzetaxv*zetax**3.0e0_DP + 12.0e0_DP * dzetaxv*zetax - 12.0e0_DP * dzetaxv)    + &
        &   6.0e0_DP * dzetaxv * (1.0e0_DP - zetax)**2.0e0_DP                                       * &
        &   (zetax**4.0e0_DP + 6.0e0_DP * zetax**2.0e0_DP - 12.0e0_DP * zetax))                     / &
        &   (4.0e0_DP  * (1.0e0_DP - zetax)**6.0e0_DP)

        dk2v = (-6.0e0_DP * dzetaxv * zetax * 8.0e0_DP * (1.0e0_DP - zetax)**2.0e0_DP               - &
        &   3.0e0_DP * zetax**2.0e0_DP * 16.0e0_DP * dzetaxv * (1.0e0_DP - zetax))                  / &
        &   (64.0e0_DP * (1.0e0_DP - zetax)**4.0e0_DP)    

        dk3v = ((6.0e0_DP * (1.0e0_DP - zetax)**3.0e0_DP)                                           * &
        &   (-4.0e0_DP*dzetaxv*zetax**3.0e0_DP + 6.0e0_DP * dzetaxv * zetax + 3.0e0_DP * dzetaxv)   - &
        &   (-zetax**4.0e0_DP + 3.0e0_DP * zetax**2.0e0_DP + 3.0e0_DP * zetax)                      * &
        &   (-18.0e0_DP * dzetaxv * (1.0e0_DP - zetax)**2.0e0_DP))                                  / &
        &   (36.0e0_DP * (1.0e0_DP - zetax)**6.0e0_DP)
   
        l_ghs = dexp(k0 + k1 * l_x0 + k2 * l_x0**2.0e0_DP + k3 * l_x0**3.0e0_DP)

        dl_ghsv = (dk0v + dk1v * l_x0 + dk2v * l_x0**2.0e0_DP + dk3v * l_x0**3.0e0_DP)              * &
        &   dexp(k0 + k1 * l_x0 + k2 * l_x0**2.0e0_DP + k3 * l_x0**3.0e0_DP)

        l_la = lachij(gri)
        l_lr = lrchij(gri)
       
        l_C = l_lr / (l_lr - l_la) * (l_lr / l_la) ** (l_la / (l_lr - l_la))       
        l_eps = epschij(gri, gri)

        l_Ba  = bchain(l_la, l_x0, l_d, l_eps)
        l_Br  = bchain(l_lr, l_x0, l_d, l_eps)
        l_B2r = bchain(2.0e0_DP * l_lr, l_x0, l_d, l_eps) 
        l_Bar = bchain(l_la + l_lr, l_x0, l_d, l_eps) 
        l_B2a = bchain(2.0e0_DP * l_la, l_x0, l_d, l_eps)

        dl_Bav  = dbchainv(l_la, l_x0, l_d, l_eps)
        dl_Brv  = dbchainv(l_lr, l_x0, l_d, l_eps)
        dl_B2rv = dbchainv(2.0e0_DP * l_lr, l_x0, l_d, l_eps) 
        dl_Barv = dbchainv(l_la + l_lr, l_x0, l_d, l_eps) 
        dl_B2av = dbchainv(2.0e0_DP * l_la, l_x0, l_d, l_eps) 

        l_a1a  = as1chain(l_la, l_d, l_eps)        
        l_a1r  = as1chain(l_lr, l_d, l_eps)
        l_a12r = as1chain(2.0e0_DP * l_lr, l_d, l_eps)
        l_a1ar = as1chain(l_la + l_lr, l_d, l_eps)
        l_a12a = as1chain(2.0e0_DP * l_la, l_d, l_eps)
        
        dl_a1av  = das1chainv(l_la, l_d, l_eps)        
        dl_a1rv  = das1chainv(l_lr, l_d, l_eps)
        dl_a12rv = das1chainv(2.0e0_DP * l_lr, l_d, l_eps)
        dl_a1arv = das1chainv(l_la + l_lr, l_d, l_eps)
        dl_a12av = das1chainv(2.0e0_DP * l_la, l_d, l_eps)

        l_G1 = 1.0e0_DP / (2.0e0_DP * PI * l_eps * l_d**3.0e0_DP) * (3.0e0_DP *     &
        &   as1grad(l_C, l_x0 ,l_d ,l_la ,l_lr ,l_eps) -                            &
        &   l_C * l_la * l_x0**l_la * (l_a1a + l_Ba) / rhos   +                     &
        &   l_C * l_lr * l_x0**l_lr * (l_a1r + l_Br) / rhos )

        dl_G1v = 1.0e0_DP / (2.0e0_DP * PI * l_eps * l_d**3.0e0_DP) * (3.0e0_DP     * &
        &   das1gradv(l_C, l_x0 ,l_d ,l_la ,l_lr ,l_eps)                            - &
        &   l_C * l_la * l_x0**l_la                                                 * & 
        &   (rhos*(dl_a1av + dl_Bav) - drhosv*(l_a1a + l_Ba))/rhos**2.0e0_DP        + &
        &   l_C * l_lr * l_x0**l_lr                                                 * &
        &   (rhos*(dl_a1rv + dl_Brv) - drhosv*(l_a1r + l_Br))/rhos**2.0e0_DP) 

        l_alpha = l_C * (1.0e0_DP / (l_la - 3.0e0_DP) - 1.0e0_DP / (l_lr - 3.0e0_DP))
        l_theta = dexp(beta_K * l_eps) - 1.0e0_DP
        
        l_zetaxs    = 0.0e0_DP
        dl_zetaxsv  = 0.0e0_DP

        do i_lngii1=1, nstypes
        do i_lngii2=1, nstypes
            l_zetaxs = l_zetaxs + Seg_array(i_lngii1)%xsi * Seg_array(i_lngii2)%xsi &
            &   * sig(i_lngii1, i_lngii2)**3.0e0_DP
        end do
        end do
        
        l_zetaxs   = l_zetaxs * PI * rhos / 6.0e0_DP * NA
        dl_zetaxsv = l_zetaxs *  drhosv / rhos

        l_gam = phik(7, 0) * (-1.0e0_DP * dtanh(phik(7, 1) * (phik(7, 2) - l_alpha)) + 1.0e0_DP) * &
        &   l_zetaxs * l_theta * dexp(phik(7, 3) * l_zetaxs + phik(7, 4) * l_zetaxs**2.0e0_DP)
        
        dl_gamv = phik(7, 0) * (-1.0e0_DP * dtanh(phik(7, 1) * (phik(7, 2) - l_alpha)) + 1.0e0_DP)   * &
        &   l_theta * dexp(phik(7, 3) * l_zetaxs + phik(7, 4) * l_zetaxs**2.0e0_DP)                  * &
        &   (dl_zetaxsv + l_zetaxs * (phik(7, 3) * dl_zetaxsv + 2.0e0_DP * phik(7, 4) * dl_zetaxsv  * &
        &   l_zetaxs))

        g2mca = 1.0e0_DP / (2.0e0_DP * PI * l_eps**2.00e0_DP * l_d**3.0e0_DP) * ( 3.0e0_DP *            &
        &   as2grad(rhos, l_C, l_x0, l_d, l_la,l_lr,l_eps) -                                            &
        &   l_eps * khs * l_C**2.0e0_DP * l_lr * l_x0**(2.0e0_DP * l_lr) * (l_a12r + l_B2r ) / rhos +   &
        &   l_eps * khs * l_C**2.0e0_DP * (l_la + l_lr) * l_x0**(l_la + l_lr)                       *   &    
        &   (l_a1ar + l_Bar) / rhos - l_eps * khs * l_C**2.0e0_DP * l_la * l_x0**(2.0e0_DP * l_la)  *   &
        &   (l_a12a + l_B2a) / rhos)
        
        dg2mcav = 1.0e0_DP / (2.0e0_DP * PI * l_eps**2.00e0_DP * l_d**3.0e0_DP)                         * &
        &   ( 3.0e0_DP * das2gradv(rhos, l_C, l_x0, l_d, l_la,l_lr,l_eps)                               - &

        &   l_eps * l_C**2.0e0_DP * l_lr * l_x0**(2.0e0_DP * l_lr)                              * ( &
        &   khs * (rhos * (dl_a12rv + dl_B2rv ) - drhosv * (l_a12r + l_B2r )) / rhos**2.0e0_DP  +   &
        &   dkhsv * (l_a12r + l_B2r ) / rhos)                                                   +   &  

        &   l_eps * l_C**2.0e0_DP * (l_la + l_lr) * l_x0**(l_la + l_lr)                         * ( &
        &   khs * (rhos * (dl_a1arv + dl_Barv ) - drhosv * (l_a1ar + l_Bar )) / rhos**2.0e0_DP  +   &
        &   dkhsv * (l_a1ar + l_Bar ) / rhos)                                                   -   &  
        
        &   l_eps * l_C**2.0e0_DP * l_la * l_x0**(2.0e0_DP * l_la)                              * ( &
        &   khs * (rhos * (dl_a12av + dl_B2av ) - drhosv * (l_a12a + l_B2a )) / rhos**2.0e0_DP  +   &
        &   dkhsv * (l_a12a + l_B2a ) / rhos)) 

        l_G2 = (1.0e0_DP + l_gam) * g2mca
        
        dl_G2v = (1.0e0_DP + l_gam) * dg2mcav + dl_gamv * g2mca

        gii = l_ghs * dexp(beta_K * l_eps * l_G1 / l_ghs + (beta_K * l_eps)**2.0e0_DP * l_G2 / l_ghs)     
        
        dgiiv = dexp(beta_K * l_eps * l_G1 / l_ghs + (beta_K * l_eps)**2.0e0_DP * l_G2 / l_ghs)      * &
        &   (dl_ghsv + l_ghs                                                                        * &
        &   (beta_K * l_eps * (dl_G1v*l_ghs - l_G1*dl_ghsv)/l_ghs**2.0e0_DP                         + &
        &   (beta_K * l_eps)**2.0e0_DP * (dl_G2v*l_ghs - l_G2*dl_ghsv)/l_ghs**2.0e0_DP))           
     
        lngii_out = dgiiv / gii

        return
    end function dlngiiv
!*************************************************************************************************** 
    function bchain( lam, xo, l_dd, l_e ) result(bchain_out)
        implicit none
        
        real(kind=DP), intent(in)   ::  lam, xo, l_dd, l_e
        real(kind=DP)               ::  ilij, jlij
        real(kind=DP)               ::  bchain_out
        
        ilij = -(xo**(3.0e0_DP - lam) - 1.0e0_DP) / (lam - 3.0e0_DP)  
        jlij = -(xo**(4.0e0_DP - lam) * (lam - 3.0e0_DP) - xo**(3.0e0_DP - lam) *  &
        &       (lam - 4.0e0_DP) - 1.0e0_DP) / ((lam - 3.0e0_DP) * (lam - 4.0e0_DP))

        bchain_out = TWOPI * rhos * l_dd**3.0e0_DP * l_e *                          &
        &   (   ((1.0e0_DP - zetax / 2.0e0_DP) / (1.0e0_DP - zetax)**3.0e0_DP) *    &
        &   ilij - (9.0e0_DP * zetax * (1.0e0_DP + zetax) / (2.0e0_DP*(1.0e0_DP -   &
        &   zetax)**3.0e0_DP)) * jlij  ) 

        return
    end function bchain 
!*************************************************************************************************** 
    function dbchainn( lam, xo, l_dd, l_e ) result(bchain_out)
        implicit none
        
        real(kind=DP), intent(in)   ::  lam, xo, l_dd, l_e
        real(kind=DP)               ::  ilij, jlij
        real(kind=DP)               ::  bchain_out
        real(kind=DP)               ::  lb_1,lb_1p,lb_2,lb_2p,lb_3,lb_3p
        
        ilij = -(xo**(3.0e0_DP - lam) - 1.0e0_DP) / (lam - 3.0e0_DP)  
        jlij = -(xo**(4.0e0_DP - lam) * (lam - 3.0e0_DP) - xo**(3.0e0_DP - lam) *  &
        &       (lam - 4.0e0_DP) - 1.0e0_DP) / ((lam - 3.0e0_DP) * (lam - 4.0e0_DP))
        
        lb_1  = TWOPI * rhos * l_dd**3.0e0_DP * l_e
        lb_1p = TWOPI * drhosn * l_dd**3.0e0_DP * l_e
        
        lb_2  = ((1.0e0_DP - zetax / 2.0e0_DP) / (1.0e0_DP - zetax)**3.0e0_DP) * ilij
        lb_2p = ilij * (-1.0e0_DP*dzetaxn / 2.0e0_DP * (1.0e0_DP - zetax)**3.0e0_DP + 3.0e0_DP*dzetaxn   * &
        &   (1.0e0_DP - zetax)**2.0e0_DP *(1.0e0_DP - zetax / 2.0e0_DP))                                / &
        &   (1.0e0_DP - zetax)**6.0e0_DP
        
        lb_3  = (9.0e0_DP * zetax * (1.0e0_DP + zetax) / (2.0e0_DP*(1.0e0_DP - zetax)**3.0e0_DP)) * jlij
        lb_3p = jlij * ((2.0e0_DP*(1.0e0_DP - zetax)**3.0e0_DP) * 9.0e0_DP * (dzetaxn * (1.0e0_DP + zetax) + &
        &   zetax * dzetaxn) + 54.0e0_DP*dzetaxn*(1.0e0_DP - zetax)**2.0e0_DP * zetax * (1.0e0_DP + zetax)) / &
        &   (4.0e0_DP*(1.0e0_DP - zetax)**6.0e0_DP)
        
        bchain_out = lb_1p * (lb_2 - lb_3) + lb_1 * (lb_2p - lb_3p) 

        return
    end function dbchainn
!***************************************************************************************************
    function dbchainv( lam, xo, l_dd, l_e ) result(bchain_out)
        implicit none
        
        real(kind=DP), intent(in)   ::  lam, xo, l_dd, l_e
        real(kind=DP)               ::  ilij, jlij
        real(kind=DP)               ::  bchain_out
        real(kind=DP)               ::  lb_1,lb_1p,lb_2,lb_2p,lb_3,lb_3p
        
        ilij = -(xo**(3.0e0_DP - lam) - 1.0e0_DP) / (lam - 3.0e0_DP)  
        jlij = -(xo**(4.0e0_DP - lam) * (lam - 3.0e0_DP) - xo**(3.0e0_DP - lam) *  &
        &       (lam - 4.0e0_DP) - 1.0e0_DP) / ((lam - 3.0e0_DP) * (lam - 4.0e0_DP))
        
        lb_1  = TWOPI * rhos * l_dd**3.0e0_DP * l_e
        lb_1p = TWOPI * drhosv * l_dd**3.0e0_DP * l_e
        
        lb_2  = ((1.0e0_DP - zetax / 2.0e0_DP) / (1.0e0_DP - zetax)**3.0e0_DP) * ilij
        lb_2p = ilij * (-dzetaxv / 2.0e0_DP * (1.0e0_DP - zetax)**3.0e0_DP + 3.0e0_DP*dzetaxv   * &
        &   (1.0e0_DP - zetax)**2.0e0_DP *(1.0e0_DP - zetax / 2.0e0_DP))                                / &
        &   (1.0e0_DP - zetax)**6.0e0_DP
        
        lb_3  = (9.0e0_DP * zetax * (1.0e0_DP + zetax) / (2.0e0_DP*(1.0e0_DP - zetax)**3.0e0_DP)) * jlij
        lb_3p = jlij * ((2.0e0_DP*(1.0e0_DP - zetax)**3.0e0_DP) * 9.0e0_DP * (dzetaxv * (1.0e0_DP + zetax) + &
        &   zetax * dzetaxv) + 54.0e0_DP*dzetaxv*(1.0e0_DP - zetax)**2.0e0_DP * zetax * (1.0e0_DP + zetax)) / &
        &   (4.0e0_DP*(1.0e0_DP - zetax)**6.0e0_DP)
        
        bchain_out = lb_1p * (lb_2 - lb_3) + lb_1 * (lb_2p - lb_3p) 

        return
    end function dbchainv
!***************************************************************************************************
    function as1chain( lam, l_dd, l_e ) result(as1chain_out)
        implicit none

        real(kind=DP), intent(in)    ::  lam, l_dd, l_e
        real(kind=DP)                ::  zxf
        real(kind=DP)                ::  as1chain_out
        
        zxf = zeff(lam)
      
        as1chain_out = -2.0e0_DP * rhos * (PI * l_e) * (l_dd**3.0e0_DP) / (lam - 3.0e0_DP)   &
        &  * ((1.0e0_DP - zxf/2.0e0_DP) / (1.0e0_DP - zxf)**3.0e0_DP)

        return
    end function as1chain
!***************************************************************************************************
    pure function das1chainn( lam, l_dd, l_e ) result(as1chain_out)
        implicit none

        real(kind=DP), intent(in)    ::  lam, l_dd, l_e
        real(kind=DP)                ::  zxf, dzxfn
        real(kind=DP)                ::  as1chain_out
        real(kind=DP)                ::  las1_1,las1_1p,las1_2,las1_2p
        
        zxf = zeff(lam)
        dzxfn = dzeffn(lam)
        
        las1_1  = -2.0e0_DP * rhos * (PI * l_e) * (l_dd**3.0e0_DP) / (lam - 3.0e0_DP)
        las1_1p = -2.0e0_DP * drhosn * (PI * l_e) * (l_dd**3.0e0_DP) / (lam - 3.0e0_DP)
        las1_2  = (1.0e0_DP - zxf/2.0e0_DP) / (1.0e0_DP - zxf)**3.0e0_DP
        las1_2p =  (-dzxfn/2.0e0_DP * (1.0e0_DP - zxf)**3.0e0_DP + 3.0e0_DP * dzxfn * (1.0e0_DP - zxf)**2.0e0_DP * &
        &   (1.0e0_DP - zxf/2.0e0_DP)) / (1.0e0_DP - zxf)**6.0e0_DP
        
        as1chain_out = las1_1 * las1_2p + las1_1p * las1_2
        
        return
    end function das1chainn
!***************************************************************************************************
    function das1chainv( lam, l_dd, l_e ) result(as1chain_out)
        implicit none

        real(kind=DP), intent(in)    ::  lam, l_dd, l_e
        real(kind=DP)                ::  zxf, dzxfv
        real(kind=DP)                ::  as1chain_out
        real(kind=DP)                ::  las1_1,las1_1p,las1_2,las1_2p
        
        zxf = zeff(lam)
        dzxfv = dzeffv(lam)
      
        las1_1  = -2.0e0_DP * rhos * (PI * l_e) * (l_dd**3.0e0_DP) / (lam - 3.0e0_DP)
        las1_1p = -2.0e0_DP * drhosv * (PI * l_e) * (l_dd**3.0e0_DP) / (lam - 3.0e0_DP)
        las1_2  = (1.0e0_DP - zxf/2.0e0_DP) / (1.0e0_DP - zxf)**3.0e0_DP
        las1_2p =  (-dzxfv/2.0e0_DP * (1.0e0_DP - zxf)**3.0e0_DP + 3.0e0_DP * dzxfv * (1.0e0_DP - zxf)**2.0e0_DP * &
        &   (1.0e0_DP - zxf/2.0e0_DP)) / (1.0e0_DP - zxf)**6.0e0_DP
        
        as1chain_out = las1_1 * las1_2p + las1_1p * las1_2
        
        return
    end function das1chainv
!***************************************************************************************************
    function as1grad( cc, xx, dd, lla, llr, leps) result(as1grad_out) 
        implicit none
        
        real(kind=DP),intent(in)    ::  cc, xx, dd, lla, llr, leps
        real(kind=DP)               ::  as1grad_out
        real(kind=DP)               ::  l_da1, l_dr1, l_const1a, l_const1r, l_za, l_zr, l_zap, l_zrp
        real(kind=DP)               ::  l_dBa, l_dBr, l_Ja, l_Jr, l_Ia, l_Ir, l_const2

        l_const1a = PI * leps*dd**3.0e0_DP / (lla - 3.0e0_DP)
        l_const1r = PI * leps*dd**3.0e0_DP / (llr - 3.0e0_DP)
       
        l_za  = zeff(lla)
        l_zr  = zeff(llr)
        l_zap = zeffprime(lla)
        l_zrp = zeffprime(llr)
               
        l_da1 = -2.0e0_DP * l_const1a * ( (1.0e0_DP - l_za/2.0e0_DP) / (1.0e0_DP - l_za)**3.0e0_DP + rhos *  &
        &   (2.5e0_DP * l_zap - l_zap * l_za) / (1.0e0_DP - l_za)**4.0e0_DP)

        l_dr1 = -2.0e0_DP * l_const1r * ( (1.0e0_DP - l_zr/2.0e0_DP) / (1.0e0_DP - l_zr)**3.0e0_DP + rhos *  &
        &   (2.5e0_DP * l_zrp - l_zrp * l_zr) / (1.0e0_DP - l_zr)**4.0e0_DP)
        
        l_const2 = 2.0e0_DP * PI * dd**3.0e0_DP * leps     
        
        l_Ia = -(xx**(3.0e0_DP - lla) - 1.0e0_DP) / (lla - 3.0e0_DP)  
        l_Ir = -(xx**(3.0e0_DP - llr) - 1.0e0_DP) / (llr - 3.0e0_DP)
   
        l_Ja = -(xx**(4.0e0_DP - lla) * (lla - 3.0e0_DP) - xx**(3.0e0_DP - lla) * (lla - 4.0e0_DP) - 1.0e0_DP) / &
        &   ((lla - 3.0e0_DP) * (lla - 4.0e0_DP))
        l_Jr = -(xx**(4.0e0_DP - llr) * (llr - 3.0e0_DP) - xx**(3.0e0_DP - llr) * (llr - 4.0e0_DP) - 1.0e0_DP) / &
        &   ((llr - 3.0e0_DP) * (llr - 4.0e0_DP))   
        
        l_dBa = l_const2 * (l_Ia * (((1.0e0_DP - zetax / 2.0e0_DP) / (1.0e0_DP - zetax)**3.0e0_DP )                 +   &
                &   zetax * (2.5e0_DP - zetax) / (1.0e0_DP - zetax)**4.0e0_DP)                                      -   &   
                &   l_Ja * (zetax/( 2.0d0 * ( 1.0e0_DP - zetax )**3.0e0_DP ) * (( 9.0e0_DP * (1.0e0_DP + zetax))    +   &
                &   (( 9.0e0_DP  * ( 1.0e0_DP + 2.0e0_DP * zetax ) * 2.0e0_DP * ( 1.0e0_DP - zetax )  )             +   &
                &   ( 9.e0_DP * ( zetax + zetax**2.0e0_DP ) * 6.0e0_DP ) ) / ( 2.0e0_DP * (1.0e0_DP - zetax) ))  ))
        
        l_dBr = l_const2 * (l_Ir * (((1.0e0_DP - zetax / 2.0e0_DP) / (1.0e0_DP - zetax)**3.0e0_DP )                 +   &
                &   zetax * (2.5e0_DP - zetax) / (1.0e0_DP - zetax)**4.0e0_DP)                                      -   &   
                &   l_Jr * (zetax/( 2.0d0 * ( 1.0e0_DP - zetax )**3.0e0_DP ) * (( 9.0e0_DP * (1.0e0_DP + zetax))    +   &
                &   (( 9.0e0_DP  * ( 1.0e0_DP + 2.0e0_DP * zetax ) * 2.0e0_DP * ( 1.0e0_DP - zetax )  )             +   &
                &   ( 9.e0_DP * ( zetax + zetax**2.0e0_DP ) * 6.0e0_DP ) ) / ( 2.0e0_DP * (1.0e0_DP - zetax) ))  ))
    
        as1grad_out = cc * ( xx**lla * (l_da1 + l_dBa) - xx**llr * (l_dr1 + l_dBr) )

        return
    end function as1grad
!***************************************************************************************************
    function das1gradn( cc, xx, dd, lla, llr, leps) result(as1grad_out) 
        implicit none
        
        real(kind=DP),intent(in)    ::  cc, xx, dd, lla, llr, leps
        real(kind=DP)               ::  as1grad_out
        real(kind=DP)               ::  l_da1, l_dr1, l_const1a, l_const1r, l_za, l_zr, l_zap, l_zrp
        real(kind=DP)               ::  l_dBa, l_dBr, l_Ja, l_Jr, l_Ia, l_Ir, l_const2
        !differentials
        real(kind=DP)               ::  dl_zan, dl_zrn, dl_zapn, dl_zrpn
        real(kind=DP)               ::  dl_da1n, dl_dr1n, dl_dBan, dl_dBrn
        real(kind=DP)               ::  temp_1,temp_2,temp_3,temp_4,temp_5,temp_6,temp_7
        real(kind=DP)               ::  temp_1p,temp_2p,temp_3p,temp_4p,temp_5p,temp_6p,temp_7p

        l_const1a = PI * leps*dd**3.0e0_DP / (lla - 3.0e0_DP)
        l_const1r = PI * leps*dd**3.0e0_DP / (llr - 3.0e0_DP)
       
        l_za  = zeff(lla)
        l_zr  = zeff(llr)
        l_zap = zeffprime(lla)
        l_zrp = zeffprime(llr)

        dl_zan  = dzeffn(lla)
        dl_zrn  = dzeffn(llr)
        dl_zapn = dzeffprimen(lla)
        dl_zrpn = dzeffprimen(llr)
             
        l_da1 = -2.0e0_DP * l_const1a * ( (1.0e0_DP - l_za/2.0e0_DP) / (1.0e0_DP - l_za)**3.0e0_DP + rhos *  &
        &   (2.5e0_DP * l_zap - l_zap * l_za) / (1.0e0_DP - l_za)**4.0e0_DP)

        l_dr1 = -2.0e0_DP * l_const1r * ( (1.0e0_DP - l_zr/2.0e0_DP) / (1.0e0_DP - l_zr)**3.0e0_DP + rhos *  &
        &   (2.5e0_DP * l_zrp - l_zrp * l_zr) / (1.0e0_DP - l_zr)**4.0e0_DP)

        dl_da1n = -2.0e0_DP * l_const1a                                                                     * &
        &   ( ((1.0e0_DP - l_za)**3.0e0_DP*(-dl_zan/2.0e0_DP)                                               + &
        &   3.0e0_DP * dl_zan * (1.0e0_DP - l_za)**2.0e0_DP * (1.0e0_DP - l_za/2.0e0_DP) )                  / &
        &   (1.0e0_DP - l_za)**6.0e0_DP                                                                     + &
        &   drhosn * (2.5e0_DP * l_zap - l_zap * l_za) / (1.0e0_DP - l_za)**4.0e0_DP                        + &
        &   rhos * ((1.0e0_DP - l_za)**4.0e0_DP*(2.5e0_DP * dl_zapn - (dl_zapn*l_za + l_zap*dl_zan))        + &
        &   4.0e0_DP * dl_zan                                                                              * &
        &   (1.0e0_DP - l_za)**3.0e0_DP * (2.5e0_DP * l_zap - l_zap * l_za)) / (1.0e0_DP - l_za)**8.0e0_DP)

        dl_dr1n = -2.0e0_DP * l_const1r                                                                     * &
        &   ( ((1.0e0_DP - l_zr)**3.0e0_DP*(-dl_zrn/2.0e0_DP)                                               + &
        &   3.0e0_DP * dl_zrn * (1.0e0_DP - l_zr)**2.0e0_DP * (1.0e0_DP - l_zr/2.0e0_DP) )                  / &
        &   (1.0e0_DP - l_zr)**6.0e0_DP                                                                     + &
        &   drhosn * (2.5e0_DP * l_zrp - l_zrp * l_zr) / (1.0e0_DP - l_zr)**4.0e0_DP                        + &
        &   rhos * ((1.0e0_DP - l_zr)**4.0e0_DP*(2.5e0_DP * dl_zrpn - (dl_zrpn*l_zr + l_zrp*dl_zrn))        + &
        &   4.0e0_DP * dl_zrn                                                                              * &
        &   (1.0e0_DP - l_zr)**3.0e0_DP * (2.5e0_DP * l_zrp - l_zrp * l_zr)) / (1.0e0_DP - l_zr)**8.0e0_DP)

        l_const2 = 2.0e0_DP * PI * dd**3.0e0_DP * leps     
        
        l_Ia = -(xx**(3.0e0_DP - lla) - 1.0e0_DP) / (lla - 3.0e0_DP)  
        l_Ir = -(xx**(3.0e0_DP - llr) - 1.0e0_DP) / (llr - 3.0e0_DP)
   
        l_Ja = -(xx**(4.0e0_DP - lla) * (lla - 3.0e0_DP) - xx**(3.0e0_DP - lla) * (lla - 4.0e0_DP) - 1.0e0_DP) / &
        &   ((lla - 3.0e0_DP) * (lla - 4.0e0_DP))
        l_Jr = -(xx**(4.0e0_DP - llr) * (llr - 3.0e0_DP) - xx**(3.0e0_DP - llr) * (llr - 4.0e0_DP) - 1.0e0_DP) / &
        &   ((llr - 3.0e0_DP) * (llr - 4.0e0_DP))   

        l_dBa = l_const2 * (l_Ia * (((1.0e0_DP - zetax / 2.0e0_DP) / (1.0e0_DP - zetax)**3.0e0_DP )                 +   &
                &   zetax * (2.5e0_DP - zetax) / (1.0e0_DP - zetax)**4.0e0_DP)                                      -   &   
                &   l_Ja * (zetax/( 2.0d0 * ( 1.0e0_DP - zetax )**3.0e0_DP ) * (( 9.0e0_DP * (1.0e0_DP + zetax))    +   &
                &   (( 9.0e0_DP  * ( 1.0e0_DP + 2.0e0_DP * zetax ) * 2.0e0_DP * ( 1.0e0_DP - zetax )  )             +   &
                &   ( 9.e0_DP * ( zetax + zetax**2.0e0_DP ) * 6.0e0_DP ) ) / ( 2.0e0_DP * (1.0e0_DP - zetax) ))  ))
        
        l_dBr = l_const2 * (l_Ir * (((1.0e0_DP - zetax / 2.0e0_DP) / (1.0e0_DP - zetax)**3.0e0_DP )                 +   &
                &   zetax * (2.5e0_DP - zetax) / (1.0e0_DP - zetax)**4.0e0_DP)                                      -   &   
                &   l_Jr * (zetax/( 2.0d0 * ( 1.0e0_DP - zetax )**3.0e0_DP ) * (( 9.0e0_DP * (1.0e0_DP + zetax))    +   &
                &   (( 9.0e0_DP  * ( 1.0e0_DP + 2.0e0_DP * zetax ) * 2.0e0_DP * ( 1.0e0_DP - zetax )  )             +   &
                &   ( 9.e0_DP * ( zetax + zetax**2.0e0_DP ) * 6.0e0_DP ) ) / ( 2.0e0_DP * (1.0e0_DP - zetax) ))  ))
        
        temp_1  = (1.0e0_DP - zetax / 2.0e0_DP) / (1.0e0_DP - zetax)**3.0e0_DP   
        temp_1p = ((-dzetaxn / 2.0e0_DP)*(1.0e0_DP - zetax)**3.0e0_DP + 3.0e0_DP*dzetaxn*(1.0e0_DP              - &
        &   zetax)**2.0e0_DP)*(1.0e0_DP - zetax / 2.0e0_DP) / (1.0e0_DP - zetax)**6.0e0_DP  
        
        temp_2  = zetax * (2.5e0_DP - zetax) / (1.0e0_DP - zetax)**4.0e0_DP
        temp_2p = dzetaxn * (2.5e0_DP - zetax) / (1.0e0_DP - zetax)**4.0e0_DP                                   + &
        &   zetax * (-dzetaxn*(1.0e0_DP - zetax)**4.0e0_DP + 4.0e0_DP*dzetaxn*(1.0e0_DP - zetax)**3.0e0_DP      * &
        &   (2.5e0_DP - zetax)) /  (1.0e0_DP - zetax)**8.0e0_DP
        
        temp_3  = zetax/( 2.0d0 * ( 1.0e0_DP - zetax )**3.0e0_DP )
        temp_3p = (dzetaxn * ( 2.0d0 * ( 1.0e0_DP - zetax )**3.0e0_DP ) + 6.0d0 * dzetaxn * zetax               * &
        &   ( 1.0e0_DP - zetax )**2.0e0_DP) / ( 4.0d0 * ( 1.0e0_DP - zetax )**6.0e0_DP )

        temp_4  = 9.0e0_DP * (1.0e0_DP + zetax) 
        temp_4p = 9.0e0_DP * dzetaxn
 
        temp_5  = 18.0e0_DP * ( 1.0e0_DP + 2.0e0_DP * zetax ) * ( 1.0e0_DP - zetax )
        temp_5p = 18.0e0_DP * (( 2.0e0_DP * dzetaxn ) * ( 1.0e0_DP - zetax )                                    - &
        &    dzetaxn * ( 1.0e0_DP + 2.0e0_DP * zetax ))
        
        temp_6  = 54.0e0_DP * ( zetax + zetax**2.0e0_DP ) 
        temp_6p = 54.0e0_DP * ( dzetaxn + 2.0e0_DP * dzetaxn * zetax ) 
        
        temp_7  = 2.0e0_DP * (1.0e0_DP - zetax)
        temp_7p = -2.0e0_DP * dzetaxn
        
        dl_dBan = l_const2 * (l_Ia * (temp_1p + temp_2p)                                                        - &
        &   l_Ja * (temp_3p*(temp_4 +(temp_5+temp_6)/temp_7)                                                    + &
        &   temp_3*(temp_4p + ((temp_5p+temp_6p)*temp_7 - temp_7p*(temp_5+temp_6))/temp_7**2.0e0_DP)))
        
        dl_dBrn = l_const2 * (l_Ir * (temp_1p + temp_2p)                                                        - &
        &   l_Jr * (temp_3p*(temp_4 +(temp_5+temp_6)/temp_7)                                                    + &
        &   temp_3*(temp_4p + ((temp_5p+temp_6p)*temp_7 - temp_7p*(temp_5+temp_6))/temp_7**2.0e0_DP)))         
      
        as1grad_out = cc * ( xx**lla * (dl_da1n + dl_dBan) - xx**llr * (dl_dr1n + dl_dBrn) )

        return
    end function das1gradn
!***************************************************************************************************
    function das1gradv( cc, xx, dd, lla, llr, leps) result(as1grad_out) 
        implicit none
        
        real(kind=DP),intent(in)    ::  cc, xx, dd, lla, llr, leps
        real(kind=DP)               ::  as1grad_out
        real(kind=DP)               ::  l_da1, l_dr1, l_const1a, l_const1r, l_za, l_zr, l_zap, l_zrp
        real(kind=DP)               ::  l_dBa, l_dBr, l_Ja, l_Jr, l_Ia, l_Ir, l_const2
        !differentials
        real(kind=DP)               ::  dl_zav, dl_zrv, dl_zapv, dl_zrpv
        real(kind=DP)               ::  dl_da1v, dl_dr1v, dl_dBav, dl_dBrv
        real(kind=DP)               ::  temp_1,temp_2,temp_3,temp_4,temp_5,temp_6,temp_7
        real(kind=DP)               ::  temp_1p,temp_2p,temp_3p,temp_4p,temp_5p,temp_6p,temp_7p

        l_const1a = PI * leps*dd**3.0e0_DP / (lla - 3.0e0_DP)
        l_const1r = PI * leps*dd**3.0e0_DP / (llr - 3.0e0_DP)
       
        l_za  = zeff(lla)
        l_zr  = zeff(llr)
        l_zap = zeffprime(lla)
        l_zrp = zeffprime(llr)

        dl_zav  = dzeffv(lla)
        dl_zrv  = dzeffv(llr)
        dl_zapv = dzeffprimev(lla)
        dl_zrpv = dzeffprimev(llr)
     
        l_da1 = -2.0e0_DP * l_const1a * ( (1.0e0_DP - l_za/2.0e0_DP) / (1.0e0_DP - l_za)**3.0e0_DP + rhos *  &
        &   (2.5e0_DP * l_zap - l_zap * l_za) / (1.0e0_DP - l_za)**4.0e0_DP)

        l_dr1 = -2.0e0_DP * l_const1r * ( (1.0e0_DP - l_zr/2.0e0_DP) / (1.0e0_DP - l_zr)**3.0e0_DP + rhos *  &
        &   (2.5e0_DP * l_zrp - l_zrp * l_zr) / (1.0e0_DP - l_zr)**4.0e0_DP)

        dl_da1v = -2.0e0_DP * l_const1a                                                                     * &
        &   ( ((1.0e0_DP - l_za)**3.0e0_DP*(-dl_zav/2.0e0_DP)                                               + &
        &   3.0e0_DP * dl_zav * (1.0e0_DP - l_za)**2.0e0_DP * (1.0e0_DP - l_za/2.0e0_DP) )                  / &
        &   (1.0e0_DP - l_za)**6.0e0_DP                                                                     + &
        &   drhosv * (2.5e0_DP * l_zap - l_zap * l_za) / (1.0e0_DP - l_za)**4.0e0_DP                        + &
        &   rhos * ((1.0e0_DP - l_za)**4.0e0_DP*(2.5e0_DP * dl_zapv - (dl_zapv*l_za + l_zap*dl_zav))        + &
        &   4.0e0_DP * dl_zav  * &!ALTERED THESE dl_zav from dl_zapv                                                                            
        &   (1.0e0_DP - l_za)**3.0e0_DP * (2.5e0_DP * l_zap - l_zap * l_za)) / (1.0e0_DP - l_za)**8.0e0_DP)

        dl_dr1v = -2.0e0_DP * l_const1r                                                                     * &
        &   ( ((1.0e0_DP - l_zr)**3.0e0_DP*(-dl_zrv/2.0e0_DP)                                               + &
        &   3.0e0_DP * dl_zrv * (1.0e0_DP - l_zr)**2.0e0_DP * (1.0e0_DP - l_zr/2.0e0_DP) )                  / &
        &   (1.0e0_DP - l_zr)**6.0e0_DP                                                                     + &
        &   drhosv * (2.5e0_DP * l_zrp - l_zrp * l_zr) / (1.0e0_DP - l_zr)**4.0e0_DP                        + &
        &   rhos * ((1.0e0_DP - l_zr)**4.0e0_DP*(2.5e0_DP * dl_zrpv - (dl_zrpv*l_zr + l_zrp*dl_zrv))        + &
        &   4.0e0_DP * dl_zrv  * &!ALTERED THESE dl_zrv from dl_zrpv  
        &   (1.0e0_DP - l_zr)**3.0e0_DP * (2.5e0_DP * l_zrp - l_zrp * l_zr)) / (1.0e0_DP - l_zr)**8.0e0_DP)

        l_const2 = 2.0e0_DP * PI * dd**3.0e0_DP * leps     
        
        l_Ia = -(xx**(3.0e0_DP - lla) - 1.0e0_DP) / (lla - 3.0e0_DP)  
        l_Ir = -(xx**(3.0e0_DP - llr) - 1.0e0_DP) / (llr - 3.0e0_DP)
   
        l_Ja = -(xx**(4.0e0_DP - lla) * (lla - 3.0e0_DP) - xx**(3.0e0_DP - lla) * (lla - 4.0e0_DP) - 1.0e0_DP) / &
        &   ((lla - 3.0e0_DP) * (lla - 4.0e0_DP))
        l_Jr = -(xx**(4.0e0_DP - llr) * (llr - 3.0e0_DP) - xx**(3.0e0_DP - llr) * (llr - 4.0e0_DP) - 1.0e0_DP) / &
        &   ((llr - 3.0e0_DP) * (llr - 4.0e0_DP))   

        l_dBa = l_const2 * (l_Ia * (((1.0e0_DP - zetax / 2.0e0_DP) / (1.0e0_DP - zetax)**3.0e0_DP )                 +   &
                &   zetax * (2.5e0_DP - zetax) / (1.0e0_DP - zetax)**4.0e0_DP)                                      -   &   
                &   l_Ja * (zetax/( 2.0d0 * ( 1.0e0_DP - zetax )**3.0e0_DP ) * (( 9.0e0_DP * (1.0e0_DP + zetax))    +   &
                &   (( 9.0e0_DP  * ( 1.0e0_DP + 2.0e0_DP * zetax ) * 2.0e0_DP * ( 1.0e0_DP - zetax )  )             +   &
                &   ( 9.e0_DP * ( zetax + zetax**2.0e0_DP ) * 6.0e0_DP ) ) / ( 2.0e0_DP * (1.0e0_DP - zetax) ))  ))
        
        l_dBr = l_const2 * (l_Ir * (((1.0e0_DP - zetax / 2.0e0_DP) / (1.0e0_DP - zetax)**3.0e0_DP )                 +   &
                &   zetax * (2.5e0_DP - zetax) / (1.0e0_DP - zetax)**4.0e0_DP)                                      -   &   
                &   l_Jr * (zetax/( 2.0d0 * ( 1.0e0_DP - zetax )**3.0e0_DP ) * (( 9.0e0_DP * (1.0e0_DP + zetax))    +   &
                &   (( 9.0e0_DP  * ( 1.0e0_DP + 2.0e0_DP * zetax ) * 2.0e0_DP * ( 1.0e0_DP - zetax )  )             +   &
                &   ( 9.e0_DP * ( zetax + zetax**2.0e0_DP ) * 6.0e0_DP ) ) / ( 2.0e0_DP * (1.0e0_DP - zetax) ))  ))
        
        temp_1  = (1.0e0_DP - zetax / 2.0e0_DP) / (1.0e0_DP - zetax)**3.0e0_DP   
        temp_1p = ((-dzetaxv / 2.0e0_DP)*(1.0e0_DP - zetax)**3.0e0_DP + 3.0e0_DP*dzetaxv*(1.0e0_DP              - &
        &   zetax)**2.0e0_DP)*(1.0e0_DP - zetax / 2.0e0_DP) / (1.0e0_DP - zetax)**6.0e0_DP 
        
        temp_2  = zetax * (2.5e0_DP - zetax) / (1.0e0_DP - zetax)**4.0e0_DP
        temp_2p = dzetaxv * (2.5e0_DP - zetax) / (1.0e0_DP - zetax)**4.0e0_DP                                   + &
        &   zetax * (-dzetaxv*(1.0e0_DP - zetax)**4.0e0_DP + 4.0e0_DP*dzetaxv*(1.0e0_DP - zetax)**3.0e0_DP      * &
        &   (2.5e0_DP - zetax)) /  (1.0e0_DP - zetax)**8.0e0_DP
        
        temp_3  = zetax/( 2.0d0 * ( 1.0e0_DP - zetax )**3.0e0_DP )
        temp_3p = (dzetaxv * ( 2.0d0 * ( 1.0e0_DP - zetax )**3.0e0_DP ) + 6.0d0 * dzetaxv * zetax               * &
        &   ( 1.0e0_DP - zetax )**2.0e0_DP) / ( 4.0d0 * ( 1.0e0_DP - zetax )**6.0e0_DP )

        temp_4  = 9.0e0_DP * (1.0e0_DP + zetax) 
        temp_4p = 9.0e0_DP * dzetaxv
 
        temp_5  = 18.0e0_DP * ( 1.0e0_DP + 2.0e0_DP * zetax ) * ( 1.0e0_DP - zetax )
        temp_5p = 18.0e0_DP * (( 2.0e0_DP * dzetaxv ) * ( 1.0e0_DP - zetax )                                    - &
        &    dzetaxv * ( 1.0e0_DP + 2.0e0_DP * zetax ))
        
        temp_6  = 54.0e0_DP * ( zetax + zetax**2.0e0_DP ) 
        temp_6p = 54.0e0_DP * ( dzetaxv + 2.0e0_DP * dzetaxv * zetax ) 
        
        temp_7  = 2.0e0_DP * (1.0e0_DP - zetax)
        temp_7p = -2.0e0_DP * dzetaxv
        
        dl_dBav = l_const2 * (l_Ia * (temp_1p + temp_2p)                                                        - &
        &   l_Ja * (temp_3p*(temp_4 +(temp_5+temp_6)/temp_7)                                                    + &
        &   temp_3*(temp_4p + ((temp_5p+temp_6p)*temp_7 - temp_7p*(temp_5+temp_6))/temp_7**2.0e0_DP)))
        
        dl_dBrv = l_const2 * (l_Ir * (temp_1p + temp_2p)                                                        - &
        &   l_Jr * (temp_3p*(temp_4 +(temp_5+temp_6)/temp_7)                                                    + &
        &   temp_3*(temp_4p + ((temp_5p+temp_6p)*temp_7 - temp_7p*(temp_5+temp_6))/temp_7**2.0e0_DP)))      
      
        as1grad_out = cc * ( xx**lla * (dl_da1v + dl_dBav) - xx**llr * (dl_dr1v + dl_dBrv) )

        return
    end function das1gradv
!***************************************************************************************************
     pure function bchain2( lam, xo, l_dd, l_e, l_r ) result(bchain2_out)
        implicit none
        
        real(kind=DP), intent(in)   ::  lam, xo, l_dd, l_e, l_r
        real(kind=DP)               ::  ilij, jlij
        real(kind=DP)               ::  bchain2_out
        
        ilij = -(xo**(3.0e0_DP - lam) - 1.0e0_DP) / (lam - 3.0e0_DP)  
        jlij = -(xo**(4.0e0_DP - lam) * (lam - 3.0e0_DP) - xo**(3.0e0_DP - lam) * (lam - 4.0e0_DP) - 1.0e0_DP) &
        &   / ((lam - 3.0e0_DP) * (lam - 4.0e0_DP))

        bchain2_out = TWOPI * l_r * l_dd**3.0e0_DP * l_e * &
        &   (   ((1.0e0_DP - zetax / 2.0e0_DP) / (1.0e0_DP - zetax)**3.0e0_DP) * ilij - (9.0e0_DP * zetax *     &
        &   (1.0e0_DP + zetax) / (2.0e0_DP * (1.0e0_DP - zetax)**3.0e0_DP)) * jlij  ) 

        return
    end function bchain2
!***************************************************************************************************   
     pure function dbchain2n( lam, xo, l_dd, l_e, l_r, dl_r ) result(bchain2_out)
        implicit none
        
        real(kind=DP), intent(in)   ::  lam, xo, l_dd, l_e, l_r, dl_r
        real(kind=DP)               ::  ilij, jlij
        real(kind=DP)               ::  bchain2_out
        real(kind=DP)               ::  lb_1,lb_1p,lb_2,lb_2p,lb_3,lb_3p
        
        ilij = -(xo**(3.0e0_DP - lam) - 1.0e0_DP) / (lam - 3.0e0_DP)  
        jlij = -(xo**(4.0e0_DP - lam) * (lam - 3.0e0_DP) - xo**(3.0e0_DP - lam) * (lam - 4.0e0_DP) - 1.0e0_DP) &
        &   / ((lam - 3.0e0_DP) * (lam - 4.0e0_DP))

        lb_1  = TWOPI * l_r * l_dd**3.0e0_DP * l_e 
        lb_1p = TWOPI * dl_r * l_dd**3.0e0_DP * l_e
        
        lb_2  = ((1.0e0_DP - zetax / 2.0e0_DP) / (1.0e0_DP - zetax)**3.0e0_DP) * ilij
        lb_2p = ilij * (-dzetaxn / 2.0e0_DP * (1.0e0_DP - zetax)**3.0e0_DP + 3.0e0_DP * dzetaxn             * &
        &   (1.0e0_DP - zetax)**2.0e0_DP * (1.0e0_DP - zetax / 2.0e0_DP)) / (1.0e0_DP - zetax)**6.0e0_DP
        
        lb_3  = (9.0e0_DP * zetax * (1.0e0_DP + zetax) / (2.0e0_DP * (1.0e0_DP - zetax)**3.0e0_DP)) * jlij
        lb_3p = -jlij * ( 2.0e0_DP * (1.0e0_DP - zetax)**3.0e0_DP * (9.0e0_DP * dzetaxn * (1.0e0_DP + zetax) + &
        &   54.0e0_DP * zetax * (1.0e0_DP + zetax) * dzetaxn * (1.0e0_DP - zetax)**2.0e0_DP     )    )      / &
        &   (4.0e0_DP *  (1.0e0_DP - zetax)**6.0e0_DP )

        bchain2_out = lb_1p * (lb_2 - lb_3) + lb_1 * (lb_2p - lb_3p) 

        return
    end function dbchain2n
!***************************************************************************************************     
     pure function dbchain2v( lam, xo, l_dd, l_e, l_r, dl_r ) result(bchain2_out)
        implicit none
        
        real(kind=DP), intent(in)   ::  lam, xo, l_dd, l_e, l_r, dl_r
        real(kind=DP)               ::  ilij, jlij
        real(kind=DP)               ::  bchain2_out
        real(kind=DP)               ::  lb_1,lb_1p,lb_2,lb_2p,lb_3,lb_3p
        
        ilij = -(xo**(3.0e0_DP - lam) - 1.0e0_DP) / (lam - 3.0e0_DP)  
        jlij = -(xo**(4.0e0_DP - lam) * (lam - 3.0e0_DP) - xo**(3.0e0_DP - lam) * (lam - 4.0e0_DP) - 1.0e0_DP) &
        &   / ((lam - 3.0e0_DP) * (lam - 4.0e0_DP))

        lb_1  = TWOPI * l_r * l_dd**3.0e0_DP * l_e 
        lb_1p = TWOPI * dl_r * l_dd**3.0e0_DP * l_e
        
        lb_2  = ((1.0e0_DP - zetax / 2.0e0_DP) / (1.0e0_DP - zetax)**3.0e0_DP) * ilij
        lb_2p = ilij * (-dzetaxv / 2.0e0_DP * (1.0e0_DP - zetax)**3.0e0_DP + 3.0e0_DP * dzetaxv             * &
        &   (1.0e0_DP - zetax)**2.0e0_DP * (1.0e0_DP - zetax / 2.0e0_DP)) / (1.0e0_DP - zetax)**6.0e0_DP
        
        lb_3  = (9.0e0_DP * zetax * (1.0e0_DP + zetax) / (2.0e0_DP * (1.0e0_DP - zetax)**3.0e0_DP)) * jlij
        lb_3  = -jlij * ( 2.0e0_DP * (1.0e0_DP - zetax)**3.0e0_DP * (9.0e0_DP * dzetaxn * (1.0e0_DP + zetax) + &
        &   54.0e0_DP * zetax * (1.0e0_DP + zetax) * dzetaxv * (1.0e0_DP - zetax)**2.0e0_DP     )    )      / &
        &   (4.0e0_DP *  (1.0e0_DP - zetax)**6.0e0_DP )

        bchain2_out = lb_1p * (lb_2 - lb_3) + lb_1 * (lb_2p - lb_3p) 

        return
    end function dbchain2v
!***************************************************************************************************   
     function as1chain2( lam, l_dd, l_e, l_r ) result(as1chain2_out)
        implicit none

        real(kind=DP), intent(in)   ::  lam, l_dd, l_e, l_r
        real(kind=DP)               ::  zxf
        real(kind=DP)               ::  as1chain2_out
     
        zxf = zeff(lam)

        as1chain2_out = -2.0e0_DP * l_r * (PI * l_e) * (l_dd**3.0e0_DP) / (lam - 3.0e0_DP) *   &
            &  ((1.0e0_DP - zxf / 2.0e0_DP) / (1.0e0_DP - zxf)**3.0e0_DP)

        return
    end function as1chain2
!***************************************************************************************************
     function das1chain2n( lam, l_dd, l_e, l_r, dl_r ) result(as1chain2_out)
        implicit none

        real(kind=DP), intent(in)   ::  lam, l_dd, l_e, l_r, dl_r
        real(kind=DP)               ::  zxf, dzxf
        real(kind=DP)               ::  as1chain2_out
        real(kind=DP)               ::  las1_1,las1_1p,las1_2,las1_2p,las1_3,las1_3p
        
        zxf  = zeff(lam)
        dzxf = dzeffn(lam)
        
        las1_1  = -2.0e0_DP * l_r * (PI * l_e) * (l_dd**3.0e0_DP) / (lam - 3.0e0_DP)
        las1_1p = -2.0e0_DP * dl_r * (PI * l_e) * (l_dd**3.0e0_DP) / (lam - 3.0e0_DP)
        
        las1_2  = 1.0e0_DP - zxf / 2.0e0_DP
        las1_2p = -dzxf / 2.0e0_DP
        
        las1_3  = (1.0e0_DP - zxf)**3.0e0_DP
        las1_3p = -3.0e0_DP * dzxf * (1.0e0_DP - zxf)**2.0e0_DP
        
         
        as1chain2_out = las1_1p * las1_2/las1_3 + las1_1 * (las1_3*las1_2p - las1_2*las1_3p)/las1_3**2.0e0_DP

        return
    end function das1chain2n
!***************************************************************************************************
    function das1chain2v( lam, l_dd, l_e, l_r, dl_r ) result(as1chain2_out)
        implicit none

        real(kind=DP), intent(in)   ::  lam, l_dd, l_e, l_r, dl_r
        real(kind=DP)               ::  zxf, dzxf
        real(kind=DP)               ::  as1chain2_out
        real(kind=DP)               ::  las1_1,las1_1p,las1_2,las1_2p,las1_3,las1_3p
        
        zxf  = zeff(lam)
        dzxf = dzeffv(lam)
       
        las1_1  = -2.0e0_DP * l_r * (PI * l_e) * (l_dd**3.0e0_DP) / (lam - 3.0e0_DP)
        las1_1p = -2.0e0_DP * dl_r * (PI * l_e) * (l_dd**3.0e0_DP) / (lam - 3.0e0_DP)
        
        las1_2  = 1.0e0_DP - zxf / 2.0e0_DP
        las1_2p = -dzxf / 2.0e0_DP
        
        las1_3  = (1.0e0_DP - zxf)**3.0e0_DP
        las1_3p = -3.0e0_DP * dzxf * (1.0e0_DP - zxf)**2.0e0_DP
        
         
        as1chain2_out = las1_1p * las1_2/las1_3 + las1_1 * (las1_3*las1_2p - las1_2*las1_3p)/las1_3**2.0e0_DP
        
        return
    end function das1chain2v
!***************************************************************************************************
    function as2grad( l_rhos, cc, xx, dd, lla, llr, leps ) result(as2grad_out)
        implicit none
        
        real(kind=DP), intent(in) ::  l_rhos,cc,xx,dd,lla,llr,leps
        
        real(kind=DP)       ::  dkhs, d_zetax, l_denom
        real(kind=DP)       ::  l_constaa, l_zaap, l_constar, l_zarp, l_constrr, l_zrrp, l_const2
        real(kind=DP)       ::  l_zaa, l_zar, l_zrr
        real(kind=DP)       ::  l_aa, l_ar, l_rr, l_daa, l_dar, l_drr, l_Baa, l_Bar, l_Brr, l_dBaa, l_dBar, l_dBrr
        real(kind=DP)       ::  diff_v, diff_dv
        real(kind=DP)       ::  l_Iaa, l_Iar, l_Irr, l_Jaa, l_Jar, l_Jrr
        real(kind=DP)       ::  as2grad_out
     
        d_zetax = zetax / l_rhos
    
        l_denom = (1.0e0_DP + 4.0e0_DP * zetax + 4.0e0_DP * zetax**2.0e0_DP - 4.0e0_DP * &
        &   zetax**3.0e0_DP + zetax**4.0e0_DP)
      
        dkhs = -4.0e0_DP * d_zetax * (1.0e0_DP - zetax)**3.0e0_DP   *                    &
        &   (l_denom + (1.0e0_DP + 2.0e0_DP * zetax - 3.0e0_DP * zetax**2.0e0_DP +       &
        &   zetax**3.0e0_DP) * (1.0e0_DP - zetax)) / l_denom**2.0e0_DP

        l_aa    =   as1chain( 2.0e0_DP * lla, dd, leps )
        l_ar    =   as1chain( llr + lla, dd, leps )
        l_rr    =   as1chain( 2.0e0_DP * llr, dd, leps )
       
        l_constaa = PI * leps * dd**3.0e0_DP / (2.0e0_DP * lla - 3.0e0_DP)
        l_zaap    = zeffprime( 2.0e0_DP * lla )
        l_constar = PI * leps * dd**3.0e0_DP / (llr + lla - 3.0e0_DP)
        l_zarp    = zeffprime( llr + lla )
        l_constrr = PI * leps * dd**3.0e0_DP / (2.0e0_DP * llr - 3.0e0_DP)
        l_zrrp    = zeffprime( 2.0e0_DP * llr )
        l_zaa     = zeff( 2.0e0_DP * lla )
        l_zar     = zeff( llr + lla )
        l_zrr     = zeff( 2.0e0_DP * llr )

        l_daa   =   -2.0e0_DP * l_constaa * ( (1.0e0_DP - l_zaa / 2.0e0_DP) / (1.0e0_DP - l_zaa)**3.0e0_DP +  &
        &    l_rhos * (2.5e0_DP * l_zaap - l_zaap * l_zaa) / (1.0e0_DP - l_zaa)**4.0e0_DP)
        l_dar   =   -2.0e0_DP * l_constar * ( (1.0e0_DP - l_zar / 2.0e0_DP) / (1.0e0_DP - l_zar)**3.0e0_DP +  &
        &    l_rhos * (2.5e0_DP * l_zarp - l_zarp * l_zar) / (1.0e0_DP - l_zar)**4.0e0_DP)
        l_drr   =   -2.0e0_DP * l_constrr * ( (1.0e0_DP - l_zrr / 2.0e0_DP) / (1.0e0_DP - l_zrr)**3.0e0_DP +  &
        &    l_rhos * (2.5e0_DP * l_zrrp - l_zrrp * l_zrr) / (1.0e0_DP - l_zrr)**4.0e0_DP)
        
        l_Baa   =   bchain(2.0e0_DP * lla, xx,dd, leps)
        l_Bar   =   bchain(llr + lla, xx, dd, leps)
        l_Brr   =   bchain(2.0e0_DP * llr, xx, dd, leps)

        l_Iaa = -(xx**(3.0e0_DP - 2.0e0_DP * lla) - 1.0e0_DP) / (2.0e0_DP * lla - 3.0e0_DP)  
        l_Iar = -(xx**(3.0e0_DP - llr - lla) - 1.0e0_DP) / (llr + lla - 3.0e0_DP) 
        l_Irr = -(xx**(3.0e0_DP - 2.0e0_DP * llr) - 1.0e0_DP) / (2.0e0_DP * llr - 3.0e0_DP)
        
        l_Jaa = -(xx**(4.0e0_DP - 2.0e0_DP * lla) * (2.0e0_DP * lla - 3.0e0_DP) - xx**(3.0e0_DP - 2.0e0_DP * lla) &
        &   * (2.0e0_DP * lla - 4.0e0_DP) - 1.0e0_DP) / ((2.0e0_DP * lla - 3.0e0_DP) * (2.0e0_DP * lla - 4.0e0_DP))
        l_Jar = -(xx**(4.0e0_DP - lla - llr) * (lla + llr - 3.0e0_DP) - xx**(3.0e0_DP - lla - llr) * (lla + llr - &
        &   4.0e0_DP) - 1.0e0_DP) / ((lla + llr - 3.0e0_DP) * (lla + llr - 4.0e0_DP))
        l_Jrr = -(xx**(4.0e0_DP - 2.0e0_DP * llr) * (2.0e0_DP * llr - 3.0e0_DP) - xx**(3.0e0_DP - 2.0e0_DP * llr) * &
        &   (2.0e0_DP * llr - 4.0e0_DP) - 1.0e0_DP) / ((2.0e0_DP * llr - 3.0e0_DP) * (2.0e0_DP * llr - 4.0e0_DP))   
        
        l_const2 = 2.0e0_DP * PI * dd**3.0e0_DP * leps 
       
        l_dBaa  =   l_const2 * ( (1.0e0_DP - zetax / 2.0e0_DP) / (1.0e0_DP - zetax)**3.0e0_DP * l_Iaa       -   &
        &        9.0e0_DP * zetax * (1.0e0_DP + zetax) / (2.0e0_DP * (1.0e0_DP - zetax)**3.0e0_DP) * l_Jaa  +   &
        &        l_rhos * (((-d_zetax / 2.0e0_DP * (1.0e0_DP - zetax) + 3.0e0_DP * d_zetax * (1.0e0_DP      -   &
        &        zetax /2.0e0_DP)) / (1.0e0_DP - zetax)**4.0e0_DP) * l_Iaa - (((9.0e0_DP * d_zetax          *   &
        &        (1.0e0_DP + zetax) + 9.0e0_DP * zetax * d_zetax) * 2.0e0_DP * (1.0e0_DP - zetax)           +   &
        &        6.0e0_DP * d_zetax * 9.0e0_DP * zetax * (1.0e0_DP + zetax)) / (4.0e0_DP * (1.0e0_DP        -   &
        &        zetax)**4.0e0_DP)) * l_Jaa )) 
        
        l_dBar  =  l_const2 * ( (1.0e0_DP - zetax / 2.0e0_DP) / (1.0e0_DP - zetax)**3.0e0_DP * l_Iar       -   &
        &        9.0e0_DP * zetax * (1.0e0_DP + zetax) / (2.0e0_DP * (1.0e0_DP - zetax)**3.0e0_DP) * l_Jar  +   &
        &        l_rhos * (((-d_zetax / 2.0e0_DP * (1.0e0_DP - zetax) + 3.0e0_DP * d_zetax * (1.0e0_DP      -   &
        &        zetax /2.0e0_DP)) / (1.0e0_DP - zetax)**4.0e0_DP) * l_Iar - (((9.0e0_DP * d_zetax          *   &
        &        (1.0e0_DP + zetax) + 9.0e0_DP * zetax * d_zetax) * 2.0e0_DP * (1.0e0_DP - zetax)           +   &
        &        6.0e0_DP * d_zetax * 9.0e0_DP * zetax * (1.0e0_DP + zetax)) / (4.0e0_DP * (1.0e0_DP        -   &
        &        zetax)**4.0e0_DP)) * l_Jar )) 
     
        l_dBrr  =  l_const2 * ( (1.0e0_DP - zetax / 2.0e0_DP) / (1.0e0_DP - zetax)**3.0e0_DP * l_Irr       -   &
        &        9.0e0_DP * zetax * (1.0e0_DP + zetax) / (2.0e0_DP * (1.0e0_DP - zetax)**3.0e0_DP) * l_Jrr  +   &
        &        l_rhos * (((-d_zetax / 2.0e0_DP * (1.0e0_DP - zetax) + 3.0e0_DP * d_zetax * (1.0e0_DP      -   &
        &        zetax /2.0e0_DP)) / (1.0e0_DP - zetax)**4.0e0_DP) * l_Irr - (((9.0e0_DP * d_zetax          *   &
        &        (1.0e0_DP + zetax) + 9.0e0_DP * zetax * d_zetax) * 2.0e0_DP * (1.0e0_DP - zetax)           +   &
        &        6.0e0_DP * d_zetax * 9.0e0_DP * zetax * (1.0e0_DP + zetax)) / (4.0e0_DP * (1.0e0_DP        -   &
        &        zetax)**4.0e0_DP)) * l_Jrr ))             

        diff_v  = xx**(2.0e0_DP * lla) * (l_aa + l_Baa) - 2.0e0_DP * xx**(lla + llr) * (l_ar + l_Bar) + &
        &       xx**(2.0e0_DP * llr) * (l_rr + l_Brr) 
        
        diff_dv = xx**(2.0e0_DP * lla) * (l_daa + l_dBaa) - 2.0e0_DP * xx**(lla + llr) * (l_dar + l_dBar) + &
        &       xx**(2.0e0_DP * llr) * (l_drr + l_dBrr)
       
        as2grad_out = 0.5e0_DP * leps * cc**2.0e0_DP *(diff_v * dkhs + diff_dv * khs)
       
        return
    end function as2grad
!***************************************************************************************************
    function das2gradn( l_rhos, cc, xx, dd, lla, llr, leps ) result(as2grad_out)
        implicit none
        
        real(kind=DP), intent(in) ::  l_rhos,cc,xx,dd,lla,llr,leps
        
        real(kind=DP)       ::  dkhs, d_zetax, l_denom
        real(kind=DP)       ::  l_constaa, l_zaap, l_constar, l_zarp, l_constrr, l_zrrp, l_const2
        real(kind=DP)       ::  l_zaa, l_zar, l_zrr
        real(kind=DP)       ::  l_aa, l_ar, l_rr, l_daa, l_dar, l_drr, l_Baa, l_Bar, l_Brr, l_dBaa, l_dBar, l_dBrr
        real(kind=DP)       ::  diff_v, diff_dv
        real(kind=DP)       ::  l_Iaa, l_Iar, l_Irr, l_Jaa, l_Jar, l_Jrr
        real(kind=DP)       ::  as2grad_out
        !differentials
        real(kind=DP)       ::  dd_zetaxn, dl_denomn, ddkhsn
        real(kind=DP)       ::  temp_1,temp_2,temp_3,temp_4,temp_5,temp_6,temp_7
        real(kind=DP)       ::  temp_1p,temp_2p,temp_3p,temp_4p,temp_5p,temp_6p,temp_7p
        real(kind=DP)       ::  dl_aan, dl_arn, dl_rrn
        real(kind=DP)       ::  dl_zaapn, dl_zarpn, dl_zrrpn
        real(kind=DP)       ::  dl_zaan, dl_zarn, dl_zrrn
        real(kind=DP)       ::  dl_daan, dl_darn, dl_drrn
        real(kind=DP)       ::  dl_Baan, dl_Barn, dl_Brrn, dl_dBaan, dl_dBarn, dl_dBrrn
        real(kind=DP)       ::  ddiff_dvn, ddiff_vn
        real(kind=DP)       ::  temp_as2u, temp_as2v, temp_as2up, temp_as2vp
        real(kind=DP)       ::  temp_as2_a, temp_as2_ap, temp_as2_b, temp_as2_bp, temp_as2_cp
        
        d_zetax = zetax / l_rhos
        dd_zetaxn = (dzetaxn*l_rhos - drhosn*zetax) / l_rhos**2.0e0_DP
        
        l_denom = 1.0e0_DP + 4.0e0_DP * zetax + 4.0e0_DP * zetax**2.0e0_DP - 4.0e0_DP * &
        &   zetax**3.0e0_DP + zetax**4.0e0_DP
        
        dl_denomn = 4.0e0_DP * dzetaxn + 8.0e0_DP * dzetaxn * zetax - 12.0e0_DP * &
        &   dzetaxn * zetax**2.0e0_DP + 4.0e0_DP * dzetaxn * zetax**3.0e0_DP
        
        dkhs = -4.0e0_DP * d_zetax * (1.0e0_DP - zetax)**3.0e0_DP   *                    &
        &   (l_denom + (1.0e0_DP + 2.0e0_DP * zetax - 3.0e0_DP * zetax**2.0e0_DP +       &
        &   zetax**3.0e0_DP) * (1.0e0_DP - zetax)) / l_denom**2.0e0_DP

        temp_1  = -4.0e0_DP * d_zetax
        temp_1p =  -4.0e0_DP * dd_zetaxn
        
        temp_2  = (1.0e0_DP - zetax)**3.0e0_DP
        temp_2p = -3.0e0_DP * dzetaxn * (1.0e0_DP - zetax)**2.0e0_DP
        
        temp_3  = l_denom + (1.0e0_DP + 2.0e0_DP * zetax - 3.0e0_DP * zetax**2.0e0_DP      + &
        &   zetax**3.0e0_DP) * (1.0e0_DP - zetax)
        temp_3p = dl_denomn + (1.0e0_DP + 2.0e0_DP * zetax - 3.0e0_DP * zetax**2.0e0_DP    + &
        &   zetax**3.0e0_DP) * (-dzetaxn)                                                  + &
        &   (2.0e0_DP * dzetaxn - 6.0e0_DP * dzetaxn * zetax                               + &
        &   3.0e0_DP * dzetaxn * zetax**2.0e0_DP) * (1.0e0_DP - zetax)
        
        temp_4  = l_denom**2.0e0_DP
        temp_4p = 2.0e0_DP * dl_denomn * l_denom
        
        ddkhsn = temp_1p*temp_2*temp_3/temp_4 + temp_1*temp_2p*temp_3/temp_4              + &
        &   temp_1*temp_2 * (temp_4*temp_3p - temp_4p*temp_3)/temp_4**2.0e0_DP

        l_aa    =   as1chain( 2.0e0_DP * lla, dd, leps )
        l_ar    =   as1chain( llr + lla, dd, leps )
        l_rr    =   as1chain( 2.0e0_DP * llr, dd, leps )
        
        l_constaa = PI * leps * dd**3.0e0_DP / (2.0e0_DP * lla - 3.0e0_DP)
        l_zaap    = zeffprime( 2.0e0_DP * lla )
        l_constar = PI * leps * dd**3.0e0_DP / (llr + lla - 3.0e0_DP)
        l_zarp    = zeffprime( llr + lla )
        l_constrr = PI * leps * dd**3.0e0_DP / (2.0e0_DP * llr - 3.0e0_DP)
        l_zrrp    = zeffprime( 2.0e0_DP * llr )
        l_zaa     = zeff( 2.0e0_DP * lla )
        l_zar     = zeff( llr + lla )
        l_zrr     = zeff( 2.0e0_DP * llr )

        dl_aan    =   das1chainn( 2.0e0_DP * lla, dd, leps )
        dl_arn    =   das1chainn( llr + lla, dd, leps )
        dl_rrn    =   das1chainn( 2.0e0_DP * llr, dd, leps )
        
        dl_zaapn    = dzeffprimen( 2.0e0_DP * lla )
        dl_zarpn    = dzeffprimen( llr + lla )
        dl_zrrpn    = dzeffprimen( 2.0e0_DP * llr )
        
        dl_zaan     = dzeffn( 2.0e0_DP * lla )
        dl_zarn     = dzeffn( llr + lla )
        dl_zrrn     = dzeffn( 2.0e0_DP * llr )

        l_daa   =   -2.0e0_DP * l_constaa * ( (1.0e0_DP - l_zaa / 2.0e0_DP) / (1.0e0_DP - l_zaa)**3.0e0_DP +  &
        &    l_rhos * (2.5e0_DP * l_zaap - l_zaap * l_zaa) / (1.0e0_DP - l_zaa)**4.0e0_DP)
        l_dar   =   -2.0e0_DP * l_constar * ( (1.0e0_DP - l_zar / 2.0e0_DP) / (1.0e0_DP - l_zar)**3.0e0_DP +  &
        &    l_rhos * (2.5e0_DP * l_zarp - l_zarp * l_zar) / (1.0e0_DP - l_zar)**4.0e0_DP)
        l_drr   =   -2.0e0_DP * l_constrr * ( (1.0e0_DP - l_zrr / 2.0e0_DP) / (1.0e0_DP - l_zrr)**3.0e0_DP +  &
        &    l_rhos * (2.5e0_DP * l_zrrp - l_zrrp * l_zrr) / (1.0e0_DP - l_zrr)**4.0e0_DP)
     
        temp_1  = 1.0e0_DP - l_zaa / 2.0e0_DP
        temp_1p = -dl_zaan / 2.0e0_DP
        
        temp_2  = (1.0e0_DP - l_zaa)**3.0e0_DP
        temp_2p = -3.0e0_DP * dl_zaan * (1.0e0_DP - l_zaa)**2.0e0_DP
        
        temp_3  = l_rhos * (2.5e0_DP * l_zaap - l_zaap * l_zaa)
        temp_3p = drhosn * (2.5e0_DP * l_zaap - l_zaap * l_zaa) + l_rhos * (2.5e0_DP * dl_zaapn - (dl_zaapn     * &
        &    l_zaa + l_zaap * dl_zaan))
        
        temp_4  = (1.0e0_DP - l_zaa)**4.0e0_DP
        temp_4p = -4.0e0_DP * dl_zaan * (1.0e0_DP - l_zaa)**3.0e0_DP
        
        dl_daan = -2.0e0_DP * l_constaa * ( (temp_2*temp_1p - temp_2p*temp_1)/temp_2**2.0e0_DP                  + &
        &   (temp_4*temp_3p - temp_4p*temp_3)/temp_4**2.0e0_DP)
               
        temp_1  = 1.0e0_DP - l_zar / 2.0e0_DP
        temp_1p = -dl_zarn / 2.0e0_DP
        
        temp_2  = (1.0e0_DP - l_zar)**3.0e0_DP
        temp_2p = -3.0e0_DP * dl_zarn * (1.0e0_DP - l_zar)**2.0e0_DP
        
        temp_3  = l_rhos * (2.5e0_DP * l_zarp - l_zarp * l_zar)
        temp_3p = drhosn * (2.5e0_DP * l_zarp - l_zarp * l_zar) + l_rhos * (2.5e0_DP * dl_zarpn - (dl_zarpn     * &
        &    l_zar + l_zarp * dl_zarn))
        
        temp_4  = (1.0e0_DP - l_zar)**4.0e0_DP
        temp_4p = -4.0e0_DP * dl_zarn * (1.0e0_DP - l_zar)**3.0e0_DP
        
        dl_darn = -2.0e0_DP * l_constar * ( (temp_2*temp_1p - temp_2p*temp_1)/temp_2**2.0e0_DP                  + &
        &   (temp_4*temp_3p - temp_4p*temp_3)/temp_4**2.0e0_DP)                
        
        temp_1  = 1.0e0_DP - l_zrr / 2.0e0_DP
        temp_1p = -dl_zrrn / 2.0e0_DP
        
        temp_2  = (1.0e0_DP - l_zrr)**3.0e0_DP
        temp_2p = -3.0e0_DP * dl_zrrn * (1.0e0_DP - l_zrr)**2.0e0_DP
        
        temp_3  = l_rhos * (2.5e0_DP * l_zrrp - l_zrrp * l_zrr)
        temp_3p = drhosn * (2.5e0_DP * l_zrrp - l_zrrp * l_zrr) + l_rhos * (2.5e0_DP * dl_zrrpn - (dl_zrrpn     * &
        &    l_zrr + l_zrrp * dl_zrrn))
        
        temp_4  = (1.0e0_DP - l_zrr)**4.0e0_DP
        temp_4p = -4.0e0_DP * dl_zrrn * (1.0e0_DP - l_zrr)**3.0e0_DP
        
        dl_drrn = -2.0e0_DP * l_constrr * ( (temp_2*temp_1p - temp_2p*temp_1)/temp_2**2.0e0_DP                  + &
        &   (temp_4*temp_3p - temp_4p*temp_3)/temp_4**2.0e0_DP) 

        l_Baa   =   bchain(2.0e0_DP * lla, xx,dd, leps)
        l_Bar   =   bchain(llr + lla, xx, dd, leps)
        l_Brr   =   bchain(2.0e0_DP * llr, xx, dd, leps)
        
        dl_Baan   =   dbchainn(2.0e0_DP * lla, xx,dd, leps)
        dl_Barn   =   dbchainn(llr + lla, xx, dd, leps)
        dl_Brrn   =   dbchainn(2.0e0_DP * llr, xx, dd, leps)

        l_Iaa = -(xx**(3.0e0_DP - 2.0e0_DP * lla) - 1.0e0_DP) / (2.0e0_DP * lla - 3.0e0_DP)  
        l_Iar = -(xx**(3.0e0_DP - llr - lla) - 1.0e0_DP) / (llr + lla - 3.0e0_DP) 
        l_Irr = -(xx**(3.0e0_DP - 2.0e0_DP * llr) - 1.0e0_DP) / (2.0e0_DP * llr - 3.0e0_DP)
        
        l_Jaa = -(xx**(4.0e0_DP - 2.0e0_DP * lla) * (2.0e0_DP * lla - 3.0e0_DP) - xx**(3.0e0_DP - 2.0e0_DP * lla) &
        &   * (2.0e0_DP * lla - 4.0e0_DP) - 1.0e0_DP) / ((2.0e0_DP * lla - 3.0e0_DP) * (2.0e0_DP * lla - 4.0e0_DP))
        l_Jar = -(xx**(4.0e0_DP - lla - llr) * (lla + llr - 3.0e0_DP) - xx**(3.0e0_DP - lla - llr) * (lla + llr - &
        &   4.0e0_DP) - 1.0e0_DP) / ((lla + llr - 3.0e0_DP) * (lla + llr - 4.0e0_DP))
        l_Jrr = -(xx**(4.0e0_DP - 2.0e0_DP * llr) * (2.0e0_DP * llr - 3.0e0_DP) - xx**(3.0e0_DP - 2.0e0_DP * llr) * &
        &   (2.0e0_DP * llr - 4.0e0_DP) - 1.0e0_DP) / ((2.0e0_DP * llr - 3.0e0_DP) * (2.0e0_DP * llr - 4.0e0_DP))   
        
        l_const2 = 2.0e0_DP * PI * dd**3.0e0_DP * leps 

        l_dBaa  =   l_const2 * ( (1.0e0_DP - zetax / 2.0e0_DP) / (1.0e0_DP - zetax)**3.0e0_DP * l_Iaa       -   &
        &        9.0e0_DP * zetax * (1.0e0_DP + zetax) / (2.0e0_DP * (1.0e0_DP - zetax)**3.0e0_DP) * l_Jaa  +   &
        &        l_rhos * (((-d_zetax / 2.0e0_DP * (1.0e0_DP - zetax) + 3.0e0_DP * d_zetax * (1.0e0_DP      -   &
        &        zetax /2.0e0_DP)) / (1.0e0_DP - zetax)**4.0e0_DP) * l_Iaa - (((9.0e0_DP * d_zetax          *   &
        &        (1.0e0_DP + zetax) + 9.0e0_DP * zetax * d_zetax) * 2.0e0_DP * (1.0e0_DP - zetax)           +   &
        &        6.0e0_DP * d_zetax * 9.0e0_DP * zetax * (1.0e0_DP + zetax)) / (4.0e0_DP * (1.0e0_DP        -   &
        &        zetax)**4.0e0_DP)) * l_Jaa )) 

        l_dBar  =  l_const2 * ( (1.0e0_DP - zetax / 2.0e0_DP) / (1.0e0_DP - zetax)**3.0e0_DP * l_Iar       -   &
        &        9.0e0_DP * zetax * (1.0e0_DP + zetax) / (2.0e0_DP * (1.0e0_DP - zetax)**3.0e0_DP) * l_Jar  +   &
        &        l_rhos * (((-d_zetax / 2.0e0_DP * (1.0e0_DP - zetax) + 3.0e0_DP * d_zetax * (1.0e0_DP      -   &
        &        zetax /2.0e0_DP)) / (1.0e0_DP - zetax)**4.0e0_DP) * l_Iar - (((9.0e0_DP * d_zetax          *   &
        &        (1.0e0_DP + zetax) + 9.0e0_DP * zetax * d_zetax) * 2.0e0_DP * (1.0e0_DP - zetax)           +   &
        &        6.0e0_DP * d_zetax * 9.0e0_DP * zetax * (1.0e0_DP + zetax)) / (4.0e0_DP * (1.0e0_DP        -   &
        &        zetax)**4.0e0_DP)) * l_Jar )) 
     
        l_dBrr  =  l_const2 * ( (1.0e0_DP - zetax / 2.0e0_DP) / (1.0e0_DP - zetax)**3.0e0_DP * l_Irr       -   &
        &        9.0e0_DP * zetax * (1.0e0_DP + zetax) / (2.0e0_DP * (1.0e0_DP - zetax)**3.0e0_DP) * l_Jrr  +   &
        &        l_rhos * (((-d_zetax / 2.0e0_DP * (1.0e0_DP - zetax) + 3.0e0_DP * d_zetax * (1.0e0_DP      -   &
        &        zetax /2.0e0_DP)) / (1.0e0_DP - zetax)**4.0e0_DP) * l_Irr - (((9.0e0_DP * d_zetax          *   &
        &        (1.0e0_DP + zetax) + 9.0e0_DP * zetax * d_zetax) * 2.0e0_DP * (1.0e0_DP - zetax)           +   &
        &        6.0e0_DP * d_zetax * 9.0e0_DP * zetax * (1.0e0_DP + zetax)) / (4.0e0_DP * (1.0e0_DP        -   &
        &        zetax)**4.0e0_DP)) * l_Jrr )) 
      
        temp_1  = (1.0e0_DP - zetax / 2.0e0_DP) / (1.0e0_DP - zetax)**3.0e0_DP 
        temp_1p = ((1.0e0_DP - zetax)**3.0e0_DP * (-dzetaxn / 2.0e0_DP)                         + &
        &   3.0e0_DP * dzetaxn * (1.0e0_DP - zetax / 2.0e0_DP) * (1.0e0_DP - zetax)**2.0e0_DP)  / &
        &   (1.0e0_DP - zetax)**6.0e0_DP 
        
        temp_2  = -9.0e0_DP * zetax * (1.0e0_DP + zetax) / (2.0e0_DP * (1.0e0_DP - zetax)**3.0e0_DP) 
        temp_2p = ( (2.0e0_DP * (1.0e0_DP - zetax)**3.0e0_DP) * (-9.0e0_DP * dzetaxn * (1.0e0_DP        + &
        &   zetax) -9.0e0_DP * zetax * dzetaxn) - (9.0e0_DP * zetax * (1.0e0_DP + zetax))               * &
        &   (6.0e0_DP * dzetaxn * (1.0e0_DP - zetax)**2.0e0_DP)) / (4.0e0_DP * (1.0e0_DP - zetax)**6.0e0_DP)
        
        temp_3 = (-d_zetax / 2.0e0_DP * (1.0e0_DP - zetax) + 3.0e0_DP * d_zetax * (1.0e0_DP      -   &
        &        zetax /2.0e0_DP)) / (1.0e0_DP - zetax)**4.0e0_DP
        temp_as2u  = -d_zetax / 2.0e0_DP * (1.0e0_DP - zetax) + 3.0e0_DP * d_zetax * (1.0e0_DP - zetax /2.0e0_DP)
        temp_as2up = -dd_zetaxn/2.0e0_DP * (1.0e0_DP - zetax) + d_zetax * dzetaxn / 2.0e0_DP  + &
        &   3.0e0_DP*dd_zetaxn * (1.0e0_DP - zetax/2.0e0_DP) - 1.5e0_DP*d_zetax * dzetaxn
        temp_as2v  = (1.0e0_DP - zetax)**4.0e0_DP
        temp_as2vp = -4.0e0_DP * dzetaxn*(1.0e0_DP - zetax)**3.0e0_DP     
        temp_3p = (temp_as2v*temp_as2up - temp_as2u*temp_as2vp) / temp_as2v**2.0e0_DP 
        
        temp_4 = ((9.0e0_DP * d_zetax          *   &
        &        (1.0e0_DP + zetax) + 9.0e0_DP * zetax * d_zetax) * 2.0e0_DP * (1.0e0_DP - zetax)           +   &
        &        6.0e0_DP * d_zetax * 9.0e0_DP * zetax * (1.0e0_DP + zetax)) / (4.0e0_DP * (1.0e0_DP        -   &
        &        zetax)**4.0e0_DP)
        temp_as2u  = (9.0e0_DP * d_zetax * (1.0e0_DP + zetax) + 9.0e0_DP * zetax * d_zetax) * 2.0e0_DP      *   &
        &   (1.0e0_DP - zetax) + 6.0e0_DP * d_zetax * 9.0e0_DP * zetax * (1.0e0_DP + zetax)
        
        temp_as2_a  = 9.0e0_DP * d_zetax * (1.0e0_DP + zetax) + 9.0e0_DP * zetax * d_zetax
        temp_as2_ap = 9.0e0_DP*dd_zetaxn * (1.0e0_DP + zetax) +  9.0e0_DP * d_zetax * dzetaxn + &
        &  9.0e0_DP * dzetaxn * d_zetax + 9.0e0_DP * zetax * dd_zetaxn
        temp_as2_b  = 2.0e0_DP * (1.0e0_DP - zetax)
        temp_as2_bp = -2.0e0_DP*dzetaxn
        temp_as2_cp = 54.0e0_DP*(dd_zetaxn * zetax * (1.0e0_DP + zetax) + d_zetax * dzetaxn *   &
        &   (1.0e0_DP + zetax) + d_zetax * zetax * dzetaxn)
        
        temp_as2up = temp_as2_ap*temp_as2_b + temp_as2_a*temp_as2_bp + temp_as2_cp
        
        temp_as2v  = 4.0e0_DP * (1.0e0_DP - zetax)**4.0e0_DP
        temp_as2vp = -16.0e0_DP * dzetaxn * (1.0e0_DP - zetax)**3.0e0_DP
        temp_4p = (temp_as2v*temp_as2up - temp_as2u*temp_as2vp) / temp_as2v**2.0e0_DP

        dl_dBaan = l_const2 * (temp_1p*l_Iaa + temp_2p*l_Jaa + l_rhos * (temp_3p * l_Iaa - temp_4p * l_Jaa) + &
        &           drhosn * (temp_3 * l_Iaa - temp_4 * l_Jaa))
        dl_dBarn = l_const2 * (temp_1p*l_Iar + temp_2p*l_Jar + l_rhos * (temp_3p * l_Iar - temp_4p * l_Jar) + &
        &           drhosn * (temp_3 * l_Iar - temp_4 * l_Jar))
        dl_dBrrn = l_const2 * (temp_1p*l_Irr + temp_2p*l_Jrr + l_rhos * (temp_3p * l_Irr - temp_4p * l_Jrr) + &
        &           drhosn * (temp_3 * l_Irr - temp_4 * l_Jrr))     

        diff_v  = xx**(2.0e0_DP * lla) * (l_aa + l_Baa) - 2.0e0_DP * xx**(lla + llr) * (l_ar + l_Bar) + &
        &       xx**(2.0e0_DP * llr) * (l_rr + l_Brr) 
        diff_dv = xx**(2.0e0_DP * lla) * (l_daa + l_dBaa) - 2.0e0_DP * xx**(lla + llr) * (l_dar + l_dBar) + &
        &       xx**(2.0e0_DP * llr) * (l_drr + l_dBrr)
        
        ddiff_vn  = xx**(2.0e0_DP * lla) * (dl_aan + dl_Baan) - 2.0e0_DP * xx**(lla + llr) * (dl_arn + dl_Barn) + &
        &       xx**(2.0e0_DP * llr) * (dl_rrn + dl_Brrn) 
        ddiff_dvn = xx**(2.0e0_DP * lla) * (dl_daan + dl_dBaan) - 2.0e0_DP*xx**(lla + llr) * (dl_darn + dl_dBarn) + &
        &       xx**(2.0e0_DP * llr) * (dl_drrn + dl_dBrrn)
   
        as2grad_out = 0.5e0_DP * leps * cc**2.0e0_DP *(ddiff_vn * dkhs + diff_v * ddkhsn  + &
        &   ddiff_dvn * khs + diff_dv * dkhsn)
   
        return
    end function das2gradn
!***************************************************************************************************
    function das2gradv( l_rhos, cc, xx, dd, lla, llr, leps ) result(as2grad_out)
        implicit none
        
        real(kind=DP), intent(in) ::  l_rhos,cc,xx,dd,lla,llr,leps
        
        real(kind=DP)       ::  dkhs, d_zetax, l_denom
        real(kind=DP)       ::  l_constaa, l_zaap, l_constar, l_zarp, l_constrr, l_zrrp, l_const2
        real(kind=DP)       ::  l_zaa, l_zar, l_zrr
        real(kind=DP)       ::  l_aa, l_ar, l_rr, l_daa, l_dar, l_drr, l_Baa, l_Bar, l_Brr, l_dBaa, l_dBar, l_dBrr
        real(kind=DP)       ::  diff_v, diff_dv
        real(kind=DP)       ::  l_Iaa, l_Iar, l_Irr, l_Jaa, l_Jar, l_Jrr
        real(kind=DP)       ::  as2grad_out
        !differentials
        real(kind=DP)       ::  dd_zetaxv, dl_denomv, ddkhsv
        real(kind=DP)       ::  temp_1,temp_2,temp_3,temp_4,temp_5,temp_6,temp_7
        real(kind=DP)       ::  temp_1p,temp_2p,temp_3p,temp_4p,temp_5p,temp_6p,temp_7p
        real(kind=DP)       ::  dl_aav, dl_arv, dl_rrv
        real(kind=DP)       ::  dl_zaapv, dl_zarpv, dl_zrrpv
        real(kind=DP)       ::  dl_zaav, dl_zarv, dl_zrrv
        real(kind=DP)       ::  dl_daav, dl_darv, dl_drrv
        real(kind=DP)       ::  dl_Baav, dl_Barv, dl_Brrv, dl_dBaav, dl_dBarv, dl_dBrrv
        real(kind=DP)       ::  ddiff_dvv, ddiff_vv
        real(kind=DP)       ::  temp_as2u, temp_as2v, temp_as2up, temp_as2vp
        real(kind=DP)       ::  temp_as2_a, temp_as2_ap, temp_as2_b, temp_as2_bp, temp_as2_cp
        
        d_zetax = zetax / l_rhos
        dd_zetaxv = (drhosv*zetax - dzetaxv*l_rhos) / l_rhos**2.0e0_DP
      
        l_denom = 1.0e0_DP + 4.0e0_DP * zetax + 4.0e0_DP * zetax**2.0e0_DP - 4.0e0_DP * &
        &   zetax**3.0e0_DP + zetax**4.0e0_DP
        
        dl_denomv = 4.0e0_DP * dzetaxv + 8.0e0_DP * dzetaxv * zetax - 12.0e0_DP * &
        &   dzetaxv * zetax**2.0e0_DP + 4.0e0_DP * dzetaxv * zetax**3.0e0_DP
      
        dkhs = -4.0e0_DP * d_zetax * (1.0e0_DP - zetax)**3.0e0_DP   *                    &
        &   (l_denom + (1.0e0_DP + 2.0e0_DP * zetax - 3.0e0_DP * zetax**2.0e0_DP +       &
        &   zetax**3.0e0_DP) * (1.0e0_DP - zetax)) / l_denom**2.0e0_DP

        temp_1  = -4.0e0_DP * d_zetax
        temp_1p =  -4.0e0_DP * dd_zetaxv
        
        temp_2  = (1.0e0_DP - zetax)**3.0e0_DP
        temp_2p = -3.0e0_DP * dzetaxv * (1.0e0_DP - zetax)**2.0e0_DP
        
        temp_3  = l_denom + (1.0e0_DP + 2.0e0_DP * zetax - 3.0e0_DP * zetax**2.0e0_DP      + &
        &   zetax**3.0e0_DP) * (1.0e0_DP - zetax)
        temp_3p = dl_denomv + (1.0e0_DP + 2.0e0_DP * zetax - 3.0e0_DP * zetax**2.0e0_DP    + &
        &   zetax**3.0e0_DP) * (-dzetaxv)                                                  + &
        &   (2.0e0_DP * dzetaxv - 6.0e0_DP * dzetaxv * zetax                               + &
        &   3.0e0_DP * dzetaxv * zetax**2.0e0_DP) * (1.0e0_DP - zetax)
        
        temp_4  = l_denom**2.0e0_DP
        temp_4p = 2.0e0_DP * dl_denomv * l_denom
        
        ddkhsv = temp_1p*temp_2*temp_3/temp_4 + temp_1*temp_2p*temp_3/temp_4              + &
        &   temp_1*temp_2 * (temp_4*temp_3p - temp_4p*temp_3)/temp_4**2.0e0_DP

        l_aa    =   as1chain( 2.0e0_DP * lla, dd, leps )
        l_ar    =   as1chain( llr + lla, dd, leps )
        l_rr    =   as1chain( 2.0e0_DP * llr, dd, leps )
        
        l_constaa = PI * leps * dd**3.0e0_DP / (2.0e0_DP * lla - 3.0e0_DP)
        l_zaap    = zeffprime( 2.0e0_DP * lla )
        l_constar = PI * leps * dd**3.0e0_DP / (llr + lla - 3.0e0_DP)
        l_zarp    = zeffprime( llr + lla )
        l_constrr = PI * leps * dd**3.0e0_DP / (2.0e0_DP * llr - 3.0e0_DP)
        l_zrrp    = zeffprime( 2.0e0_DP * llr )
        l_zaa     = zeff( 2.0e0_DP * lla )
        l_zar     = zeff( llr + lla )
        l_zrr     = zeff( 2.0e0_DP * llr )

        dl_aav    =   das1chainv( 2.0e0_DP * lla, dd, leps )
        dl_arv    =   das1chainv( llr + lla, dd, leps )
        dl_rrv    =   das1chainv( 2.0e0_DP * llr, dd, leps )
   
        dl_zaapv    = dzeffprimev( 2.0e0_DP * lla )
        dl_zarpv    = dzeffprimev( llr + lla )
        dl_zrrpv    = dzeffprimev( 2.0e0_DP * llr )
        
        dl_zaav     = dzeffv( 2.0e0_DP * lla )
        dl_zarv     = dzeffv( llr + lla )
        dl_zrrv     = dzeffv( 2.0e0_DP * llr )

        l_daa   =   -2.0e0_DP * l_constaa * ( (1.0e0_DP - l_zaa / 2.0e0_DP) / (1.0e0_DP - l_zaa)**3.0e0_DP +  &
        &    l_rhos * (2.5e0_DP * l_zaap - l_zaap * l_zaa) / (1.0e0_DP - l_zaa)**4.0e0_DP)
        l_dar   =   -2.0e0_DP * l_constar * ( (1.0e0_DP - l_zar / 2.0e0_DP) / (1.0e0_DP - l_zar)**3.0e0_DP +  &
        &    l_rhos * (2.5e0_DP * l_zarp - l_zarp * l_zar) / (1.0e0_DP - l_zar)**4.0e0_DP)
        l_drr   =   -2.0e0_DP * l_constrr * ( (1.0e0_DP - l_zrr / 2.0e0_DP) / (1.0e0_DP - l_zrr)**3.0e0_DP +  &
        &    l_rhos * (2.5e0_DP * l_zrrp - l_zrrp * l_zrr) / (1.0e0_DP - l_zrr)**4.0e0_DP)
        
        temp_1  = 1.0e0_DP - l_zaa / 2.0e0_DP
        temp_1p = -dl_zaav / 2.0e0_DP
        
        temp_2  = (1.0e0_DP - l_zaa)**3.0e0_DP
        temp_2p = -3.0e0_DP * dl_zaav * (1.0e0_DP - l_zaa)**2.0e0_DP
        
        temp_3  = l_rhos * (2.5e0_DP * l_zaap - l_zaap * l_zaa)
        temp_3p = drhosv * (2.5e0_DP * l_zaap - l_zaap * l_zaa) + l_rhos * (2.5e0_DP * dl_zaapv - (dl_zaapv     * &
        &    l_zaa + l_zaap * dl_zaav))
        
        temp_4  = (1.0e0_DP - l_zaa)**4.0e0_DP
        temp_4p = -4.0e0_DP * dl_zaav * (1.0e0_DP - l_zaa)**3.0e0_DP
        
        dl_daav = -2.0e0_DP * l_constaa * ( (temp_2*temp_1p - temp_2p*temp_1)/temp_2**2.0e0_DP                  + &
        &   (temp_4*temp_3p - temp_4p*temp_3)/temp_4**2.0e0_DP)
             
        temp_1  = 1.0e0_DP - l_zar / 2.0e0_DP
        temp_1p = -dl_zarv / 2.0e0_DP
        
        temp_2  = (1.0e0_DP - l_zar)**3.0e0_DP
        temp_2p = -3.0e0_DP * dl_zarv * (1.0e0_DP - l_zar)**2.0e0_DP
        
        temp_3  = l_rhos * (2.5e0_DP * l_zarp - l_zarp * l_zar)
        temp_3p = drhosv * (2.5e0_DP * l_zarp - l_zarp * l_zar) + l_rhos * (2.5e0_DP * dl_zarpv - (dl_zarpv     * &
        &    l_zar + l_zarp * dl_zarv))
        
        temp_4  = (1.0e0_DP - l_zar)**4.0e0_DP
        temp_4p = -4.0e0_DP * dl_zarv * (1.0e0_DP - l_zar)**3.0e0_DP
        
        dl_darv = -2.0e0_DP * l_constar * ( (temp_2*temp_1p - temp_2p*temp_1)/temp_2**2.0e0_DP                  + &
        &   (temp_4*temp_3p - temp_4p*temp_3)/temp_4**2.0e0_DP)        
          
        temp_1  = 1.0e0_DP - l_zrr / 2.0e0_DP
        temp_1p = -dl_zrrv / 2.0e0_DP
        
        temp_2  = (1.0e0_DP - l_zrr)**3.0e0_DP
        temp_2p = -3.0e0_DP * dl_zrrv * (1.0e0_DP - l_zrr)**2.0e0_DP
        
        temp_3  = l_rhos * (2.5e0_DP * l_zrrp - l_zrrp * l_zrr)
        temp_3p = drhosv * (2.5e0_DP * l_zrrp - l_zrrp * l_zrr) + l_rhos * (2.5e0_DP * dl_zrrpv - (dl_zrrpv     * &
        &    l_zrr + l_zrrp * dl_zrrv))
        
        temp_4  = (1.0e0_DP - l_zrr)**4.0e0_DP
        temp_4p = -4.0e0_DP * dl_zrrv * (1.0e0_DP - l_zrr)**3.0e0_DP
        
        dl_drrv = -2.0e0_DP * l_constrr * ( (temp_2*temp_1p - temp_2p*temp_1)/temp_2**2.0e0_DP                  + &
        &   (temp_4*temp_3p - temp_4p*temp_3)/temp_4**2.0e0_DP) 

        l_Baa   =   bchain(2.0e0_DP * lla, xx,dd, leps)
        l_Bar   =   bchain(llr + lla, xx, dd, leps)
        l_Brr   =   bchain(2.0e0_DP * llr, xx, dd, leps)
        
        dl_Baav   =   dbchainv(2.0e0_DP * lla, xx,dd, leps)
        dl_Barv   =   dbchainv(llr + lla, xx, dd, leps)
        dl_Brrv   =   dbchainv(2.0e0_DP * llr, xx, dd, leps)

        l_Iaa = -(xx**(3.0e0_DP - 2.0e0_DP * lla) - 1.0e0_DP) / (2.0e0_DP * lla - 3.0e0_DP)  
        l_Iar = -(xx**(3.0e0_DP - llr - lla) - 1.0e0_DP) / (llr + lla - 3.0e0_DP) 
        l_Irr = -(xx**(3.0e0_DP - 2.0e0_DP * llr) - 1.0e0_DP) / (2.0e0_DP * llr - 3.0e0_DP)
        
        l_Jaa = -(xx**(4.0e0_DP - 2.0e0_DP * lla) * (2.0e0_DP * lla - 3.0e0_DP) - xx**(3.0e0_DP - 2.0e0_DP * lla) &
        &   * (2.0e0_DP * lla - 4.0e0_DP) - 1.0e0_DP) / ((2.0e0_DP * lla - 3.0e0_DP) * (2.0e0_DP * lla - 4.0e0_DP))
        l_Jar = -(xx**(4.0e0_DP - lla - llr) * (lla + llr - 3.0e0_DP) - xx**(3.0e0_DP - lla - llr) * (lla + llr - &
        &   4.0e0_DP) - 1.0e0_DP) / ((lla + llr - 3.0e0_DP) * (lla + llr - 4.0e0_DP))
        l_Jrr = -(xx**(4.0e0_DP - 2.0e0_DP * llr) * (2.0e0_DP * llr - 3.0e0_DP) - xx**(3.0e0_DP - 2.0e0_DP * llr) * &
        &   (2.0e0_DP * llr - 4.0e0_DP) - 1.0e0_DP) / ((2.0e0_DP * llr - 3.0e0_DP) * (2.0e0_DP * llr - 4.0e0_DP))   
        
        l_const2 = 2.0e0_DP * PI * dd**3.0e0_DP * leps 

        l_dBaa  =   l_const2 * ( (1.0e0_DP - zetax / 2.0e0_DP) / (1.0e0_DP - zetax)**3.0e0_DP * l_Iaa       -   &
        &        9.0e0_DP * zetax * (1.0e0_DP + zetax) / (2.0e0_DP * (1.0e0_DP - zetax)**3.0e0_DP) * l_Jaa  +   &
        &        l_rhos * (((-d_zetax / 2.0e0_DP * (1.0e0_DP - zetax) + 3.0e0_DP * d_zetax * (1.0e0_DP      -   &
        &        zetax /2.0e0_DP)) / (1.0e0_DP - zetax)**4.0e0_DP) * l_Iaa - (((9.0e0_DP * d_zetax          *   &
        &        (1.0e0_DP + zetax) + 9.0e0_DP * zetax * d_zetax) * 2.0e0_DP * (1.0e0_DP - zetax)           +   &
        &        6.0e0_DP * d_zetax * 9.0e0_DP * zetax * (1.0e0_DP + zetax)) / (4.0e0_DP * (1.0e0_DP        -   &
        &        zetax)**4.0e0_DP)) * l_Jaa )) 

        l_dBar  =  l_const2 * ( (1.0e0_DP - zetax / 2.0e0_DP) / (1.0e0_DP - zetax)**3.0e0_DP * l_Iar       -   &
        &        9.0e0_DP * zetax * (1.0e0_DP + zetax) / (2.0e0_DP * (1.0e0_DP - zetax)**3.0e0_DP) * l_Jar  +   &
        &        l_rhos * (((-d_zetax / 2.0e0_DP * (1.0e0_DP - zetax) + 3.0e0_DP * d_zetax * (1.0e0_DP      -   &
        &        zetax /2.0e0_DP)) / (1.0e0_DP - zetax)**4.0e0_DP) * l_Iar - (((9.0e0_DP * d_zetax          *   &
        &        (1.0e0_DP + zetax) + 9.0e0_DP * zetax * d_zetax) * 2.0e0_DP * (1.0e0_DP - zetax)           +   &
        &        6.0e0_DP * d_zetax * 9.0e0_DP * zetax * (1.0e0_DP + zetax)) / (4.0e0_DP * (1.0e0_DP        -   &
        &        zetax)**4.0e0_DP)) * l_Jar )) 
     
        l_dBrr  =  l_const2 * ( (1.0e0_DP - zetax / 2.0e0_DP) / (1.0e0_DP - zetax)**3.0e0_DP * l_Irr       -   &
        &        9.0e0_DP * zetax * (1.0e0_DP + zetax) / (2.0e0_DP * (1.0e0_DP - zetax)**3.0e0_DP) * l_Jrr  +   &
        &        l_rhos * (((-d_zetax / 2.0e0_DP * (1.0e0_DP - zetax) + 3.0e0_DP * d_zetax * (1.0e0_DP      -   &
        &        zetax /2.0e0_DP)) / (1.0e0_DP - zetax)**4.0e0_DP) * l_Irr - (((9.0e0_DP * d_zetax          *   &
        &        (1.0e0_DP + zetax) + 9.0e0_DP * zetax * d_zetax) * 2.0e0_DP * (1.0e0_DP - zetax)           +   &
        &        6.0e0_DP * d_zetax * 9.0e0_DP * zetax * (1.0e0_DP + zetax)) / (4.0e0_DP * (1.0e0_DP        -   &
        &        zetax)**4.0e0_DP)) * l_Jrr )) 
        
        temp_1  = (1.0e0_DP - zetax / 2.0e0_DP) / (1.0e0_DP - zetax)**3.0e0_DP 
        temp_1p = ((1.0e0_DP - zetax)**3.0e0_DP * (-dzetaxv / 2.0e0_DP)                         + &
        &   3.0e0_DP * dzetaxv * (1.0e0_DP - zetax / 2.0e0_DP) * (1.0e0_DP - zetax)**2.0e0_DP)  / &
        &   (1.0e0_DP - zetax)**6.0e0_DP 
        
        temp_2  = -9.0e0_DP * zetax * (1.0e0_DP + zetax) / (2.0e0_DP * (1.0e0_DP - zetax)**3.0e0_DP) 
        temp_2p = ( (2.0e0_DP * (1.0e0_DP - zetax)**3.0e0_DP) * (-9.0e0_DP * dzetaxv * (1.0e0_DP        + &
        &   zetax) -9.0e0_DP * zetax * dzetaxv) - (9.0e0_DP * zetax * (1.0e0_DP + zetax))               * &
        &   (6.0e0_DP * dzetaxv * (1.0e0_DP - zetax)**2.0e0_DP)) / (4.0e0_DP * (1.0e0_DP - zetax)**6.0e0_DP)
        
        temp_3 = (-d_zetax / 2.0e0_DP * (1.0e0_DP - zetax) + 3.0e0_DP * d_zetax * (1.0e0_DP      -   &
        &        zetax /2.0e0_DP)) / (1.0e0_DP - zetax)**4.0e0_DP
        temp_as2u  = -d_zetax / 2.0e0_DP * (1.0e0_DP - zetax) + 3.0e0_DP * d_zetax * (1.0e0_DP - zetax /2.0e0_DP)
        temp_as2up = -dd_zetaxv/2.0e0_DP * (1.0e0_DP - zetax) + d_zetax * dzetaxv / 2.0e0_DP  + &
        &   3.0e0_DP*dd_zetaxv * (1.0e0_DP - zetax/2.0e0_DP) - 1.5e0_DP*d_zetax * dzetaxv
        temp_as2v  = (1.0e0_DP - zetax)**4.0e0_DP
        temp_as2vp = -4.0e0_DP * dzetaxv*(1.0e0_DP - zetax)**3.0e0_DP     
        temp_3p = (temp_as2v*temp_as2up - temp_as2u*temp_as2vp) / temp_as2v**2.0e0_DP 

        temp_4 = ((9.0e0_DP * d_zetax          *   &
        &        (1.0e0_DP + zetax) + 9.0e0_DP * zetax * d_zetax) * 2.0e0_DP * (1.0e0_DP - zetax)           +   &
        &        6.0e0_DP * d_zetax * 9.0e0_DP * zetax * (1.0e0_DP + zetax)) / (4.0e0_DP * (1.0e0_DP        -   &
        &        zetax)**4.0e0_DP)
        temp_as2u  = (9.0e0_DP * d_zetax * (1.0e0_DP + zetax) + 9.0e0_DP * zetax * d_zetax) * 2.0e0_DP      *   &
        &   (1.0e0_DP - zetax) + 6.0e0_DP * d_zetax * 9.0e0_DP * zetax * (1.0e0_DP + zetax)

        temp_as2_a  = 9.0e0_DP * d_zetax * (1.0e0_DP + zetax) + 9.0e0_DP * zetax * d_zetax
        temp_as2_ap = 9.0e0_DP*dd_zetaxv * (1.0e0_DP + zetax) +  9.0e0_DP * d_zetax * dzetaxv + &
        &  9.0e0_DP * dzetaxv * d_zetax + 9.0e0_DP * zetax * dd_zetaxv
        temp_as2_b  = 2.0e0_DP * (1.0e0_DP - zetax)
        temp_as2_bp = -2.0e0_DP*dzetaxv
        temp_as2_cp = 54.0e0_DP*(dd_zetaxv * zetax * (1.0e0_DP + zetax) + d_zetax * dzetaxv *   &
        &   (1.0e0_DP + zetax) + d_zetax * zetax * dzetaxv)

        temp_as2up = temp_as2_ap*temp_as2_b + temp_as2_a*temp_as2_bp + temp_as2_cp
        
        temp_as2v  = 4.0e0_DP * (1.0e0_DP - zetax)**4.0e0_DP
        temp_as2vp = -16.0e0_DP * dzetaxv * (1.0e0_DP - zetax)**3.0e0_DP
        temp_4p = (temp_as2v*temp_as2up - temp_as2u*temp_as2vp) / temp_as2v**2.0e0_DP
        
        dl_dBaav = l_const2 * (temp_1p*l_Iaa + temp_2p*l_Jaa + l_rhos * (temp_3p * l_Iaa - temp_4p * l_Jaa) + &
        &           drhosv * (temp_3 * l_Iaa - temp_4 * l_Jaa))
        dl_dBarv = l_const2 * (temp_1p*l_Iar + temp_2p*l_Jar + l_rhos * (temp_3p * l_Iar - temp_4p * l_Jar) + &
        &           drhosv * (temp_3 * l_Iar - temp_4 * l_Jar))
        dl_dBrrv = l_const2 * (temp_1p*l_Irr + temp_2p*l_Jrr + l_rhos * (temp_3p * l_Irr - temp_4p * l_Jrr) + &
        &           drhosv * (temp_3 * l_Irr - temp_4 * l_Jrr))     

        diff_v  = xx**(2.0e0_DP * lla) * (l_aa + l_Baa) - 2.0e0_DP * xx**(lla + llr) * (l_ar + l_Bar) + &
        &       xx**(2.0e0_DP * llr) * (l_rr + l_Brr) 
        diff_dv = xx**(2.0e0_DP * lla) * (l_daa + l_dBaa) - 2.0e0_DP * xx**(lla + llr) * (l_dar + l_dBar) + &
        &       xx**(2.0e0_DP * llr) * (l_drr + l_dBrr)
        
        ddiff_vv  = xx**(2.0e0_DP * lla) * (dl_aav + dl_Baav) - 2.0e0_DP * xx**(lla + llr) * (dl_arv + dl_Barv) + &
        &       xx**(2.0e0_DP * llr) * (dl_rrv + dl_Brrv) 
        ddiff_dvv = xx**(2.0e0_DP * lla) * (dl_daav + dl_dBaav) - 2.0e0_DP*xx**(lla + llr) * (dl_darv + dl_dBarv) + &
        &       xx**(2.0e0_DP * llr) * (dl_drrv + dl_dBrrv)
   
        as2grad_out = 0.5e0_DP * leps * cc**2.0e0_DP *(ddiff_vv * dkhs + diff_v * ddkhsv  + &
        &   ddiff_dvv * khs + diff_dv * dkhsv)

        return
    end function das2gradv
!***************************************************************************************************
end module Chain
