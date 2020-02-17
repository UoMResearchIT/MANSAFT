!***************************************************************************************************
!***************************************************************************************************
!   SAFT Module to:
!       1. Calculate all zeff terms
!***************************************************************************************************
!
!***************************************************************************************************
module Z_eff
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
    pure function zeff( lam_in ) result(zeff_out)
        implicit none
        
        integer                         ::  i_zeff
        real(kind=DP),intent(in)        ::  lam_in
        real(kind=DP)                   ::  zeff_out
        real(kind=DP)                   ::  ca(1:4)
        
        do i_zeff=1, 4
            ca(i_zeff) = c_array(i_zeff, 1) * 1.0e0_DP + c_array(i_zeff, 2) * 1.0e0_DP /  &
            &   lam_in + c_array(i_zeff, 3) * 1.0e0_DP / lam_in**2.0e0_DP                 &
            &   + c_array(i_zeff, 4) * 1.0e0_DP / lam_in**3.0e0_DP
        end do        

        zeff_out = ca(1) * zetax + ca(2) * zetax**2.0e0_DP + ca(3) * zetax**3.0e0_DP +  &
        &   ca(4) * zetax**4.0e0_DP      

        return
    end function zeff  
!***************************************************************************************************    
    pure function dzeffn( lam_in ) result(zeff_out)
        implicit none
        
        integer                         ::  i_zeff
        real(kind=DP),intent(in)        ::  lam_in
        real(kind=DP)                   ::  zeff_out
        real(kind=DP)                   ::  ca(1:4)
        
        do i_zeff=1, 4
            ca(i_zeff) = c_array(i_zeff, 1) * 1.0e0_DP + c_array(i_zeff, 2) * 1.0e0_DP /  &
            &   lam_in + c_array(i_zeff, 3) * 1.0e0_DP / lam_in**2.0e0_DP                 &
            &   + c_array(i_zeff, 4) * 1.0e0_DP / lam_in**3.0e0_DP
        end do        

        zeff_out = ca(1) * dzetaxn + ca(2) * 2.0e0_DP*dzetaxn*zetax + ca(3) * 3.0e0_DP*dzetaxn*zetax**2.0e0_DP +  &
        &   ca(4) * 4.0e0_DP*dzetaxn*zetax**3.0e0_DP      

        return
    end function dzeffn
!*************************************************************************************************** 
    pure function dzeffv( lam_in ) result(zeff_out)
        implicit none
        
        integer                         ::  i_zeff
        real(kind=DP),intent(in)        ::  lam_in
        real(kind=DP)                   ::  zeff_out
        real(kind=DP)                   ::  ca(1:4)
        
        do i_zeff=1, 4
            ca(i_zeff) = c_array(i_zeff, 1) * 1.0e0_DP + c_array(i_zeff, 2) * 1.0e0_DP /  &
            &   lam_in + c_array(i_zeff, 3) * 1.0e0_DP / lam_in**2.0e0_DP                 &
            &   + c_array(i_zeff, 4) * 1.0e0_DP / lam_in**3.0e0_DP
        end do        

        zeff_out = ca(1) * dzetaxv + ca(2) * 2.0e0_DP*dzetaxv*zetax + ca(3) * 3.0e0_DP*dzetaxv*zetax**2.0e0_DP +  &
        &   ca(4) * 4.0e0_DP*dzetaxv*zetax**3.0e0_DP       

        return
    end function dzeffv  
!***************************************************************************************************
    pure function zeffprime( lamin ) result(zeffprime_out)
        implicit none
        
        integer                     ::  i_zeffp
        real(kind=DP), intent(in)   ::  lamin
        real(kind=DP)               ::  ca(1:4), l_zetax
        real(kind=DP)               ::  zeffprime_out
        
        l_zetax = zetax / rhos
      
        do i_zeffp=1, 4
            ca(i_zeffp) = c_array(i_zeffp, 1) * 1.0e0_DP + c_array(i_zeffp, 2) * 1.0e0_DP / lamin +     &
            &   c_array(i_zeffp, 3) * 1.0e0_DP / lamin**2.0e0_DP                                        &
            &   + c_array(i_zeffp, 4) * 1.0e0_DP / lamin**3.0e0_DP
        end do        

        zeffprime_out = ca(1) * l_zetax + 2.0e0_DP * zetax * ca(2) * l_zetax           +       &
        &   3.0e0_DP * zetax**2.0e0_DP * ca(3) * l_zetax + 4.0e0_DP * zetax**3.0e0_DP *       &
        &   ca(4) * l_zetax    

        return
    end function zeffprime  
!***************************************************************************************************
    pure function dzeffprimen( lamin ) result(zeffprime_out)
        implicit none
        
        integer                     ::  i_zeffp
        real(kind=DP), intent(in)   ::  lamin
        real(kind=DP)               ::  ca(1:4), l_zetax, dl_zetaxn
        real(kind=DP)               ::  zeffprime_out
        
        l_zetax   = zetax / rhos
        dl_zetaxn = (rhos * dzetaxn - drhosn * zetax) / rhos**2.0e0_DP
      
        do i_zeffp=1, 4
            ca(i_zeffp) = c_array(i_zeffp, 1) * 1.0e0_DP + c_array(i_zeffp, 2) * 1.0e0_DP / lamin +     &
            &   c_array(i_zeffp, 3) * 1.0e0_DP / lamin**2.0e0_DP                                        &
            &   + c_array(i_zeffp, 4) * 1.0e0_DP / lamin**3.0e0_DP
        end do        

        zeffprime_out = ca(1) * dl_zetaxn + 2.0e0_DP * ca(2) * (zetax*dl_zetaxn + dzetaxn*l_zetax)          + &
        &   3.0e0_DP * ca(3) * (zetax**2.0e0_DP*dl_zetaxn + 2.0e0_DP*zetax*dzetaxn*l_zetax)                 + &
        &   4.0e0_DP * ca(4) * (zetax**3.0e0_DP*dl_zetaxn + 3.0e0_DP*zetax**2.0e0_DP*dzetaxn*l_zetax)  

        return
    end function dzeffprimen 
!***************************************************************************************************
    pure function dzeffprimev( lamin ) result(zeffprime_out)
        implicit none
        
        integer                     ::  i_zeffp
        real(kind=DP), intent(in)   ::  lamin
        real(kind=DP)               ::  ca(1:4), l_zetax, dl_zetaxv
        real(kind=DP)               ::  zeffprime_out
        
        l_zetax   = zetax / rhos
        dl_zetaxv = (rhos * dzetaxv - drhosv * zetax) / rhos**2.0e0_DP
      
        do i_zeffp=1, 4
            ca(i_zeffp) = c_array(i_zeffp, 1) * 1.0e0_DP + c_array(i_zeffp, 2) * 1.0e0_DP / lamin +     &
            &   c_array(i_zeffp, 3) * 1.0e0_DP / lamin**2.0e0_DP                                        &
            &   + c_array(i_zeffp, 4) * 1.0e0_DP / lamin**3.0e0_DP
        end do        

        zeffprime_out = ca(1) * dl_zetaxv + 2.0e0_DP * ca(2) * (zetax*dl_zetaxv + dzetaxv*l_zetax)          + &
        &   3.0e0_DP * ca(3) * (zetax**2.0e0_DP*dl_zetaxv + 2.0e0_DP*zetax*dzetaxv*l_zetax)                 + &
        &   4.0e0_DP * ca(4) * (zetax**3.0e0_DP*dl_zetaxv + 3.0e0_DP*zetax**2.0e0_DP*dzetaxv*l_zetax)  

        return
    end function dzeffprimev 
!***************************************************************************************************
!***************************************************************************************************
end module Z_eff
