!***************************************************************************************************
!***************************************************************************************************
!   SAFT Module to:
!       1. Calculate A ideal
!***************************************************************************************************
!***************************************************************************************************
module Ideal
!***************************************************************************************************
!Modules
!=======    
    use Types           ! Definitions of types and double precision
    use Global          ! Important global parameters 
!***************************************************************************************************
    implicit none
!***************************************************************************************************
    public  ::  A_ideal, A_ideal_dv, A_ideal_dn
!***************************************************************************************************
    contains
!***************************************************************************************************
    function A_ideal() result(a_res)     ![1]  A1      [2] 5
        implicit none
        
        integer         ::  i_a
        real(kind=DP)   ::  a_res
       
        a_res = -1.0e0_DP
        
        do i_a=1,nctypes  
            if(Comp_array(i_a)%xi/=0.0e0_DP) a_res = a_res + Comp_array(i_a)%xi * &
            & dlog(Comp_array(i_a)%rho * NA * Comp_array(i_a)%db3)   !NOTE [1] A1 misprint - xi outside LN                       
        end do

        return
    end function A_ideal
!***************************************************************************************************   
    function A_ideal_dn( i_part ) result(a_res)
        implicit none
        
        integer                 ::  i_a
        integer, intent(in)     ::  i_part
        real(kind=DP)           ::  a_res

        a_res = 0.0e0_DP

        do i_a=1,nctypes
            if(Comp_array(i_a)%xi/=0.0e0_DP) a_res = a_res + Comp_array(i_a)%dxin*dlog(Comp_array(i_a)%nm/v * Comp_array(i_a)%db3)     
        end do
        
        a_res = nsum *( a_res +  Comp_array(i_part)%xi * 1.0e0_DP/Comp_array(i_part)%nm ) +  A_ideal()   
 
        return
    end function A_ideal_dn       
!***************************************************************************************************         
    function A_ideal_dv() result(a_res)     ![1]  A1      [2] 5
        implicit none
        
        integer         ::  i_a
        real(kind=DP)   ::  a_res
       
        a_res = 0.0e0_DP
        
        do i_a=1,nctypes  
            if(Comp_array(i_a)%xi/=0.0e0_DP) a_res = a_res + Comp_array(i_a)%xi / v                                 
        end do    

        return
    end function A_ideal_dv
!***************************************************************************************************  
end module Ideal
