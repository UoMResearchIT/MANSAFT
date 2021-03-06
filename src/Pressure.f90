!***************************************************************************************************
!   SAFT Module to:
!       1. Calculate pressure
!
!       Updated March 2017 to use analytical derivatives
!
!***************************************************************************************************
!
!***************************************************************************************************
module Pressure
!***************************************************************************************************
!Modules
!=======    
    use Types           ! Definitions of types and double precision
    use Global          ! Important global parameters
    use Setup           ! To setup the system
    use Ideal           ! Calculate ideal A
    use Mono            ! Calculate mono A
    use Chain           ! Calculate chain A
    use Assoc           ! Calculate assoc A
    use Ion             ! Calculate ion A   
    use Diff            ! Differentials setup
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
end module Pressure
!***************************************************************************************************
