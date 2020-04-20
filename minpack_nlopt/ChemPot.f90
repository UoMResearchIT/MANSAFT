!***************************************************************************************************
!   SAFT Module to:
!       1. Calculate chemical potential
!
!       Updated March 2017 to use analytical integrals
!
!***************************************************************************************************
!
!***************************************************************************************************
module ChemPot
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
    use Diff            ! Differentials module  
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
end module ChemPot
!***************************************************************************************************
