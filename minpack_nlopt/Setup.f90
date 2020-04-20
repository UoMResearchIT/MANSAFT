!***************************************************************************************************
!   SAFT Module to:
!       1. SETUP COMMON DATA
!***************************************************************************************************
!   Created by SJ Halstead Nov 2015
!       Update Jan 2016
!           -   Making more efficient / functions
!       Update Mar 2016
!           -   Removed units nm -> m for clarity
!           -   Made all arrays to deallocate if previously allocated
!       Update June 2016
!           -   Tidied and efficiency approved
!***************************************************************************************************
!
!***************************************************************************************************
module Setup
!***************************************************************************************************
!MODULES
use Global     !CONTAINS GLOBAL VARIABLES
use Types
!***************************************************************************************************
    implicit none
    
    integer             ::  i_set1, i_set2, i_set3, i_set4
!***************************************************************************************************
    public  ::  Sys_setup
    private ::  Debroglie_set, Nm_set, Rho_set, Xsi_set, Cij_set,    &
    &   Dij_set, Zetal_set, Zetax_set, Alpha_set, F_set, Chain_set
!***************************************************************************************************
contains
!***************************************************************************************************
    subroutine Sys_setup( is_mu_in )
        implicit none
        
        logical, optional, intent(in)   :: is_mu_in
        logical                         :: is_mu
        
        is_mu = .false.
        if(present(is_mu_in)) is_mu = is_mu_in
        
        call Debroglie_set()         !thermal debroglie wavelengths and beta     
        if(is_mu) then
            call Xi_set()            !determine xi of each molecule
        else
            call Nm_set()                !determine nm of each molecule       
        end if
        call Rho_set()               !Number density and spherical segment density         
        call Xsi_set()               !mole fractions of components
        call Cij_set()               !Cij values  and on
        call Dij_set()               !Dij effective diameters 
        call Zetal_set               !Moments of number density        
        call Zetax_set()             !Packing fractions   
        call Alpha_set()             !for mono
        call F_set()                 !for mono
        call Chain_set()             !for chains
        
        return
    end subroutine Sys_setup
!***************************************************************************************************
    subroutine Debroglie_set() !DEGROGLIE WAVELENGTH
        implicit none
        
        real(kind=DP)       ::  db_local
        
        !Set beta values, depending on T (2nd is because energy values in units of K)
        beta    = 1/(KB * t)
        beta_K  = 1/t
        
        do i_set1=1,nctypes
            db_local = H / dsqrt( TWOPI * Comp_array(i_set1)%m * KB * t)    
            Comp_array(i_set1)%db3 = db_local**3.e0_DP                     
        end do

        return
    end subroutine Debroglie_set
!*************************************************************************************************** 
    subroutine  Xi_set() !Calculates mole fraction - for Mu calc only
        implicit none   
        
        real(kind=DP)       ::  nm_tot

        nm_tot = 0.e0_DP
                       
        do i_set1=1, nctypes
            nm_tot = nm_tot + Comp_array(i_set1)%nm
        end do
        
        do i_set1=1, nctypes
            Comp_array(i_set1)%xi = Comp_array(i_set1)%nm / nm_tot
        end do
 
        nsum=nm_tot
    
        return
    end subroutine Xi_set
!***************************************************************************************************  
    subroutine  Nm_set() !Equates number of molecules to mole fraction
        implicit none   
                       
        do i_set1=1,nctypes    
            Comp_array(i_set1)%nm = Comp_array(i_set1)%xi * NA
        end do
        
        return
    end subroutine Nm_set
!*************************************************************************************************** 
    subroutine  Rho_set()     
        implicit none

        real(kind=DP)    ::  trho, tsum
      
        tsum = sum(Comp_array(:)%nm)

        !Density moles per m^3
        rho = 0.0e0_DP       
   
        do i_set1=1, nctypes
            trho = Comp_array(i_set1)%nm / v / NA

            rho = rho + trho
            Comp_array(i_set1)%rho = trho 
        end do
  
        !Density of spherical segments moles per m^3
        rhos = 0.e0_DP 

        do i_set1=1,nctypes
        do i_set2=1,nstypes
            rhos = rhos + Comp_array(i_set1)%xi * Comp_array(i_set1)%comp(i_set2) &  
            &   * Seg_array(i_set2)%sf * Seg_array(i_set2)%nseg          
        end do
        end do

        rhos = rhos * rho   !moles per m^3

        return
    end subroutine Rho_set
!***************************************************************************************************
    subroutine Xsi_set()     ![1] A8  
        implicit none
        
        real(kind=DP)    ::  xsi_tot, xsi_sum
        
        xsi_tot=0.0e0_DP
        
        do i_set1=1, nctypes
        do i_set2=1, nstypes
            xsi_tot = xsi_tot + Comp_array(i_set1)%comp(i_set2) * Seg_array(i_set2)%nseg *Comp_array(i_set1)%xi &
            &   * Seg_array(i_set2)%sf
        end do
        end do
        
        do i_set1=1, nstypes
        
            xsi_sum=0.0e0_DP
            
            do i_set2=1, nctypes
                xsi_sum = xsi_sum +  (Comp_array(i_set2)%comp(i_set1) * Seg_array(i_set1)%sf *      &  !Including shape factor!
                &                  Seg_array(i_set1)%nseg * Comp_array(i_set2)%xi) / xsi_tot
            end do
            
        !mole fraction for each segments
            Seg_array(i_set1)%xsi = xsi_sum           
        end do   
  
        return
    end subroutine Xsi_set
!***************************************************************************************************
    subroutine Cij_set() ![1] A3 Mie correction scaling factor
        implicit none
        
        if(allocated(cij)) deallocate(cij)
        allocate(cij(1:nstypes, 1:nstypes))
        
        do i_set1=1, nstypes
        do i_set2=1, nstypes
            cij(i_set1, i_set2) = ( lr(i_set1, i_set2) / (lr(i_set1, i_set2) - la(i_set1, i_set2)) ) * (lr(i_set1, i_set2) /   &
            &   la(i_set1, i_set2)) ** (la(i_set1, i_set2) / (lr(i_set1, i_set2) - la(i_set1, i_set2)))
        end do
        end do
      
        return 
    end subroutine Cij_set    
!***************************************************************************************************
    subroutine Dij_set()     ![1] A9 effective diameter of Mie segments 
        implicit none

        real(kind=DP)           ::  lg_sum,mie_x,ba2
        real(kind=DP)           ::  t_eps,t_sig,t_cij,t_lr,t_la
        
        if(allocated(dij)) deallocate(dij)
        allocate(dij(1:nstypes, 1:nstypes))

        !Gaussian-Legendre sum
        do i_set1=1, nstypes
            t_eps = eps(i_set1, i_set1)
            t_sig = sig(i_set1, i_set1)
            
            t_cij = cij(i_set1, i_set1)
            t_lr  = lr(i_set1, i_set1) 
            t_la  = la(i_set1, i_set1)
            
            ba2 = t_sig / 2.0e0_DP     !G-L limit
            lg_sum = 0.0e0_DP
     
            do i_set2=1, N_GL
                  mie_x = ba2 * (X_GL(i_set2) + 1.0e0_DP)              
                  lg_sum = lg_sum + W_GL(i_set2) * (1.0e0_DP - &
                  & dexp( -beta_K * Mieii( mie_x, t_sig, t_eps, t_cij, t_lr, t_la)))
            end do

            dij(i_set1, i_set1) = ba2 * lg_sum      !In metres               
        end do
  
        do i_set1=1, nstypes
        do i_set2=1, nstypes      
            dij(i_set1, i_set2) = ( dij(i_set1, i_set1) + dij(i_set2, i_set2)) / 2.0e0_DP    !1 A46           
        end do
        end do

        return
    end subroutine Dij_set
!***************************************************************************************************
    subroutine Zetal_set()     ![1]  A7 Moments of number density
        implicit none
        
        real(kind=DP)           ::  t_dii

        zl = 0.0e0_DP
        
        do i_set1=1,nstypes
            t_dii = dij(i_set1,i_set1)     
           
            zl(0) = zl(0) + Seg_array(i_set1)%xsi
            zl(1) = zl(1) + Seg_array(i_set1)%xsi * t_dii
            zl(2) = zl(2) + Seg_array(i_set1)%xsi * t_dii**2.0e0_DP
            zl(3) = zl(3) + Seg_array(i_set1)%xsi * t_dii**3.0e0_DP          
        end do

        zl(0) = zl(0) * PI * rhos / 6.0e0_DP * NA              
        zl(1) = zl(1) * PI * rhos / 6.0e0_DP * NA              
        zl(2) = zl(2) * PI * rhos / 6.0e0_DP * NA            
        zl(3) = zl(3) * PI * rhos / 6.0e0_DP * NA            

        return
    end subroutine Zetal_set
!***************************************************************************************************
    subroutine Zetax_set()         ![1] A13 packing fractions
        implicit none

        zetax = 0.0e0_DP
        
        do i_set1=1, nstypes
        do i_set2=1, nstypes
            zetax = zetax + Seg_array(i_set1)%xsi * Seg_array(i_set2)%xsi * dij(i_set1,i_set2)**3.0e0_DP
        end do
        end do

        zetax_sum = zetax      

        zetax = zetax * PI * rhos / 6.0e0_DP * NA
      
        KHS = ((1.d0 - zetax)**4.0e0_DP) / (1.d0 + 4.0e0_DP * zetax + 4.0e0_DP * zetax**2 - 4.0e0_DP *  &
        &       zetax**3.0e0_DP + zetax**4.0e0_DP)  ![1] A21 isothermal compressibility

        zetabar = 0.0e0_DP

        do i_set1=1, nstypes
        do i_set2=1, nstypes
            zetabar = zetabar + Seg_array(i_set1)%xsi * Seg_array(i_set2)%xsi * sig(i_set1,i_set2)**3.0e0_DP 
        end do
        end do
        !Averaged packing fractions
        zetabar = NA * zetabar * PI * rhos / 6.0e0_DP

    end subroutine Zetax_set         
!***************************************************************************************************
    subroutine Alpha_set()  ![1] A24
        implicit none
        
        if(allocated(alpha)) deallocate(alpha)
        allocate(alpha(1:nstypes, 1:nstypes))
        
        do i_set1=1, nstypes
        do i_set2=1, nstypes
            alpha(i_set1, i_set2) = cij(i_set1, i_set2) * (1.0e0_DP /( la(i_set1, i_set2)   &
            &   - 3.0e0_DP) - 1.0e0_DP / (lr(i_set1, i_set2) - 3.0e0_DP)) 
        end do
        end do
        
        return
    end subroutine  Alpha_set
!**************************************************************************************************
    subroutine F_set()    ![1] A26
        implicit none
        
        double precision    ::  alpha_l,f_numer,f_denom
           
        if(allocated(fkij)) deallocate(fkij)           
        allocate(fkij(1:6, 1:nstypes, 1:nstypes))
        
        do i_set1=1, nstypes
        
            do i_set2=1, nstypes
                alpha_l = alpha(i_set1, i_set2)
                
                do i_set3=1, 6    
                
                    f_numer = 0.0e0_DP
                    f_denom = 0.0e0_DP
                    
                do i_set4 = 0, 3
                    f_numer = f_numer + phik(i_set3, i_set4) * alpha_l**i_set4
                end do
                    
                    do i_set4=4, 6  
                        f_denom = f_denom + phik(i_set3, i_set4) * alpha_l**(i_set4-3)
                    end do
                    
                    fkij(i_set3, i_set1, i_set2) = f_numer / (1.0e0_DP + f_denom)
                
                end do
            end do
        end do
         
        return
    end subroutine F_set
!**************************************************************************************************    
    subroutine Chain_set() !Set up values for chain calculation
        implicit none

        real(kind=DP)    ::  zki(1:nctypes, 1:nstypes), zsum
        
        zki = 0.0e0_DP ![2] 13
        
        do i_set1=1, nctypes       
            zsum = 0.0e0_DP
        
            do i_set2=1, nstypes
                zsum = zsum + Comp_array(i_set1)%comp(i_set2) * Seg_array(i_set2)%nseg * Seg_array(i_set2)%sf
            end do
            
            do i_set2=1, nstypes
                zki(i_set1, i_set2) =  Comp_array(i_set1)%comp(i_set2) * Seg_array(i_set2)%nseg *       &
                    &        Seg_array(i_set2)%sf / zsum   
            end do
        end do
        
        if(allocated(sig3ij)) deallocate(sig3ij,d3ij,epschij,lachij,lrchij)  
        allocate(sig3ij(1:nctypes, 1:nctypes), d3ij(1:nctypes) ,                  &
        &       epschij(1:nctypes, 1:nctypes), lachij(1:nctypes), lrchij(1:nctypes))
        
        !All averaged parameters for molecule
        sig3ij  = 0.0e0_DP     ![2] 14  
        d3ij    = 0.0e0_DP     ![2] 15
        epschij = 0.0e0_DP     ![2] 16
        lachij  = 0.0e0_DP     ![2] 17
        lrchij  = 0.0e0_DP
        sigx3   = 0.0e0_DP     ![2] 25 
        
        do i_set1=1, nctypes   !Form molecular averages
            do i_set2=1, nstypes
            do i_set3=1, nstypes
                sig3ij(i_set1, i_set1)= sig3ij(i_set1, i_set1) + zki(i_set1, i_set2) *      &
                &                       zki(i_set1, i_set3) * sig(i_set2, i_set3)**3.0e0_DP

                d3ij(i_set1)          = d3ij(i_set1) + zki(i_set1, i_set2) * zki(i_set1, i_set3) *  &
                &                       dij(i_set2, i_set3)**3.0e0_DP            
                epschij(i_set1, i_set1) = epschij(i_set1,i_set1) + zki(i_set1, i_set2) * zki(i_set1, i_set3) *  &
                &                         eps(i_set2, i_set3)
                lachij(i_set1)          = lachij(i_set1) + zki(i_set1, i_set2) * zki(i_set1, i_set3) *      &
                &                         la(i_set2, i_set3)
                lrchij(i_set1)          = lrchij(i_set1) + zki(i_set1, i_set2) * zki(i_set1, i_set3) *      &
                &                         lr(i_set2, i_set3)              
            end do
            end do           
        end do

        do i_set1=1, nctypes   !Estimate intermolecular averages (for Aassoc)
        do i_set2=1, nctypes
            if(i_set1/=i_set2) then
                sig3ij(i_set1, i_set2)   =   (sig3ij(i_set2, i_set2) + sig3ij(i_set1, i_set1)) /2.0e0_DP
                epschij(i_set1, i_set2)  =   dsqrt(sig3ij(i_set2, i_set2) * sig3ij(i_set1, i_set1)) / &
                    &   sig3ij(i_set1, i_set2)  * dsqrt(epschij(i_set2, i_set2) * epschij(i_set1, i_set1))
            end if
        end do
        end do
    
        do i_set1=1, nstypes
        do i_set2=1, nstypes
            sigx3 = sigx3 + Seg_array(i_set1)%xsi * Seg_array(i_set2)%xsi * sig(i_set1, i_set2)**3.0e0_DP
        end do
        end do

        return
    end subroutine Chain_set 
!***************************************************************************************************
    function Mieii(xm_in, sm_in, em_in, cm_in, lrm_in, lam_in) result(mie_out) !A2 for self-interaction  
        implicit none
        
        real(kind=DP),intent(in)    ::  xm_in,sm_in,em_in,cm_in,lrm_in,lam_in
        real(kind=DP)               ::  mie_out
        real(kind=DP)               ::  t_mie

        t_mie = sm_in / xm_in        
        mie_out = cm_in * em_in * (t_mie**lrm_in - t_mie**lam_in)

        return
    end function MIEII
!***************************************************************************************************

!***************************************************************************************************
end module Setup
!***************************************************************************************************