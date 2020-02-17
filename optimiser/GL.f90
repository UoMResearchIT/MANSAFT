!***************************************************************************************************
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
!***************************************************************************************************
!   SAFT Module to:
!       1. Calculate GaussianLegednre quadrature abscissas and weights
!***************************************************************************************************
!
!   This module is slightly modified from:
!       http://web.aeromech.usyd.edu.au//wwwcomp/subrout.html
!
!       Gaussm3.f90
!
!***************************************************************************************************
module GL
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
	subroutine GL_arrays( )

        implicit none
        integer                     :: i, j, m
        real(kind=DP)               :: p1, p2, p3, pp, z, z1
        real(kind=DP), parameter    :: EPS=3.0e-15_DP       

	    m = (N_GL + 1) / 2
	    
	    do i=1, m				
     		z = dcos( PI * (i-0.25e0_DP) / (N_GL + 0.5e0_DP) )

100     	p1 = 1.0e0_DP
        	p2 = 0.0e0_DP

        	do j = 1, N_GL
                p3 = p2
                p2 = p1
                p1 = ((2.0e0_DP * j - 1.0e0_DP) * z * p2 - (j - 1.0e0_DP) * p3) / j
        	enddo

        	pp = N_GL * (z * p1 - p2) / (z * z - 1.0e0_DP)
        	z1 = z
        	z  = z1 - p1 / pp            

        	if (dabs(z - z1) > EPS) goto  100

      	    X_GL(i) =  - z                    	
            X_GL(N_GL+1-i) =  + z               
            W_GL(i) = 2.0d0/((1.0d0-z*z)*pp*pp) 
      	    W_GL(N_GL+1-i) = W_GL(i)               
        end do     
    end subroutine GL_arrays
!***************************************************************************************************
end module GL
