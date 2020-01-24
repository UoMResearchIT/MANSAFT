!   SAFT Module to:
!       1. READ OPTIMISATION INPUT
!***************************************************************************************************
!   Created by SJ Halstead Sept 2016
!
!
!***************************************************************************************************
!
!***************************************************************************************************
module Input_opt_mod
!***************************************************************************************************
!MODULES
use Global_mod      ! CONTAINS GLOBAL VARIABLES
use Types_mod       ! Definition of types and numerical precision
use Zig_mod         ! RNG
!***************************************************************************************************
    implicit none
!***************************************************************************************************
contains
!***************************************************************************************************
    subroutine Read_opt( n_param, max_param, min_param, in_vals )
!--------------------------------------------------------------------------
        implicit none
!--------------------------------------------------------------------------        
!Variables
        integer, intent(out)                    ::  n_param
        real(kind=DP), allocatable, intent(out) ::  max_param(:), min_param(:)    
        real(kind=DP), allocatable, intent(out) ::  in_vals(:,:) 
        
        integer                                 ::  cinter1, cinter2
        character(len=80)                       ::  line_l
        character(len=50)                       ::  fileopt 

        integer                                 ::  ierr
!seed
        integer                     ::  seed
        integer, dimension(1:8)     ::  info            
!--------------------------------------------------------------------------
!Open optimise input        
        call getarg(2,fileopt)
        
        open(11, file=fileopt, status="old", action="read", iostat=ierr)
    
        if(ierr/=0) then
            print*, "Failed to open ",fileopt
            stop
        end if
!--------------------------------------------------------------------------        
!Nparams to fit
        read(11,*) n_param
        
        allocate( max_param(1:n_param), min_param(1:n_param), &
        &   param_key(1:n_param), param_index(1:n_param,1:2), param_index2(1:n_param,1:2))
!--------------------------------------------------------------------------       
!Read params
!key:
! 1 sig ii; 2 eps ii; 3 lr ii; 4 la ii; 5 sf; 6 nseg; 7 eps ij; 8 lr ij; 9 eHB; 10 kHB  
        do cinter1=1, n_param       
            read(11,'(a80)') line_l
            read(line_l,*) cinter2
        
            if((cinter2==1).or.(cinter2==2).or.(cinter2==3).or.(cinter2==4).or.&
            (cinter2==5).or.(cinter2==6)) then
                read(line_l,*) param_key(cinter1), param_index(cinter1,1), &
                &   min_param(cinter1), max_param(cinter1)
            else if((cinter2==7).or.(cinter2==8))then
                read(line_l,*) param_key(cinter1), param_index(cinter1,1:2),&
                &   min_param(cinter1), max_param(cinter1)
            else if((cinter2==9).or.(cinter2==10))then
                read(line_l,*) param_key(cinter1), param_index(cinter1,1:2),&
                &   param_index2(cinter1,1:2), min_param(cinter1), max_param(cinter1)
            else
                stop "parameter not recognised"
            end if
        end do

        close(11)
!--------------------------------------------------------------------------
!Set rng
        call date_and_time( values=info )
        seed=mod(((info(3) + info(7) + info(1)) * info(6) * info(5) +   &
        &   info(8) + info(2)), info(8) * info(8))      
        call zigset(seed) !SEED ZIG RNG
!--------------------------------------------------------------------------
!Set simplex input
        allocate( in_vals(1:n_param + 1, 1:n_param))

        do cinter1=1, n_param+1
        do cinter2=1, n_param
            in_vals(cinter1, cinter2) = (max_param(cinter2) -        &
            &       min_param(cinter2)) * uni() + min_param(cinter2)
        end do
        end do              
!--------------------------------------------------------------------------  
        return
    end subroutine Read_opt        
!--------------------------------------------------------------------------    

!***************************************************************************************************
end module Input_opt_mod
!***************************************************************************************************