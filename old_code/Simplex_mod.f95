!***************************************************************************************************
!   Module to:
!       Optimise a set of parameters using SIMPLEX method
!***************************************************************************************************
!   Created by SJ Halstead July 2016
!   Updates:
!       1. Sept 2016 - rearranged modules for FORTRAN standards
!***************************************************************************************************
!
!***************************************************************************************************
module Simplex_mod
!***************************************************************************************************
!MODULES        !NO EXTERNAL VALUES
use Zig_mod 
use Error_mod   ! calculate simplex errors
!***************************************************************************************************
    implicit none
!***************************************************************************************************
contains
!***************************************************************************************************
!***************************************************************************************************
    subroutine Simp(nsimp, simpin, simperr, simpout, maxit, maxparam, minparam)
        implicit none
        !-------------------------------------------------------------------------------------------
        !Input / Output
        integer,intent(in)              ::  nsimp                       !Number of paramters
        integer,intent(in)              ::  maxit                       !Maximum number of iterations
        real(kind=DP)                   ::  simpin(1:nsimp+1, 1:nsimp)  !Initial parameter guesses   
        real(kind=DP),intent(out)       ::  simperr                     !Final error
        real(kind=DP),intent(out)       ::  simpout(1:nsimp)            !Optimised parameters
        real(kind=DP)                   ::  simpfun                     !Function used for optimisation
        real(kind=DP),intent(in)        ::  maxparam(1:nsimp), minparam(1:nsimp)
        !-------------------------------------------------------------------------------------------
        !Convergence criteria
        real(kind=DP), parameter     ::  ERR_CRIT=1.D-4
        real(kind=DP)                ::  err_cyc(1:nsimp*2), err_vold(1:100)
        !-------------------------------------------------------------------------------------------
        !Local variables
        real(kind=DP)                ::  error_array(1:nsimp+1)             !error outputs for each point
        real(kind=DP)                ::  point_array(1:nsimp+1, 1:nsimp)    !ordered array, best to worst
        real(kind=DP)                ::  mid_array(1:nsimp)
        real(kind=DP)                ::  reflect_array(1:nsimp), reflect_error
        real(kind=DP)                ::  extend_array(1:nsimp), extend_error, error_sum,error_old
        real(kind=DP)                ::  contract_array(1:nsimp), contract_error, simp_prev(1:nsimp)
        integer                      ::  i_simp1,i_simp2,mainit,nbounce
        logical                      ::  collapsed, possible, cycle, real_value, isetup
        !-------------------------------------------------------------------------------------------
        !Check initial values are possible
        point_array=simpin
     
        print*, "Initial SIMPLEX values:"
        print*, "======================"
        print*,
     
        do i_simp1=1, nsimp+1
        
        do i_simp2=1, nsimp  
            if((point_array(i_simp1,i_simp2) > maxparam(i_simp2)).or.&
            &   (point_array(i_simp1,i_simp2) < minparam(i_simp2))) then
                print*, "Initial input values outside of limits, specfically:"
                print*, "    Point ", i_simp1, " value ", i_simp2
                print*, "    Value is: ",point_array(i_simp1, i_simp2)
                print*, "    Max and Min are: ",maxparam(i_simp2), minparam(i_simp2)
                stop
            end if
            
        end do
        
        print*,i_simp1,point_array(i_simp1,:)
        
        end do          
        
        print*,
        print*,
        print*, "Calculating initial errors now"
        print*,
        print*,
        
        call flush()
        !-------------------------------------------------------------------------------------------
        !Calculate initial errors
        do i_simp1=1, nsimp+1         
            if((properties%type=='opt').or.(properties%type=='OPT').or.(properties%type=='Opt')) then
                error_array(i_simp1) = Opt_err(nsimp, point_array(i_simp1, 1:nsimp))
            else if((properties%type=='eopt').or.(properties%type=='EOPT')) then                
                error_array(i_simp1) = EOpt_err(nsimp, point_array(i_simp1, 1:nsimp))
            else if((properties%type=='ilopt').or.(properties%type=='ILOPT')) then
                error_array(i_simp1) = ILOpt_err(nsimp, point_array(i_simp1, 1:nsimp),isetup)
            end if
            
            !Check not NaN
            real_value=.true.
            if(error_array(i_simp1)/=error_array(i_simp1)) real_value=.false.
                
            do while(.not.real_value)
                print*, "initial value NaN: ",i_simp1,point_array(i_simp1, 1:nsimp)
            
                call Bounce(point_array(i_simp1, 1:nsimp), nsimp, maxparam, minparam)
                
                if((properties%type=='opt').or.(properties%type=='OPT').or.(properties%type=='Opt')) then
                    error_array(i_simp1) = Opt_err(nsimp, point_array(i_simp1, 1:nsimp))
                else if((properties%type=='eopt').or.(properties%type=='EOPT')) then                
                    error_array(i_simp1) = EOpt_err(nsimp, point_array(i_simp1, 1:nsimp))
                else if((properties%type=='ilopt').or.(properties%type=='ILOPT')) then
                    error_array(i_simp1) = ILOpt_err(nsimp, point_array(i_simp1, 1:nsimp),isetup)
                end if
                
                if(error_array(i_simp1)==error_array(i_simp1)) real_value=.true.
            end do    

            print*, i_simp1, error_array(i_simp1)         
            call flush()

        end do
      
        call Order(nsimp, error_array, point_array)
        error_sum = sum(error_array)
        mainit = 0
        nbounce = 0
        err_cyc = 1.e10_DP
        cycle = .false.

        do       
            !Calculate mid point 
            mid_array=0.e0_DP
            
            do i_simp1=1, 2    !Take best 2 points
            do i_simp2=1, nsimp
                mid_array(i_simp2) = mid_array(i_simp2) + point_array(i_simp1,i_simp2)
            end do    
            end do
       
            mid_array = mid_array / 2.0e0_DP

            !REFLECT
            possible=.true.
            do i_simp1 = 1, nsimp
                reflect_array(i_simp1) = 2.0e0_DP * mid_array(i_simp1) - point_array(nsimp+1, i_simp1)                   
                if((reflect_array(i_simp1) > maxparam(i_simp1)).or.(reflect_array(i_simp1) <    &
                &   minparam(i_simp1))) possible=.false.
            end do
      
            if(possible) then
                if((properties%type=='opt').or.(properties%type=='OPT').or.(properties%type=='Opt')) then
                    reflect_error = Opt_err(nsimp, reflect_array(1:nsimp))
                else if((properties%type=='eopt').or.(properties%type=='EOPT')) then                
                    reflect_error = EOpt_err(nsimp, reflect_array(1:nsimp))
                else if((properties%type=='ilopt').or.(properties%type=='ILOPT')) then
                    reflect_error = ILOpt_err(nsimp, reflect_array(1:nsimp))
                end if
            end if
       
            if((reflect_error <= error_array(nsimp+1)).and.(possible)) then
                !Extend
                possible = .TRUE.
                do i_simp1=1, nsimp
                    extend_array(i_simp1) = reflect_array(i_simp1) + mid_array(i_simp1) - point_array(nsimp+1,i_simp1)                    
                    if((extend_array(i_simp1) > maxparam(i_simp1)).or.(extend_array(i_simp1) <  &
                    &   minparam(i_simp1))) possible = .false.    
                end do

                if(possible) then
                    if((properties%type=='opt').or.(properties%type=='OPT').or.(properties%type=='Opt')) then
                        extend_error = Opt_err(nsimp, extend_array(1:nsimp) )
                    else if((properties%type=='eopt').or.(properties%type=='EOPT')) then                
                        extend_error = EOpt_err(nsimp, extend_array(1:nsimp))
                    else if((properties%type=='ilopt').or.(properties%type=='ILOPT')) then
                        extend_error = ILOpt_err(nsimp, extend_array(1:nsimp))
                    end if
                end if

                if((extend_error <= reflect_error).and.(possible)) then    !Extend is best
                    do i_simp1 = 1, nsimp                   
                        point_array(nsimp+1, i_simp1) = extend_array(i_simp1)                        
                    end do    
                       
                    error_array(nsimp+1) = extend_error                    
                else                                    !Reflect is best
              
                    do i_simp1 = 1, nsimp
                        point_array(nsimp+1, i_simp1) = reflect_array(i_simp1)                     
                    end do   
                     
                    error_array(nsimp+1) = reflect_error                    
                end if
            else
                !CONTRACT (reflect and extend not good)
                do i_simp1=1,nsimp            
                    contract_array(i_simp1) = (mid_array(i_simp1) + point_array(nsimp+1, i_simp1)) / 2.0e0_DP     
                end do
        
                if((properties%type=='opt').or.(properties%type=='OPT').or.(properties%type=='Opt')) then
                    contract_error = Opt_err(nsimp, contract_array(1:nsimp))
                else if((properties%type=='eopt').or.(properties%type=='EOPT')) then                
                    contract_error = EOpt_err(nsimp, contract_array(1:nsimp))
                else if((properties%type=='ilopt').or.(properties%type=='ILOPT')) then
                    contract_error = ILOpt_err(nsimp, contract_array(1:nsimp))
                end if

                do i_simp1=1, nsimp               
                    point_array(nsimp+1, i_simp1) = contract_array(i_simp1)
                end do
              
                error_array(nsimp+1) = contract_error                
            end if
         !-------------------------------------------------------------------------------------------   
            !Check new value not NaN
            real_value=.true.
            if(error_array(nsimp+1)/=error_array(nsimp+1)) real_value=.false.
                
            do while(.not.real_value)
                call Bounce(point_array(nsimp+1, 1:nsimp), nsimp, maxparam, minparam)
                error_array(nsimp+1) = Opt_err(nsimp, point_array(nsimp+1, 1:nsimp))
                if(error_array(nsimp+1)==error_array(nsimp+1)) real_value=.true.
            end do    

            call order(nsimp, error_array, point_array) !Reorder arrays, best to worst 
        !-------------------------------------------------------------------------------------------  
            error_old = error_sum
            error_sum = sum(error_array)
        !-------------------------------------------------------------------------------------------
        !CHECK CONVERGENCE
            if(error_array(1) < ERR_CRIT) then
                do i_simp1 = 1, nsimp
                    simpout(i_simp1) = point_array(1, i_simp1)    
                end do
                
                simperr         = error_array(1)
                print*,  "Converged ",mainit               
                exit
            end if
               
            !CHECK SIMPLEX HASN'T COLLAPSED
            collapsed = .true.
            
            do i_simp1 = 1, nsimp
                if(dabs(point_array(1, i_simp1) - point_array(nsimp+1, i_simp1) / dabs(point_array(1, i_simp1))) &
                &   > ERR_CRIT) collapsed = .false.
            end do
            
            if(collapsed) then
                do i_simp1 = 1, nsimp
                    simpout(i_simp1) = point_array(1, i_simp1)    
                end do
                
                simperr         = error_array(1)
                print*,  "Collapsed ",mainit
                
                exit
            end if
            
            !GONE AS FAR AS POSSIBLE
            if(dabs(error_sum - error_old) / error_sum < ERR_CRIT) then          
                do i_simp1 = 1, nsimp
                    simpout(i_simp1) = point_array(1, i_simp1)    
                end do
                
                simperr         = error_array(1)
                print*,  "Gone as far as possible ",mainit
                exit
            end if
            
            if(mainit>100) then          
                if(error_array(1) == err_vold(1)) then                
                    do i_simp1 = 1, nsimp
                        simpout(i_simp1) = point_array(1, i_simp1)                     
                    end do
                    
                    simperr         = error_array(1)

                    print*,  "Error not decreases for 100 iterations ",mainit
                    exit
                end if
            end if
            
            !TOO MANY BOUNCES?
            if(nbounce > 50) then            
                do i_simp1=1,nsimp
                    simpout(i_simp1) = point_array(1, i_simp1)
                end do

                simperr         = error_array(1)
                print*,  "Bounced enough ",mainit
                exit
            end if
            
            !CHECK IF TOO MANY ITERATIONS
            mainit = mainit + 1             
       
            if(mainit > maxit) stop "Maximum iterations reached"
            
            !CHECK CYCLIC ERRORS
            if(mainit > nsimp * 2) then
            do i_simp1 = 1, nsimp * 2
       
                if(dabs(error_array(nsimp+1) / err_cyc(i_simp1) - 1.0e0_DP) <= 1.0e-3_DP) then 
                    call Bounce(point_array(nsimp+1, 1:nsimp), nsimp, maxparam, minparam)
                    
                    if((properties%type=='opt').or.(properties%type=='OPT').or.(properties%type=='Opt')) then
                        error_array(nsimp+1) = Opt_err(nsimp, point_array(nsimp+1, 1:nsimp))
                    else if((properties%type=='eopt').or.(properties%type=='EOPT')) then                
                        error_array(nsimp+1) = EOpt_err(nsimp, point_array(nsimp+1, 1:nsimp))
                    else if((properties%type=='ilopt').or.(properties%type=='ILOPT')) then
                        error_array(nsimp+1) = ILOpt_err(nsimp, point_array(nsimp+1, 1:nsimp))
                    end if
    
                    call order(nsimp, error_array, point_array)
                    nbounce = nbounce + 1

                    print*, 
                    print*,  " ***   BOUNCE  *** "
                    print*, 
                    go to 90
                end if
            end do
              
            !CHECK CYCLIC ERRORS
            cycle = .false.
            do i_simp1=1, nsimp * 2          
                if(dabs(error_array(nsimp+1)-err_cyc(i_simp1))/error_array(nsimp+1) <= 1.0e-3_DP) &
                &   cycle = .true.         
            end do
            
            if(cycle) then  
                    call Bounce(point_array(nsimp+1, 1:nsimp), nsimp, maxparam, minparam)
                    
                    if((properties%type=='opt').or.(properties%type=='OPT').or.(properties%type=='Opt')) then
                        error_array(nsimp+1) = Opt_err(nsimp, point_array(nsimp+1, 1:nsimp))
                    else if((properties%type=='eopt').or.(properties%type=='EOPT')) then                
                        error_array(nsimp+1) = EOpt_err(nsimp, point_array(nsimp+1, 1:nsimp))
                    else if((properties%type=='ilopt').or.(properties%type=='ILOPT')) then
                        error_array(nsimp+1) = ILOpt_err(nsimp, point_array(nsimp+1, 1:nsimp))
                    end if
                    
                    call Order(nsimp, error_array, point_array)
                    nbounce = nbounce + 1
                    
                    print*,  " ***   BOUNCE  *** "
                    go to 90      
                end if
            
            end if
          
            !Set cycle errors            
 90         do i_simp1=2, nsimp * 2
                err_cyc(i_simp1-1) = err_cyc(i_simp1)
            end do
            
            do i_simp1=2, 100
                err_vold(i_simp1-1) = err_vold(i_simp1) 
            end do
            
            err_cyc(nsimp * 2) = error_array(nsimp+1)
            err_vold(100)      = error_array(1)

            print*, mainit," errors: ", error_array
            print*, "Best point: ",point_array(1,:)
            print*, "#####################################################################"
            call flush()         
        end do
              
        return
    end subroutine Simp
!***************************************************************************************************
    subroutine Order(nsimp_l, error_array_l, point_array_l)
        implicit none

        integer, intent(in)             ::  nsimp_l                         !Number of paramters
        real(kind=DP)                   ::  error_array_l(1:nsimp_l+1)        !error outputs for each point
        real(kind=DP)                   ::  point_array_l(1:nsimp_l+1, 1:nsimp_l)!ordered array, best to worst
        real(kind=DP)                   ::  temp_array(1:nsimp_l), temp_err  !temporary data for bubble sort
        integer                         ::  corder, corder2
        logical                         ::  ordered
  
        ordered = .false.
        
        do while(.not.ordered)
            ordered = .true.
        do corder=1, nsimp_l
            if(error_array_l(corder+1) < error_array_l(corder)) then
                !Temporary store values
                temp_err = error_array_l(corder)    
            
                do corder2 = 1, nsimp_l
                    temp_array(corder2) = point_array_l(corder, corder2)
                end do
            
                !Move values
                error_array_l(corder) = error_array_l(corder+1)
                
                do corder2=1, nsimp_l
                    point_array_l(corder, corder2) = point_array_l(corder+1, corder2)
                end do
            
                !Move temp values
                error_array_l(corder+1) = temp_err
                
                do corder2=1, nsimp_l
                    point_array_l(corder+1,corder2) = temp_array(corder2)
                end do
                
                ordered = .false.
            end if
        end do
        end do

        return
    end subroutine Order
!***************************************************************************************************
   subroutine Bounce(point_array_l, nsimp_l, maxin, minin)
        implicit none
        
        integer, intent(in)             ::  nsimp_l
        real(kind=DP), intent(inout)    ::  point_array_l(1:nsimp_l)
        integer                         ::  cbounce
        real(kind=DP), intent(in)       ::  maxin(:),minin(:)
        
        do cbounce=1, nsimp_l
            point_array_l(cbounce) = (maxin(cbounce) - minin(cbounce)) * uni() + minin(cbounce)
        end do
        
        return
    end subroutine
!***************************************************************************************************
end module Simplex_mod
!***************************************************************************************************