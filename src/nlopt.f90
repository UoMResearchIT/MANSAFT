Module nlopt

	use Types           ! Definitions of types and double precision
	use Global          ! Important global parameters
	use Pressure
	use Input
	use Input_opt   ! Read optimisation parameters
	use Vol
	use solver
	use ChemPot
!     .. Implicit None Statement ..
!    Implicit None
!     .. Accessibility Statements ..
      
    Contains
	
	subroutine opt_vle(output_opt)
	  implicit none
!	  external						:: myfunc
	  real(kind=DP), allocatable		:: lb(:), ub(:), x(:), output_opt(:)
	  real(kind=DP)					:: minf
	  integer 						:: opt, ires, i
!	  integer					:: opt_num    !This in Global mod
	  include 'nlopt.f'
	  
!	  opt_num = 4												!Set the number of optimized parameters
	  allocate (lb(1:opt_num),ub(1:opt_num),x(1:opt_num))  
	  
	  call nlo_create(opt, NLOPT_LN_BOBYQA, opt_num)			!Set the algorithm
	  
	  do i=1, opt_num
			lb(i) = min_num(i) 									!Set lower bounds
			ub(i) = max_num(i) 									!Set upper bounds 
			x(i) = init_values(i) 								!Set initial guesses
	  end do
	  
	  call nlo_set_lower_bounds(ires, opt, lb)					!Set lower bounds
	    
	  call nlo_set_upper_bounds(ires, opt, ub)					!Set upper bounds
	  
	  call nlo_set_min_objective(ires, opt, myfunc, 0)			!Set minimal or maximal objiective functions
	  
	  call nlo_set_xtol_rel(ires, opt, 1.D-4)					!Set convergence conditions
	  
	  ! initial guesses											!Set initial guesses
	  !x(1) = 3.0
	  
	  call nlo_optimize(ires, opt, x, minf)
      if (ires.lt.0) then
         write(*,*) 'nlopt failed!'
      else
!         write(*,*) 'found min at ', x(1:opt_num)
         write(*,102) 'min f = ', minf
		 do i=1, opt_num
		   output_opt(i) = x(i)
		 end do
      end if
	  
102 Format ((2X,A10),(2X,F6.4))  
	end subroutine opt_vle
	
	subroutine myfunc(val, n_in, x, grad, need_gradient, f_data)
	  !implicit none
	  real(kind=DP)					:: val, x(n_in), grad(n_in), mu_sum
	  integer 						:: n_in, need_gradient, i, j
	  real(kind=DP),allocatable     :: output_vle(:), mu_vle(:)
	  real(Kind=DP)				::  p_oj, y_oj(2), p_squ, y_squ(2)
	  
	  do i = 1, opt_num
			if(param_key(i)==1) then               
				!min_bound(i)=min_num(i) * ANG
				!max_bound(i)=max_num(i) * ANG
				sig(param_index(i,1),param_index(i,1)) = x(i)* ANG/1000  !1000 is to make nag get converged
        !print*, Seg_array(param_index(i,1))%sig, x(i)
            else if(param_key(i)==2) then
				!min_bound(i)=min_num(i)
				!max_bound(i)=max_num(i)
				eps(param_index(i,1),param_index(i,1)) = x(i)/10     !10 is to make nag get converged 
            else if(param_key(i)==3) then
				!min_bound(i)=min_num(i)
				!max_bound(i)=max_num(i)
				lr(param_index(i,1),param_index(i,1)) = x(i)/100   !100 is to make nag get converged
            else if(param_key(i)==4) then
				!min_bound(i)=min_num(i)
				!max_bound(i)=max_num(i)
				la(param_index(i,1),param_index(i,1)) = x(i)
            else if(param_key(i)==5) then
				!min_bound(i)=min_num(i)
				!max_bound(i)=max_num(i)
				Seg_array(param_index(i,1))%sf = x(i)/10000 !10000 is to make nag get converged
            else if(param_key(i)==6) then
				!min_bound(i)=min_num(i)
				!max_bound(i)=max_num(i)
				Seg_array(param_index(i,1))%nseg = x(i)
			else if(param_key(i)==7)then
				!min_bound(i)=min_num(i)
				!max_bound(i)=max_num(i)
				eps(param_index(i,1), param_index(i,2))=x(i)/10  !10 is to make nag get converged
				eps(param_index(i,2), param_index(i,1))=x(i)/10 !10 is to make nag get converged				
            else if(param_key(i)==8)then
				!min_bound(i)=min_num(i)
				!max_bound(i)=max_num(i)
				lr(param_index(i,1), param_index(i,2))=x(i)/100 !100 is to make nag get converged
				lr(param_index(i,2), param_index(i,1))=x(i)/100 !100 is to make nag get converged
            else if(param_key(i)==9)then
				!min_bound(i)=min_num(i)
				!max_bound(i)=max_num(i)
				ehb(param_index(i,1),param_index(i,2), &
				& param_index2(i,1),param_index2(i,2)) = x(i)
				ehb(param_index2(i,1),param_index2(i,2), &
				& param_index(i,1),param_index(i,2)) = x(i)
            else if(param_key(i)==10)then
				!min_bound(i)=min_num(i) * ANG**3.0_DP * NA
				!max_bound(i)=max_num(i) * ANG**3.0_DP * NA
				khb(param_index(i,1),param_index(i,2), &
				& param_index2(i,1),param_index2(i,2)) = x(i)* ANG**3.0_DP * NA/10 !10 is to make nag get converged
				khb(param_index2(i,1),param_index2(i,2), &
				& param_index(i,1),param_index(i,2)) = x(i)* ANG**3.0_DP * NA/10 !10 is to make nag get converged
			end if
	  end do
		
		p_oj = 0
		y_oj(1)=0
		y_oj(2)=0
		p_squ=0
		y_squ(1)=0
		y_squ(2)=0
		


	  if (properties%opt_l(6)) then
	  allocate(output_vle(1:nctypes+1), mu_vle(1:nctypes))
		do i = 1, properties%nmu
          Comp_array(1:nctypes)%xi = properties%ximu(i,1:nctypes)
		  t = properties%t_opt(i)
		  p = properties%p_opt(i)*1000
		  comp_x(1:nctypes) = Comp_array(1:nctypes)%xi
		  P = 1.5*P
		  v = Vol_dens_g( )
		  !mu_1 = Mu_res(1)
		  !mu_2 = Mu_res(2)
		  x_init(1) = v
		  x_init(2) = 8.314*t/(properties%p_opt(i)*1000)
		  mu_sum = 0
		  do j=1, nctypes
 		    mu_vle(j) = dexp(Mu_res(j)/(8.314*t))*comp_x(j)
			mu_sum = mu_sum + mu_vle(j)
		  end do
		  do j=1, nctypes-1
		    x_init(j+2) =  mu_vle(j)/mu_sum
		  end do

		  call solve_nle(output_vle)
		  
			p_oj=p_oj+ABS(((P_L+P_V)/2000-properties%p_opt(i))/properties%p_opt(i))
			p_squ=p_squ+(((P_L+P_V)/2000-properties%p_opt(i))/properties%p_opt(i))**2
			y_oj(1)=y_oj(1)+ABS((output_vle(3)-properties%yimu(i,1))/properties%yimu(i,1))
			y_squ(1)=y_squ(1)+((output_vle(3)-properties%yimu(i,1))/properties%yimu(i,1))**2
			y_oj(2)=y_oj(2)+ABS((1.0_DP-output_vle(3)-properties%yimu(i,2))/properties%yimu(i,2))
			y_squ(2)=y_squ(2)+((1.0_DP-output_vle(3)-properties%yimu(i,2))/properties%yimu(i,2))**2
		end do
			p_oj=p_oj/properties%nmu
			p_squ=p_squ/properties%nmu
			y_oj(1)=y_oj(1)/properties%nmu
			y_squ(1)=y_squ(1)/properties%nmu
			y_oj(2)=y_oj(2)/properties%nmu
			y_squ(2)=y_squ(2)/properties%nmu
	  end if		
        
		
		val = p_squ + y_squ(1) + y_squ(2)
   !print*, "f          x(i)"
      write(*, 101)"f=",val, "p_oj=", p_oj, "y_oj1=", y_oj(1), "y_oj2=", y_oj(2),"x(i)=",x(1:opt_num)
101 Format ((2X,A2),(2X,F6.4),(2X,A5),(2X,F6.4),(2X,A6),(2X,F6.4),(2X,A6),&
          &(2X,F6.4),(2X,A5),10(2X,F10.4))   
	end subroutine myfunc


	
end module  
	 
	
