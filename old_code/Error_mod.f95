!***************************************************************************************************
!   SAFT Module to:
!       1. Calculate optimisation error
!***************************************************************************************************
!
!***************************************************************************************************
module Error_mod
!***************************************************************************************************
!Modules
!=======
    use Types_mod       ! Definitions of types and double precision
    use Global_mod      ! Important global parameters
    use Press_mod       ! Pressure function
    use Vol_mod         ! Volume function
    use Mu_mod          ! Mu function
    use Therm_mod       ! Heat capacity, includes enthalpy function
    use Pure_phase_mod  ! Calculate equilibrium pressure
!***************************************************************************************************
contains
!***************************************************************************************************
    function Opt_err( nval_l, val_l ) result(error_f )
        implicit none
        
        integer, intent(in)             ::  nval_l
        integer                         ::  i_err, i_opterr1, i_opterr2, seg1_opt
        real(kind=DP), intent(in)       ::  val_l(1:nval_l)
        real(kind=DP)                   ::  error_f
        real(kind=DP)                   ::  esum1, esum2, esum3, esum4, esum5
        real(kind=DP), parameter        ::  EFACT1=1.0e0_DP, EFACT2=1.0e0_DP, EFACT3=1.0e0_DP, &
        &   EFACT4=1.0e0_DP, EFACT5=1.0e0_DP !These are weights
        
!key:
! 1 sig ii; 2 eps ii; 3 lr ii; 4 la ii; 5 sf; 6 nseg; 7 eps ij; 8 lr ij; 9 eHB; 10 kHB        
        !single beads
        do i_err=1, nval_l
            seg1_opt = param_index(i_err,1)
        
            if(param_key(i_err)==1) then
                Seg_array(seg1_opt)%sig = val_l(i_err) * ANG
            else if(param_key(i_err)==2) then
                Seg_array(seg1_opt)%eps = val_l(i_err) 
            else if(param_key(i_err)==3) then
                Seg_array(seg1_opt)%lr = val_l(i_err) 
            else if(param_key(i_err)==4) then
                Seg_array(seg1_opt)%la = val_l(i_err)
            else if(param_key(i_err)==5) then
                Seg_array(seg1_opt)%sf = val_l(i_err) 
            else if(param_key(i_err)==6) then
                Seg_array(seg1_opt)%nseg = val_l(i_err)
            end if
        end do
        
        !apply combining rules (as in Input.f95)
        do i_opterr2=1,nstypes      !Full 2D array
        do i_opterr1=i_opterr2,nstypes          ![1] A45
            sig(i_opterr1, i_opterr2)   =  (Seg_array(i_opterr1)%sig + Seg_array(i_opterr2)%sig) / 2.0_DP 
       
            lr(i_opterr1, i_opterr2) = dsqrt((Seg_array(i_opterr1)%lr - 3.0_DP) * &
            &   (Seg_array(i_opterr2)%lr - 3.0_DP)) + 3.0_DP     ![1] A56
            
            eps(i_opterr1, i_opterr2)  =  dsqrt(Seg_array(i_opterr1)%sig**3.0_DP * Seg_array(i_opterr2)%sig**3.0) &
            &   / sig(i_opterr1,i_opterr2)**3.0_DP &
            &   * dsqrt(Seg_array(i_opterr1)%eps * Seg_array(i_opterr2)%eps)                           ![1] A51 
            
            !Reflect 
            eps(i_opterr2, i_opterr1) = eps(i_opterr1, i_opterr2)
            sig(i_opterr2, i_opterr1) = sig(i_opterr1, i_opterr2)
            lr(i_opterr2, i_opterr1)  = lr(i_opterr1, i_opterr2)
        end do
        end do  
      
        !reapply cross terms from input file
        do i_opterr1=1,un_neps
            eps(un_ieps(i_opterr1),un_jeps(i_opterr1)) = un_eps(i_opterr1) 
            eps(un_jeps(i_opterr1),un_ieps(i_opterr1)) = un_eps(i_opterr1) 
        end do
        
        do i_opterr1=1,un_nlr
            lr(un_ilr(i_opterr1),un_jlr(i_opterr1)) = un_lr(i_opterr1) 
            lr(un_jlr(i_opterr1),un_ilr(i_opterr1)) = un_lr(i_opterr1) 
        end do
       
!key:
! 1 sig ii; 2 eps ii; 3 lr ii; 4 la ii; 5 sf; 6 nseg; 7 eps ij; 8 lr ij; 9 eHB; 10 kHB         
        !apply other optimised 2 index terms
        do i_err=1, nval_l
            if(param_key(i_err)==7) then
                eps(param_index(i_err,1), param_index(i_err,2)) = val_l(i_err)
            else if(param_key(i_err)==8) then
                lr(param_index(i_err,1), param_index(i_err,2)) = val_l(i_err)
            else if(param_key(i_err)==9) then
                ehb(param_index(i_err,1),param_index2(i_err,1),&
                &   param_index(i_err,2),param_index2(i_err,2)) = val_l(i_err)
            else if(param_key(i_err)==10) then
                khb(param_index(i_err,1),param_index2(i_err,1),&
                &   param_index(i_err,2),param_index2(i_err,2)) = val_l(i_err) * ANG**3.0_DP * NA
            end if
        end do

!phew! calculate the error....  
        esum1 = 0.0e0_DP
        esum2 = 0.0e0_DP
        esum3 = 0.0e0_DP
        esum4 = 0.0e0_DP
        esum5 = 0.0e0_DP

        if(properties%opt_l(1))  call Liquid_errv( esum1 )  
        if(properties%opt_l(2))  call Liquid_errc( esum2 ) 
        if(properties%opt_l(3))  call Liquid_errh( esum3 ) 
        if(properties%opt_l(5).or.properties%opt_l(5))  call Mu_err( esum4, esum5 )        
        error_f  = EFACT1*esum1 + EFACT2*esum2 + EFACT3*esum3 + EFACT4*esum4 + EFACT5*esum5

        call flush( )  
        
        return
    end function Opt_err
!***************************************************************************************************
    function EOpt_err( nval_l, val_l ) result(error_f )
        implicit none
        
        integer, intent(in)             ::  nval_l
        integer                         ::  i_opterr1, i_opterr2
        real(kind=DP), intent(in)       ::  val_l(1:nval_l)
        real(kind=DP)                   ::  error_f, vl, vv, p, lgam, dh1, dh2
        real(kind=DP)                   ::  esum1, esum2
        real(kind=DP), parameter        ::  EFACT1=1.0e0_DP, EFACT2=1.0e0_DP, EFACT3=1.0e0_DP, &
                                        &   EFACT4=1.0e0_DP !These are weights

        stop "Opt_err method not updated as yet"  
        
        return
    end function EOpt_err
!***************************************************************************************************
    function ILOpt_err( nval_l, val_l, is_initial ) result( error_f )
        implicit none
        
        integer, intent(in)         ::  nval_l
        logical, optional           ::  is_initial
        real(kind=DP), intent(in)   ::  val_l(1:nval_l)
        integer                     ::  i_err, i_opterr1, i_opterr2, seg1_opt
        real(kind=DP)               ::  esum1, error_f

        if(present(is_initial)) then
            is_initial=.true.
            print*, "testing: ",val_l(:)
        end if
!key:
! 1 sig ii; 2 eps ii; 3 lr ii; 4 la ii; 5 sf; 6 nseg; 7 eps ij; 8 lr ij; 9 eHB; 10 kHB        
        !single beads
        do i_err=1, nval_l
            seg1_opt = param_index(i_err,1)
        
            if(param_key(i_err)==1) then
                Seg_array(seg1_opt)%sig = val_l(i_err) * ANG
            else if(param_key(i_err)==2) then
                Seg_array(seg1_opt)%eps = val_l(i_err) 
            else if(param_key(i_err)==3) then
                Seg_array(seg1_opt)%lr = val_l(i_err) 
            else if(param_key(i_err)==4) then
                Seg_array(seg1_opt)%la = val_l(i_err)
            else if(param_key(i_err)==5) then
                Seg_array(seg1_opt)%sf = val_l(i_err) 
            else if(param_key(i_err)==6) then
                Seg_array(seg1_opt)%nseg = val_l(i_err)
            end if
        end do
        
        !apply combining rules (as in Input.f95)
        do i_opterr2=1,nstypes      !Full 2D array
        do i_opterr1=i_opterr2,nstypes          ![1] A45
            sig(i_opterr1, i_opterr2)   =  (Seg_array(i_opterr1)%sig + Seg_array(i_opterr2)%sig) / 2.0_DP 
       
            lr(i_opterr1, i_opterr2) = dsqrt((Seg_array(i_opterr1)%lr - 3.0_DP) * &
            &   (Seg_array(i_opterr2)%lr - 3.0_DP)) + 3.0_DP     ![1] A56
            
            eps(i_opterr1, i_opterr2)  =  dsqrt(Seg_array(i_opterr1)%sig**3.0_DP * Seg_array(i_opterr2)%sig**3.0) &
            &   / sig(i_opterr1,i_opterr2)**3.0_DP &
            &   * dsqrt(Seg_array(i_opterr1)%eps * Seg_array(i_opterr2)%eps)                           ![1] A51 
            
            !Reflect 
            eps(i_opterr2, i_opterr1) = eps(i_opterr1, i_opterr2)
            sig(i_opterr2, i_opterr1) = sig(i_opterr1, i_opterr2)
            lr(i_opterr2, i_opterr1)  = lr(i_opterr1, i_opterr2)
        end do
        end do  
      
        !reapply cross terms from input file
        do i_opterr1=1,un_neps
            eps(un_ieps(i_opterr1),un_jeps(i_opterr1)) = un_eps(i_opterr1) 
            eps(un_jeps(i_opterr1),un_ieps(i_opterr1)) = un_eps(i_opterr1) 
        end do
        
        do i_opterr1=1,un_nlr
            lr(un_ilr(i_opterr1),un_jlr(i_opterr1)) = un_lr(i_opterr1) 
            lr(un_jlr(i_opterr1),un_ilr(i_opterr1)) = un_lr(i_opterr1) 
        end do
       
!key:
! 1 sig ii; 2 eps ii; 3 lr ii; 4 la ii; 5 sf; 6 nseg; 7 eps ij; 8 lr ij; 9 eHB; 10 kHB         
        !apply other optimised 2 index terms
        do i_err=1, nval_l
            if(param_key(i_err)==7) then
                eps(param_index(i_err,1), param_index(i_err,2)) = val_l(i_err)
            else if(param_key(i_err)==8) then
                lr(param_index(i_err,1), param_index(i_err,2)) = val_l(i_err)
            else if(param_key(i_err)==9) then
                ehb(param_index(i_err,1),param_index2(i_err,1),&
                &   param_index(i_err,2),param_index2(i_err,2)) = val_l(i_err)
            else if(param_key(i_err)==10) then
                khb(param_index(i_err,1),param_index2(i_err,1),&
                &   param_index(i_err,2),param_index2(i_err,2)) = val_l(i_err) * ANG**3.0_DP * NA
            end if
        end do

!phew! calculate the error....                
        call Liquid_errv( esum1 )       
        error_f  = esum1

        call flush( )  

        return
    end function ILOpt_err
!***************************************************************************************************
    subroutine Mu_err( sum1, sum2 )
        implicit none
   
        real(kind=DP), intent(out)  ::    sum1, sum2
        integer                     ::    mu_seg
        real(kind=DP)               ::    aa(1:4), xi_old, mu1, mu2, xi_array(1:4), vl_e, vv_e
        real(kind=DP), parameter    ::    DN=1.0e-6_DP
        integer                     ::    cmu, i_mu
  
        sum1=0.0e0_DP
        sum2=0.0e0_DP

        do mu_seg=1,nstypes
        do i_mu=1, properties%nmu           
            t=properties%t_opt(i_mu)
            Comp_array(1:nctypes)%xi = properties%ximu(i_mu,1:nctypes)
                 
            if(properties%opt_l(5)) then      
                v=properties%v_l(i_mu) 
                mu1=Mu( mu_seg )
                
                v=properties%v_v(i_mu)
                mu2=Mu( mu_seg )
             
                sum1 = sum1 + dabs((mu1 - mu2) / mu2)
            end if
        end do
        end do
        
        do i_mu=1, properties%nmu        
            t=properties%t_opt(i_mu)
            Comp_array(1:nctypes)%xi = properties%ximu(i_mu,1:nctypes)
            
            if(properties%opt_l(4)) then
                call Pure_phase( vl_e, vv_e )
                v = vv_e
                sum2 = sum2 + ((properties%p_opt(i_mu)-PRESS())/properties%p_opt(i_mu))**2.0e0_DP          
            end if   
        end do
        
        
        sum1 = sum1/REAL(properties%nmu)/nstypes
        sum2 = sum2/REAL(properties%nmu)     
      
        return
    end subroutine Mu_err
!***************************************************************************************************
    subroutine Liquid_errv( sum3 )!liq vol
        implicit none

        real(kind=DP),intent(out)   ::  sum3
        integer                     ::  c_liq
        integer, parameter          ::  ZERO=0
        logical                     ::  v_reasonable
        
        sum3=0.0e0_DP

        do c_liq=1, properties%nliqv     
            Comp_array(1:nctypes)%xi = properties%xiv(c_liq,1:nctypes)
            t = properties%t_liqv(c_liq)
            p = properties%p_liqv(c_liq) 

            v = Vol_dens_g(v_reasonable)   
            
            if(.not.v_reasonable) then
                sum3=sum3/ZERO
                go to 99
            end if
                       
            sum3 = sum3 + dabs((properties%v_liq(c_liq) - v)/properties%v_liq(c_liq))         
        end do    

        sum3 = sum3 / real(properties%nliqv)

 99     return
    end subroutine Liquid_errv
!***************************************************************************************************
    subroutine Liquid_errc( sum4 )!liq Cp
        implicit none

        real(kind=DP),intent(out)   ::  sum4
        integer                     ::  c_liq
        real(kind=DP)               ::  cp_l

        sum4=0.0e0_DP

        do c_liq=1, properties%nliqc 
            Comp_array(1:nctypes)%xi = properties%xic(c_liq,1:nctypes)
            t = properties%t_liqc(c_liq)
            p = properties%p_liqc(c_liq) 

            cp_l = Cp_ex( )
            sum4 = sum4 + dabs((properties%cp(c_liq) - cp_l) / properties%cp(c_liq))
!print*,"cp ",cp_l,properties%cp(c_liq),sum4        
        end do    

        sum4 = sum4 / real(properties%nliqc)

        return
    end subroutine Liquid_errc
!***************************************************************************************************
    subroutine Liquid_errh( sum5 )
        implicit none

        real(kind=DP),intent(out)   ::  sum5
        integer                     ::  c_liq
        real(kind=DP)               ::  h_val, h1, h2

        sum5=0.0e0_DP

        do c_liq=1, properties%nliqh
            t = properties%t_liqh(c_liq)
            p = properties%p_liqh(c_liq) 
            
            Comp_array(1)%xi = 1.0e0_DP
            Comp_array(2)%xi = 0.0e0_DP
            h1 = Enth( )
            
            Comp_array(1)%xi = 0.0e0_DP
            Comp_array(2)%xi = 1.0e0_DP
            h2 = Enth( )
            
            Comp_array(1:nctypes)%xi = properties%xih(c_liq,1:nctypes)
            h_val = Enth( ) - Comp_array(1)%xi*h1 - Comp_array(2)%xi*h2        
            sum5 = sum5 + ((properties%h(c_liq) - h_val)/properties%h(c_liq))**2.0e0_DP
        end do    

        sum5 = sum5 / real(properties%nliqh)

        return
    end subroutine Liquid_errh
!***************************************************************************************************

!***************************************************************************************************
!***************************************************************************************************
!***************************************************************************************************
end module Error_mod
!***************************************************************************************************
