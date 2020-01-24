!***************************************************************************************************
!   SAFT Module to:
!       1. Calculate pure phase equilibrium
!       
!       Recast for efficiency April 2017
!
!***************************************************************************************************
!
!***************************************************************************************************
module Pure_phase_mod
!***************************************************************************************************
!Modules
!=======    
    use Types_mod           ! Definitions of types and double precision
    use Global_mod          ! Important global parameters
    use Press_mod           ! Calculate pressue
    use Mu_mod              ! Calculate chemical potential
!***************************************************************************************************
    implicit none
!***************************************************************************************************
    contains
!***************************************************************************************************    
    subroutine Pure_phase( vl, vv )    
!***************************************************************************************************     
        real(kind=DP), intent(out)      ::  vl, vv
        integer, parameter              ::  nvpoints=5500, nvstep=500
        real(kind=DP)                   ::  dv_local, dp_local, v_local(1:nvpoints), p_local(1:nvpoints)

        real(kind=DP), parameter        ::  p_max=100000000.0e0_DP !100 MPa        
        real(kind=DP)                   ::  v_accept(1:5), mu_accept(1:5), lmu_temp
        
        real(kind=DP)                   ::  v_big, v_small, mu_big, mu_small, con1, con2
        real(kind=DP)                   ::  vold_big, vold_small, muold_big, muold_small, pold
        
        real(kind=DP)                   ::  grad_local, grad1, grad2
        
        integer                         ::  i_ph, i_ph2, i_ph3  
        integer                         ::  cvp, accept1, accept2   
        
        logical                         ::  converged         
!***************************************************************************************************
        dv_local=1.0e-7_DP
        i_ph=0
        
        do i_ph2=1,11
        do i_ph3=1, nvstep
            i_ph=i_ph+1
            v = v + dv_local
            v_local(i_ph) = v
            p_local(i_ph) = Press()                                                          
        end do 
            dv_local = dv_local * 10.0e0_DP
        end do         
!---------------------------------------------------------------------------------
        open(31, file="temp.data", status="unknown", action="write")
        
        p=0.0e0_DP
        dp_local = 1.0e-3_DP
        i_ph2=0
        
        do while(p < p_max)       
            cvp=0
            p = p+dp_local
           
            !Find volumes corresponding to desired P                
            do i_ph=1, nvpoints   

                if(((p >= p_local(i_ph)).and.(p <= p_local(i_ph+1))).or.              &
                &(p <= p_local(i_ph)).and.(p >= p_local(i_ph+1))) then
            !Intepolate             
                    grad_local = (p_local(i_ph+1) - p_local(i_ph)) / (v_local(i_ph + 1) - v_local(i_ph))  
                    v = (p - p_local(i_ph)) / grad_local + v_local(i_ph)
                  
                    if(v==v) then !filter out spurious volumes
                           
                        lmu_temp = Mu(1)!filter out spurious mu
                        if((lmu_temp==lmu_temp).and.(lmu_temp<0.0e0_DP)) then
                            cvp = cvp + 1
                            v_accept(cvp) = v             
                            mu_accept(cvp) = lmu_temp 
                        end if
                    end if
                end if
  
            end do                   
!---------------------------------------------------------------------------------        
            if(cvp>1) then           
                accept1=1
                accept2=2
          
                !find lowest Mu
                do i_ph=3, cvp
                    if(mu_accept(i_ph)<mu_accept(accept1)) then
                        if(mu_accept(accept1)<mu_accept(accept2)) then
                            accept2=i_ph
                        else
                            accept1=i_ph
                        end if
                    else if(mu_accept(i_ph)<accept2) then
                        accept2=i_ph
                    end if
                end do
         
                !determine large / small
                if(v_accept(accept1)<v_accept(accept2)) then
                    v_small  = v_accept(accept1)
                    mu_small = mu_accept(accept1)
                    v_big    = v_accept(accept2)
                    mu_big   = mu_accept(accept2)
                else
                    v_small  = v_accept(accept2)
                    mu_small = mu_accept(accept2)
                    v_big    = v_accept(accept1)
                    mu_big   = mu_accept(accept1)                
                end if

                write(31,*) p, v_big, mu_big, v_small, mu_small                               
            end if
              
            i_ph2 = i_ph2+1
            if(mod(i_ph2,1000)==0) dp_local=dp_local*10.0e0_DP  

        end do
        
        close(31)   
!---------------------------------------------------------------------------------
!Analyse
        open(31, file="temp.data", status="old", action="read")
        
        i_ph = 0
        converged = .false.

        do
            read(31, *, end=99) p, v_big, mu_big, v_small, mu_small
            
            if(i_ph /= 0) then
                
               if((( mu_big > muold_small ).and.( muold_big < mu_small )).or.&
                    &(( mu_big < muold_small ).and.( muold_big > mu_small ))) then 
        
                    grad1   =   (mu_big - muold_big) / (p - pold)
                    con1    =    mu_big - grad1 * p
                    grad2   =   (mu_small - muold_small) / (p - pold)
                    con2    =    mu_small - grad2 * p
        
                    if((grad1>0.0e0_DP).and.(grad2>0.0e0_DP)) then
                        p   =   (con2 - con1) / (grad1 - grad2)
                        
                        grad1   =   (v_small - vold_small) / (p - pold)   
                        grad2   =   (v_big - vold_big) / (p - pold)

                        v_small = grad1 * (p - pold) + vold_small
                        v_big   = grad2 * (p - pold) + vold_big
                           
                        converged = .true.
                        go to 99
                    end if
                end if
            end if
            
            pold            =   p
            vold_big        =   v_big
            muold_big       =   mu_big
            vold_small      =   v_small
            muold_small     =   mu_small
            
            i_ph = i_ph + 1
            
        end do

 99     close(31, status="delete")

        if(converged) then
            vl = v_small
            vv = v_big  
        else
            print*,"damn"
            stop          
        end if   
!***************************************************************************************************    
        return
    end subroutine Pure_phase
!*************************************************************************************************** 
    subroutine E_vle( l_vl, l_vv, l_p ) 
        implicit none
        
        real(kind=DP), intent(out)  ::  l_vl, l_vv, l_p
        integer, parameter  ::  nvap=600, nliq=3000
        real(kind=DP)       ::  ldv, lrho
        real(kind=DP)       ::  xi_old(1:nctypes)
        real(kind=DP)       ::  l_vvap(1:nvap), l_pvap(1:nvap), l_mvap(1:nvap) 
        real(kind=DP)       ::  l_vliq(1:nliq), l_pliq(1:nliq), l_mliq(1:nliq)
        real(kind=DP)       ::  l_m1, l_m2, l_c1, l_c2
        logical             ::  t_ion_switch
        integer             ::  i_evle, i_evle2, c_evle, ndos

        !vapour
        xi_old(1:nctypes) = Comp_array(1:nctypes)%xi
        Comp_array(1:nctypes)%xi = 0.0e0_DP
        Comp_array(1)%xi = 1.0e0_DP
        t_ion_switch = ion_switch
        ion_switch = .false.
          
        lrho=0.0e0_DP
        c_evle=0
        
        do i_evle2=1,5
            if(i_evle2==1) then
                ldv=1.00e-3_DP
                ndos=100
            else if(i_evle2==2) then
                ldv=1.00e-2_DP
                ndos=100
            else if(i_evle2==3) then
                ldv=1.00e-1_DP
                ndos=100
            else if(i_evle2==4) then
                ldv=1.00e0_DP
                ndos=100
            else if(i_evle2==5) then
                ldv=1.00e1_DP
                ndos=200
            end if 
        do i_evle=1,ndos      
            lrho = lrho+ldv
            c_evle=c_evle+1
            
            l_vvap(c_evle) = 1.0e0_DP / lrho
            v              = l_vvap(c_evle)
            l_pvap(c_evle) = Press( )
            l_mvap(c_evle) = Mu(1)                      
        end do
        end do
      
        !liquid
        Comp_array(1:nctypes)%xi = xi_old(1:nctypes)
        ion_switch = t_ion_switch
 
        lrho=4.0e4_DP        
        ldv=1.00e1_DP
   
        do i_evle=1,nliq   
            lrho = lrho + ldv                     
            l_vliq(i_evle) = 1.0e0_DP / lrho
            v              = l_vliq(i_evle) 
            l_pliq(i_evle) = Press( )
            l_mliq(i_evle) = Mu_neut(1)                         
        end do
     
        !scan
        do i_evle  = 1, nliq-1
        do i_evle2 = 1, nvap-1     
            if(((l_pliq(i_evle)<l_pvap(i_evle2)).and.(l_pliq(i_evle+1)>l_pvap(i_evle2+1))).or.&
            &((l_pliq(i_evle)>l_pvap(i_evle2)).and.(l_pliq(i_evle+1)<l_pvap(i_evle2+1)))) then
          
            if(((l_mliq(i_evle)<l_mvap(i_evle2)).and.(l_mliq(i_evle+1)>l_mvap(i_evle2+1))).or.&
            &((l_mliq(i_evle)>l_mvap(i_evle2)).and.(l_mliq(i_evle+1)<l_mvap(i_evle2+1)))) then
                              
                !linear interpolation
                l_m1 = (l_mliq(i_evle+1) - l_mliq(i_evle)) / (l_pliq(i_evle+1) - l_pliq(i_evle))
                l_c1 = l_mliq(i_evle) - l_m1 * l_pliq(i_evle) 
               
                l_m2 = (l_mvap(i_evle2+1) - l_mvap(i_evle2)) / (l_pvap(i_evle2+1) - l_pvap(i_evle2))
                l_c2 = l_mvap(i_evle2) - l_m2 * l_pvap(i_evle2) 
              
                p = (l_c2 - l_c1) / (l_m1 - l_m2)
                l_p = p
           
                !volumes
                l_m1 = (l_pliq(i_evle+1) - l_pliq(i_evle))/ (l_vliq(i_evle+1) - l_vliq(i_evle))
                l_m2 = (l_pvap(i_evle2+1) - l_pvap(i_evle2))/ (l_vvap(i_evle2+1) - l_vvap(i_evle2))  
                 
                l_vl = (l_p - l_pliq(i_evle)) / l_m1 + l_vliq(i_evle)
                l_vv = (l_p - l_pvap(i_evle2)) / l_m2 + l_vvap(i_evle2) 
                            
                go to 99
               
            end if 
            end if        
        end do
        end do

 99     return
    end subroutine E_vle
!***************************************************************************************************   

!***************************************************************************************************   

!***************************************************************************************************   

!***************************************************************************************************       
end module Pure_phase_mod
!***************************************************************************************************