!***************************************************************************************************
!   SAFT Module to:
!       1. READ INPUT AND ASSIGN ARRAYS
!***************************************************************************************************
!   Created by SJ Halstead Oct 2015
!       1. Tidied 23/6/2016
!***************************************************************************************************
!
!***************************************************************************************************
module Input
!***************************************************************************************************
!MODULES
use Global     !CONTAINS GLOBAL VARIABLES
use Types
use GL         !Calculate G-L dabscissas ans weights
!***************************************************************************************************
    implicit none
!***************************************************************************************************
contains
!***************************************************************************************************
    subroutine Read_input( testinput )
!--------------------------------------------------------------------------
        implicit none
!--------------------------------------------------------------------------        
!Variables        
        !local do counters
        integer         ::  i1,i2,i3,i4,ierr
        integer         ::  j1,j2,j3,j4
        !input file
        character(len=80)               ::  filein
        character(len=80)               ::  line
        character(len=100)               ::  longline
        character(len=*),optional       ::  testinput
        !temp values
        real(kind=DP)   ::  doub,doub2
        !gl input
        integer                         ::  cfin        
        character(len=80)               ::  path,path1        
!--------------------------------------------------------------------------        
!Read file        
        if(present(testinput)) then
            filein = testinput
        else  
            call getarg(1,filein)        
        end if
   
        open(11, file=filein, status="old", action="read", iostat=ierr)
     
        if(ierr/=0) then
            print*, "Failed to open ",filein
            stop
        end if
    
        read(11,'(a80)') line
        read(11,*) 
      
        
        print*, "Reading file:"
        write(*,'(a5, a80)') "     ", line       
!--------------------------------------------------------------------------             
!Segments      
        chain_switch = .false.       
        assoc_switch = .false.
        ion_switch   = .false.
  
        read(11,*) nstypes       

        allocate(Seg_array(1:nstypes), stat=ierr)
      
        if(ierr/=0) then
            print*, "Failed to allocate segment array, nstypes = ",nstypes
            stop
        end if
      
        do i1=1, nstypes
            read(11,*, iostat=ierr)                                                 &
            &   Seg_array(i1)%lab, Seg_array(i1)%sig, Seg_array(i1)%eps,            &
            &   Seg_array(i1)%lr, Seg_array(i1)%la, Seg_array(i1)%sf,               &
            &   Seg_array(i1)%nseg, Seg_array(i1)%m, Seg_array(i1)%q,               & 
            &   Seg_array(i1)%sigb, Seg_array(i1)%nassoc                                

            !Set ion_switch
            if(.not.ion_switch) then
                if(dabs(Seg_array(i1)%q)>1.0e-6_DP) ion_switch=.true.
            end if
            
            if(ierr/=0) then
                print*, "Problem reading component ",i1
                stop
            end if
            
            Seg_array(i1)%charged=.true. 
            if(dabs(Seg_array(i1)%q) <=1.0e-6_DP) Seg_array(i1)%charged=.false.  

            !Number of associations types
            i3 = Seg_array(i1)%nassoc
            allocate(Seg_array(i1)%nassoctyp(1:i3), stat=ierr)
            
            if(ierr/=0) then
                print*, "Failed to allocate associations array for component ",i1
                stop
            end if
            
            !Number of sites for each association type on a segment
            do i2=1,i3
                read(11,*) Seg_array(i1)%nassoctyp(i2)               
            end do
            
            !Adjust units to SI
            Seg_array(i1)%sig  = Seg_array(i1)%sig  * ANG
            Seg_array(i1)%sigb = Seg_array(i1)%sigb * ANG
            Seg_array(i1)%m    = Seg_array(i1)%m    * AMU
        end do
        
        read(11,*)      
!--------------------------------------------------------------------------      
!Form interaction arrays        
        allocate(eps(1:nstypes, 1:nstypes), stat=ierr)
        if(ierr/=0) then
            print*, "Failed to allocate eps array "
            stop
        end if
            
        allocate(sig(1:nstypes, 1:nstypes), stat=ierr)
        if(ierr/=0) then
            print*, "Failed to allocate sig array "
            stop
        end if
        
        allocate(la(1:nstypes, 1:nstypes), stat=ierr)
        if(ierr/=0) then
            print*, "Failed to allocate lambda association array "
            stop
        end if
        
        allocate(lr(1:nstypes, 1:nstypes), stat=ierr)
        if(ierr/=0) then
            print*, "Failed to allocate lambda repulsion array "
            stop
        end if
       
        do i2=1,nstypes                                                                 !Full 2D array
        do i1=i2,nstypes
            sig(i1, i2)   =  (Seg_array(i1)%sig + Seg_array(i2)%sig) / 2.0_DP           ![1] A45
            eps(i1, i2)  =  dsqrt(Seg_array(i1)%sig**3.0_DP * Seg_array(i2)%sig**3.0) / sig(i1,i2)**3.0_DP &
            &   * dsqrt(Seg_array(i1)%eps * Seg_array(i2)%eps)                           ![1] A51
            
            la(i1, i2) = dsqrt((Seg_array(i1)%la - 3.0_DP) * (Seg_array(i2)%la - 3.0_DP)) + 3.0_DP     ![1] A56
            lr(i1, i2) = dsqrt((Seg_array(i1)%lr - 3.0_DP) * (Seg_array(i2)%lr - 3.0_DP)) + 3.0_DP     ![1] A56
            
            !Reflect 
            sig(i2, i1) = sig(i1, i2)
            eps(i2, i1) = eps(i1, i2)
            la(i2, i1)  = la(i1, i2)
            lr(i2, i1)  = lr(i1, i2)
        end do
        end do
!--------------------------------------------------------------------------     
!Read unlikes
        !unlike with specific eps       
        read(11,*) i2
        un_neps = i2
        allocate(un_eps(1:un_neps), un_ieps(1:un_neps), un_jeps(1:un_neps) )
      
        do i1=1,i2
            read(11,*) j1, j2, doub
          
            eps(j1,j2) = doub
            eps(j2,j1) = doub
            
            un_eps(i1)  = doub
            un_ieps(i1) = j1
            un_jeps(i1) = j2
        end do
       
        read(11,*)
        read(11,*) i2
        un_nlr = i2
        allocate(un_lr(1:un_nlr), un_ilr(1:un_nlr), un_jlr(1:un_nlr) )
       
        do i1=1,i2
            read(11,*) j1, j2, doub
         
            lr(j1,j2) = doub
            lr(j2,j1) = doub
            
            un_lr(i1)  = doub
            un_ilr(i1) = j1
            un_jlr(i1) = j2
        end do
              
        read(11,*)   
!--------------------------------------------------------------------------      
!Read associations
        allocate(ehb(1:nstypes, 1:3, 1:nstypes, 1:3), stat=ierr)
        if(ierr/=0) then
            print*, "Failed to allocate ehb association array "
            stop
        end if
        
        allocate(khb(1:nstypes, 1:3, 1:nstypes, 1:3), stat=ierr)
        if(ierr/=0) then
            print*, "Failed to allocate khb association array "
            stop
        end if
        
        ehb=0.0_DP   !energy term
        khb=0.0_DP   !volume term
        
        read(11,*) i3
        
        !Set assoc switch
        if(i3/=0) assoc_switch=.true.

        do i1=1,i3
            read(11,*) j1,j2,j3,j4,doub,doub2        
        
            ehb(j1,j2,j3,j4) = doub
            khb(j1,j2,j3,j4) = doub2 * ANG**3.0_DP * NA  !Scale to m^3 mol^-1
            
            ehb(j3,j4,j1,j2) = doub
            khb(j3,j4,j1,j2) = doub2 * ANG**3.0_DP * NA  !Scale to m^3 mol^-1
        end do                                             
!--------------------------------------------------------------------------     
!Read molecule descriptions
        read(11,*) nctypes   !number of molecules / components
        allocate(Comp_array(1:nctypes),stat=ierr)
       
        if(ierr/=0) then
            print*, "Failed to allocate comp array "
            stop
        end if    

        do i1=1,nctypes
            !Molecule description
            read(11,*, iostat=ierr)  Comp_array(i1)%lab, Comp_array(i1)%xi, Comp_array(i1)%m,     &
            &                        Comp_array(i1)%dv, Comp_array(i1)%dt,                        &   
            &                        Comp_array(i1)%ms   !name, mole fraction, mass, number of segments

            if(ierr/=0) then
                print*, "Problem reading component array, component ",i1
                stop
            end if
        
            Comp_array(i1)%solv = .true.
            if((dabs(Comp_array(i1)%dv)<1.0e-6_DP).and.(dabs(Comp_array(i1)%dv)<1.0e-6_DP))       &
            &   Comp_array(i1)%solv=.false.
        
            Comp_array(i1)%m = Comp_array(i1)%m * AMU  !to kg per particle
            
            allocate(Comp_array(i1)%comp(1:nstypes), stat=ierr)
            if(ierr/=0) then
                print*, "Failed to allocate comp molecule list array "
                stop
            end if
        
            Comp_array(i1)%comp = 0
            
            !List number of each type of segment in molecule
            do i2 = 1, Comp_array(i1)%ms
                read(11,*) j1, j2
                Comp_array(i1)%comp(j1) = j2                              
            end do
            !set no. segments to total no.
            Comp_array(i1)%ms = sum(Comp_array(i1)%comp(1:Comp_array(i1)%ms))
            !Set chain switch         
            if(.not.chain_switch) then
                if(Comp_array(i1)%ms>1) chain_switch=.true.
                
                do i2=1, nstypes
                    if(Comp_array(i1)%comp(i2)/=0) then
                        if(Seg_array(i2)%nseg>1.0e0_DP) then
                            chain_switch=.true.
                        end if
                    end if
                end do
                
            end if

        end do
    
        read(11,*)     
!--------------------------------------------------------------------------     
!Calculation type
        read(11,*) properties%type, properties%n

        if(allocated(properties%v)) deallocate(properties%v)
        allocate(properties%v(1:properties%n), stat=ierr)
        if(ierr/=0) then
            print*, "Failed to allocate property v array "
            stop
        end if 
        
        if(allocated(properties%t)) deallocate(properties%t)
        allocate(properties%t(1:properties%n), stat=ierr)
        if(ierr/=0) then
            print*, "Failed to allocate property t array "
            stop
        end if 
        
        if(allocated(properties%p)) deallocate(properties%p)
        allocate(properties%p(1:properties%n), stat=ierr)
        if(ierr/=0) then
            print*, "Failed to allocate property p array "
            stop
        end if 
        
        if(allocated(properties%phase)) deallocate(properties%phase)
        allocate(properties%phase(1:properties%n), stat=ierr)
        if(ierr/=0) then
            print*, "Failed to allocate property phase array "
            stop
        end if 
        
        if(allocated(properties%xi)) deallocate(properties%xi)
        allocate(properties%xi(1:properties%n, 1:nctypes), stat=ierr)
        if(ierr/=0) then
            print*, "Failed to allocate property xi array "
            stop
        end if 
        
        select case(properties%type)
            !Pressure
            case('p','P')  
                do i1=1,properties%n
                    read(11,*) properties%t(i1), properties%v(i1)
                end do
            !Pressure for varying conc
            case('px','PX')  
                do i1=1,properties%n
                    read(11,*) properties%t(i1), properties%v(i1), properties%xi(i1, 1:nctypes)
                end do
            !Volume    
            case('v','V')   
                do i1=1,properties%n
                    read(11,*) properties%t(i1), properties%p(i1)
                end do
            !Volume for varying mole fractions   
            case('vx','VX')   
                do i1=1,properties%n                    
                    read(11,*) properties%t(i1), properties%p(i1), properties%xi(i1, 1:nctypes)
                end do
            !Pure phase equilibrium
            case('pha1','PHA1','Pha1')
                if(nctypes>1) stop "Method  PHA1  is only for single component (pure) phase equilibrium"
                
                do i1=1,properties%n
                    read(11,*) properties%t(i1)
                end do
            !Binary VLE
            case('vle','VLE')
                if(nctypes/=2) print*, "Method  VLE  is only for binary liquid-vapour equilibrium"
               
                do i1=1,properties%n
                    read(11,*) properties%t(i1), properties%p(i1),properties%xi(i1, 1)                  
                end do
               
            !Activity coefficient
            case('act','ACT','Act')
                read(11,*) t, p
                read(11,*) properties%cation, properties%cat_stoich
                read(11,*) properties%anion, properties%an_stoich
                     
                do i1=1,properties%n
                    read(11,*) properties%xi(i1, 1:nctypes)
                end do   
            !Cp heat capacity at constant pressure
            case('cp','Cp','CP')  
                read(11,*) properties%A, properties%B, properties%C, properties%D, properties%E
            
                do i1=1,properties%n
                    read(11,*) properties%xi(i1, 1:nctypes),properties%t(i1), properties%p(i1)
                end do 
            !Themodynamic state data
            case('Therm','therm','THERM')  
                read(11,*) properties%A, properties%B, properties%C, properties%D, properties%E
                
                do i1=1,properties%n
                    read(11,*) properties%xi(i1, 1:nctypes),properties%t(i1), properties%p(i1)
                end do 
            !Gibbs free energy
            case('Gval','gval','GVAL')  
                do i1=1,properties%n
                    read(11,*) properties%xi(i1, 1:nctypes),properties%t(i1), properties%p(i1)
                end do
            !Excess chemical potential
            case('Mu','mu','MU')  
                do i1=1,properties%n
                    read(11,*) properties%xi(i1, 1:nctypes), properties%v(i1), properties%t(i1)
                end do                 
!********************************************************************************************************************************
            !special - temporary testng input choice
            case('spec')
!                do i1=1,properties%n
!                    read(11,*) properties%t(i1),properties%p(i1)!,properties%xi(i1,1:nctypes)
!                end do
!********************************************************************************************************************************                                               
            !Optimiser
            case('OPT','Opt','opt')   
                read(11,*) properties%opt(:) 
                properties%opt_l(:) = .false.
                
                if((properties%opt(1)=='y').or.(properties%opt(1)=='Y')) properties%opt_l(1)=.true. !Liq V
                if((properties%opt(2)=='y').or.(properties%opt(2)=='Y')) properties%opt_l(2)=.true. !Liq Cp
                if((properties%opt(3)=='y').or.(properties%opt(3)=='Y')) properties%opt_l(3)=.true. !Liq H
                if((properties%opt(4)=='y').or.(properties%opt(4)=='Y')) properties%opt_l(4)=.true. !VLE Psat
                if((properties%opt(5)=='y').or.(properties%opt(5)=='Y')) properties%opt_l(5)=.true. !VLE V
				if((properties%opt(6)=='y').or.(properties%opt(6)=='Y')) properties%opt_l(6)=.true. !VLE binary  yichun
                       
                read(11,*)
                
                !Liquid volumes
                if(properties%opt_l(1)) then
                    read(11,*) properties%nliqv
                    
                    allocate(properties%t_liqv(1:properties%nliqv), stat=ierr)               
                    allocate(properties%p_liqv(1:properties%nliqv), stat=ierr)
                    allocate(properties%v_liq(1:properties%nliqv), stat=ierr)
                    allocate(properties%xiv(1:properties%nliqv,1:nctypes), stat=ierr)               
                    
                    do i1=1, properties%nliqv
                        read(11,*) properties%xiv(i1,1:nctypes),properties%t_liqv(i1),properties%p_liqv(i1), &
                        &   properties%v_liq(i1)
                    end do
                    
                    read(11,*)
                end if
                
                !Liquid heat capacities
                if(properties%opt_l(2)) then
                    read(11,*) properties%nliqc
                    read(11,*) properties%A, properties%B, properties%C, properties%D, properties%E
                    
                    allocate(properties%t_liqc(1:properties%nliqc), stat=ierr)               
                    allocate(properties%p_liqc(1:properties%nliqc), stat=ierr)
                    allocate(properties%xic(1:properties%nliqc,1:nctypes), stat=ierr)               
                    allocate(properties%cp(1:properties%nliqc), stat=ierr)
                    
                    do i1=1, properties%nliqc
                        read(11,*) properties%xic(i1,1:nctypes),properties%t_liqc(i1),properties%p_liqc(i1), &
                        &   properties%cp(i1)
                    end do
                    
                    read(11,*)
                end if
                
                !Liquid total enthalpies
                if(properties%opt_l(3)) then
                    read(11,*) properties%nliqh
                    
                    allocate(properties%t_liqh(1:properties%nliqh), stat=ierr)               
                    allocate(properties%p_liqh(1:properties%nliqh), stat=ierr)
                    allocate(properties%xih(1:properties%nliqh,1:nctypes), stat=ierr)               
                    allocate(properties%h(1:properties%nliqh), stat=ierr)
                    
                    do i1=1, properties%nliqh
                        read(11,*) properties%xih(i1,1:nctypes),properties%t_liqh(i1),properties%p_liqh(i1), &
                        &   properties%h(i1)
                    end do
                    
                    read(11,*)
                end if
                 
                if((properties%opt_l(4)).or.(properties%opt_l(5))) then
                    if(nctypes/=1) stop "Can only optimise pure VLE data at present!"
                
                    read(11,*) properties%nmu
                    allocate(properties%v_l(1:properties%nmu))
                    allocate(properties%v_v(1:properties%nmu))
                    allocate(properties%p_opt(1:properties%nmu))
                    allocate(properties%t_opt(1:properties%nmu))
                    allocate(properties%ximu(1:properties%nmu,1:nctypes))
                    
                    properties%v_l(:)   = 0.0e0_DP
                    properties%v_v(:)   = 0.0e0_DP
                    properties%p_opt(:) = 0.0e0_DP
                    properties%t_opt(:) = 0.0e0_DP
                    
                    do i1=1, properties%nmu
                        read(11,'(a100)') longline  
                        
                        if((properties%opt_l(4)).and.(properties%opt_l(5))) then
                            read(longline,*) properties%ximu(i1,1:nctypes), properties%t_opt(i1), properties%p_opt(i1), &
                                &   properties%v_l(i1), properties%v_v(i1)
                        else if(properties%opt_l(5)) then
                            read(longline,*)  properties%ximu(i1,1:nctypes), properties%t_opt(i1), properties%p_opt(i1),  &
                            &   properties%v_l(i1), properties%v_v(i1)
                        else if(properties%opt_l(4)) then    
                            read(longline,*)  properties%ximu(i1,1:nctypes), properties%t_opt(i1),  &
                            &   properties%v_l(i1), properties%v_v(i1)
                        else
                            print*, "Why am I in this loop?"
                        end if
                    end do
                end if
				
				!The followings are added by yichun
				if(properties%opt_l(6)) then
                    read(11,*) properties%nmu
                    
                    allocate(properties%t_opt(1:properties%nmu), stat=ierr)               
                    allocate(properties%p_opt(1:properties%nmu), stat=ierr)
                    allocate(properties%ximu(1:properties%nmu,1:nctypes), stat=ierr)               
                    allocate(properties%yimu(1:properties%nmu,1:nctypes), stat=ierr)
                    
                    do i1=1, properties%nmu
                        read(11,*) properties%t_opt(i1),properties%p_opt(i1),properties%ximu(i1,1:nctypes), &
                        &   properties%yimu(i1,1:nctypes)
                    end do
                    
                end if
!********************************************************************************************************************************                                               
            !Electrolyte optimiser
            case('eopt','EOPT')   
                
                read(11,*) properties%nliqv_e
                
                allocate(properties%xi_liqve(1:properties%nliqv_e,1:nctypes),    &
                &   properties%t_liqve(1:properties%nliqv_e),  &
                &   properties%p_liqve(1:properties%nliqv_e),  &
                &   properties%v_liqve(1:properties%nliqv_e))
                
                do i1=1, properties%nliqv_e
                    read(11,*) properties%t_liqve(i1),properties%p_liqve(i1),properties%xi_liqve(i1,1:nctypes), &
                    &   properties%v_liqve(i1)
                end do
                
                read(11,*)
                read(11,*) properties%ndhmix
                
                allocate(properties%xidhmixi(1:properties%ndhmix,1:nctypes), &
                &   properties%xidhmixf(1:properties%ndhmix,1:nctypes), &
                &   properties%tdhmix(1:properties%ndhmix), properties%pdhmix(1:properties%ndhmix), &
                &   properties%dhmix(1:properties%ndhmix))
                
                do i1=1, properties%ndhmix
                    read(11,*) properties%tdhmix(i1), properties%pdhmix(i1),  &
                    &   properties%xidhmixi(i1,1:nctypes), properties%xidhmixf(i1,1:nctypes), properties%dhmix(i1)
                end do
                
!********************************************************************************************************************************
            !Ionic liquid optimiser
            case('ilopt','ILOPT')   
                
                read(11,*) properties%nliqv
                
                allocate(properties%xiv(1:properties%nliqv,1:nctypes),    &
                &   properties%t_liqv(1:properties%nliqv),  &
                &   properties%p_liqv(1:properties%nliqv),  &
                &   properties%v_liq(1:properties%nliqv))
                
                do i1=1, properties%nliqv
                    read(11,*) properties%t_liqv(i1),properties%p_liqv(i1),properties%xiv(i1,1:nctypes), &
                    &   properties%v_liq(i1)
                end do            
                
                !Liquid heat capacities
                read(11,*)
                read(11,*) properties%nliqc
                    
                allocate(properties%t_liqc(1:properties%nliqc), stat=ierr)               
                allocate(properties%p_liqc(1:properties%nliqc), stat=ierr)
                allocate(properties%xic(1:properties%nliqc,1:nctypes), stat=ierr)               
                allocate(properties%cp(1:properties%nliqc), stat=ierr)
                
                do i1=1, properties%nliqc
                    read(11,*) properties%t_liqc(i1),properties%p_liqc(i1),properties%xic(i1,1:nctypes), &
                    &   properties%cp(i1)
                end do
!********************************************************************************************************************************             
            case default
                print*, "Unrecognised property ", properties%type
                stop
        end select
!*************************************************************************************************** 
!Read in / setup up constant arrays
        !coefficients in [1] A17 (effective packing fraction)
        c_array(1, 1) = 0.81096e0_DP
        c_array(1, 2) = 1.7888e0_DP
        c_array(1, 3) = -37.578e0_DP
        c_array(1, 4) = 92.284e0_DP
        
        c_array(2, 1) = 1.0205e0_DP
        c_array(2, 2) = -19.341e0_DP
        c_array(2, 3) = 151.26e0_DP
        c_array(2, 4) = -463.50e0_DP
        
        c_array(3, 1) = -1.9057e0_DP
        c_array(3, 2) = 22.845e0_DP
        c_array(3, 3) = -228.14e0_DP
        c_array(3, 4) = 973.92e0_DP
        
        c_array(4, 1) = 1.0885e0_DP
        c_array(4, 2) = -6.1962e0_DP
        c_array(4, 3) = 106.98e0_DP
        c_array(4, 4) = -677.64e0_DP
        
        !coefficients in [1] A26 (part of monomer contribution)
        phik      = 0.e0_DP
        
        phik(1, 0) = 7.5365557e0_DP
        phik(1, 1) = -37.60463e0_DP
        phik(1, 2) = 71.745953e0_DP
        phik(1, 3) = -46.83552e0_DP
        phik(1, 4) = -2.467982e0_DP
        phik(1, 5) = -0.50272e0_DP
        phik(1, 6) = 8.0956883e0_DP
        
        phik(2, 0) = -359.44e0_DP
        phik(2, 1) = 1825.6e0_DP
        phik(2, 2) = -3168.0e0_DP
        phik(2, 3) = 1884.2e0_DP
        phik(2, 4) = -0.82376e0_DP
        phik(2, 5) = -3.1935e0_DP
        phik(2, 6) = 3.7090e0_DP
        
        phik(3, 0) = 1550.9e0_DP
        phik(3, 1) = -5070.1e0_DP
        phik(3, 2) = 6534.6e0_DP
        phik(3, 3) = -3288.7e0_DP
        phik(3, 4) = -2.7171e0_DP
        phik(3, 5) = 2.0883e0_DP
        phik(3, 6) = 0.0000e0_DP
        
        phik(4, 0) = -1.19932e0_DP
        phik(4, 1) = 9.063632e0_DP
        phik(4, 2) = -17.9482e0_DP
        phik(4, 3) = 11.34027e0_DP
        phik(4, 4) = 20.52142e0_DP
        phik(4, 5) = -56.6377e0_DP
        phik(4, 6) = 40.53683e0_DP
        
        phik(5, 0) = -1911.28e0_DP
        phik(5, 1) = 21390.175e0_DP
        phik(5, 2) = -51320.7e0_DP
        phik(5, 3) = 37064.54e0_DP
        phik(5, 4) = 1103.742e0_DP
        phik(5, 5) = -3264.61e0_DP
        phik(5, 6) = 2556.181e0_DP
        
        phik(6, 0) = 9236.9e0_DP
        phik(6, 1) = -129430e0_DP
        phik(6, 2) = 357230e0_DP
        phik(6, 3) = -315530e0_DP
        phik(6, 4) = 1390.2e0_DP
        phik(6, 5) = -4518.2e0_DP
        phik(6, 6) = 4241.6e0_DP
        
        phik(7, 0) = 10.e0_DP
        phik(7, 1) = 10.e0_DP
        phik(7, 2) = 0.57e0_DP
        phik(7, 3) = -6.7e0_DP
        phik(7, 4) = -8.e0_DP
        phik(7, 5) = 0.e0_DP
        phik(7, 6) = 0.e0_DP
        
        !coefficients in [2] 24  (part of association,  parameters originally from [3])
        c_x=0.e0_DP
        
        c_x(0, 0) = 7.56425183020431e-02_DP
        c_x(0, 1) = -1.28667137050961e-01_DP
        c_x(0, 2) = 1.28350632316055e-01_DP
        c_x(0, 3) = -7.25321780970292e-02_DP
        c_x(0, 4) = 2.57782547511452e-02_DP
        c_x(0, 5) = -6.01170055221687e-03_DP
        c_x(0, 6) = 9.33363147191978e-04_DP
        c_x(0, 7) = -9.55607377143667e-05_DP
        c_x(0, 8) = 6.19576039900837e-06_DP
        c_x(0, 9) = -2.30466608213628e-07_DP
        c_x(0, 10) = 3.74605718435540e-09_DP
        c_x(1, 0) = 1.34228218276565e-01_DP
        c_x(1, 1) = -1.82682168504886e-01_DP
        c_x(1, 2) = 7.71662412959262e-02_DP
        c_x(1, 3) = -7.17458641164565e-04_DP
        c_x(1, 4) = -8.72427344283170e-03_DP
        c_x(1, 5) = 2.97971836051287e-03_DP
        c_x(1, 6) = -4.84863997651451e-04_DP
        c_x(1, 7) = 4.35262491516424e-05_DP
        c_x(1, 8) = -2.07789181640066e-06_DP
        c_x(1, 9) = 4.13749349344802e-08_DP
        c_x(2, 0) = -5.65116428942893e-01_DP
        c_x(2, 1) = 1.00930692226792e+00_DP
        c_x(2, 2) = -6.60166945915607e-01_DP
        c_x(2, 3) = 2.14492212294301e-01_DP
        c_x(2, 4) = -3.88462990166792e-02_DP
        c_x(2, 5) = 4.06016982985030e-03_DP
        c_x(2, 6) = -2.39515566373142e-04_DP
        c_x(2, 7) = 7.25488368831468e-06_DP
        c_x(2, 8) = -8.58904640281928e-08_DP
        c_x(3, 0) = -3.87336382687019e-01_DP
        c_x(3, 1) = -2.11614570109503e-01_DP
        c_x(3, 2) = 4.50442894490509e-01_DP
        c_x(3, 3) = -1.76931752538907e-01_DP
        c_x(3, 4) = 3.17171522104923e-02_DP
        c_x(3, 5) = -2.91368915845693e-03_DP
        c_x(3, 6) = 1.30193710011706e-04_DP
        c_x(3, 7) = -2.14505500786531e-06_DP
        c_x(4, 0) = 2.13713180911797e+00_DP
        c_x(4, 1) = -2.02798460133021e+00_DP
        c_x(4, 2) = 3.36709255682693e-01_DP
        c_x(4, 3) = 1.18106507393722e-03_DP
        c_x(4, 4) = -6.00058423301506e-03_DP
        c_x(4, 5) = 6.26343952584415e-04_DP
        c_x(4, 6) = -2.03636395699819e-05_DP
        c_x(5, 0) = -3.00527494795524e-01_DP
        c_x(5, 1) = 2.89920714512243e+00_DP
        c_x(5, 2) = -5.67134839686498e-01_DP
        c_x(5, 3) = 5.18085125423494e-02_DP
        c_x(5, 4) = -2.39326776760414e-03_DP
        c_x(5, 5) = 4.15107362643844e-05_DP
        c_x(6, 0) = -6.21028065719194e+00_DP
        c_x(6, 1) = -1.92883360342573e+00_DP
        c_x(6, 2) = 2.84109761066570e-01_DP
        c_x(6, 3) = -1.57606767372364e-02_DP
        c_x(6, 4) = 3.68599073256615e-04_DP
        c_x(7, 0) = 1.16083532818029e+01_DP
        c_x(7, 1) = 7.42215544511197e-01_DP
        c_x(7, 2) = -8.23976531246117e-02_DP
        c_x(7, 3) = 1.86167650098254e-03_DP
        c_x(8, 0) = -1.02632535542427e+01_DP
        c_x(8, 1) = -1.25035689035085e-01_DP
        c_x(8, 2) = 1.14299144831867e-02_DP
        c_x(9, 0) = 4.65297446837297e+00_DP
        c_x(9, 1) = -1.92518067137033e-03_DP
        c_x(10, 0) = -8.67296219639940e-01_DP
        
!Gaussian Legendre setup for [1] A9 integral
        call GL_arrays( )
!***************************************************************************************************

!*************************************************************************************************** 
!Housekeeping
        close(11)   
!***************************************************************************************************  
        return
    end subroutine Read_input        
!***************************************************************************************************
end module Input
