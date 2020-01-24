!   SAFT Module to:
!       1. Calculate phase equilibrium
!
!       April 2018
!
!***************************************************************************************************
!
!***************************************************************************************************
module Phase_simplex_mod
!***************************************************************************************************
!Modules
!=======  
    use Types_mod           ! Definitions of types and double precision
    use Global_mod          ! Important global parameters
    use Mu_mod              ! Calculate chemical potential
    use Press_mod           ! Calculate pressure for given volume
    use Zig_mod
    use Vol_mod
!***************************************************************************************************
    implicit none  
    
    !liquid / vapour definitions
    real(kind=DP)       ::  vliqs(1:2), vvaps(1:2), xilim(1:2), yilim(1:2)  
    
    !volume interpolation arrays
    integer                 ::  nvint
    real(kind=DP), allocatable      ::  xi_arr(:), vl_arr(:)
    real(kind=DP), allocatable      ::  yi_arr(:), vv_arr(:)
    real(kind=DP)           ::  delxi, delyi
    
    !simplex
    real(kind=DP)           ::  simp_bin(1:3,1:2), simp_err(1:3) 
    logical                 ::  failed
 
 contains   
!***************************************************************************************************    
 subroutine Bin_simp()
!***************************************************************************************************    
    implicit none    
    
    !counters
    integer                 ::  ibin, cbin1, nsimpit
    
    !intial setup
    real(kind=DP)           ::  vext(1:2), lxi1, lxi2, xitrans, vfact, vtemp(1:2)
    logical                 ::  liqsmallx
    
    !RNG
    integer                 ::  lseed, dandt(1:8)

    !simplex
    real(kind=DP)           ::  ref_err, ext_err, con_err
    real(kind=DP)           ::  ref_val(1:2), ext_val(1:2), con_val(1:2), mid_val(1:2)
    logical                 ::  logref, logext
!------------------------------------------------------------------------------------------
!Main do loop 
    do ibin = 1, properties%n 
        t = properties%t(ibin)
        p = properties%p(ibin)
        failed=.false.
!------------------------------------------------------------------------------------------
!Determine volume extremes at T & P
        Comp_array(1)%xi = 0.001e0_DP
        Comp_array(2)%xi = 1.0e0_DP - Comp_array(1)%xi
        vext(1) = Vol_dens_g()
        
        Comp_array(1)%xi = 0.999e0_DP
        Comp_array(2)%xi = 1.0e0_DP - Comp_array(1)%xi
        vext(2) = Vol_dens_g()
        
        if(vext(1)<vext(2)/2.0e0_DP) then
            vliqs(1) = vext(1)
            vvaps(2) = vext(2)
            
            xilim(1) = 0.005e0_DP
            yilim(2) = 0.995e0_DP
            liqsmallx = .true.
            
            vfact=vvaps(2)/vliqs(1)/5.0e0_DP
            if(vfact<2.0e0_DP) vfact=2.0e0_DP
        else if(vext(2)<vext(1)/2.0e0_DP) then
            vliqs(1) = vext(2)
            vvaps(2) = vext(1)
            
            xilim(2) = 0.995e0_DP
            yilim(1) = 0.005e0_DP
            liqsmallx = .false.
            
            vfact=vvaps(1)/vliqs(2)/5.0e0_DP
            if(vfact<2.0e0_DP) vfact=2.0e0_DP
        else
            print*, t, p, " volumes imply heterogeneous"
            go to 99
        end if

!Determine transition point to 1.0e-4_DP
        lxi1 = 0.005e0_DP
        lxi2 = 0.995e0_DP

        xitrans=1.0e0_DP

        do while(xitrans>1.0e-3_DP)
            Comp_array(1)%xi = (lxi1+lxi2)/2.0e0_DP
            Comp_array(2)%xi = 1.0e0_DP - Comp_array(1)%xi
            v = Vol_dens_g()      
           
            if(v/vliqs(1)>vfact)then !2.0e0_DP) then
                lxi2 = Comp_array(1)%xi
            else
                lxi1 = Comp_array(1)%xi
            end if
            
            xitrans = abs(lxi2-lxi1)         
        end do
        
        Comp_array(1)%xi = lxi1
        Comp_array(2)%xi = 1.0e0_DP - Comp_array(1)%xi
        vtemp(1) = Vol_dens_g()
        
        Comp_array(1)%xi = lxi2
        Comp_array(2)%xi = 1.0e0_DP - Comp_array(1)%xi
        vtemp(2) = Vol_dens_g()
        
        if(vtemp(1)<vtemp(2)) then  !1 is liquid
            if(liqsmallx) then
                vliqs(2) = vtemp(1)
                vvaps(1) = vtemp(2)
                
                xilim(2) = lxi1
                yilim(1) = lxi2
            else
                vliqs(1) = vtemp(1)
                vvaps(2) = vtemp(2)
                
                xilim(1) = lxi1
                yilim(2) = lxi2
            end if
        else    !2 is liquid
            if(liqsmallx) then
                vliqs(2) = vtemp(2)
                vvaps(1) = vtemp(1)
                
                xilim(2) = lxi1
                yilim(1) = lxi2
            else
                vliqs(1) = vtemp(2)
                vvaps(2) = vtemp(1)
                
                xilim(1) = lxi1
                yilim(2) = lxi2
            end if
        
        end if               
!------------------------------------------------------------------------------------------
!make V interpolation arrays
        nvint = 20

 10     call V_interp()     
!------------------------------------------------------------------------------------------
!Call grid for best single point
        call Gridder()      
!------------------------------------------------------------------------------------------
!Initiate RNG
        if(ibin==1) then
            call date_and_time(values=dandt)
            lseed = (dandt(1)*dandt(8) + dandt(2)*dandt(7) + dandt(3)*dandt(6)*dandt(5))*dandt(8)   
            call zigset(lseed)
        end if          
!------------------------------------------------------------------------------------------    
!Initiate simplex values 2 & 3
        do cbin1=2,3
            simp_bin(cbin1,1) = uni() * (xilim(2)-xilim(1)) + xilim(1)        
            simp_bin(cbin1,2) = uni() * (yilim(2)-yilim(1)) + yilim(1)
            
            simp_err(cbin1) = F_err(simp_bin(cbin1,1:2))
            if(failed) go to 99
        end do
!------------------------------------------------------------------------------------------    
        call Bin_order()
!------------------------------------------------------------------------------------------
!Simplex
        nsimpit = 0        
        write(*,'(i8,a)',advance="no") nsimpit, char(13)

        do
!mid value
            mid_val(1) = (simp_bin(1,1) + simp_bin(2,1) ) / 2.0e0_DP
            mid_val(2) = (simp_bin(1,2) + simp_bin(2,2) ) / 2.0e0_DP
!reflect
            ref_val(1) = 2.0e0_DP * mid_val(1) - simp_bin(3,1) 
            ref_val(2) = 2.0e0_DP * mid_val(2) - simp_bin(3,2)          
            
            logref=.true.
            if(ref_val(1)<xilim(1)) logref=.false.
            if(ref_val(1)>xilim(2)) logref=.false.
            if(ref_val(2)<yilim(1)) logref=.false.
            if(ref_val(2)>yilim(2)) logref=.false.
            
            if(logref) ref_err = F_err(ref_val(1:2))
            if(failed) go to 99
            
            if((ref_err<simp_err(3)).and.(logref)) then
!extend        
                ext_val(1) = 3.0e0_DP * mid_val(1) - 2.0e0_DP * simp_bin(3,1) 
                ext_val(2) = 3.0e0_DP * mid_val(2) - 2.0e0_DP * simp_bin(3,2)                               
            
                logext=.true.
                if(ext_val(1)<xilim(1)) logext=.false.
                if(ext_val(1)>xilim(2)) logext=.false.
                if(ext_val(2)<yilim(1)) logext=.false.
                if(ext_val(2)>yilim(2)) logext=.false.
                
                if(logext) ext_err = F_err(ext_val(1:2))
                if(failed) go to 99
!update            
                if((ref_err<ext_err).or.(.not.logext)) then
                    simp_bin(3,1:2) = ref_val(1:2)
                    simp_err(3)     = ref_err
                else
                    simp_bin(3,1:2) = ext_val(1:2)
                    simp_err(3)     = ext_err
                end if
                
                call Bin_order()
            else
!contract                            
                con_val(1) = (mid_val(1) + simp_bin(3,1))/2.0e0_DP 
                con_val(2) = (mid_val(2) + simp_bin(3,2))/2.0e0_DP    
                
                con_err = F_err(con_val(1:2))
                if(failed) go to 99
!update            
                if((ref_err<con_err).and.(logref)) then
                    simp_bin(3,1:2) = ref_val(1:2)
                    simp_err(3)     = ref_err
                else
                    simp_bin(3,1:2) = con_val(1:2)
                    simp_err(3)     = con_err
                end if
                
                call Bin_order()  
             end if   
!------------------------------------------------------------------------------------------                
            nsimpit = nsimpit+1
            write(*,'(i8,a)',advance="no") nsimpit, char(13)
!------------------------------------------------------------------------------------------               
!check convergence
            if(simp_err(1)<1.0e0_DP) then
                write(*,'(i8,f10.3,e10.3,5e10.3)') nsimpit,t,p,simp_bin(1,1:2),simp_err(1)
                go to 99
            end if
!------------------------------------------------------------------------------------------
!bounce
            if(mod(nsimpit,50)==0) then
                call Bouncey()
            end if
!------------------------------------------------------------------------------------------
!check maxiterations
            if(nsimpit>300) then
                write(*,'(i8,f10.3,e10.3,5e10.3)') nsimpit,t,p,simp_bin(1,1:2),simp_err(1)
                go to 99                
            end if          
!------------------------------------------------------------------------------------------                                
        end do !simplex loop            
!------------------------------------------------------------------------------------------                       
 99 end do !property loop
!------------------------------------------------------------------------------------------
    return
end subroutine Bin_simp
!***************************************************************************************************

!***************************************************************************************************

!***************************************************************************************************
function F_err(valin) result(valout)
    implicit none

    real(kind=DP), intent(in)   ::  valin(1:2)  
    real(kind=DP)               ::  valout, lmu(1:2), vmu(1:2)
    real(kind=DP)               ::  lvliq, lvvap
    integer                     ::  c1ferr, cixi, ciyi

    do c1ferr=1,nvint-1
        if((xi_arr(c1ferr)<valin(1)).and.(xi_arr(c1ferr+1)>=valin(1))) then
            cixi=c1ferr
         
            Comp_array(1)%xi = valin(1)
            Comp_array(2)%xi = 1.0e0_DP - Comp_array(1)%xi
       
            v=Vol_dens_g2(vl_arr(cixi),vl_arr(cixi+1))
            if(v>1.0e5_DP) go to 11
            
            lmu(1) = Mu(1)
            lmu(2) = Mu(2)
            
            go to 10
        end if
    end do
    
    print*,"Correpsonding xi value not found in F_err...",valin(1), xilim(1:2)
    failed=.true.

 11 Comp_array(1)%xi = valin(1)
    Comp_array(2)%xi = 1.0e0_DP - Comp_array(1)%xi
    v = Vol_dens_g()
    lmu(1) = Mu(1)
    lmu(2) = Mu(2)
            
 10 do c1ferr=1,nvint-1
        if((yi_arr(c1ferr)<valin(2)).and.(yi_arr(c1ferr+1)>=valin(2))) then
            ciyi=c1ferr       
          
            Comp_array(1)%xi = valin(2)
            Comp_array(2)%xi = 1.0e0_DP - Comp_array(1)%xi
         
            v=Vol_dens_g2(vv_arr(ciyi),vv_arr(ciyi+1))
            if(v>1.0e5_DP) go to 21
            
            vmu(1) = Mu(1)
            vmu(2) = Mu(2)
            
            go to 20
        end if
    end do
    
    print*,"Correpsonding yi value not found in F_err...",valin(2), yilim(1:2)
    failed=.true.
    
 21 Comp_array(1)%xi = valin(1)
    Comp_array(2)%xi = 1.0e0_DP - Comp_array(1)%xi
    v = Vol_dens_g()
    vmu(1) = Mu(1)
    vmu(2) = Mu(2)
    
 20 valout = abs(lmu(1)-vmu(1)) + abs(lmu(2)-vmu(2))

    return
end function F_err
!***************************************************************************************************
subroutine Gridder()
    implicit none
    
    integer                 ::  igrid, ilgrid, ivgrid
    real(kind=DP)           ::  grid_err 
    real(kind=DP)           ::  lmu(1:nvint,1:2), vmu(1:nvint,1:2)
   
    do igrid=1, nvint
        Comp_array(1)%xi = xi_arr(igrid)
        Comp_array(2)%xi = 1.0e0_DP - Comp_array(1)%xi
        v = vl_arr(igrid)
        lmu(igrid,1) = Mu(1)
        lmu(igrid,2) = Mu(2)
        
        Comp_array(1)%xi = yi_arr(igrid)
        Comp_array(2)%xi = 1.0e0_DP - Comp_array(1)%xi
        v = vv_arr(igrid)
        vmu(igrid,1) = Mu(1)
        vmu(igrid,2) = Mu(2)       
    end do
   
    !find 3 best options
    simp_err(1)=1.0e10_DP
    simp_err(2)=2.0e10_DP
    
    do ilgrid=1,nvint
    do ivgrid=1,nvint
    
        grid_err = abs(lmu(ilgrid,1)-vmu(ivgrid,1)) + abs(lmu(ilgrid,2)-vmu(ivgrid,2))

        !test
        if(grid_err<simp_err(1)) then
            simp_err(1)   = grid_err
            simp_bin(1,1) = xi_arr(ilgrid)
            simp_bin(1,2) = yi_arr(ivgrid)       
        end if

    end do
    end do

    return
end subroutine Gridder 
!***************************************************************************************************
subroutine bin_order()
    implicit none
    
    integer         ::  c1_order
    logical         ::  lordered
    real(kind=DP)   ::  tr_order
    
    lordered=.false.
    
    do while(.not.lordered)
        lordered=.true.
    do c1_order=1,2
        if(simp_err(c1_order)>simp_err(c1_order+1)) then
            lordered=.false.
            
            tr_order             = simp_err(c1_order)    
            simp_err(c1_order)   = simp_err(c1_order+1)
            simp_err(c1_order+1) = tr_order
           
            tr_order               = simp_bin(c1_order,1)  
            simp_bin(c1_order,1)   = simp_bin(c1_order+1,1) 
            simp_bin(c1_order+1,1) = tr_order
            
            tr_order               = simp_bin(c1_order,2)   
            simp_bin(c1_order,2)   = simp_bin(c1_order+1,2) 
            simp_bin(c1_order+1,2) = tr_order
        end if
    end do
    end do
    
    return
end subroutine bin_order
!***************************************************************************************************
subroutine V_interp()
    implicit none
    
    integer ::  c1vint
    
    if(allocated(xi_arr)) deallocate(xi_arr, vl_arr, yi_arr, vv_arr)
 
    allocate(xi_arr(1:nvint), vl_arr(1:nvint), yi_arr(1:nvint), vv_arr(1:nvint))

    delxi = (xilim(2) - xilim(1))/real(nvint-1) 
    delyi = (yilim(2) - yilim(1))/real(nvint-1)

    do c1vint=1,nvint
        Comp_array(1)%xi = xilim(1) + real(c1vint-1)*delxi
        Comp_array(2)%xi = 1.0e0_DP - Comp_array(1)%xi
        xi_arr(c1vint)    = Comp_array(1)%xi
        vl_arr(c1vint)    = Vol_dens_g()
        
        Comp_array(1)%xi = yilim(1) + real(c1vint-1)*delyi
        Comp_array(2)%xi = 1.0e0_DP - Comp_array(1)%xi
        yi_arr(c1vint)    = Comp_array(1)%xi
        vv_arr(c1vint)    = Vol_dens_g()
    end do  

    return
end subroutine V_interp
!***************************************************************************************************
subroutine Bouncey()
    implicit none
        
    simp_bin(2,1) = uni() * (xilim(2)-xilim(1)) + xilim(1)    
    simp_bin(2,2) = uni() * (yilim(2)-yilim(1)) + yilim(1)
    simp_err(2) = F_err(simp_bin(2,1:2))
    
    simp_bin(3,1) = uni() * (xilim(2)-xilim(1)) + xilim(1)    
    simp_bin(3,2) = uni() * (yilim(2)-yilim(1)) + yilim(1)
    simp_err(3) = F_err(simp_bin(3,1:2))
    
    call bin_order() 
    
    return
end subroutine Bouncey
!***************************************************************************************************
end module Phase_simplex_mod
!*************************************************************************************************** 
!***************************************************************************************************