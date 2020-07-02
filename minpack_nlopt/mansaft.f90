program pha1
    use Types           ! Definitions of types and double precision
    use Global          ! Important global parameters
    Use Pressure
    use Input
    use Vol
    use ChemPot
    use solver
    use pure_phase
      implicit none
      integer             ::  i,j
      real(kind=DP)      :: mu_sum, output_pha(2)
      real(kind=DP),allocatable           :: output_vle(:), mu_vle(:)
      call Read_input( )
      allocate(Mu_L(1:nctypes), Mu_V(1:nctypes), x_init(1:nctypes+1))
        select case(properties%type)
          case('vle','VLE')

                print*, "V-L Phase equilibrium"
                print*, "====================="
                print*, " "
                print*, "  T(K)     P(kPa)     V_L          V_V       y(i)    "  
                print*, " "
                
                    allocate(output_vle(1:nctypes+1), comp_x(1:nctypes), mu_vle(1:nctypes))

                    
                    do i = 1, properties%n                    
                          t = properties%t(i)
                          p = properties%p(i)*1000
                          !x_1 = properties%xi(i, 1)
                          Comp_array(1:nctypes)%xi = properties%xi(i,1:nctypes)
                          comp_x(1:nctypes) = Comp_array(1:nctypes)%xi
                          P = 1.5*P
                          v = Vol_dens_g( )
                          !mu_1 = Mu_res(1)
                          !mu_2 = Mu_res(2)
                          x_init(1) = v
                          x_init(2) = 8.314*t/(properties%p(i)*1000)
                          mu_sum = 0
                          do j=1, nctypes
                            mu_vle(j) = dexp(Mu_res(j)/(8.314*t))*comp_x(j)
                            mu_sum = mu_sum + mu_vle(j)
                          end do
                          do j=1, nctypes-1
                            x_init(j+2) = mu_vle(j)/mu_sum
                          end do                         
                          call solve_nle(output_vle)
                          write(*,100) t, (P_L+P_V)/2000, output_vle(1), output_vle(2),output_vle(3:nctypes+1)
                    end do
				
			  case('pha1','PHA1','Pha1')

                print*, "Phase equilibrium"
                print*, "====================="
                print*, " "
                print*, "  T(K)     P(kPa)    V_L          V_V        "  
                print*, " "
		
		            
                    
                    do i = 1, properties%n
                          t = properties%t(i)
                          p = properties%p(i)*1000
                          Comp_array(1)%xi = 1.0_DP
                          P = 2*P
                          v = Vol_dens_g( )
                          x_init(1) = v
                          x_init(2) = 8.314*t/(properties%p(i)*1000)
                          ! print*, mu_1, mu_2, x_init_1, x_init_2, x_init_3
                        
                          call solve_pha1(output_pha )
                          write(*,200) t, (P_L+P_V)/2000,output_pha(1), output_pha(2)
                    end do
                    
        end select           
      stop
100   Format (2X, F7.2, 2X, F8.2, 2X, E11.4, 2X, E11.4, 2X, F5.3)
200   Format (2X, F7.2, 2X, F8.2, 2X, E11.4, 2X, E11.4) 

end
