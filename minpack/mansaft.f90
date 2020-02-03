program pha1
    use Types           ! Definitions of types and double precision
    use Global          ! Important global parameters
    Use Pressure
    use Input
    use Vol
    use ChemPot
    use solver
      implicit none
      integer             ::  i
      real(kind=DP)      :: values(3), mu_1, mu_2
      call Read_input( )
        select case(properties%type)
          case('vle','VLE')

                print*, "V-L Phase equilibrium"
                print*, "====================="
                print*, " "
                print*, "  T(K)     P(kPa)    x1     y1     V_L          V_V        "  
                print*, " "
                
                    do i = 1, properties%n                    
                          t = properties%t(i)
                          p = properties%p(i)*1000
                          x_1 = properties%xi(i, 1)
                          Comp_array(1)%xi = x_1
                          Comp_array(2)%xi = 1.0 - x_1
                          P = 1.5*P
                          v = Vol_dens_g( )
                          mu_1 = Mu_res(1)
                          mu_2 = Mu_res(2)
                          x_init_1 = v
                          x_init_2 = 8.314*t/(properties%p(i)*1000)
                          x_init_3 = exp(mu_1/(8.314*t))*x_1/(exp(mu_1/(8.314*t))*x_1+exp(mu_2/(8.314*t))*(1-x_1))
                                                        
                          call solve_nle(values)
                          write(*,100) t, (P_L+P_V)/2000, x_1, values(3),values(1), values(2) 
                    end do
          end select           
      stop
100 Format ('  ',F7.2,'  ',F8.2,'  ',F5.3,'  ',F5.3,'  'E11.4,'  ',E11.4)  
end
