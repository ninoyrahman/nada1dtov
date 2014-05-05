  SUBROUTINE OppSnyder
  USE tmp_mdparam
  USE tmp_mdgrid
  USE tmp_mdmetric
  USE MD_Boundary

  IMPLICIT NONE

  INTEGER :: i
  REAL(KIND=double) :: Rcap0, r0, M,rho0
 REAL(KIND=double), PARAMETER :: PID=3.14159265358979d0


  M= 1.d0
  Rcap0 = 5.d0*M
  r0 = Rcap0*(1.d0-M/Rcap0 + sqrt(1.d0-2.d0*M/Rcap0))*0.5d0
  rho0 = M*3.d0/(4.d0*PID*Rcap0**3)

    DO i=imin-ghzx,imax+ghzx
    
       if (xh(i).le.r0) then 
          
          psi(i) = sqrt((1.d0+sqrt(1.d0-2.d0*M/Rcap0)*r0*Rcap0*Rcap0)/(2.d0*r0**3+M*xh(i)*xh(i)))
       else  
          psi(i) = 1.d0+M/(2.d0*xh(i))

       end if

      xsi(i) = log(phi(i))
      xsi2(i) = exp(-2.d0*xsi(i))
      small_a(i) = 1.d0
      small_b(i) = 1.d0
      betaux(i) = 0.d0
      A_a(i) = 0.d0
      trK(i) = 0.d0
      delta_x(i) = 0.d0
      capB(i) = 0.d0

      alp(i) = 1.d0
      rho_matter(i) = rho0!rhoenthalpy!rho(i)*(1. + eps(i))
      S_a(i) = 0.d0
      S_b(i) = 0.d0
      j_r(i)  = 0.d0

      lorentz(i) = 1.d0
      RHSdeltar(i) = 0.d0

      grr(i) = psi(i)**4 	
      gtt(i) = grr(i)*xh(i)**2
      gpp(i) = grr(i)*xh(i)**2
      
      det = grr(i)*gtt(i)*gpp(i)
      sqdetg(i)= sqrt(det)

      p(i) = 0.d0
      eps(i) = 0.d0
      rho(i) =rho0

      h(i) = 1.d0+eps(i)+p(i)/rho(i)

      dens(i) = rho(i)*sqdetg(i)      
      tau(i) = sqdetg(i)*(rho(i)*h(i) -p(i))-dens(i)
      srval(i) = 0.d0
      

   ENDDO

  call Boundary_zaxis_func(xsi2) 
  call Boundary_zaxis_func(xsi) 
  call Boundary_zaxis_func(small_a) 
  call Boundary_zaxis_func(small_b) 
  call Boundary_zaxis_func(trK) 
  call Boundary_zaxis_func(A_a) 
  call Boundary_zaxis_asymfunc(delta_x) 
  call Boundary_zaxis_func(alp)
  call Boundary_zaxis_asymfunc(betaux) 

  call Boundary_zaxis_func(rho) 
  call Boundary_zaxis_func(eps) 
  call Boundary_zaxis_func(p) 
  call Boundary_zaxis_func(lorentz) 
  call Boundary_zaxis_asymfunc(velr) 
  call Boundary_zaxis_asymfunc(vr) 
  



END SUBROUTINE OppSnyder

















