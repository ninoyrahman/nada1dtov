  SUBROUTINE Gauge_dyn
  USE tmp_mdparam
  USE tmp_mdgrid
  USE tmp_mdmetric
  USE MD_Boundary

  IMPLICIT NONE

  INTEGER :: i





    DO i=imin-ghzx,imax+ghzx

      small_a(i) = 1.d0

      small_b(i) = 1.d0

      betaux(i) = 0.d0

      A_a(i) = 0.d0

      alp(i) = 1.d0 + (0.01d0*xh(i)*xh(i)*(exp(-(xh(i)-5.0)**2/1.) + exp(-(xh(i)+5.0)**2/1.)))/ & 
             & (1.d0 + xh(i)*xh(i))

      xsi(i) = 0.d0

      trK(i) = 0.d0

      delta_x(i) = 0.d0


      rho_matter(i) = 0.d0

      S_a(i) = 0.d0

      S_b(i) = 0.d0

      j_r(i)  = 0.d0

      capB(i) = 0.d0

      RHSdeltar(i) = 0.d0

      xsi2(i) = exp(-2.d0*xsi(i))


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



  write(*,*) 'INITIAL LAPSE------>>>>',alp(1)


END SUBROUTINE Gauge_dyn

















