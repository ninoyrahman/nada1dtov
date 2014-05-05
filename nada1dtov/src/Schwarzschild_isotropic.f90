  SUBROUTINE Schwarzschild_isotropic
  USE tmp_mdparam
  USE tmp_mdgrid
  USE tmp_mdmetric
  USE MD_Boundary

  IMPLICIT NONE

  INTEGER :: i




    DO i=imin-ghzx,imax+ghzx
    
       grr(i) = (1.d0 + 1.d0/(2.d0*xh(i)))**4

      small_a(i) = 1.d0

      small_b(i) = 1.d0

      betaux(i) = 0.d0

      A_a(i) = 0.d0

      psi(i) = (1.d0 + 1.d0/(2.d0*xh(i)))

      alp(i) = psi(i)**(-2) !puncture gauge

!      alp(i) = 1.d0

      xsi(i) = log(1.d0 + 1.d0/(2.d0*xh(i)))

      trK(i) = 0.d0

      delta_x(i) = 0.d0

      rho_matter(i) = 0.d0

      S_a(i) = 0.d0

      S_b(i) = 0.d0

      j_r(i)  = 0.d0

      capB(i) = 0.d0

      RHSdeltar(i) = 0.d0

      xsi2(i) = exp(-2.d0*xsi(i))

! quasi isotropic coordinates:

!      alp(i) = (1.d0 - 1.d0/(2.d0*xh(i)))/(1.d0 + 1.d0/(2.d0*xh(i)))
 
!!$      if (xh(i).lt.0.2d0) then 
!!$
!!$      alp(i) = 0.005d0 
!!$
!!$      xsi(i) = xsi(i) - xsi(i)*0.1
!!$      xsi2(i) = exp(-2.d0*xsi(i))

!!$      small_a(i) = 1.d0
!!$
!!$      small_b(i) = 1.d0
!!$
!!$      betaux(i) = 0.d0
!!$
!!$      A_a(i) = 0.d0
!!$
!!$      psi(i) = (1.d0 + 1.d0/(2.d0*0.51d0))
!!$
!!$      grr(i) = (1.d0 + 1.d0/(2.d0*0.51d0))**4
!!$
!!$!      xsi(i) = log(1.d0 + 1.d0/(2.d0*0.51d0))
!!$
!!$      trK(i) = 0.d0
!!$
!!$      delta_x(i) = 0.d0
!!$
!!$      rho_matter(i) = 0.d0
!!$
!!$      S_a(i) = 0.d0
!!$
!!$      S_b(i) = 0.d0
!!$
!!$      j_r(i)  = 0.d0
!!$
!!$      capB(i) = 0.d0
!!$
!!$      RHSdeltar(i) = 0.d0
!!$
!!$      xsi2(i) = exp(-2.d0*xsi(i))


 !      else if ((xh(i).gt.0.2d0).and.(xh(i).lt.0.5d0)) then 

!!$     if (xh(i).lt.0.5d0) then 
!!$
!!$      xsi(i) = xsi(i) !- xsi(i)*0.1
!!$      xsi2(i) = exp(-2.d0*xsi(i))
!!$      small_a(i) = 1.d0
!!$
!!$      small_b(i) = 1.d0
!!$
!!$      betaux(i) = 0.d0
!!$
!!$      A_a(i) = 0.d0
!!$
!!$!      psi(i) = (1.d0 + 1.d0/(2.d0*0.51d0))
!!$!      grr(i) = (1.d0 + 1.d0/(2.d0*0.51d0))**4
!!$
!!$       psi(i) = (1.d0 + 1.d0/(2.d0*xh(i)))
!!$       grr(i) = (1.d0 + 1.d0/(2.d0*xh(i)))**4
!!$
!!$!      xsi(i) = log(1.d0 + 1.d0/(2.d0*0.51d0))
!!$
!!$      trK(i) = 0.d0
!!$
!!$      delta_x(i) = 0.d0
!!$
!!$      rho_matter(i) = 0.d0
!!$
!!$      S_a(i) = 0.d0
!!$
!!$      S_b(i) = 0.d0
!!$
!!$      j_r(i)  = 0.d0
!!$
!!$      capB(i) = 0.d0
!!$
!!$      RHSdeltar(i) = 0.d0
!!$
!!$      xsi2(i) = exp(-2.d0*xsi(i))
!!$
!!$       alp(i) = 0.005d0 
!!$!      psi(i) = (1.d0 + 1.d0/(2.d0*xh(i)))
!!$
!!$!      alp(i) = psi(i)**(-2) !puncture gauge
!!$
!!$
!!$      end if

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


  END SUBROUTINE Schwarzschild_isotropic

















