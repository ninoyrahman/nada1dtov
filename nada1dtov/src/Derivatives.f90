  SUBROUTINE Derivatives
  USE tmp_mdparam
  USE tmp_mdgrid
  USE tmp_mdmetric

  IMPLICIT NONE


  INTEGER :: i,k,q
  integer :: coun1,coun2,coun_rate
  REAL(KIND=double) ::timelapse


  REAL(KIND=double):: det,invdet
  REAL(KIND=double), PARAMETER :: onethird = 1.d0/3.d0


!!$ call SYSTEM_CLOCK(coun1,coun_rate)
!!$   call SYSTEM_CLOCK(coun2,coun_rate)
!!$   timelapse = dfloat(coun2-coun1)/dfloat(coun_rate)
!!$   write(*,*) 'time step1', timelapse
!!$  


      xsi2 = exp(-2.d0*xsi)
  call Derivatives_func_12(xsi2,dxdxsi2,dxxdxsi2)


!  call Derivatives_func_12(xsi,dxdxsi,dxxdxsi)

  dxdxsi = -dxdxsi2/(2.d0*xsi2)

  dxxdxsi = -( (dxxdxsi2*xsi2 - dxdxsi2**2)/(2.d0*xsi2*xsi2))  

  dxdxsi_x2 = -dxdxsi2*xsi2/(2.d0)

  dxxdxsi_x2 = -(dxxdxsi2*xsi2 - dxdxsi2**2)/(2.d0)  




  call Derivatives_func_12(alp,dxdalp,dxxdalp)

  call Derivatives_func_12(betaux,dxdbetaux,dxxdbetaux)

  dxdbetaux_x = (dxdbetaux*xh-betaux)/(xh*xh)


  call Derivatives_func_12(small_a,dxdsmall_a,dxxdsmall_a)

  call Derivatives_func_12(small_b,dxdsmall_b,dxxdsmall_b)



  call Derivatives_func_1(delta_x,dxddelta_x)

  call Derivatives_func_1(trK,dxdtrK)

  call Derivatives_func_1(A_a,dxdA_a)

!!$  do i =1,imax
!!$  write(*,*) grr(i) - exp(4.d0*xsi(i))*small_a(i)
!!$  enddo
!!$  stop

  grr = exp(4.d0*xsi)*small_a
  gtt = exp(4.d0*xsi)*small_b*xh**2
  gpp = exp(4.d0*xsi)*small_b*xh**2

  urr = 1.d0/grr
  utt = 1.d0/gtt
  upp = 1.d0/gpp


  krr  = exp(4.d0*xsi)*small_a*A_a+grr*trK/3.d0
  ktt  = -0.5d0*exp(4.d0*xsi)*small_b*xh*xh*A_a+gtt*trK/3.d0
  kpp  = ktt

  sqdetg = sqrt(grr*gtt*gpp)


  call Derivatives_func_1(grr,dxdgrr)
  call Derivatives_func_1(gtt,dxdgtt)
  call Derivatives_func_1(gpp,dxdgpp)



  do i=imin,imax

  dxdab2(i) = small_b(i)*small_b(i)*dxdsmall_a(i)+ small_a(i)*2.d0*small_b(i)*dxdsmall_b(i)



 !laplacian of the lapse:

  laplacian_alp(i) = (1.d0/(small_a(i)*exp(4.d0*xsi(i))))*(dxxdalp(i) - dxdalp(i)*(dxdsmall_a(i)/(2.d0*small_a(i))- & 
                & dxdsmall_b(i)/small_b(i) - 2.d0*dxdxsi(i) - 2.d0/xh(i)))


!  laplacian_alp(i) = xsi2(i)**2/(small_a(i))*(dxxdalp(i) - dxdalp(i)*(dxdsmall_a(i)/(2.d0*small_a(i))- & 
!                & dxdsmall_b(i)/small_b(i) - 2.d0*dxdxsi(i) - 2.d0/xh(i)))


!  laplacian_alp(i) = xsi2(i)**2/(small_a(i))*(dxxdalp(i) - dxdalp(i)*(dxdsmall_a(i)/(2.d0*small_a(i))- & 
!                & dxdsmall_b(i)/small_b(i) - 2.d0/xh(i))) + 2.d0*dxdalp(i)*dxdxsi_x2(i)/small_a(i)



 !Conformal divergence of the shift:

! Confdiv_beta(i) = dxdbetaux(i) + betaux(i)*(dxdab2(i)/(2.d0*small_a(i)*small_b(i)*small_b(i)) +2.d0/xh(i))

 Confdiv_beta(i) = dxdbetaux(i) + betaux(i)*(dxdsmall_a(i)/(2.d0*small_a(i)) + dxdsmall_b(i)/small_b(i) +2.d0/xh(i))



! dxdConfdiv_beta(i) = dxxdbetaux(i) + dxdbetaux(i)*(dxdsmall_a(i)/(2.d0*small_a(i))+dxdsmall_b(i)/small_b(i) + 2.d0/xh(i)) + & 
!               &    betaux(i)*( (dxxdsmall_a(i) - dxdsmall_a(i))/(2.d0*small_a(i)) + & 
!               &    (small_b(i)*dxxdsmall_b(i)-dxdsmall_b(i)*dxdsmall_b(i))/(small_b(i)*small_b(i)) - 2.d0/(xh(i)*xh(i)))

 dxdConfdiv_beta(i) = dxxdbetaux(i) + dxdbetaux(i)*(dxdsmall_a(i)/(2.d0*small_a(i))+dxdsmall_b(i)/small_b(i) + 2.d0/xh(i)) + & 
               &    betaux(i)*( (small_a(i)*dxxdsmall_a(i) - dxdsmall_a(i)**2)/(2.d0*small_a(i)**2) + & 
               &    (small_b(i)*dxxdsmall_b(i)-dxdsmall_b(i)*dxdsmall_b(i))/(small_b(i)*small_b(i)) - 2.d0/(xh(i)*xh(i)))


 !Cov_alp:

! Cov_alp(i) = xsi2(i)**2/(small_a(i))*(dxxdalp(i) - dxdalp(i)*(dxdsmall_a(i)/(2.d0*small_a(i))+ & 
!                &  2.d0*dxdxsi(i) ))

! Cov_alp(i) = (1.d0/(small_a(i)*exp(4.d0*xsi(i))))*(dxxdalp(i) - dxdalp(i)*(dxdsmall_a(i)/(2.d0*small_a(i))+ & 
!                &  2.d0*dxdxsi(i) ))


 Cov_alp(i) = xsi2(i)**2/(small_a(i))*(dxxdalp(i) - dxdalp(i)*(dxdsmall_a(i)/(2.d0*small_a(i)))) & 
            & - 2.D0*dxdalp(i)*dxdxsi_x2(i)/small_a(i)


 end do

! write(*,*) laplacian_alp(1), laplacian_alp(2) 

!!$
!!$  OPEN(UNIT=50, FILE='data/laplacianalp.dat',STATUS='UNKNOWN',ACTION= 'WRITE',POSITION='append')
!!$
!!$100 FORMAT (2E22.11)
!!$110 FORMAT (A9,E22.11)
!!$
!!$ !  call Derivatives_func_12(xsi,dxdxsi,dxxdxsi)
!!$
!!$     write(50,*)''
!!$     write(50,110)'"time = ',time
!!$    DO i=imin,imax
!!$       write(50,*) xh(i),laplacian_alp(i),(dxxdalp(i) - dxdalp(i)*(- 2.d0*dxdxsi(i) ))!- 2.d0/xh(i)))
!!$
!!$    ENDDO 
!!$
!!$
!!$
!!$    stop

END SUBROUTINE Derivatives


