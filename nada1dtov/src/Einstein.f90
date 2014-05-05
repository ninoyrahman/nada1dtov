  SUBROUTINE Einstein(RHS)
  USE tmp_mdparam
  USE tmp_mdgrid
  USE tmp_mdmetric


  IMPLICIT NONE
 
  REAL(KIND=double), DIMENSION(1:12,1-ghzx:nx+ghzx),INTENT(INOUT) :: RHS


  REAL(KIND=double), PARAMETER :: onesixth = 1.d0/6.d0
  REAL(KIND=double), PARAMETER :: onethird = 1.d0/3.d0 

  INTEGER :: i
  REAL(KIND=double) :: fac,KO_disp
  REAL(KIND=double), DIMENSION(1-ghzx:nx+ghzx) :: Rrr,R



 
!        adpx(i,k) = betaux(i,k)*inv12dx* (phi(i+3,k)-6.d0*phi(i+2,k)+18.d0*phi(i+1,k)- & 
!                    & 10.d0*phi(i,k)-3.d0*phi(i-1,k))



    do i=imin,imax


       if (betaux(i).ge.0.d0) then 

!Eq. 6.7: xsi

       KO_disp = eps1/16.*(xsi(i-2)-4.d0*xsi(i-1)+6.d0*xsi(i)-4.d0*xsi(i+1)+xsi(i+2))

       RHS(1,i) = betaux(i)*inv12dx*(xsi(i+3)-6.d0*xsi(i+2)+18.d0*xsi(i+1)- & 
                    & 10.d0*xsi(i)-3.d0*xsi(i-1)) + onesixth*Confdiv_beta(i) - onesixth*alp(i)*trK(i) + KO_disp



!Eq. 6.7: xsi2


!!$
!!$       KO_disp = eps1/16.*(xsi2(i-2)-4.d0*xsi2(i-1)+6.d0*xsi2(i)-4.d0*xsi2(i+1)+xsi2(i+2))
!!$
!!$
!!$       RHS(1,i) = betaux(i)*inv12dx*(xsi2(i+3)-6.d0*xsi2(i+2)+18.d0*xsi2(i+1)- & 
!!$                    & 10.d0*xsi2(i)-3.d0*xsi2(i-1)) + onethird*xsi2(i)*(-Confdiv_beta(i) + alp(i)*trK(i)) + KO_disp
!!$



!Eq. 6.8: small_a

       KO_disp = eps1/16.*(small_a(i-2)-4.d0*small_a(i-1)+6.d0*small_a(i)-4.d0*small_a(i+1)+small_a(i+2))

!       RHS(2,i) = betaux(i)*dxdsmall_a(i)+2.d0*small_a(i)*dxdbetaux(i)- & 
!                & twothird*small_a(i)*Confdiv_beta(i)-2.d0*alp(i)*small_a(i)*A_a(i) + KO_disp

       RHS(2,i) = betaux(i)*inv12dx*(small_a(i+3)-6.d0*small_a(i+2)+18.d0*small_a(i+1)- & 
                    & 10.d0*small_a(i)-3.d0*small_a(i-1)) +2.d0*small_a(i)*dxdbetaux(i)- & 
                & twothird*small_a(i)*Confdiv_beta(i)-2.d0*alp(i)*small_a(i)*A_a(i) + KO_disp
 

!Eq. 6.9: small_b

       KO_disp = eps1/16.*(small_b(i-2)-4.d0*small_b(i-1)+6.d0*small_b(i)-4.d0*small_b(i+1)+small_b(i+2))

!       RHS(3,i) = betaux(i)*dxdsmall_b(i)+2.d0*small_b(i)*betaux(i)/xh(i)- & 
!                & twothird*small_b(i)*Confdiv_beta(i)+alp(i)*small_b(i)*A_a(i) + KO_disp

       RHS(3,i) = betaux(i)*inv12dx*(small_b(i+3)-6.d0*small_b(i+2)+18.d0*small_b(i+1)- & 
                    & 10.d0*small_b(i)-3.d0*small_b(i-1))+2.d0*small_b(i)*betaux(i)/xh(i)- & 
                & twothird*small_b(i)*Confdiv_beta(i)+alp(i)*small_b(i)*A_a(i) + KO_disp
 


!Eq. 6.13: trK

       KO_disp = eps1/16.*(trK(i-2)-4.d0*trK(i-1)+6.d0*trK(i)-4.d0*trK(i+1)+trK(i+2))

!!$       RHS(4,i) = betaux(i)*dxdtrK(i)-laplacian_alp(i)+alp(i)*(A_a(i)*A_a(i)+0.5*A_a(i)*A_a(i)+onethird*trK(i)*trK(i)) & 
!!$!                & + 4.d0*3.141592653589793238462d0*alp(i)*(rho_matter(i) + SS(i))
!!$                & + 4.d0*3.141592653589793238462d0*alp(i)*(rho_matter(i) + S_a(i) + 2.d0*S_b(i)) + KO_disp

       RHS(4,i) = betaux(i)*inv12dx*(trK(i+3)-6.d0*trK(i+2)+18.d0*trK(i+1)- & 
                & 10.d0*trK(i)-3.d0*trK(i-1)) -laplacian_alp(i)+alp(i)*(A_a(i)*A_a(i)+0.5*A_a(i)*A_a(i)+onethird*trK(i)*trK(i))+ & 
                & 4.d0*3.14159265358979323d0*alp(i)*(rho_matter(i) + S_a(i) + 2.d0*S_b(i)) + KO_disp

 
!Eq. 6.18:
!!$
!!$       fac = -1.d0/(small_a(i)*exp(4.d0*xsi(i)))

       fac = -xsi2(i)**2/(small_a(i))

!!$       Rrr(i) = fac*(dxxdsmall_a(i)/(2.d0*small_a(i))-small_a(i)*dxddelta_x(i)-0.75d0*(dxdsmall_a(i)/small_a(i))**2 + & 
!!$              & 0.5d0*(dxdsmall_b(i)/small_b(i))**2 - 0.5d0*delta_x(i)*dxdsmall_a(i)+dxdsmall_a(i)/(xh(i)*small_b(i)) + & 
!!$              & (2.d0*(1.d0-small_a(i)/small_b(i))/(xh(i)*xh(i)))*(1.d0+xh(i)*dxdsmall_b(i)/small_b(i)) & 
!!$              & + 4.d0*dxxdxsi(i) - 2.d0*dxdxsi(i)*(dxdsmall_a(i)/small_a(i)- dxdsmall_b(i)/small_b(i) - 2.d0/xh(i)))


       Rrr(i) = fac*(dxxdsmall_a(i)/(2.d0*small_a(i))-small_a(i)*dxddelta_x(i)-0.75d0*(dxdsmall_a(i)/small_a(i))**2 + & 
              & 0.5d0*(dxdsmall_b(i)/small_b(i))**2 - 0.5d0*delta_x(i)*dxdsmall_a(i)+dxdsmall_a(i)/(xh(i)*small_b(i)) + & 
              & (2.d0*(1.d0-small_a(i)/small_b(i))/(xh(i)*xh(i)))*(1.d0+xh(i)*dxdsmall_b(i)/small_b(i))) - (1.d0/small_a(i))* & 
              & (4.d0*dxxdxsi_x2(i) - 2.d0*dxdxsi_x2(i)*(dxdsmall_a(i)/small_a(i)- dxdsmall_b(i)/small_b(i) - 2.d0/xh(i)))



!Eq. 6.19:

!!$       R(i) = fac*(dxxdsmall_a(i)/(2.d0*small_a(i))+dxxdsmall_b(i)/small_b(i)-small_a(i)*dxddelta_x(i)- & 
!!$            & (dxdsmall_a(i)/small_a(i))*(dxdsmall_a(i)/small_a(i))+0.5d0*(dxdsmall_b(i)/small_b(i))* & 
!!$            & (dxdsmall_b(i)/small_b(i))+2.d0/(xh(i)*small_b(i))*(3.d0-small_a(i)/small_b(i))*dxdsmall_b(i)+ & 
!!$            & 4.d0/(xh(i)*xh(i))*(1.d0-small_a(i)/small_b(i))+8.d0*(dxxdxsi(i)+dxdxsi(i)**2)- & 
!!$            & 8.d0*dxdxsi(i)*(dxdsmall_a(i)/(2.d0*small_a(i))-dxdsmall_b(i)/(small_b(i))-2.d0/xh(i)))


       R(i) = fac*(dxxdsmall_a(i)/(2.d0*small_a(i))+dxxdsmall_b(i)/small_b(i)-small_a(i)*dxddelta_x(i)- & 
            & (dxdsmall_a(i)/small_a(i))*(dxdsmall_a(i)/small_a(i))+0.5d0*(dxdsmall_b(i)/small_b(i))* & 
            & (dxdsmall_b(i)/small_b(i))+2.d0/(xh(i)*small_b(i))*(3.d0-small_a(i)/small_b(i))*dxdsmall_b(i)+ & 
            & 4.d0/(xh(i)*xh(i))*(1.d0-small_a(i)/small_b(i))) - (1.d0/small_a(i))*8.d0*(dxxdxsi_x2(i)+dxdxsi(i)*dxdxsi_x2(i)) + & 
            & (1.d0/small_a(i))*8.d0*dxdxsi_x2(i)*(dxdsmall_a(i)/(2.d0*small_a(i))-dxdsmall_b(i)/(small_b(i))-2.d0/xh(i))



!!$       write(*,*)'R',R(i),Rrr(i)
!!$
!!$       stop

!!$  OPEN(UNIT=50, FILE='data/laplacianalp.dat',STATUS='UNKNOWN',ACTION= 'WRITE',POSITION='append')
!!$
!!$100 FORMAT (2E22.11)
!!$110 FORMAT (A9,E22.11)
!!$
!!$
!!$     write(50,*)''
!!$     write(50,110)'"time = ',time
!!$    DO i=imin,imax
!!$       write(50,100) xh(i),laplacian_alp(i)
!!$    ENDDO 
!!$


!Eq. 6.16: A_a

       KO_disp = eps1/16.*(A_a(i-2)-4.d0*A_a(i-1)+6.d0*A_a(i)-4.d0*A_a(i+1)+A_a(i+2))

!!$       RHS(5,i) = betaux(i)*dxdA_a(i)-(Cov_alp(i)-onethird*laplacian_alp(i))+alp(i)*(Rrr(i)-onethird*R(i)) & 
!!$                & +alp(i)*trK(i)*A_a(i) & 
!!$                & -16.d0*3.141592653589793238462d0*alp(i)*(S_a(i) - S_b(i))+ KO_disp


       RHS(5,i) = betaux(i)*inv12dx*(A_a(i+3)-6.d0*A_a(i+2)+18.d0*A_a(i+1)- & 
                    & 10.d0*A_a(i)-3.d0*A_a(i-1)) -(Cov_alp(i)-onethird*laplacian_alp(i))+alp(i)*(Rrr(i)-onethird*R(i)) & 
                & +alp(i)*trK(i)*A_a(i) & 
                & -16.d0*3.141592653589793238462d0*alp(i)*(S_a(i) - S_b(i))+ KO_disp


!       RHS(5,i) = 0.d0 
!Eq. 6.20: delta_x

       KO_disp = eps1/16.*(delta_x(i-2)-4.d0*delta_x(i-1)+6.d0*delta_x(i)-4.d0*delta_x(i+1)+delta_x(i+2))



!OK without matter source term
!!$
!!$       RHS(6,i) = betaux(i)*inv12dx*(delta_x(i+3)-6.d0*delta_x(i+2)+18.d0*delta_x(i+1)- & 
!!$                    & 10.d0*delta_x(i)-3.d0*delta_x(i-1))  -delta_x(i)*dxdbetaux(i)+dxxdbetaux(i)/small_a(i) + & 
!!$                & 2.d0*dxdbetaux_x(i)/small_b(i)+ onethird*(dxdConfdiv_beta(i)/small_a(i) + & 
!!$                & 2.d0*delta_x(i)*Confdiv_beta(i))- 2.d0/small_a(i)*(A_a(i)*dxdalp(i)+alp(i)*dxdA_a(i)) + & 
!!$                & 2.d0*alp(i)*(A_a(i)*delta_x(i)-2.d0/(xh(i)*small_b(i))*(A_a(i)+0.5d0*A_a(i))) + & 
!!$                & 2.d0*alp(i)/small_a(i)*(dxdA_a(i)-twothird*dxdtrK(i)+6.d0*A_a(i)*dxdxsi(i) + &
!!$                & (A_a(i)+0.5d0*A_a(i))*(2.d0/xh(i)+dxdsmall_b(i)/small_b(i))) + KO_disp


       RHS(6,i) = betaux(i)*inv12dx*(delta_x(i+3)-6.d0*delta_x(i+2)+18.d0*delta_x(i+1)- & 
                    & 10.d0*delta_x(i)-3.d0*delta_x(i-1))  -delta_x(i)*dxdbetaux(i)+dxxdbetaux(i)/small_a(i) + & 
                & 2.d0*dxdbetaux_x(i)/small_b(i)+ onethird*(dxdConfdiv_beta(i)/small_a(i) + & 
                & 2.d0*delta_x(i)*Confdiv_beta(i))- 2.d0/small_a(i)*(A_a(i)*dxdalp(i)+alp(i)*dxdA_a(i)) + & 
                & 2.d0*alp(i)*(A_a(i)*delta_x(i)-2.d0/(xh(i)*small_b(i))*(A_a(i)+0.5d0*A_a(i))) + & 
                & 2.d0*alp(i)/small_a(i)*(dxdA_a(i)-twothird*dxdtrK(i)+6.d0*A_a(i)*dxdxsi(i) + &
                & (A_a(i)+0.5d0*A_a(i))*(2.d0/xh(i)+dxdsmall_b(i)/small_b(i))- & 
                & 8.d0*3.141592653589793238462d0*j_r(i)) + KO_disp


!       RHS(6,i) = 0.d0

!Lapse

       KO_disp = eps1/16.*(alp(i-2)-4.d0*alp(i-1)+6.d0*alp(i)-4.d0*alp(i+1)+alp(i+2))

       RHS(7,i) = -2.d0*alp(i)*trK(i) + KO_disp


!       RHS(7,i) = -2.d0*alp(i)*psi(i)*trK(i) + KO_disp

!       RHS(7,i) = -2.d0*alp(i)*trK(i) + betaux(i)*inv12dx*(alp(i+3)-6.d0*alp(i+2)+18.d0*alp(i+1)- & 
!                    & 10.d0*alp(i)-3.d0*alp(i-1)) + KO_disp

!       RHS(7,i) = -alp(i)*alp(i)*trK(i) !+ KO_disp !gauge dynamics test

!Shift:

       KO_disp = eps1/16.*(betaux(i-2)-4.d0*betaux(i-1)+6.d0*betaux(i)-4.d0*betaux(i+1)+betaux(i+2))


!        RHS(8,i) = 0.d0 !gauge dynamics test

      RHS(8,i) = capB(i) + KO_disp


      RHS(9,i) = 0.75d0*RHSdeltar(i) !-  0.1d0/1.6d0*capB(i)

!      RHS(9,i) = RHSdeltar(i)





 else if (betaux(i).lt.0.d0) then 


!        adxsix(i,k) = betaux(i,k)*inv12dx* (-xsi(i-3,k)+6.d0*xsi(i-2,k)-18.d0*xsi(i-1,k)+ & 
!                    & 10.d0*xsi(i,k)+3.d0*xsi(i+1,k))


!Eq. 6.7: xsi


       KO_disp = eps1/16.*(xsi(i-2)-4.d0*xsi(i-1)+6.d0*xsi(i)-4.d0*xsi(i+1)+xsi(i+2))

       RHS(1,i) = betaux(i)*inv12dx*(-xsi(i-3)+6.d0*xsi(i-2)-18.d0*xsi(i-1)+ & 
                    & 10.d0*xsi(i)+3.d0*xsi(i+1)) + onesixth*Confdiv_beta(i) - onesixth*alp(i)*trK(i) + KO_disp

!!$

!Eq. 6.7: xsi2



!!$       KO_disp = eps1/16.*(xsi2(i-2)-4.d0*xsi2(i-1)+6.d0*xsi2(i)-4.d0*xsi2(i+1)+xsi2(i+2))
!!$
!!$!       RHS(1,i) = betaux(i)*dxdxsi2(i) + onethird*xsi2(i)*(-Confdiv_beta(i) + alp(i)*trK(i)) + KO_disp
!!$
!!$
!!$       RHS(1,i) = betaux(i)*inv12dx*(-xsi2(i-3)+6.d0*xsi2(i-2)-18.d0*xsi2(i-1)+ & 
!!$                    & 10.d0*xsi2(i)+3.d0*xsi2(i+1)) + onethird*xsi2(i)*(-Confdiv_beta(i) + alp(i)*trK(i)) + KO_disp
!!$



!Eq. 6.8: small_a

       KO_disp = eps1/16.*(small_a(i-2)-4.d0*small_a(i-1)+6.d0*small_a(i)-4.d0*small_a(i+1)+small_a(i+2))

!       RHS(2,i) = betaux(i)*dxdsmall_a(i)+2.d0*small_a(i)*dxdbetaux(i)- & 
!                & twothird*small_a(i)*Confdiv_beta(i)-2.d0*alp(i)*small_a(i)*A_a(i) + KO_disp

       RHS(2,i) = betaux(i)*inv12dx*(-small_a(i-3)+6.d0*small_a(i-2)-18.d0*small_a(i-1)+ & 
                    & 10.d0*small_a(i)+3.d0*small_a(i+1)) +2.d0*small_a(i)*dxdbetaux(i)- & 
                & twothird*small_a(i)*Confdiv_beta(i)-2.d0*alp(i)*small_a(i)*A_a(i) + KO_disp
 

!Eq. 6.9: small_b

       KO_disp = eps1/16.*(small_b(i-2)-4.d0*small_b(i-1)+6.d0*small_b(i)-4.d0*small_b(i+1)+small_b(i+2))

!       RHS(3,i) = betaux(i)*dxdsmall_b(i)+2.d0*small_b(i)*betaux(i)/xh(i)- & 
!                & twothird*small_b(i)*Confdiv_beta(i)+alp(i)*small_b(i)*A_a(i) + KO_disp

       RHS(3,i) = betaux(i)*inv12dx*(-small_b(i-3)+6.d0*small_b(i-2)-18.d0*small_b(i-1)+ & 
                    & 10.d0*small_b(i)+3.d0*small_b(i+1))+2.d0*small_b(i)*betaux(i)/xh(i)- & 
                & twothird*small_b(i)*Confdiv_beta(i)+alp(i)*small_b(i)*A_a(i) + KO_disp
 


!Eq. 6.13: trK

       KO_disp = eps1/16.*(trK(i-2)-4.d0*trK(i-1)+6.d0*trK(i)-4.d0*trK(i+1)+trK(i+2))

!!$       RHS(4,i) = betaux(i)*dxdtrK(i)-laplacian_alp(i)+alp(i)*(A_a(i)*A_a(i)+0.5*A_a(i)*A_a(i)+onethird*trK(i)*trK(i)) & 
!!$!                & + 4.d0*3.141592653589793238462d0*alp(i)*(rho_matter(i) + SS(i))
!!$                & + 4.d0*3.141592653589793238462d0*alp(i)*(rho_matter(i) + S_a(i) + 2.d0*S_b(i)) + KO_disp

       RHS(4,i) = betaux(i)*inv12dx*(-trK(i-3)+6.d0*trK(i-2)-18.d0*trK(i-1)+ & 
             & 10.d0*trK(i)+3.d0*trK(i+1)) -laplacian_alp(i)+alp(i)*(A_a(i)*A_a(i)+0.5*A_a(i)*A_a(i)+onethird*trK(i)*trK(i))+ & 
             & 4.d0*3.14159265358979323d0*alp(i)*(rho_matter(i) + S_a(i) + 2.d0*S_b(i)) + KO_disp

 
!Eq. 6.18:

!!$       fac = -1.d0/(small_a(i)*exp(4.d0*xsi(i)))

       fac = -xsi2(i)**2/(small_a(i))
!!$
!!$       Rrr(i) = fac*(dxxdsmall_a(i)/(2.d0*small_a(i))-small_a(i)*dxddelta_x(i)-0.75d0*(dxdsmall_a(i)/small_a(i))**2 + & 
!!$              & 0.5d0*(dxdsmall_b(i)/small_b(i))**2 - 0.5d0*delta_x(i)*dxdsmall_a(i)+dxdsmall_a(i)/(xh(i)*small_b(i)) + & 
!!$              & (2.d0*(1.d0-small_a(i)/small_b(i))/(xh(i)*xh(i)))*(1.d0+xh(i)*dxdsmall_b(i)/small_b(i)) & 
!!$              & + 4.d0*dxxdxsi(i) - 2.d0*dxdxsi(i)*(dxdsmall_a(i)/small_a(i)- dxdsmall_b(i)/small_b(i) - 2.d0/xh(i)))


       Rrr(i) = fac*(dxxdsmall_a(i)/(2.d0*small_a(i))-small_a(i)*dxddelta_x(i)-0.75d0*(dxdsmall_a(i)/small_a(i))**2 + & 
              & 0.5d0*(dxdsmall_b(i)/small_b(i))**2 - 0.5d0*delta_x(i)*dxdsmall_a(i)+dxdsmall_a(i)/(xh(i)*small_b(i)) + & 
              & (2.d0*(1.d0-small_a(i)/small_b(i))/(xh(i)*xh(i)))*(1.d0+xh(i)*dxdsmall_b(i)/small_b(i))) - (1.d0/small_a(i))* & 
              & (4.d0*dxxdxsi_x2(i) - 2.d0*dxdxsi_x2(i)*(dxdsmall_a(i)/small_a(i)- dxdsmall_b(i)/small_b(i) - 2.d0/xh(i)))



!Eq. 6.19:

!!$       R(i) = fac*(dxxdsmall_a(i)/(2.d0*small_a(i))+dxxdsmall_b(i)/small_b(i)-small_a(i)*dxddelta_x(i)- & 
!!$            & (dxdsmall_a(i)/small_a(i))*(dxdsmall_a(i)/small_a(i))+0.5d0*(dxdsmall_b(i)/small_b(i))* & 
!!$            & (dxdsmall_b(i)/small_b(i))+2.d0/(xh(i)*small_b(i))*(3.d0-small_a(i)/small_b(i))*dxdsmall_b(i)+ & 
!!$            & 4.d0/(xh(i)*xh(i))*(1.d0-small_a(i)/small_b(i))+8.d0*(dxxdxsi(i)+dxdxsi(i)**2)- & 
!!$            & 8.d0*dxdxsi(i)*(dxdsmall_a(i)/(2.d0*small_a(i))-dxdsmall_b(i)/(small_b(i))-2.d0/xh(i)))


       R(i) = fac*(dxxdsmall_a(i)/(2.d0*small_a(i))+dxxdsmall_b(i)/small_b(i)-small_a(i)*dxddelta_x(i)- & 
            & (dxdsmall_a(i)/small_a(i))*(dxdsmall_a(i)/small_a(i))+0.5d0*(dxdsmall_b(i)/small_b(i))* & 
            & (dxdsmall_b(i)/small_b(i))+2.d0/(xh(i)*small_b(i))*(3.d0-small_a(i)/small_b(i))*dxdsmall_b(i)+ & 
            & 4.d0/(xh(i)*xh(i))*(1.d0-small_a(i)/small_b(i))) - (1.d0/small_a(i))*8.d0*(dxxdxsi_x2(i)+dxdxsi(i)*dxdxsi_x2(i)) + & 
            & (1.d0/small_a(i))*8.d0*dxdxsi_x2(i)*(dxdsmall_a(i)/(2.d0*small_a(i))-dxdsmall_b(i)/(small_b(i))-2.d0/xh(i))
!!$



!Eq. 6.16: A_a

       KO_disp = eps1/16.*(A_a(i-2)-4.d0*A_a(i-1)+6.d0*A_a(i)-4.d0*A_a(i+1)+A_a(i+2))

!!$       RHS(5,i) = betaux(i)*dxdA_a(i)-(Cov_alp(i)-onethird*laplacian_alp(i))+alp(i)*(Rrr(i)-onethird*R(i)) & 
!!$                & +alp(i)*trK(i)*A_a(i) & 
!!$                & -16.d0*3.141592653589793238462d0*alp(i)*(S_a(i) - S_b(i))+ KO_disp


       RHS(5,i) = betaux(i)*inv12dx*(-A_a(i-3)+6.d0*A_a(i-2)-18.d0*A_a(i-1)+ & 
                    & 10.d0*A_a(i)+3.d0*A_a(i+1)) -(Cov_alp(i)-onethird*laplacian_alp(i))+alp(i)*(Rrr(i)-onethird*R(i)) & 
                & +alp(i)*trK(i)*A_a(i) & 
                & -16.d0*3.141592653589793238462d0*alp(i)*(S_a(i) - S_b(i))+ KO_disp

 
!       RHS(5,i) = 0.d0

!Eq. 6.20: delta_x

       KO_disp = eps1/16.*(delta_x(i-2)-4.d0*delta_x(i-1)+6.d0*delta_x(i)-4.d0*delta_x(i+1)+delta_x(i+2))


!!$       RHS(6,i) = betaux(i)*dxddelta_x(i)-delta_x(i)*dxdbetaux(i)+dxxdbetaux(i)/small_a(i) + & 
!!$                & 2.d0*dxdbetaux_x(i)/small_b(i)+ onethird*(dxdConfdiv_beta(i)/small_a(i) + & 
!!$                & 2.d0*delta_x(i)*Confdiv_beta(i))- 2.d0/small_a(i)*(A_a(i)*dxdalp(i)+alp(i)*dxdA_a(i)) + & 
!!$                & 2.d0*alp(i)*(A_a(i)*delta_x(i)-2.d0/(xh(i)*small_b(i))*(A_a(i)+0.5d0*A_a(i))) + & 
!!$                & 2.d0*alp(i)/small_a(i)*(dxdA_a(i)-twothird*dxdtrK(i)+6.d0*A_a(i)*dxdxsi(i) + &
!!$                & (A_a(i)+0.5d0*A_a(i))*(2.d0/xh(i)+dxdsmall_b(i)/small_b(i))) + KO_disp
!!$

!OK without matter source term
!!$
!!$       RHS(6,i) = betaux(i)*inv12dx*(delta_x(i+3)-6.d0*delta_x(i+2)+18.d0*delta_x(i+1)- & 
!!$                    & 10.d0*delta_x(i)-3.d0*delta_x(i-1))  -delta_x(i)*dxdbetaux(i)+dxxdbetaux(i)/small_a(i) + & 
!!$                & 2.d0*dxdbetaux_x(i)/small_b(i)+ onethird*(dxdConfdiv_beta(i)/small_a(i) + & 
!!$                & 2.d0*delta_x(i)*Confdiv_beta(i))- 2.d0/small_a(i)*(A_a(i)*dxdalp(i)+alp(i)*dxdA_a(i)) + & 
!!$                & 2.d0*alp(i)*(A_a(i)*delta_x(i)-2.d0/(xh(i)*small_b(i))*(A_a(i)+0.5d0*A_a(i))) + & 
!!$                & 2.d0*alp(i)/small_a(i)*(dxdA_a(i)-twothird*dxdtrK(i)+6.d0*A_a(i)*dxdxsi(i) + &
!!$                & (A_a(i)+0.5d0*A_a(i))*(2.d0/xh(i)+dxdsmall_b(i)/small_b(i))) + KO_disp


       RHS(6,i) = betaux(i)*inv12dx*(-delta_x(i-3)+6.d0*delta_x(i-2)-18.d0*delta_x(i-1)+ & 
                    & 10.d0*delta_x(i)+3.d0*delta_x(i+1))  -delta_x(i)*dxdbetaux(i)+dxxdbetaux(i)/small_a(i) + & 
                & 2.d0*dxdbetaux_x(i)/small_b(i)+ onethird*(dxdConfdiv_beta(i)/small_a(i) + & 
                & 2.d0*delta_x(i)*Confdiv_beta(i))- 2.d0/small_a(i)*(A_a(i)*dxdalp(i)+alp(i)*dxdA_a(i)) + & 
                & 2.d0*alp(i)*(A_a(i)*delta_x(i)-2.d0/(xh(i)*small_b(i))*(A_a(i)+0.5d0*A_a(i))) + & 
                & 2.d0*alp(i)/small_a(i)*(dxdA_a(i)-twothird*dxdtrK(i)+6.d0*A_a(i)*dxdxsi(i) + &
                & (A_a(i)+0.5d0*A_a(i))*(2.d0/xh(i)+dxdsmall_b(i)/small_b(i))- & 
                & 8.d0*3.141592653589793238462d0*j_r(i)) + KO_disp


!       RHS(6,i) = 0.d0

!Lapse

       KO_disp = eps1/16.*(alp(i-2)-4.d0*alp(i-1)+6.d0*alp(i)-4.d0*alp(i+1)+alp(i+2))

       RHS(7,i) = -2.d0*alp(i)*trK(i) + KO_disp

!       RHS(7,i) = -2.d0*alp(i)*psi(i)*trK(i) + KO_disp


!!$       RHS(7,i) = -2.d0*alp(i)*trK(i) +  betaux(i)*inv12dx*(-alp(i-3)+6.d0*alp(i-2)-18.d0*alp(i-1)+ & 
!!$                    & 10.d0*alp(i)+3.d0*alp(i+1)) + KO_disp

!       RHS(7,i) = -alp(i)*alp(i)*trK(i) !+ KO_disp !gauge dynamics test

!Shift:

       KO_disp = eps1/16.*(betaux(i-2)-4.d0*betaux(i-1)+6.d0*betaux(i)-4.d0*betaux(i+1)+betaux(i+2))


!        RHS(8,i) = 0.d0 !gauge dynamics test


      RHS(8,i) = capB(i) + KO_disp



      RHS(9,i) = 0.75d0*RHSdeltar(i) !- 0.1d0/1.6d0*capB(i)
       
!      RHS(9,i) = RHSdeltar(i)


   end if 








 end do




!!$  OPEN(UNIT=50, FILE='data/Ricci.dat',STATUS='UNKNOWN',ACTION= 'WRITE',POSITION='append')
!!$
!!$100 FORMAT (2E22.11)
!!$110 FORMAT (A9,E22.11)
!!$
!!$
!!$     write(50,*)''
!!$     write(50,110)'"time = ',time
!!$    DO i=imin,imax
!!$       write(50,*) xh(i), exp(4.d0*xsi(i))*Rrr(i),R(i)
!!$
!!$    ENDDO 
!!$
!!$

!!$
!!$
!!$
!!$  OPEN(UNIT=51, FILE='data/RHS_K.dat',STATUS='UNKNOWN',ACTION= 'WRITE',POSITION='append')
!!$
!!$
!!$
!!$     write(51,*)''
!!$     write(51,110)'"time = ',time
!!$    DO i=imin,imax 
!!$       write(51,*) xh(i), RHS(4,i)-4.d0*3.14159265358979323d0*alp(i)*(rho_matter(i) + S_a(i) + 2.d0*S_b(i))
!!$
!!$    ENDDO 
!!$


!    stop





END SUBROUTINE Einstein
