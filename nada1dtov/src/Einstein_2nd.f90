  SUBROUTINE Einstein_2nd(RHS2,RHS3)
  USE tmp_mdparam
  USE tmp_mdgrid
  USE tmp_mdmetric


  IMPLICIT NONE

  REAL(KIND=double), DIMENSION(1-ghzx:nx+ghzx),INTENT(INOUT) :: RHS2,RHS3

  REAL(KIND=double), PARAMETER :: onesixth = 1.d0/6.d0
  REAL(KIND=double), PARAMETER :: onethird = 1.d0/3.d0 

  INTEGER :: i
  REAL(KIND=double) :: fac
  REAL(KIND=double), DIMENSION(1-ghzx:nx+ghzx) :: Rrr,R

  REAL(KIND=double) :: KO_disp

    do i=imin,imax

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




!!$       RHS2(i) = alp(i)*(Rrr(i)-onethird*R(i))
!!$
!!$


       KO_disp = eps1/16.*(A_a(i-2)-4.d0*A_a(i-1)+6.d0*A_a(i)-4.d0*A_a(i+1)+A_a(i+2))

!!$
!!$       RHS3(i) = betaux(i)*dxdA_a(i)-(Cov_alp(i)-onethird*laplacian_alp(i)) & 
!!$                & +alp(i)*trK(i)*A_a(i) !+ KO_disp
!!$


       RHS2(i) = -(Cov_alp(i)-onethird*laplacian_alp(i))+alp(i)*(Rrr(i)-onethird*R(i))


!!$       RHS3(i) = betaux(i)*dxdA_a(i) +alp(i)*trK(i)*A_a(i) & 
!!$               &  -16.d0*3.141592653589793238462d0*alp(i)*(S_a(i) - S_b(i))+ KO_disp



       if (betaux(i).ge.0.d0) then 

       RHS3(i) = betaux(i)*inv12dx*(A_a(i+3)-6.d0*A_a(i+2)+18.d0*A_a(i+1)- & 
                    & 10.d0*A_a(i)-3.d0*A_a(i-1))+alp(i)*trK(i)*A_a(i) & 
               &  -16.d0*3.141592653589793238462d0*alp(i)*(S_a(i) - S_b(i))+ KO_disp


       else 


       RHS3(i) = betaux(i)*inv12dx*(-A_a(i-3)+6.d0*A_a(i-2)-18.d0*A_a(i-1)+ & 
                    & 10.d0*A_a(i)+3.d0*A_a(i+1))+alp(i)*trK(i)*A_a(i) & 
               &  -16.d0*3.141592653589793238462d0*alp(i)*(S_a(i) - S_b(i))+ KO_disp


       end if 


!!$       RHS3(i) = 0.d0
!!$       RHS2(i) = 0.d0
 
 end do







END SUBROUTINE Einstein_2nd
