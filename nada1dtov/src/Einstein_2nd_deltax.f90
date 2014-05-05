  SUBROUTINE Einstein_2nd_deltax(RHS2,RHS3)
  USE tmp_mdparam
  USE tmp_mdgrid
  USE tmp_mdmetric


  IMPLICIT NONE

  REAL(KIND=double), DIMENSION(1-ghzx:nx+ghzx),INTENT(INOUT) :: RHS2,RHS3

  REAL(KIND=double), PARAMETER :: onesixth = 1.d0/6.d0
  REAL(KIND=double), PARAMETER :: onethird = 1.d0/3.d0 

  INTEGER :: i
  REAL(KIND=double) :: fac

  REAL(KIND=double) :: KO_disp

    do i=imin,imax


!!$       RHS2(i) = - 2.d0/small_a(i)*(A_a(i)*dxdalp(i)+alp(i)*dxdA_a(i)) + & 
!!$                & 2.d0*alp(i)*(-2.d0/(xh(i)*small_b(i))*(A_a(i)+0.5d0*A_a(i))) + & 
!!$                & 2.d0*alp(i)/small_a(i)*(dxdA_a(i)-twothird*dxdtrK(i)+6.d0*A_a(i)*dxdxsi(i) + &
!!$                & (A_a(i)+0.5d0*A_a(i))*(2.d0/xh(i)+dxdsmall_b(i)/small_b(i)) -&
!!$                & 8.d0*3.141592653589793238462d0*j_r(i)) & 
!!$                & - delta_x(i)*dxdbetaux(i)+dxxdbetaux(i)/small_a(i) + & 
!!$                & 2.d0*dxdbetaux_x(i)/small_b(i)+ onethird*(dxdConfdiv_beta(i)/small_a(i) + & 
!!$                & 2.d0*delta_x(i)*Confdiv_beta(i))

       RHS2(i) = - 2.d0/small_a(i)*(A_a(i)*dxdalp(i)+alp(i)*dxdA_a(i)) + & 
                & 2.d0*alp(i)*(-2.d0/(xh(i)*small_b(i))*(A_a(i)+0.5d0*A_a(i))) + & 
                & 2.d0*alp(i)/small_a(i)*(dxdA_a(i)-twothird*dxdtrK(i)+6.d0*A_a(i)*dxdxsi(i) + &
                & (A_a(i)+0.5d0*A_a(i))*(2.d0/xh(i)+dxdsmall_b(i)/small_b(i)) -&
                & 8.d0*3.141592653589793238462d0*j_r(i)) & 
                & - delta_x(i)*dxdbetaux(i)+dxxdbetaux(i)/small_a(i) + & 
                & 2.d0*dxdbetaux_x(i)/small_b(i)+ onethird*(dxdConfdiv_beta(i)/small_a(i)) 





       KO_disp = eps1/16.*(delta_x(i-2)-4.d0*delta_x(i-1)+6.d0*delta_x(i)-4.d0*delta_x(i+1)+delta_x(i+2))

!!$       RHS3(i) = betaux(i)*dxddelta_x(i) + & 
!!$                & 2.d0*alp(i)*(A_a(i)*delta_x(i))  + KO_disp


       if (betaux(i).ge.0.d0) then 

       RHS3(i) = betaux(i)*inv12dx*(delta_x(i+3)-6.d0*delta_x(i+2)+18.d0*delta_x(i+1)- & 
                    & 10.d0*delta_x(i)-3.d0*delta_x(i-1)) + & 
                & 2.d0*alp(i)*(A_a(i)*delta_x(i)) +   2.d0*delta_x(i)*Confdiv_beta(i)/3.d0 + KO_disp


       else 

       RHS3(i) = betaux(i)*inv12dx*(-delta_x(i-3)+6.d0*delta_x(i-2)-18.d0*delta_x(i-1)+ & 
                    & 10.d0*delta_x(i)+3.d0*delta_x(i+1)) + & 
                & 2.d0*alp(i)*(A_a(i)*delta_x(i)) +   2.d0*delta_x(i)*Confdiv_beta(i)/3.d0 + KO_disp


       end if 


!       RHS2(i) = 0.d0
!       RHS3(i) = 0.d0
 
  end do







END SUBROUTINE Einstein_2nd_deltax
