  SUBROUTINE Einstein_2nd_trK(RHS2,RHS3)
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


!       RHS2(i) =  -laplacian_alp(i)+alp(i)*(A_a(i)*A_a(i)+0.5*A_a(i)*A_a(i))

       RHS2(i) =  -laplacian_alp(i)


       KO_disp = eps1/16.*(trK(i-2)-4.d0*trK(i-1)+6.d0*trK(i)-4.d0*trK(i+1)+trK(i+2))




!!$       if (betaux(i).ge.0.d0) then 
!!$
!!$       RHS3(i) = betaux(i)*inv12dx*(trK(i+3)-6.d0*trK(i+2)+18.d0*trK(i+1)- & 
!!$                    & 10.d0*trK(i)-3.d0*trK(i-1)) + alp(i)*(onethird*trK(i)*trK(i)) + &
!!$                &  4.d0*3.141592653589793238462d0*alp(i)*(rho_matter(i) + S_a(i) + 2.d0*S_b(i)) + KO_disp
!!$
!!$       else 
!!$
!!$
!!$       RHS3(i) = betaux(i)*inv12dx*(-trK(i-3)+6.d0*trK(i-2)-18.d0*trK(i-1)+ & 
!!$                    & 10.d0*trK(i)+3.d0*trK(i+1)) + alp(i)*(onethird*trK(i)*trK(i)) + &
!!$                &  4.d0*3.141592653589793238462d0*alp(i)*(rho_matter(i) + S_a(i) + 2.d0*S_b(i)) + KO_disp
!!$
!!$
!!$       end if 



       if (betaux(i).ge.0.d0) then 

       RHS3(i) = betaux(i)*inv12dx*(trK(i+3)-6.d0*trK(i+2)+18.d0*trK(i+1)- & 
                    & 10.d0*trK(i)-3.d0*trK(i-1)) + alp(i)*(onethird*trK(i)*trK(i)) +alp(i)*(A_a(i)*A_a(i)+0.5*A_a(i)*A_a(i)) + &
                &  4.d0*3.141592653589793238462d0*alp(i)*(rho_matter(i) + S_a(i) + 2.d0*S_b(i)) + KO_disp

       else 


       RHS3(i) = betaux(i)*inv12dx*(-trK(i-3)+6.d0*trK(i-2)-18.d0*trK(i-1)+ & 
                    & 10.d0*trK(i)+3.d0*trK(i+1)) + alp(i)*(onethird*trK(i)*trK(i)) +alp(i)*(A_a(i)*A_a(i)+0.5*A_a(i)*A_a(i)) + &
                &  4.d0*3.141592653589793238462d0*alp(i)*(rho_matter(i) + S_a(i) + 2.d0*S_b(i)) + KO_disp


       end if 

 

 end do



! write(*,*) RHS2(1), RHS3(1)



END SUBROUTINE Einstein_2nd_trK
