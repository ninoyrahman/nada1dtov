  SUBROUTINE Einstein_2nd_beta(RHS2,RHS3)
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


       KO_disp = eps1/16.*(betaux(i-2)-4.d0*betaux(i-1)+6.d0*betaux(i)-4.d0*betaux(i+1)+betaux(i+2))


!       RHS2(i) =   capB(i) + KO_disp

       RHS2(i) =   0.d0 !gauge waves
 

       RHS3(i) = 0.d0


    end do




END SUBROUTINE Einstein_2nd_beta
