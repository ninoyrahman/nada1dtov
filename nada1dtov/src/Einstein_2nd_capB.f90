  SUBROUTINE Einstein_2nd_capB(RHS2,RHS3)
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





       RHS2(i) =   0.75d0*RHSdeltar(i) !-  0.1d0/1.6d0*capB(i)



       RHS3(i) = 0.d0


 

 end do



! write(*,*) RHS2(1), RHS3(1)



END SUBROUTINE Einstein_2nd_capB
