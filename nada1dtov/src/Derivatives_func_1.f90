
  SUBROUTINE Derivatives_func_1(func,dxdfunc)
  USE tmp_mdparam
  USE tmp_mdgrid

  IMPLICIT NONE

  INTEGER :: i,k

  REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx), INTENT(IN) :: func
  REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx), INTENT(OUT) :: dxdfunc



 !$OMP PARALLEL
 !$OMP DO PRIVATE(i)

    do i=imin,imax

      dxdfunc(i) = (func(i-2)-8.d0*func(i-1)+ 8.d0*func(i+1)-func(i+2))*inv12dx

    end do


 !$OMP END DO
 !$OMP END PARALLEL


  END SUBROUTINE Derivatives_func_1

