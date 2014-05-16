!
! Norm.f90
!
!  Created on: May 8, 2014
!      Author: ninoy
!

MODULE NORM
    IMPLICIT NONE

CONTAINS

    !! calculate L2 norm
    SUBROUTINE calculate_normL2(normL2)

        USE tmp_mdparam
        USE tmp_mdgrid
        USE tmp_mdhydro
        implicit none

        integer i
        REAL(KIND=double), INTENT(OUT):: normL2

        normL2 = 0.0
        call OMP_SET_NUM_THREADS(4)
        !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)
        !$OMP PARALLEL DO SCHEDULE(STATIC) REDUCTION(+:normL2)
        do i=imin,imax
!            print *, "i=", i, "rho=", rho(i)
            if(.not.isnan(rho(i)) .and. .not.isnan(rhoi(i))) then
                normL2 = normL2 + (rho(i) - rhoi(i))*(rho(i) - rhoi(i))
            endif
        enddo
        !$OMP END PARALLEL DO
        !$OMP END PARALLEL

        normL2 = SQRT(normL2/(nx+2.0*ghzx))

    END SUBROUTINE calculate_normL2

    !! calculate Infinite norm
    SUBROUTINE calculate_normInf(normInf)

        USE tmp_mdparam
        USE tmp_mdgrid
        USE tmp_mdhydro
        implicit none

        integer i
        REAL(KIND=double), INTENT(OUT):: normInf

        normInf = 0.0
        call OMP_SET_NUM_THREADS(4)
        !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)
        !$OMP PARALLEL DO SCHEDULE(STATIC) REDUCTION(MAX:normInf)
        do i=imin,imax
            if(.not.isnan(rho(i)) .and. .not.isnan(rhoi(i))) then
                if(abs(rho(i) - rhoi(i)) .gt. normInf) then
                    normInf = abs(rho(i) - rhoi(i))
                endif
            endif
        enddo
        !$OMP END PARALLEL DO
        !$OMP END PARALLEL

    END SUBROUTINE calculate_normInf



END MODULE NORM




