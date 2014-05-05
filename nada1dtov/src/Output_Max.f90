  SUBROUTINE Output_Max
  USE tmp_mdparam
  USE tmp_mdgrid
  USE tmp_mdmetric
  USE tmp_mdhydro
  USE MD_HydroSubroutines

  IMPLICIT NONE
  
  INTEGER :: i,k,ierror
  REAL(kind=double) :: rhomax,maxtmp,funtmp


  OPEN(UNIT=50, FILE='data/max_rho.dat',STATUS='UNKNOWN',ACTION= 'WRITE',POSITION='append',IOSTAT= ierror)
  OPEN(UNIT=51, FILE='data/max_betaux.dat',STATUS='UNKNOWN',ACTION= 'WRITE',POSITION='append',IOSTAT= ierror)
  OPEN(UNIT=52, FILE='data/min_alp.dat',STATUS='UNKNOWN',ACTION= 'WRITE',POSITION='append',IOSTAT= ierror)




100 FORMAT (2E22.11)
110 FORMAT (A9,E22.11)



    write(50,*) time*4.9255408d-6,maxval(rho)!/ 1.28d-3

    write(51,*) time,maxval(betaux)

    write(52,*) time,alp(1),alp(3)

    call l2norm_ham

!    do i=30,50
      CLOSE(unit=50)
      CLOSE(unit=51)
      CLOSE(unit=52)
!    ENDDO

      call AH_Finder

  ENDSUBROUTINE Output_Max











