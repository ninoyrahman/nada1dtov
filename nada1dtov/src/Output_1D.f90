  SUBROUTINE Output_1D
  USE tmp_mdparam
  USE tmp_mdgrid
  USE tmp_mdmetric
  USE tmp_mdhydro
  USE MD_HydroSubroutines

  IMPLICIT NONE
  
  INTEGER :: i,k,ierror
  REAL(kind=double) :: rhomax,maxtmp,funtmp

  OPEN(UNIT=30, FILE='data/small_a.dat',STATUS='UNKNOWN',ACTION= 'WRITE',POSITION='append',IOSTAT= ierror)
  OPEN(UNIT=31, FILE='data/small_b.dat',STATUS='UNKNOWN',ACTION= 'WRITE',POSITION='append',IOSTAT= ierror)
  OPEN(UNIT=32, FILE='data/alp.dat',STATUS='UNKNOWN',ACTION= 'WRITE',POSITION='append',IOSTAT= ierror)
  OPEN(UNIT=33, FILE='data/betaux.dat',STATUS='UNKNOWN',ACTION= 'WRITE',POSITION='append',IOSTAT= ierror)
  OPEN(UNIT=34, FILE='data/trK.dat',STATUS='UNKNOWN',ACTION= 'WRITE',POSITION='append',IOSTAT= ierror)
  OPEN(UNIT=35, FILE='data/delta_x.dat',STATUS='UNKNOWN',ACTION= 'WRITE',POSITION='append',IOSTAT= ierror)
  OPEN(UNIT=36, FILE='data/A_a.dat',STATUS='UNKNOWN',ACTION= 'WRITE',POSITION='append',IOSTAT= ierror)
  OPEN(UNIT=37, FILE='data/xsi.dat',STATUS='UNKNOWN',ACTION= 'WRITE',POSITION='append',IOSTAT= ierror)

  OPEN(UNIT=38, FILE='data/capB.dat',STATUS='UNKNOWN',ACTION= 'WRITE',POSITION='append',IOSTAT= ierror)

  OPEN(UNIT=39, FILE='data/rho.dat',STATUS='UNKNOWN',ACTION= 'WRITE',POSITION='append',IOSTAT= ierror)
  OPEN(UNIT=40, FILE='data/eps.dat',STATUS='UNKNOWN',ACTION= 'WRITE',POSITION='append',IOSTAT= ierror)
  OPEN(UNIT=41, FILE='data/p.dat',STATUS='UNKNOWN',ACTION= 'WRITE',POSITION='append',IOSTAT= ierror)
  OPEN(UNIT=42, FILE='data/velr.dat',STATUS='UNKNOWN',ACTION= 'WRITE',POSITION='append',IOSTAT= ierror)

  OPEN(UNIT=43, FILE='data/xsi2.dat',STATUS='UNKNOWN',ACTION= 'WRITE',POSITION='append',IOSTAT= ierror)

  OPEN(UNIT=44, FILE='data/lorentz.dat',STATUS='UNKNOWN',ACTION= 'WRITE',POSITION='append',IOSTAT= ierror)

  OPEN(UNIT=45, FILE='data/dens.dat',STATUS='UNKNOWN',ACTION= 'WRITE',POSITION='append',IOSTAT= ierror)
  OPEN(UNIT=46, FILE='data/sr.dat',STATUS='UNKNOWN',ACTION= 'WRITE',POSITION='append',IOSTAT= ierror)
  OPEN(UNIT=47, FILE='data/tau.dat',STATUS='UNKNOWN',ACTION= 'WRITE',POSITION='append',IOSTAT= ierror)




100 FORMAT (2E22.11)
110 FORMAT (A9,E22.11)


     write(30,*)''
     write(30,110)'"time = ',time
    DO i=imin-ghzx,imax+ghzx
       write(30,100) xh(i),small_a(i)
    ENDDO 

     write(31,*)''
     write(31,110)'"time = ',time
    DO i=imin-ghzx,imax+ghzx
       write(31,100) xh(i),small_b(i)
    ENDDO


     write(32,*)''
     write(32,110)'"time = ',time
    DO i=imin-ghzx,imax+ghzx
       write(32,100) xh(i),alp(i)
    ENDDO

     write(33,*)''
     write(33,110)'"time = ',time
    DO i=imin-ghzx,imax+ghzx
       write(33,100) xh(i),betaux(i)
    ENDDO


     write(34,*)''
     write(34,110)'"time = ',time
    DO i=imin-ghzx,imax+ghzx
       write(34,*) xh(i),trK(i)
    ENDDO

     write(35,*)''
     write(35,110)'"time = ',time
    DO i=imin-ghzx,imax+ghzx
       write(35,100) xh(i),delta_x(i)
    ENDDO


     write(36,*)''
     write(36,110)'"time = ',time
    DO i=imin-ghzx,imax+ghzx
       write(36,100) xh(i),A_a(i)
    ENDDO

     write(37,*)''
     write(37,110)'"time = ',time
    DO i=imin-ghzx,imax+ghzx
       write(37,*) xh(i),xsi(i)
    ENDDO


     write(38,*)''
     write(38,110)'"time = ',time
    DO i=imin-ghzx,imax+ghzx
       write(38,*) xh(i),capB(i)
    ENDDO


     write(39,*)''
     write(39,110)'"time = ',time
    DO i=imin-ghzx,imax+ghzx
       write(39,100) xh(i),rho(i)
    ENDDO

     write(40,*)''
     write(40,110)'"time = ',time
    DO i=imin-ghzx,imax+ghzx
       write(40,100) xh(i),eps(i)
    ENDDO

     write(41,*)''
     write(41,110)'"time = ',time
    DO i=imin-ghzx,imax+ghzx
       write(41,100) xh(i),p(i)
    ENDDO


     write(42,*)''
     write(42,110)'"time = ',time
    DO i=imin-ghzx,imax+ghzx
       write(42,100) xh(i),velr(i)
    ENDDO


     write(43,*)''
     write(43,110)'"time = ',time
    DO i=imin-ghzx,imax+ghzx
       write(43,*) xh(i),xsi2(i)
    ENDDO



     write(44,*)''
     write(44,110)'"time = ',time
    DO i=imin-ghzx,imax+ghzx
       write(44,*) xh(i),lorentz(i)
    ENDDO


     write(45,*)''
     write(45,110)'"time = ',time
    DO i=imin-ghzx,imax+ghzx
       write(45,*) xh(i),dens(i)
    ENDDO


     write(46,*)''
     write(46,110)'"time = ',time
    DO i=imin-ghzx,imax+ghzx
       write(46,*) xh(i),srval(i)
    ENDDO

     write(47,*)''
     write(47,110)'"time = ',time
    DO i=imin-ghzx,imax+ghzx
       write(47,*) xh(i),tau(i)
    ENDDO




    do i=30,47
      CLOSE(unit=i)
    ENDDO


    call Constraints

    call BH

   ENDSUBROUTINE Output_1D











