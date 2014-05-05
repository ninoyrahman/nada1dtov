
!*******************************************************************************
!                     -------------------
!                        NAME : BSSN_1D
!                     -------------------
!
!             1 DIMENSIONAL CODE SOLVING THE FULL SET OF EINSTEIN's EQUATIONS
!             IN SPHERICAL COORDINATES
!            -------------------------------------------------------------------
!   AUTHOR : PEDRO J. MONTERO MURIEL
!            -------------------------------------------------------------------
!*******************************************************************************


  PROGRAM main_BSSN_1D

  USE tmp_mdparam
  USE tmp_mdgrid
  USE tmp_mdmetric
  USE MD_Boundary
  USE MD_HydroSubroutines
  IMPLICIT NONE

 
  INTEGER :: k,ierror,i,point,q,iter1,threads,OMP_GET_NUM_THREADS
  REAL(KIND=double) :: timelapse,testlapse
  integer :: count1,count2,count_rate


100 FORMAT (2E20.9)


      WRITE(*,*) ''
      WRITE(*,*) '    |-------------------------|'
      WRITE(*,*) '    |  STARTING 1D EVOLUTION  |'
      WRITE(*,*) '    |-------------------------|'
      WRITE(*,*) ''



!  call omp_set_num_threads(8)
!$OMP PARALLEL PRIVATE(threads)
!   threads=OMP_GET_NUM_THREADS()!
!   write(*,*)'Number of threads for OpnenMP2 = ',threads
!$OMP END PARALLEL

!   call Read_parfile

   call Grid1D

!   call Init_data 
 
!   call BH_ID 
!   call Schwarzschild_isotropic


   call TOV_source

!   call Shock

!   call Gauge_dyn

   call Output_1D
   call Output_Max

 

 
  write(*,*) ' '
  write(*,*) ' *****************************************'
  write(*,*) '  Starting Time Evolution'
  write(*,*) ' *****************************************'
  write(*,*) ' '
  write(*,*) ' '


    write(*,*) 'Time:',time


   write(*,*)'call time'
   call SYSTEM_CLOCK(count1,count_rate)


    lastiter = 50000

    write(*,*) 'Total Time:',dt*lastiter

  do iter=1,lastiter

    time  = time+dt

    call Runge_kutta_2nd_BCs

    IF (MOD(iter,100).eq.0) THEN

      call Output_1D
      write(*,*) 'Output time,it:',time,iter

    END IF

    IF (MOD(iter,100).eq.0) THEN

      write(*,*) 'Output_max time,it:',time,iter

      call Output_Max


    END IF

   enddo



   call SYSTEM_CLOCK(count2,count_rate)
   timelapse = dfloat(count2-count1)/dfloat(count_rate)
   print*,'time total', timelapse
   

   write(*,*) ' *****************************************'
   write(*,*) ' Done last time update'
   write(*,*) ' *****************************************'

 END PROGRAM MAIN_BSSN_1D



