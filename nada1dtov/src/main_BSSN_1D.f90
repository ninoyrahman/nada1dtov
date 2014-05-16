
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
    USE tmp_mdhydro
    USE tmp_mdsimdata
    USE MD_Boundary
    USE MD_HydroSubroutines
    USE NORM
    USE Data_IO
    IMPLICIT NONE

    INTEGER :: k,ierror,i,point,q,iter1,threads,OMP_GET_NUM_THREADS
    REAL(KIND=double) :: timelapse,testlapse
    INTEGER :: count1,count2,count_rate
    REAL(KIND=double) :: normL2, normInf

!    REAL(KIND=double) :: maxtime
!    INTEGER :: outputtime
!    INTEGER :: hydroOnly



100 FORMAT (2E20.9)


    WRITE(*,*) ''
    WRITE(*,*) '    |-------------------------|'
    WRITE(*,*) '    |  STARTING 1D EVOLUTION  |'
    WRITE(*,*) '    |-------------------------|'
    WRITE(*,*) ''


    call read_simdata


    !  call omp_set_num_threads(8)
    !!$OMP PARALLEL PRIVATE(threads)
    !   threads=OMP_GET_NUM_THREADS()!
    !   write(*,*)'Number of threads for OpnenMP2 = ',threads
    !!$OMP END PARALLEL

    !   call Read_parfile

    call Grid1D

    !   call Init_data
 
    !   call BH_ID
    !   call Schwarzschild_isotropic


!    call TOV_source

!    call Shock

    call Shock_Tube

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

    !    hydroonly = read_hydro_only()
    !    maxtime = read_maxtime()
    !    outputtime = read_outputtime()
    lastiter = 1000000


    write(*,*) 'Total Time:',maxtime
    write(*,*) 'Step per Output:',outputtime
    write(*,*) 'Max Iteration:',lastiter
    write(*,*) 'Hydro Only:',hydroonly
    write(*,*) 'Reconstruction Scheme:', reconstruction_scheme

    do iter=1,lastiter

        time  = time+dt

        if(time .ge. maxtime) exit

        if(hydroonly .eq. 1) then
            call Runge_kutta_2nd_BCs_HydroOnly
        else
            call Runge_kutta_2nd_BCs
        endif

        IF (MOD(iter,outputtime).eq.0) THEN

            call Output_1D
            write(*,*) 'Output time,it:',time,iter

        END IF

        IF (MOD(iter,outputtime).eq.0) THEN

            write(*,*) 'Output_max time,it:',time,iter
            call Output_Max
            call calculate_normL2(normL2)
            call calculate_normInf(normInf)
            call write_norm(iter, normL2, normInf)

        END IF

    enddo


    call write_final_norm(iter, maxtime, normL2, normInf)

    call SYSTEM_CLOCK(count2,count_rate)
    timelapse = dfloat(count2-count1)/dfloat(count_rate)
    print*,'time total', timelapse


    write(*,*) ' *****************************************'
    write(*,*) ' Done last time update'
    write(*,*) ' *****************************************'



END PROGRAM MAIN_BSSN_1D



