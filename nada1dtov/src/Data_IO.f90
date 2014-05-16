!
! Data_IO.f90
!
!  Created on: May 9, 2014
!      Author: ninoy
!

MODULE Data_IO
    USE tmp_mdparam
    USE tmp_mdgrid
    USE tmp_mdsimdata
    IMPLICIT NONE

CONTAINS

    !! read and return time step size factor from "parameter.dat" file
    !! time step size = dtfact*space step size
    REAL(KIND=double) function read_dtfact()

        implicit none
        INTEGER i, status
        REAL :: tmp_dtfact
        CHARACTER(LEN=300) :: line

        open(unit=1000, FILE='parameter.dat', STATUS='UNKNOWN', ACTION= 'READ')

        do
            read(1000,*,IOSTAT=status) line
            if(IS_IOSTAT_END(status)) then
                exit
            endif
            if(line .eq. "#dtfact") then
                read(1000,*) tmp_dtfact
                exit
            endif
        enddo

        close(unit=1000)

        read_dtfact = tmp_dtfact

    end function read_dtfact


    INTEGER function read_nx()

        implicit none
        INTEGER :: tmp_nx

        open(unit=1001, FILE='parameter.dat', STATUS='UNKNOWN', ACTION= 'READ')

        read(1001,*)
        read(1001,*)
        read(1001,*)
        read(1001,*) tmp_nx

        close(unit=1001)

        read_nx = tmp_nx

    end function read_nx


    !! read and return end time of simulation from "parameter.dat" file
    REAL(KIND=double) function read_maxtime()

        implicit none
        INTEGER i, status
        REAL :: tmp_maxtime
        CHARACTER(LEN=300) :: line

        open(unit=1002, FILE='parameter.dat', STATUS='UNKNOWN', ACTION= 'READ')

        do
            read(1002,*,IOSTAT=status) line
            if(IS_IOSTAT_END(status)) then
                exit
            endif
            if(line .eq. "#maxtime") then
                read(1002,*) tmp_maxtime
                exit
            endif
        enddo
        close(unit=1002)

        read_maxtime = tmp_maxtime

    end function read_maxtime

    INTEGER function read_outputtime()

        implicit none
        INTEGER i, status
        INTEGER :: tmp_outputtime
        CHARACTER(LEN=300) :: line

        open(unit=1003, FILE='parameter.dat', STATUS='UNKNOWN', ACTION= 'READ')

        do
            read(1003,*,IOSTAT=status) line
            if(IS_IOSTAT_END(status)) then
                exit
            endif
            if(line .eq. "#outputtime") then
                read(1003,*) tmp_outputtime
                exit
            endif
        enddo
        close(unit=1003)

        read_outputtime = tmp_outputtime

    end function read_outputtime

    INTEGER function read_reconstruction_scheme()

        implicit none
        INTEGER i, status
        INTEGER :: tmp_reconstruction_scheme
        CHARACTER(LEN=300) :: line

        open(unit=1004, FILE='parameter.dat', STATUS='UNKNOWN', ACTION= 'READ')

        do
            read(1004,*,IOSTAT=status) line
            if(IS_IOSTAT_END(status)) then
                exit
            endif
            if(line .eq. "#reconstruction") then
                read(1004,*) tmp_reconstruction_scheme
                exit
            endif
        enddo
        close(unit=1004)

        read_reconstruction_scheme = tmp_reconstruction_scheme

    end function read_reconstruction_scheme

    INTEGER function read_hydro_only()

        implicit none
        INTEGER i, status
        INTEGER :: tmp_hydro_only
        CHARACTER(LEN=300) :: line

        open(unit=1005, FILE='parameter.dat', STATUS='UNKNOWN', ACTION= 'READ')

        do
            read(1005,*,IOSTAT=status) line
            if(IS_IOSTAT_END(status)) then
                exit
            endif
            if(line .eq. "#HydroOnly") then
                read(1005,*) tmp_hydro_only
                exit
            endif
        enddo
        close(unit=1005)

        read_hydro_only = tmp_hydro_only

    end function read_hydro_only


    SUBROUTINE read_simdata()

        implicit none

        maxtime = read_maxtime()
        hydroonly = read_hydro_only()
        outputtime = read_outputtime()
        reconstruction_scheme = read_reconstruction_scheme()
        dtfact =  read_dtfact()

    END SUBROUTINE read_simdata


    SUBROUTINE write_norm(iter, normL2, normInf)

        implicit none

        INTEGER, INTENT(IN) :: iter
        REAL(KIND=double), INTENT(IN) :: normL2, normInf

        open(unit=2000, FILE='data/error.dat', STATUS='UNKNOWN', ACTION= 'WRITE',POSITION='append')

        write(2000,*) iter,normL2,normInf

        close(unit=2000)

    END SUBROUTINE write_norm


    SUBROUTINE write_final_norm(maxiter, maxtime, normL2, normInf)

        implicit none

        INTEGER , INTENT(IN) :: maxiter
        REAL(KIND=double), INTENT(IN) :: maxtime, normL2, normInf
        CHARACTER(LEN=300) :: line
        CHARACTER(LEN=300) :: FMT
        INTEGER :: reconstruction_scheme

        reconstruction_scheme = read_reconstruction_scheme()
        FMT = "(A5,A1,F5.2,A1,F7.4,A1,I10,A1,F12.2,A1,E14.4E3,A1,E14.4E3)"

        open(unit=2001, FILE='final_error.dat', STATUS='UNKNOWN', ACTION= 'WRITE',POSITION='append')

        if(reconstruction_scheme .eq. 1) then
            write(2001,FMT) "PPM",",",dtfact,",",dx,",",maxiter,",",maxtime,",",normL2,",",normInf
        else
            write(2001,FMT) "MC",",",dtfact,",",dx,",",maxiter,",",maxtime,",",normL2,",",normInf
        endif

        close(unit=2001)

    END SUBROUTINE write_final_norm



END MODULE Data_IO




