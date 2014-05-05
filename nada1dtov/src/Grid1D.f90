  SUBROUTINE Grid1D
  USE tmp_mdgrid
  USE tmp_mdparam

  IMPLICIT NONE

  INTEGER :: i
  REAL(KIND=double) :: convtime
  REAL(KIND=double), DIMENSION(1:5) :: xa


    time = 0.d0
  
    convtime= 4.9255408d-6

!    dx = 0.1 !gauge waves
 
    dx = 0.2
 
!    dtfact =  0.5*1.5703926303627d-2!0.6

    dtfact =  0.6

    write(*,*) 'dtfact',dtfact

!    stop

    gridtype=1

   dt  = dtfact*dx

!Define some parameters 

  inv2dx = 1.d0 / (2.d0*dx) 

  inv12dx = 1.d0 / (12.d0*dx) 

  invdx = 1.d0 / dx

  invsqdx = 1.d0 / (dx*dx) 

  inv12sqdx = 1.d0 / (12.d0*dx*dx) 


! if(gridtype.eq.1) then 

   xh(imin-ghzx) = -(3.d0*dx/2.d0)-dx

! else if(gridtype.eq.0) then 

!   xh(imin-ghzx) = -2.*dx

! end if


  DO i=imin-ghzx+1,imax+ghzx
    xh(i) = xh(i-1)+dx
 ENDDO  


        write(*,*) ''
        write(*,*) 'GRID1D---->:nx',nx
        write(*,*) 'GRID1D---->:dx',dx,dt
        write(*,*) ''
        write(*,*) ''
        write(*,*) ''
        write(*,*) 'GRID1D---->:xmin',xh(imin)
        write(*,*) 'GRID1D---->:xmax',xh(imax)

        write(*,*) ''




 ENDSUBROUTINE Grid1D









