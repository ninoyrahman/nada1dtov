  MODULE tmp_mdgrid
  USE tmp_mdparam

!******************************************************
! Pedro Montero - 27/11/02
!------------------------------------------------------
! This model defines the grid 
!****************************************************** 


! these wil be read from the parameter file
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    INTEGER,save :: lastiter,outfreq,Space,xsiderivatives,outfreq2D
    INTEGER,save :: InitSpc,boundary,Lapse,Shift
    REAL(KIND=double),save :: dtfact,etapar
    REAL(KIND=DOUBLE),save :: rho_cent_tov,rho_atm,rho_thr
    INTEGER,save ::Azz_adjust,Slicing
    REAL(KIND=DOUBLE),save :: gamma,kp
    INTEGER,save :: irec,rsolver,imethod,xsimethod,gridtype,order,interp,driver
    INTEGER,save :: advect,checkpoint,outfreqmax,eqsymmetry
    INTEGER,save :: nuclear,nospecies
    
    REAL(KIND=DOUBLE),save :: Zcno
    REAL(KIND=DOUBLE) :: ADM_mass_init

    
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  REAL(KIND=double) :: dx,dt

  REAL(KIND=double) :: inv2dx,inv2dt
  REAL(KIND=double) :: inv12dx

  REAL(KIND=double) :: invdx

  REAL(KIND=double) :: invsqdx
  REAL(KIND=double) :: inv12sqdx


  REAL(KIND=double) :: xh(-ghzx+1:nx+ghzx)

  REAL(KIND=double) :: time

  INTEGER :: iter







  ENDMODULE tmp_mdgrid



